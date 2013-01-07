#include "MSBlast.h"
#include "auxfun.h"
#include "DeNovoSolutions.h"

//-model CID_IT_TRYP -msb_id_file MSB\pks_prots.txt -msb_keep_file MSB\keep.txt -msb_exclude_file MSB\exclude.txt -msb_query_name MSB\SPNF\SPNF -msb_map_file MSB\nrdb95_map.txt -msb_examine_results MSB/SPBOTH/SPBOTH_MSB_raw2a.txt

bool getScanStringFromTitle(const string& title, string& scanStr)
{
	const size_t len = title.length();
	if (len<=7)
		return false;

	if (title[len-1] == 'a' && title[len-2] == 't' && title[len-3]=='d' && 
		title[len-6] == '.' && title[len-4]== '.')
	{
		size_t pos = len-7;
		int numDots = 0;
		while (pos>0 && numDots<2)
		{
			--pos;
			if (title[pos] == '.')
				++numDots;
		}
		if (numDots == 2)
		{
			string scanStr = title.substr(pos+1,len-7-pos);
			return true;
		}
	}
	return false;
}


void MSBlastCollector::generateMsBlastSequences(const char* msb_name,
												const vector<string>& list_vector, 
												AllScoreModels *model,
												float min_filter_prob,
												size_t maxSequencesPerSet,
												bool indOutputCumulativeProbs)
{
	config_ = model->get_config();
	PeptideRankScorer *rank_model = model->get_rank_model_ptr(1);
	PrmGraph *prm_ptr = NULL;
	vector<PrmGraph *> prm_ptrs;

	int total_benchmark = 0;
	int correct_benchmark = 0;
	double totalScore = 0.0;

	const string dnv_name = msb_name + std::string("_dnv.txt");
	const string full_name = msb_name + std::string("_full.txt");
	ofstream dnv_stream(dnv_name.c_str());
	ofstream full_stream(full_name.c_str());

	if (! dnv_stream.good() || ! full_stream.good())
		error("problem creating msb files: ",msb_name);

	for (size_t fileIdx=0; fileIdx<list_vector.size(); fileIdx++) 
	{
		const char *spectra_file = list_vector[fileIdx].c_str();

		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(spectra_file, config_);

		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();


		cout << "Processing: " << fileIdx << "  " << spectra_file << " ("
			 << sl.getNumHeaders() << " spectra)" << endl;

		if (sl.getNumHeaders()<=0)
			continue;

		for (size_t sc=0; sc<sl.getNumHeaders(); sc++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(sc);

			AnnotatedSpectrum as;
			if (! as.readSpectrum(sa, header))
			{
				SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
				nonConstHeader->setSpectraFileIndexInList(fileIdx);
				header->printStats(dnv_stream);
				dnv_stream << "#could not read spectrum correctly..." << endl;
				cout << sc << "\t";
				header->printStats(cout, false);
				cout << " #could not read spectrum correctly..." << endl;
				continue;
			}

			SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
			nonConstHeader->setSpectraFileIndexInList(fileIdx);
			
			header->printStats(dnv_stream, false);
			cout << sc << "\t";
			header->printStats(cout, false);
			
			if (as.getNumPeaks()<5)
			{
				dnv_stream << endl << "# too few peaks..." << endl;
				cout << " # too few peaks..." << endl;
				continue;
			}

			vector<SeqPath> solutions;
			solutions.clear();
			if ( as.getCharge() >= model->get_max_score_model_charge())
			{
				dnv_stream << endl << "# Charge " << as.getCharge() << " not supported yet..." << endl << endl;
				cout << " # Charge " << as.getCharge() << " not supported yet..." << endl;
				continue;
			}

			bool perform_rerank=false;
			int rerank_size_idx = NEG_INF;
			int num_sols_to_generate_before_ranking=300;
			float spectrum_quality = 0.0;
		
			if (1)
			{
				vector<mass_t> pms_with_19;
				vector<int>    charges;
				pms_with_19.clear();
				charges.clear();		
			
				// output m/z and prob values for the different charge states
				vector<PmcSqsChargeRes> pmcsqs_res;
				model->selectPrecursorMassesAndCharges(config_, as, pms_with_19, charges, &pmcsqs_res);
				
				if (pmcsqs_res.size()>charges[0])
				{
					const float sqs = pmcsqs_res[charges[0]].sqs_prob;
					if (sqs<min_filter_prob)
					{
						dnv_stream << endl << "# low quality, skipping: " << sqs << endl << endl;
						cout	   << " # low quality, skipping: " << sqs << endl;
						continue;
					}
					else
						dnv_stream << "\tSQS: " << sqs << endl;
				}

				if (pmcsqs_res.size()>charges[0])
					spectrum_quality = pmcsqs_res[charges[0]].sqs_prob;

				if (prm_ptrs.size()<pms_with_19.size())
					prm_ptrs.resize(pms_with_19.size(),NULL);

				if (pms_with_19[0]<100.0)
				{
					dnv_stream << endl << "# Could not process spectrum..." << endl << endl;
					cout << " # Could not process spectrum..." << endl;
					continue;
				}
				
				const int num_rerank_sols_per_charge[] = {0,300,300,300,100,100,100,100,100,100,100,100,100,100};
				const int num_rerank_per_charge =  num_rerank_sols_per_charge[charges[0]];
				if (rank_model && num_sols_to_generate_before_ranking<num_rerank_per_charge)
						num_sols_to_generate_before_ranking=num_rerank_per_charge;
				
				generate_denovo_solutions_from_several_pms(
					prm_ptrs,
					model,
					&as,
					true, 
					num_sols_to_generate_before_ranking,
					MIN_MSB_DENOVO_SEQ_LENGTH,
					MAX_MSB_DENOVO_SEQ_LENGTH,
					pms_with_19,
					charges,
					solutions,
					false);

				// use charge of top scoring solution
				if (solutions.size()>1)
				{
					const int sol_charge = solutions[0].charge;
					for (size_t j=1; j<solutions.size(); j++)
					{
						if (solutions[j].charge != sol_charge)
						{
							if (j==solutions.size()-1)
							{
								solutions.pop_back();
							}
							else
							{
								solutions[j]=solutions[solutions.size()-1];
								solutions.pop_back();
								j--;
							}
						}
					}
				}
			}

			if (rank_model && solutions.size()>0)
			{
				rerank_size_idx = config_->calc_size_idx(solutions[0].charge,solutions[0].pm_with_19);
				if (rank_model->get_ind_part_model_was_initialized(solutions[0].charge,rerank_size_idx))
					perform_rerank=true;
			}

			vector<score_pair> score_idx_pairs;
			if (perform_rerank)
			{
				rank_model->scoreDenovoSequences(solutions, as, score_idx_pairs, rerank_size_idx);
				if (score_idx_pairs.size() == 0)
				{
					dnv_stream << endl << "# Could not process spectrum..." << endl << endl;
					cout << " # Could not process spectrum..." << endl;
					continue;
				}
				for (size_t i=0; i<score_idx_pairs.size(); i++)
					solutions[score_idx_pairs[i].idx].rerank_score = score_idx_pairs[i].score;
				sort(score_idx_pairs.begin(),score_idx_pairs.end());
			}
			else
			{
				score_idx_pairs.resize(solutions.size());
				for (size_t i=0; i<solutions.size(); i++)
					score_idx_pairs[i].idx=i;
			}

			// for debug
			bool had_pep = false;
			bool had_correct = false;

			if (solutions.size() == 0)
			{
				dnv_stream << endl << "# No solutions found." << endl << endl;
				cout << " # No solutions found." << endl;
			}
			else 
			{
				dnv_stream << "#Index\t";
				dnv_stream << "RnkScr\t";
				if (indOutputCumulativeProbs)
					dnv_stream << "CumProb\t";

				dnv_stream << "PnvScr\tN-Gap\tC-Gap\t[M+H]\tCharge\tSequence" << endl;

				if ( indOutputCumulativeProbs)
				{
					for (size_t i=0; i<solutions.size() && i<maxSequencesPerSet; i++)
					{
						const int idx = (perform_rerank ? score_idx_pairs[i].idx : i);
						const vector<PathPos>& positions = solutions[idx].positions;
						solutions[idx].prm_ptr->calc_amino_acid_probs(solutions[idx],i);
					}
					const int first_sol_charge = solutions[0].charge;
					const int first_sol_size_idx = config_->calc_size_idx(first_sol_charge,solutions[0].pm_with_19);
					CumulativeSeqProbModel* csp_model = (CumulativeSeqProbModel* )model->get_cumulative_seq_prob_model_ptr(0);
			
					csp_model->calc_cumulative_seq_probs(first_sol_charge, first_sol_size_idx, 
						spectrum_quality, score_idx_pairs, solutions); 
				}

				for (size_t i=0; i<solutions.size() && i<maxSequencesPerSet; i++) 
				{
					const int idx = (perform_rerank ? score_idx_pairs[i].idx : i);
					mass_t c_gap=solutions[idx].pm_with_19 - solutions[idx].c_term_mass;
					if (c_gap<24.0)
						c_gap = 0;

					dnv_stream << setprecision(3) << fixed << i << "\t";
					if (perform_rerank)
					{
						dnv_stream << score_idx_pairs[i].score << "\t";
					}
					else
						dnv_stream << -999 << "\t";

					if (indOutputCumulativeProbs)
						dnv_stream << solutions[idx].cumulative_seq_prob << "\t";

					dnv_stream << solutions[idx].path_score << "\t";
					dnv_stream << solutions[idx].n_term_mass << "\t";
					dnv_stream << c_gap << "\t";
					dnv_stream << solutions[idx].pm_with_19 << "\t";
					dnv_stream << solutions[idx].charge << "\t";
					dnv_stream << solutions[idx].seq_str;	

					if (header->getPeptideStr().length() >2)
					{
						if (solutions[idx].check_if_correct(header->getPeptideStr(),config_))
						{
							dnv_stream << " *";
							if (! had_correct)
							{
								correct_benchmark++;
								had_correct=true;
							}
						}
						had_pep=true;
					}
					dnv_stream << endl;
				}
			}
			if (had_pep) // for annotated spectra (benchmark)
				total_benchmark++;

			dnv_stream << endl;

			// create MSBlast seqeunces
			MSBSequenceSet seqSet;
			seqSet.convertSeqPathsToMSBSequences(config_, solutions, maxSequencesPerSet);

			string msbLine;
			seqSet.createMSBFullLine(msbLine);
			if (msbLine.length()>4)
				full_stream << msbLine << endl;
			if (seqSet.getMSBSequences().size()>0)
			{
				cout << "\t" << seqSet.getMSBSequences().size() <<" \t" << seqSet.getMSBSequences()[0].msbScore;
			}
			else
				cout << "\t" << 0;
			
			addToExistingSets(seqSet, maxSequencesPerSet*3);

			if (had_pep)
			{
				const string& pep = header->getPeptideStr();
				const vector<MSBSequence>& msbs = seqSet.getMSBSequences();
				float maxScore=0.0;
				for (size_t i=0; i<msbs.size(); i++)
				{
					float score = msbs[i].calcMatchScore(config_, pep);
					if (score > maxScore)
						maxScore = score;
				}
				totalScore += maxScore;
				cout << " s: " << maxScore;
			}
			if (had_correct)
				cout << " *";

			cout << endl;
		}
	}

	dnv_stream.close();
	full_stream.close();

	if (totalScore != 0.0)
		cout << "TotalScore = " << totalScore << endl;

	cout << endl << "Done..." << endl;
	cout << "Created" << dnv_name << endl;
	cout << "Created " << full_name << endl;

	/////////////////////////////////////////////////////////////////
	// this part works only if the spectra are annotated (benchmark)
	/////////////////////////////////////////////////////////////////
	if (total_benchmark>0)
	{
		cout << "#Correct spectra " << correct_benchmark << "/" << total_benchmark << " (" <<
			fixed << setprecision(3) << (double)correct_benchmark/(double)total_benchmark << ")" << endl;
	}

}


bool compMSBStringScore(const MSBString& lhs, const MSBString& rhs)
{
	return (lhs.msbScore > rhs.msbScore);
}


bool compMsbStringSet(const MSBStringSet& lhs, const MSBStringSet& rhs)
{
	if (lhs.indMarkedForKeep && ! rhs.indMarkedForKeep)
		return true;

	if (rhs.indMarkedForKeep && ! lhs.indMarkedForKeep)
		return false;

	if (rhs.indMarkedForExclusion && ! lhs.indMarkedForExclusion)
		return true;

	if (lhs.indMarkedForExclusion && ! rhs.indMarkedForExclusion)
		return false;

	if (lhs.indFromMainFullFile && ! rhs.indFromMainFullFile)
		return true;

	if (rhs.indFromMainFullFile && ! lhs.indFromMainFullFile)
		return false;

	// if they have the same "standing" go by the score
	return (lhs.strings.size()>0 && (rhs.strings.size() == 0 || lhs.strings[0].msbScore > rhs.strings[0].msbScore));
}

void MSBlastCollector::writeMsBlastQuery(const char* name, size_t maxQuerySize, float minMsbScore)
{
	for (size_t i=0; i<msbStringSets.size(); i++)
		sort(msbStringSets[i].strings.begin(), msbStringSets[i].strings.end(), compMSBStringScore);

	sort(msbStringSets.begin(), msbStringSets.end(), compMsbStringSet);

	string query_file = name + std::string("_query.txt");
	ofstream query_stream(query_file.c_str());

	if (! query_stream.good())
		error("could not open query file for writing: ",query_file.c_str());

	cout << "Creating " << query_file;
	cout.flush();

	size_t usedSize =0;
	size_t idx=0;
	while (idx<msbStringSets.size() && usedSize + 256 <maxQuerySize)
	{
		vector<MSBString>& vec = msbStringSets[idx].strings;
		if (vec[0].msbScore < minMsbScore)
			break;

		ostringstream oss;
		oss << vec[0].fileIndex << "\t" << vec[0].scanNumber << "\t" << setprecision(3) << fixed << vec[0].mz << "\t" << vec.size() << "\t" << vec[0].msbScore << "\t";
		for (size_t i=0; i<vec.size(); i++)
		{
			oss << "-" << vec[i].seqStr;
		}

		const string line = oss.str();
		size_t numChars=0;
		for (size_t i=0; i<line.length(); i++)
			if (line[i]>='A' && line[i]<='Z')
				numChars++;

		usedSize += numChars;
		query_stream << line << endl;
		idx++;
	}
	query_stream.close();
	cout << " ...Done" << endl;
	cout << "Wrote " << idx << " sets of sequences (" << usedSize << " chars)" << endl;
}


void MSBlastCollector::generateMsBlastQueryWithHistory(const char* name,  
													   map<string,int>& peptideHistory,
													   size_t maxQuerySize, float minMsbScore)
{
	for (size_t i=0; i<msbStringSets.size(); i++)
		sort(msbStringSets[i].strings.begin(), msbStringSets[i].strings.end(), compMSBStringScore);

	sort(msbStringSets.begin(), msbStringSets.end(), compMsbStringSet);

	string query_file = name + std::string("_query.txt");
	ofstream query_stream(query_file.c_str());

	if (! query_stream.good())
		error("could not open query file for writing: ",query_file.c_str());

	cout << "Creating " << query_file;
	cout.flush();

	size_t numKept=0, numExcluded=0, numLinesWritten=0;
	size_t numFromMain = 0;
	size_t usedSize =0;
	size_t idx=0;
	while (idx<msbStringSets.size() && usedSize + 256 <maxQuerySize)
	{
		MSBStringSet&       peptideSet = msbStringSets[idx++];
		vector<MSBString>& vec = peptideSet.strings;
		if (vec[0].msbScore < minMsbScore)
			break;

		ostringstream oss;
		oss << vec[0].fileIndex << "\t" << vec[0].scanNumber << "\t" << setprecision(3) << fixed << vec[0].mz << "\t" << vec[0].msbScore << "\t";

		// if marked for keep, all peptides will be in the query
		if (peptideSet.indMarkedForKeep)
		{
			oss << "\t" << vec.size() << "\t" << vec[0].msbScore << "\t";
			for (size_t i=0; i<vec.size(); i++)
			{
				peptideHistory[vec[i].seqStr]++;
				oss << "-" << vec[i].seqStr;
			}
			numKept++;
		}
		else if (peptideSet.indMarkedForExclusion)
		{
			numExcluded++;
			continue;
		}
		else // if not marked for keep, we need to choose which peptides will be included; this depends on which ones already hit, and which ones are
			 // in the 
		{	
			vector<int>  histCounts(vec.size(),0);
			vector<bool> falgs(vec.size(),false);
			int numWritten=0;
			for (size_t i=0; i<vec.size(); i++)
			{
				map<string,int>::iterator it = peptideHistory.find(vec[i].seqStr);
				if (it != peptideHistory.end())
					histCounts[i]=it->second;

				if (vec[i].indFoundInMap)
				{
					oss << "-" << vec[i].seqStr;
					peptideHistory[vec[i].seqStr]++;
					numWritten++;
				}
			}
	
			int maxSeqs=3;
			if (vec[0].msbScore>7.0)
				maxSeqs=4;
			if (vec[0].msbScore>8.0)
				maxSeqs=5;
			if (vec[0].msbScore>9.0)
				maxSeqs=20;

			for (size_t i=0; i<vec.size() && numWritten<maxSeqs; i++)
			{
				if (histCounts[i] == 0 && ! vec[i].indFoundInMap)
				{
					oss << "-" << vec[i].seqStr;
					peptideHistory[vec[i].seqStr]++;
					numWritten++;
				}
			}

			if (numWritten == 0)
			{
				for (size_t i=0; i<vec.size(); i++)
				{
					if (histCounts[i] == 1 && ! vec[i].indFoundInMap)
					{
						oss << "-" << vec[i].seqStr;
						peptideHistory[vec[i].seqStr]++;
						numWritten++;
						break;
					}
				}
			}
			if (numWritten == 0)
				continue;
		}

		const string line = oss.str();
		size_t numChars=0;
		for (size_t i=0; i<line.length(); i++)
			if (line[i]>='A' && line[i]<='Z')
				numChars++;

		usedSize += numChars;
		query_stream << line << endl;

		numLinesWritten++;
		if (peptideSet.indFromMainFullFile)
			numFromMain++;
	}
	query_stream.close();
	cout << " ...Done" << endl;
	cout << "Wrote " << numLinesWritten << " sets of sequences (" << usedSize << " chars), " << numFromMain << " where from the main results file." << endl;
	cout << "Scanned a total of " << idx << " lines, " << numKept << " were marked for keeping and " << numExcluded << " were marked for exclusion." << endl;
}


void findMinScoreAndPos(const vector<MSBString>& vec, float& score, size_t& pos)
{
	if (vec.size() == 0)
	{
		score = -9999.0;
		pos   = 0;
		return;
	}

	score = vec[0].msbScore;
	pos	  = 0;

	for (size_t i=1; i<vec.size(); i++)
	{
		if (vec[i].msbScore<score)
		{
			score = vec[i].msbScore;
			pos = i;
		}
	}
}


void MSBlastCollector::addToExistingSets(MSBSequenceSet& sequenceSet, size_t maxSequencesPerSet)
{
	const vector<MSBSequence>& sequences = sequenceSet.getMSBSequences();
	vector<MSBString> msbStrings(sequences.size());

	for (size_t i=0; i<sequences.size(); i++)
	{
		msbStrings[i].seqStr      = sequences[i].makeSeqString(config_);
		msbStrings[i].msbScore    = sequences[i].msbScore;
		msbStrings[i].fileIndex   = sequences[i].fileIndex;
		msbStrings[i].mz	      = sequences[i].mz;
		msbStrings[i].scanNumber  = sequences[i].scanNumber;
	}

	size_t setIdx = MAX_SIZE_T;
	for (size_t i=0; i<msbStrings.size(); i++)
	{
		map<string,size_t>::iterator it = coveredStrings.find(msbStrings[i].seqStr);
		if (it != coveredStrings.end())
		{
			setIdx = it->second;
			break;
		}

		// try suffixes of string if 
		if (msbStrings[i].seqStr[0] == 'B' && msbStrings[i].seqStr.length()>10)
		{
			const size_t lastStartIdx = msbStrings[i].seqStr.length()-8;
			for (size_t j=0; j<lastStartIdx; j++)
			{
				map<string,size_t>::iterator it = coveredStrings.find(msbStrings[i].seqStr.substr(j));
				if (it != coveredStrings.end())
				{
					setIdx = it->second;
					break;
				}
			}
			if (setIdx != MAX_SIZE_T)
				break;
		}
	}


	if (setIdx == MAX_SIZE_T)
	{
		setIdx = msbStringSets.size();
		msbStringSets.resize(setIdx+1);
	}

	vector<MSBString>& vec = msbStringSets[setIdx].strings;

	// add to vector
	for (size_t i=0; i<msbStrings.size(); i++)
	{
		float  minScore=0.0;
		size_t minPos=0;

		if (vec.size() >= maxSequencesPerSet)
			findMinScoreAndPos(vec, minScore, minPos);

		for (size_t i=0; i<msbStrings.size(); i++)
		{
			const MSBString& newStr = msbStrings[i];
			if (vec.size()>=maxSequencesPerSet && newStr.msbScore <= minScore)
				continue;
			size_t j;
			for (j=0; j<vec.size(); j++)
				if (vec[j] == newStr)
					break;

			if (j<vec.size())
			{
				if (vec[j].msbScore < newStr.msbScore)
					vec[j] = newStr; // replace with higher score source
				continue;
			}

			coveredStrings[newStr.seqStr]=setIdx;

			if (vec.size()<maxSequencesPerSet)
			{
				vec.push_back(newStr);
			}
			else
				vec[minPos] = newStr;

			if (vec.size() >= maxSequencesPerSet)
				findMinScoreAndPos(vec, minScore, minPos);
		}
	}
}


// this function is a bit heuristic; it makes a set and gives arbitrary
// scores (according to the highest score recorded in the full.txt file)
void MSBlastCollector::addToExistingSets(const char* line, size_t maxSequencesPerSet, float minMsbScore)
{
	// parse line
	istringstream iss(line);
	int fileIdx=-1, scan=-1, n=0;
	float mz=-1.0, score = -1.0;
	string rest;

	iss >> fileIdx >> scan >> mz >> n >> score >> rest;
	if (fileIdx<0 || scan<0 || mz < 0.0)
	{
		cout << "Bad line: " << line << endl;
		return;
	}

	if (n<=0 || score <minMsbScore)
		return;

	assert(n<1000);
	vector<MSBString> msbStrings(n);

	for (size_t i=0; i<rest.length(); i++)
		if (rest[i] == '-')
			rest[i] = ' ';

	iss.clear();
	iss.str(rest);
	for (int i=0; i<n; i++)
	{
		iss >> msbStrings[i].seqStr;
		msbStrings[i].msbScore    = score - i*0.2;
		msbStrings[i].fileIndex   = fileIdx;
		msbStrings[i].mz	      = mz;
		msbStrings[i].scanNumber  = scan;
	}

	size_t setIdx = MAX_SIZE_T;
	for (size_t i=0; i<msbStrings.size(); i++)
	{
		map<string,size_t>::iterator it = coveredStrings.find(msbStrings[i].seqStr);
		if (it != coveredStrings.end())
		{
			setIdx = it->second;
			break;
		}

		// try suffixes of string if 
		if (msbStrings[i].seqStr[0] == 'B' && msbStrings[i].seqStr.length()>10)
		{
			const size_t lastStartIdx = msbStrings[i].seqStr.length()-8;
			for (size_t j=0; j<lastStartIdx; j++)
			{
				map<string,size_t>::iterator it = coveredStrings.find(msbStrings[i].seqStr.substr(j));
				if (it != coveredStrings.end())
				{
					setIdx = it->second;
					break;
				}
			}
			if (setIdx != MAX_SIZE_T)
				break;
		}
	}


	if (setIdx == MAX_SIZE_T)
	{
		setIdx = msbStringSets.size();
		msbStringSets.resize(setIdx+1);
		msbStringSets[setIdx].indFromMainFullFile = false;
	}

	vector<MSBString>& vec = msbStringSets[setIdx].strings;

	// add to vector
	for (size_t i=0; i<msbStrings.size(); i++)
	{
		float  minScore=0.0;
		size_t minPos=0;

		if (vec.size() >= maxSequencesPerSet)
			findMinScoreAndPos(vec, minScore, minPos);

		for (size_t k=0; k<msbStrings.size(); k++)
		{
			const MSBString& newStr = msbStrings[k];
			if (vec.size()>=maxSequencesPerSet && newStr.msbScore <= minScore)
				continue;
			size_t j;
			for (j=0; j<vec.size(); j++)
				if (vec[j].seqStr == newStr.seqStr)
					break;

			if (j<vec.size())
			{
				if (vec[j].msbScore < newStr.msbScore)
					vec[j] = newStr; // replace with higher score source
				continue;
			}

			coveredStrings[newStr.seqStr]=setIdx;

			if (vec.size()<maxSequencesPerSet)
			{
				vec.push_back(newStr);
			}
			else
				vec[minPos] = newStr;

			if (vec.size() >= maxSequencesPerSet)
				findMinScoreAndPos(vec, minScore, minPos);
		}
	}
}




void MSBlastCollector::readFullFileIntoCollector(const char* file, size_t maxSequencesPerSet, float minMsbScore)
{
	ifstream ifs(file);
	if (! ifs.good())
		error("Could not open file for reading: ",file);

	char buffer[1024];
	while (! ifs.eof())
	{
		ifs.getline(buffer,1024);
		if (ifs.gcount()>7)
			addToExistingSets(buffer, maxSequencesPerSet, minMsbScore);
	}

	ifs.close();
}



// Takes a bunch of MS-Blast "_full.txt" files and creates a single one
// good for merging parallel jobs or adding to an existing query
void MSBlastCollector::mergeAndSplitQueries(const char* msb_name,
							   const vector<string>& list_vector,
							   size_t maxSequencesPerSet,
							   size_t maxQuerySize, 
							   float minMsbScore)
{
	for (size_t i=0; i<list_vector.size(); i++)
	{
		ifstream ifs(list_vector[i].c_str());
		if (! ifs.good())
		{
			cout << "Warning: could not open file " << list_vector[i] << endl;
			continue;
		}

		char buffer[1024];
		while (! ifs.eof())
		{
			ifs.getline(buffer,1024);
			if (ifs.gcount() < 7)
				continue;

			addToExistingSets(buffer);
		}
	}

	for (size_t i=0; i<msbStringSets.size(); i++)
		sort(msbStringSets[i].strings.begin(), msbStringSets[i].strings.end(), compMSBStringScore);


	cout.flush();

	size_t totalSize = 0;
	size_t usedSize = 0;
	size_t setIdx=0;
	size_t fileIdx=0;
	ofstream full_stream;
	string   full_file;
	while (setIdx<msbStringSets.size())
	{
		vector<MSBString>& vec = msbStringSets[setIdx++].strings;
		if (vec[0].msbScore<minMsbScore)
			continue;

		ostringstream oss;
		oss << vec[0].fileIndex << "\t" << vec[0].scanNumber << "\t" << setprecision(3) << fixed << vec[0].mz << "\t" << vec.size() << "\t" << vec[0].msbScore << "\t";
		for (size_t i=0; i<vec.size() && i<maxSequencesPerSet; i++)
			oss << "-" << vec[i].seqStr;

		const string& line = oss.str();
		size_t numChars=0;
		for (size_t i=0; i<line.length(); i++)
			if (line[i]>='A' && line[i]<='Z')
				numChars++;

		if (usedSize == 0)
		{
			ostringstream ossf;
			ossf << msb_name << "_query_pt_" << fileIdx << ".txt";
			full_file = ossf.str();
			fileIdx++;
			if (full_stream.is_open())
				full_stream.close();
			full_stream.open(full_file.c_str());
			if (! full_stream.good())
				error("Could not open file for writing: ",full_file.c_str());
		}

		usedSize += numChars;
		full_stream << line << endl;
		if (usedSize + 256 > maxQuerySize)
		{
			totalSize += usedSize;
			usedSize = 0;
			full_stream.close();
			cout << full_file << "\t" << setIdx << "\t" << totalSize << endl;
		}
		
	}
	
	if (full_stream.is_open())
	{
		totalSize += usedSize;
		full_stream.close();
		cout << full_file << "\t" << setIdx << "\t" << totalSize << endl;
	}

	full_stream.close();
	cout << endl << "Done... " << endl;
	cout << "Wrote " << setIdx << " sets of sequences (" << totalSize << " chars) to " << fileIdx << " query files." << endl;
}


bool isMarkedForExclusion(const MSBStringSet& msb ) { return msb.indMarkedForExclusion; }

void MSBlastCollector::markAccordingToKeywords(const MSBlastMapFile& mapFile,const set<string>& goodPeptides, const set<string>& badPeptides)
{
	size_t numKeep=0, numExclude=0, numHitMap=0;
	for (size_t i=0; i<msbStringSets.size(); i++)
	{
		MSBStringSet& msb = msbStringSets[i];
		for (size_t j=0; j<msb.strings.size(); j++)
		{
			const string& str = msb.strings[j].seqStr;
			
			if (goodPeptides.find(str) != goodPeptides.end())
			{
				msb.indMarkedForKeep = true;
			}
			else if (badPeptides.find(str) != badPeptides.end())
			{
				msb.indMarkedForExclusion = true;
			}
			else if (str.length()>9)// try substrings
			{
				const size_t lastStartIdx = str.length()-8;
				for (size_t i=0; i<lastStartIdx; i++)
					if (badPeptides.find(str.substr(i)) != badPeptides.end())
					{
						msb.indMarkedForExclusion = true;
						break;
					}
			}

			if (mapFile.checkIfPeptideInMap(str))
			{
				msb.strings[j].indFoundInMap = true;
				numHitMap++;
			}
		}
		if (msb.indMarkedForKeep && msb.indMarkedForExclusion)
			msb.indMarkedForExclusion = false;

		if (msb.indMarkedForKeep)
			numKeep++;

		if (msb.indMarkedForExclusion)
			numExclude++;
	}

	cout << "Marked " << numKeep << " sets for keeping because their protein matched good keyword" << endl;
	cout << "Marked " << numExclude << " sets for exclusion because their protein matched bad keyword (they were also removed)" << endl;
	cout << numHitMap << " peptides in query hit the map file." << endl;

	cout << "size before: " << msbStringSets.size() << endl;
	sort(msbStringSets.begin(), msbStringSets.end(), compMsbStringSet);


	while (msbStringSets.size()>0 && msbStringSets.back().indMarkedForExclusion)
		msbStringSets.pop_back();

	cout << "size after : " << msbStringSets.size() << endl;
}


void MSBlastCollector::generateMsBlastFinalQuery(const char* msb_name,
								   const vector<string>& list_vector)
{
	
}



void create_iterative_query(const char* msb_query_name,
							const char* msb_map_file,
							const char* full_file, 
							const char* seconday_file,
							const char* id_file,
							const char* keep_file,
							const char* exclude_file,
							size_t maxSequencesPerSet,
							size_t maxQuerySize, 
							float minMsbScore)
{

	// read the summary file, and keywords
	const string historyPath = msb_query_name + std::string("_pep_hist.txt");
	
	MSBlastMapFile mapFile;

	ifstream ifs(msb_map_file);
	if (ifs.good())
	{
		ifs.close();
		mapFile.readFile(msb_map_file);
	}
	else
	{
		cout << endl << "Warning did not find a map file!" << endl;
		ifs.close();
	}

	MSBlastKeywordFile keepKeywords, excludeKeywords;

	if (! (keep_file && keepKeywords.readKeywordFile(keep_file)))
		cout << "Warning: did not find keywords file of proteins to keep: " << keep_file << endl;
	if (! (exclude_file && excludeKeywords.readKeywordFile(exclude_file)))
		cout << "Warning: did not find keywords file of proteins to exclude: " << exclude_file << endl;

	MSBlastIdFile idsToKeep;
	if (! (id_file && idsToKeep.readIdFile(id_file)))
		cout << "Warning: did not find ids file with good protein ids: " << id_file << endl;

	const map<string, int>& proteinNames = mapFile.getProteinNames();
	const map<string, vector<int> >& peptides = mapFile.getPeptides();

	int maxIdx=0;
	for (map<string, int>::const_iterator it = proteinNames.begin(); it != proteinNames.end(); it++)
		if (it->second>maxIdx)
			maxIdx=it->second;

	vector<bool> proteinExludeIndicators(maxIdx+1,false);
	vector<bool> proteinKeepIndicators(maxIdx+1,false);

	// test protein names
	for (map<string, int>::const_iterator it = proteinNames.begin(); it != proteinNames.end(); it++)
	{
		string lowerName = it->first;
		for (size_t i=0; i<lowerName.length(); i++)
			lowerName[i]=tolower(lowerName[i]);

		proteinExludeIndicators[it->second] = excludeKeywords.checkForMatch(lowerName);
		proteinKeepIndicators[it->second]   = keepKeywords.checkForMatch(lowerName);
		proteinKeepIndicators[it->second]   = idsToKeep.checkForIdMatch(lowerName);
	}

	// test peptides. Any peptide that is marked fot exclusion without a specific keep mark too will be excluded
	set<string> badPeptides, goodPeptides; // specific peptides that should be kept or thrown away
	for (map<string, vector<int> >::const_iterator it=peptides.begin(); it != peptides.end(); it++)
	{
		bool foundExclude=false;
		size_t i;
		for (i=0; i<it->second.size(); i++)
			if (proteinExludeIndicators[it->second[i]])
				break;
		if (i<it->second.size())
			foundExclude = true;

		bool foundKeep = false;
		for (i=0; i<it->second.size(); i++)
			if (proteinKeepIndicators[it->second[i]])
				break;

		if (i<it->second.size())
			foundKeep = true;

		if (foundExclude && ! foundKeep)
			badPeptides.insert(it->first);
		

		if (foundKeep)
			goodPeptides.insert(it->first);		
	}
	cout << "Found " << goodPeptides.size() << " peptides marked for inclusion (belong to proteins in keep list)" << endl;
	cout << "Found " << badPeptides.size() << " peptides marked for exclusion (belong to proteins in exclusion list)" << endl;

	// read peptide history
	map<string,int> peptideHistory;
	int generationIdx=-1;
	ifstream histfs(historyPath.c_str());
	if (histfs.good())
	{
		histfs >> generationIdx;
		while (! histfs.eof())
		{
			string pep=std::string();
			int count=-1;
			histfs >> pep >> count;
			if (count>0)
				peptideHistory[pep]=count;	
		}
		cout << "Read " << peptideHistory.size() << " peptides from " << historyPath.c_str() << endl;
		histfs.close();
	}
	else
	{
		cout << "Warning: could not find " << historyPath << ".\nCreating new file!" << endl;
		histfs.close();
	}

	
	
	// read files into sequence sets (all results)
	MSBlastCollector queryCollector;

/*	queryCollector.readFullFileIntoCollector(full_file, maxSequencesPerSet*3, minMsbScore, true);
	if (seconday_file)
		queryCollector.readFullFileIntoCollector(seconday_file, maxSequencesPerSet*3, minMsbScore, false);

	queryCollector.markAccordingToKeywords(mapFile, goodPeptides, badPeptides);*/
	
	ostringstream oss;
	oss << msb_query_name << "_" << ++generationIdx;
	queryCollector.generateMsBlastQueryWithHistory(oss.str().c_str(), peptideHistory, maxQuerySize, minMsbScore);

	// write new peptide history
	ofstream ofs(historyPath.c_str());
	if (! ofs.good())
		error("Could not open peptide history for writing: ",historyPath.c_str());

	ofs << generationIdx << endl;
	for (map<string,int>::const_iterator it = peptideHistory.begin(); it != peptideHistory.end(); it++)
		ofs << it->first << "\t" << it->second << endl;
	ofs.close();
	cout << "Wrote " << peptideHistory.size() << " peptides to " << historyPath << endl;
}



struct ProteinResults {
	ProteinResults() : name(std::string()), length(0), totalScore(0), numPeptides(0) {}
	string name;
	int	   length;
	int	   totalScore;
	int	   numPeptides;
};


void examine_msb_results(const char* msb_results_file, 
						 const char* msb_map_file,
						 const char* id_file, 
						 const char* keep_file, 
						 const char* exclude_file)
{
	// read the summary gile, and keywords
	MSBlastMapFile mapFile;
	ifstream ifs(msb_map_file);
	if (ifs.good())
	{
		ifs.close();
		mapFile.readFile(msb_map_file);
	}
	else
	{
		cout << endl << "Warning did not find a map file!" << endl;
		ifs.close();
	}

	MSBlastKeywordFile keepKeywords, excludeKeywords;

	if (keep_file && ! keepKeywords.readKeywordFile(keep_file))
		cout << "Warning: did not find keywords file of proteins to keep: " << keep_file << endl;
	if (exclude_file && ! excludeKeywords.readKeywordFile(exclude_file))
		cout << "Warning: did not find keywords file of proteins to exclude: " << exclude_file << endl;

	MSBlastIdFile idsToKeep;
	if (id_file && ! idsToKeep.readIdFile(id_file))
		cout << "Warning: did not find ids file with good protein ids: " << id_file << endl;

	const map<string, int>& proteinNames = mapFile.getProteinNames();
	const map<string, vector<int> >& peptides = mapFile.getPeptides();

	int maxIdx=0;
	for (map<string, int>::const_iterator it = proteinNames.begin(); it != proteinNames.end(); it++)
		if (it->second>maxIdx)
			maxIdx=it->second;

	vector<bool> proteinExludeIndicators(maxIdx+1,false);
	vector<bool> proteinKeepIndicators(maxIdx+1,false);

	// test protein names
	for (map<string, int>::const_iterator it = proteinNames.begin(); it != proteinNames.end(); it++)
	{
		proteinExludeIndicators[it->second] = excludeKeywords.checkForMatch(it->first);
		proteinKeepIndicators[it->second]   = keepKeywords.checkForMatch(it->first);
		proteinKeepIndicators[it->second]   = idsToKeep.checkForIdMatch(it->first);
	}

	// test peptides. Any peptide that is marked fot exclusion without a specific keep mark too will be excluded
	set<string> badPeptides, goodPeptides; // specific peptides that should be kept or thrown away
	for (map<string, vector<int> >::const_iterator it=peptides.begin(); it != peptides.end(); it++)
	{
		bool foundExclude=false;
		size_t i;
		for (i=0; i<it->second.size(); i++)
			if (proteinExludeIndicators[it->second[i]])
				break;
		if (i<it->second.size())
			foundExclude = true;

		bool foundKeep = false;
		for (i=0; i<it->second.size(); i++)
			if (proteinKeepIndicators[it->second[i]])
				break;

		if (i<it->second.size())
			foundKeep = true;

		if (foundExclude && ! foundKeep)
			badPeptides.insert(it->first);
		

		if (foundKeep)
			goodPeptides.insert(it->first);		
	}
	cout << "Found " << goodPeptides.size() << " peptides marked for inclusion (belong to proteins in keep list)" << endl;
	cout << "Found " << badPeptides.size() << " peptides marked for exclusion (belong to proteins in exclusion list)" << endl;


	const size_t fileSize = getFileSize(msb_results_file);
	ifstream ifs_res(msb_results_file);
	if (! ifs_res.good())
		error("Could not open MSBlast results file for reading: ",msb_results_file);

	cout << "Parsing : " << msb_results_file << " (" << fileSize << " bytes)" << endl;

	char* buffer = new char[fileSize+1024];
	ifs_res.read(buffer, fileSize);
	ifs_res.close();
	buffer[fileSize]='\0';

	// remove terminating strings if they are somewhere in the buffer
	for (size_t i=0; i<fileSize; i++)
		if (buffer[i]=='\0')
			buffer[i]=' ';

	// find the protein results
	const char* alignmentsPtr = strstr(buffer,"Alignments:");
	if (! alignmentsPtr)
		error("Could not find \"Alignments:\" in results file");

	size_t numProteinsExamined=0;
	size_t numProteinsAdded=0;
	size_t numPeptidesExamined=0;
	size_t numPeptidesAdded=0;
	size_t numRefsAdded =0;

	// loop until all proteins are done

	vector<ProteinResults> goodResults, excludedResults, neutralResults;
	const char* ptr = alignmentsPtr;
	while (ptr)
	{
		const char* protPtr = strstr(ptr, "^ =");
		if (! protPtr)
			break;
		
		const char* lengthPtr = strstr(protPtr,"Length = ");
		if (! lengthPtr)
			error("Could not find length of protein!");

		int protLength = atoi(lengthPtr+9);
		assert(protLength>0 && protLength < 1E7);

		// parse protein name
		assert(lengthPtr - 5 > protPtr);
		string orgName = std::string(protPtr+4, lengthPtr - protPtr - 5);
		string name = std::string();
		// convert white space
		for (size_t i=0; i<orgName.length(); i++)
		{
			if (orgName[i] == '\t' || orgName[i] == '\r' || orgName[i] == '\n')
				orgName[i] = ' ';
			if (i>0 && orgName[i-1]==' ' && orgName[i] == ' ')
				continue;
			name.push_back(orgName[i]);
		}

		// remove trailing whitespace
		assert(name.length()>5);
		size_t last=name.length()-1;
		while (last>0 && name[last] == ' ')
			last--;
		if (last<name.length()-1)
			name.erase(last+1);

		const char* totalScorePtr = strstr(lengthPtr,"Total Score:");
		if (! totalScorePtr)
			error("Could not find total score");

		
		const char* nextProt = strstr(lengthPtr,"^ =");
		if (! nextProt)
		{
			nextProt = strstr(lengthPtr,"Parameters:");
			if (! nextProt)
				error("Could not find teminating \"Parameters:\"");
		}

		numProteinsExamined++;

		// parse the query hits
		vector<string> sequences;
		const char* p = lengthPtr;
		while (1)
		{
			p=strstr(p,"Query:");
			if (! p || p>=nextProt)
				break;
			
			istringstream iss(std::string(p,64));
			string dummy=std::string(), idx=std::string(), pep=std::string();
			iss >> dummy >> idx >> pep;
			if (pep.length()<3)
			{
				char* pp=const_cast<char*>(nextProt);
				*pp='\0';
				error("Problem parsing line: ",p);
			}
			p+=16;
			numPeptidesExamined++;
			sequences.push_back(pep);
		}
		ptr = lengthPtr;

		// add protein if needed
		name=name.substr(0,128);
		if (! mapFile.checkIfProteinInMap(name))
		{
			cout << "Protein name not found in map: " << name << endl;
			cout << "should run -msb_process_results before!" << endl;
			exit(0);
		}

		ProteinResults pr;
		pr.name = name;
		pr.length = protLength;
		pr.numPeptides = sequences.size();
		pr.totalScore = atoi(totalScorePtr+12);

		string lowerName = name;
		for (size_t i=0; i<name.length(); i++)
			lowerName[i] = tolower(name[i]);

		if (idsToKeep.checkForIdMatch(lowerName) || keepKeywords.checkForMatch(lowerName))
		{
			goodResults.push_back(pr);
		}
		else if (excludeKeywords.checkForMatch(lowerName))
		{
			excludedResults.push_back(pr);
		}
		else
			neutralResults.push_back(pr);
	}

	cout << endl << "Read results with " << numProteinsExamined << " identified proteins and " << numPeptidesExamined << " peptides:" << endl;
	cout << goodResults.size() << " proteins with good ids or key words in their names." << endl;
	cout << excludedResults.size() << " proteins with bad keywords in their names." << endl;
	cout << neutralResults.size()  << " proteins that did not hit any keywords." << endl;

	if (goodResults.size())
	{
		cout << endl << "Good protein hits:" << endl;
		cout		 << "-----------------" << endl;
		cout << "#\tScore\t#Peps\tLength\tName" << endl;
		for (size_t i=0; i<goodResults.size(); i++)
			cout << i+1 << "\t" << goodResults[i].totalScore << "\t" << goodResults[i].numPeptides << "\t" << goodResults[i].length << "\t"
				 << goodResults[i].name << endl;
	}

	if (neutralResults.size())
	{
		cout << endl << "Neutral protein hits:" << endl;
		cout <<         "--------------------" << endl;
		cout << "#\tScore\t#Peps\tLength\tName" << endl;
		for (size_t i=0; i<neutralResults.size(); i++)
			cout << i+1 << "\t" << neutralResults[i].totalScore << "\t" << neutralResults[i].numPeptides << "\t" << neutralResults[i].length << "\t"
				 << neutralResults[i].name << endl;
	}

}



struct ReportEntry {
	bool operator< (const ReportEntry& rhs) const 
	{
		return (score>rhs.score);
	}
	string fullMSBEntry;
	float score;
	int length;
	int numPeptides;
	string name; //
};


void writeReportForProteinSet(const char* name, const map<string, ReportEntry>& prots)
{
	vector<ReportEntry> entries;
	for (map<string, ReportEntry>::const_iterator it=prots.begin(); it != prots.end(); it++)
		entries.push_back(it->second);
	sort(entries.begin(), entries.end());

	if (entries.size()<1)
		return;

	string summaryName = std::string(name) + "_sum.txt";
	string allResultsName = std::string(name) + "_allres.txt";
	ofstream summaryStream(summaryName.c_str());
	ofstream allResultsStream(allResultsName.c_str());
	if (! summaryStream.good() || ! allResultsStream.good())
		error("Could not open result files for writing: ",name);

	summaryStream << "#Num\tScore\t#peptides Length  Protein name" << endl;
	for (size_t i=0; i<entries.size(); i++)
	{
		summaryStream << "#" << i+1 << "\t" << entries[i].score << "\t" << entries[i].numPeptides << "\t" << entries[i].length << "\t" <<
			entries[i].name << endl;

		allResultsStream << "#" << i+1 << "\t" << entries[i].fullMSBEntry << endl << endl;
	}
	summaryStream.close();
	allResultsStream.close();
	cout << "Wrote results of " << entries.size() << " protein hits to:" << endl;
	cout << "\t" << summaryName << endl;
	cout << "\t" << allResultsName << endl;
}


void create_msb_reprot( const char* list_file,
						const char* msb_query_name,
						const char* msb_map_file,
						const char* msb_id_file,
						const char* msb_keep_file,
						const char* msb_exclude_file)
{
	vector<string> paths;
	readListOfPaths(list_file, paths);
	if (paths.size()<1)
		error("Must supply file with full list of paths to MS-Blast results.");

	MSBlastMapFile mapFile;
	ifstream ifs(msb_map_file);
	if (ifs.good())
	{
		ifs.close();
		mapFile.readFile(msb_map_file);
	}
	else
		error("Must provide valid map file!");
	
	MSBlastKeywordFile keepKeywords, excludeKeywords;

	if (msb_keep_file && ! keepKeywords.readKeywordFile(msb_keep_file))
		cout << "Warning: did not find keywords file of proteins to keep: " << msb_keep_file << endl;
	if (msb_exclude_file && ! excludeKeywords.readKeywordFile(msb_exclude_file))
		cout << "Warning: did not find keywords file of proteins to exclude: " << msb_exclude_file << endl;

	MSBlastIdFile idsToKeep;
	if (msb_id_file && ! idsToKeep.readIdFile(msb_id_file))
		cout << "Warning: did not find ids file with good protein ids: " << msb_id_file << endl;

	const map<string, int>& proteinNames = mapFile.getProteinNames();
	const map<string, vector<int> >& peptides = mapFile.getPeptides();

	int maxIdx=0;
	for (map<string, int>::const_iterator it = proteinNames.begin(); it != proteinNames.end(); it++)
		if (it->second>maxIdx)
			maxIdx=it->second;

	vector<bool> proteinExludeIndicators(maxIdx+1,false);
	vector<bool> proteinKeepIndicators(maxIdx+1,false);

	// test protein names
	for (map<string, int>::const_iterator it = proteinNames.begin(); it != proteinNames.end(); it++)
	{
		proteinExludeIndicators[it->second] = excludeKeywords.checkForMatch(it->first);
		proteinKeepIndicators[it->second]   = keepKeywords.checkForMatch(it->first);
		proteinKeepIndicators[it->second]   = idsToKeep.checkForIdMatch(it->first);
	}

	map< string, ReportEntry> goodResults, excludedResults, otherResults;

	for (size_t f=0; f<paths.size(); f++)
	{
		const size_t fileSize = getFileSize(paths[f].c_str());
		ifstream ifs_res(paths[f].c_str());
		if (! ifs_res.good())
			error("Could not open MSBlast results file for reading: ",paths[f].c_str());

		cout << "Parsing : " << paths[f] << " (" << fileSize << " bytes)" << endl;

		char* buffer = new char[fileSize+1024];
		ifs_res.read(buffer, fileSize);
		ifs_res.close();
		buffer[fileSize]='\0';

		// remove terminating strings if they are somewhere in the buffer
		for (size_t i=0; i<fileSize; i++)
			if (buffer[i]=='\0')
				buffer[i]=' ';

		// find the protein results
		const char* alignmentsPtr = strstr(buffer,"Alignments:");
		if (! alignmentsPtr)
			error("Could not find \"Alignments:\" in results file");

		size_t numProteinsExamined=0;
	
		// loop until all proteins are done
		const char* ptr = alignmentsPtr;
		while (ptr)
		{
			const char* protPtr = strstr(ptr, "^ =");
			if (! protPtr)
				break;
			
			const char* lengthPtr = strstr(protPtr,"Length = ");
			if (! lengthPtr)
				error("Could not find length of protein!");

			int protLength = atoi(lengthPtr+9);
			assert(protLength>0 && protLength < 1E7);

			// parse protein name
			assert(lengthPtr - 5 > protPtr);
			string orgName = std::string(protPtr+4, lengthPtr - protPtr - 5);
			string name = std::string();
			// convert white space
			for (size_t i=0; i<orgName.length(); i++)
			{
				if (orgName[i] == '\t' || orgName[i] == '\r' || orgName[i] == '\n')
					orgName[i] = ' ';
				if (i>0 && orgName[i-1]==' ' && orgName[i] == ' ')
					continue;
				name.push_back(orgName[i]);
			}

			// remove trailing whitespace
			assert(name.length()>5);
			size_t last=name.length()-1;
			while (last>0 && name[last] == ' ')
				last--;
			if (last<name.length()-1)
				name.erase(last+1);

			const char* totalScorePtr = strstr(lengthPtr,"Total Score:");
			if (! totalScorePtr)
				error("Could not find total score");

			
			const char* nextProt = strstr(lengthPtr,"^ =");
			if (! nextProt)
			{
				nextProt = strstr(lengthPtr,"Parameters:");
				if (! nextProt)
					error("Could not find teminating \"Parameters:\"");
			}

			numProteinsExamined++;

			// parse the query hits
			vector<string> sequences;
			const char* p = lengthPtr;
			while (1)
			{
				p=strstr(p,"Query:");
				if (! p || p>=nextProt)
					break;
				
				istringstream iss(std::string(p,64));
				string dummy=std::string(), idx=std::string(), pep=std::string();
				iss >> dummy >> idx >> pep;
				if (pep.length()<3)
				{
					char* pp=const_cast<char*>(nextProt);
					*pp='\0';
					error("Problem parsing line: ",p);
				}
				p+=16;
				sequences.push_back(pep);
			}
			ptr = lengthPtr;

			// add protein if needed
			name=name.substr(0,128);
			if (! mapFile.checkIfProteinInMap(name))
			{
				cout << "Protein name not found in map: " << name << endl;
				cout << "should run -msb_process_results before!" << endl;
				exit(0);
			}

			ReportEntry newEntry;
			newEntry.fullMSBEntry = std::string(protPtr,nextProt-protPtr);
			newEntry.length = protLength;
			newEntry.numPeptides  = sequences.size();
			newEntry.score		  = atoi(totalScorePtr+12);
			newEntry.name		  = name;

			string lowerName = name;
			for (size_t i=0; i<name.length(); i++)
				lowerName[i] = tolower(name[i]);

			if (idsToKeep.checkForIdMatch(lowerName) || keepKeywords.checkForMatch(lowerName))
			{
				map<string, ReportEntry>::iterator it = goodResults.find(name);
				if ( it == goodResults.end() || newEntry.score > it->second.score)
					goodResults[name]=newEntry;
			}
			else if (excludeKeywords.checkForMatch(lowerName))
			{
				map<string, ReportEntry>::iterator it = excludedResults.find(name);
				if ( it == excludedResults.end() || newEntry.score > it->second.score)
					excludedResults[name]=newEntry;
			}
			else
			{
				map<string, ReportEntry>::iterator it = otherResults.find(name);
				if ( it == otherResults.end() || newEntry.score > it->second.score)
					otherResults[name]=newEntry;
			}
		}
		delete [] buffer;
		cout << "	found " << numProteinsExamined << " protein entries. Protein set sizes: [good="
			 << goodResults.size() << " ,excluded=" << excludedResults.size() << " ,other=" << otherResults.size() << "]" << endl << endl;
	}

	// write final reports
	string goodName = std::string(msb_query_name) + "_good";
	string excludedName = std::string(msb_query_name) + "_excluded";
	string otherName = std::string(msb_query_name) + "_other";

	writeReportForProteinSet(goodName.c_str(), goodResults);
	writeReportForProteinSet(excludedName.c_str(), excludedResults);
	writeReportForProteinSet(otherName.c_str(), otherResults);
}







