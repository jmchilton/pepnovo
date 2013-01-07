#include "PMCSQS.h"
#include "SpectraList.h"
#include "auxfun.h"


bool PMCSQS_Scorer::read_pmc_rank_models(const Config *_config, const char *file_name)
{
	config_ = _config;

	string path = config_->get_resource_dir() + "/" + file_name;
	ifstream in_stream(path.c_str(),ios::in);
	if (! in_stream.good())
	{
		cout << "Warning: couldn't open pmc rank model for reading: " << path << endl;
		return false;
	}


	char buff[512];
	int numCharges=-1;

	in_stream.getline(buff,256);
	istringstream iss1(buff);

	frag_pair_sum_offset=NEG_INF;
	bin_increment=NEG_INF;
	iss1 >> bin_increment >> this->frag_pair_sum_offset;
	if (frag_pair_sum_offset==NEG_INF || bin_increment == NEG_INF)
	{
		cout << "Error in pmc model file!" << endl;
		exit(1);
	}

	in_stream.getline(buff,256);
	istringstream iss(buff);

	iss >> numCharges;
	maximalChargeWithModels_=numCharges-1;
	
	pmc_rank_models.resize(numCharges);
	pmcMassThresholds_.resize(numCharges);
	pmcMzBiases_.resize(numCharges);


	int i;
	for (i=0; i<numCharges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);
		int num_threshes=0;
		iss >> num_threshes;
		
		pmcMassThresholds_[i].resize(num_threshes,NEG_INF);
		int j;
		for (j=0; j<num_threshes; j++)
			iss >> pmcMassThresholds_[i][j];
	}

	for (i=0; i<numCharges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);
		int num_biases=0;
		iss >> num_biases;
		
		pmcMzBiases_[i].resize(num_biases,NEG_INF);
		int j;
		for (j=0; j<num_biases; j++)
			iss >> pmcMzBiases_[i][j];
	}
	
	// read Boost models
	for (i=0; i<numCharges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);

		int num_models=-1;
		iss >> num_models;

		if (num_models<0)
		{
			cout << "Error: bad parsing of PMCR model file!" << endl;
			exit(0);
		}
		pmc_rank_models[i].resize(num_models,NULL);

		int j;
		for (j=0; j<num_models; j++)
		{
			pmc_rank_models[i][j]=new RankBoostModel;
			pmc_rank_models[i][j]->read_rankboost_model(in_stream);
		}
		
	}
	in_stream.close();


	this->ind_initialized_pmcr = true;
	return true;
}


void PMCSQS_Scorer::write_pmc_rank_models(const char *path) const
{
	ofstream out_stream(path,ios::out);
	if (! out_stream.good())
	{
		cout << "Error: couldn't open pmc model for writing: " << path << endl;
		exit(1);
	}

	out_stream << this->bin_increment << " " << this->frag_pair_sum_offset << endl;
	out_stream << this->pmcMassThresholds_.size() << endl;
	out_stream << setprecision(3) << fixed;

	int i;
	for (i=0; i<this->pmcMassThresholds_.size(); i++)
	{
		out_stream << pmcMassThresholds_[i].size();
		int j;
		for (j=0; j<pmcMassThresholds_[i].size(); j++)
			out_stream << " " << pmcMassThresholds_[i][j];
		out_stream << endl;
	}

	
	for (i=0; i<this->pmcMzBiases_.size(); i++)
	{
		out_stream << pmcMzBiases_[i].size();
		int j;
		for (j=0; j<pmcMzBiases_[i].size(); j++)
			out_stream << " " << pmcMzBiases_[i][j];
		out_stream << endl;
	}

	// write RankBoost models
	for (i=0; i<pmc_rank_models.size(); i++)
	{
		int j;
		
		if (pmc_rank_models[i].size()==1 && ! pmc_rank_models[i][0])
		{
			out_stream << 0 << endl;
			continue;
		}

		out_stream << pmc_rank_models[i].size() << endl;
		for (j=0; j<pmc_rank_models[i].size(); j++)
		{
			if (pmc_rank_models[i][j])
			{
				pmc_rank_models[i][j]->write_rankboost_model(out_stream,true);
			}
			else
			{
				cout << "Error: non intialized rank pmc model!" << endl;
				exit(1);
			}
		}
	}
	
	out_stream.close();
}


void PMCSQS_Scorer::set_pmc_mass_thresholds(int option)
{
	if (option==0)
	{
		pmcMassThresholds_ = config_->get_size_thresholds();
		int i;
		for (i=0; i<pmcMassThresholds_.size(); i++)
		{
			if (pmcMassThresholds_[i].size()>0)
			{
				if (pmcMassThresholds_[i][pmcMassThresholds_[i].size()-1]>10000)
					pmcMassThresholds_[i].pop_back();
			}
		}
	}

	if (option==1)
	{
		pmcMassThresholds_.resize(4);
		pmcMassThresholds_[1].push_back(1150.0);
		pmcMassThresholds_[1].push_back(1400.0);
 
		pmcMassThresholds_[2].push_back(1100.0);
		pmcMassThresholds_[2].push_back(1300.0);
		pmcMassThresholds_[2].push_back(1600.0);
		pmcMassThresholds_[2].push_back(1900.0);
		pmcMassThresholds_[2].push_back(2400.0);

		pmcMassThresholds_[3].push_back(1950.0);
		pmcMassThresholds_[3].push_back(2450.0);
		pmcMassThresholds_[3].push_back(3000.0);
	}

	if (option==2)
	{
		pmcMassThresholds_.resize(4);
		
		pmcMassThresholds_[2].push_back(1300.0);
	
		pmcMassThresholds_[2].push_back(1900.0);
	}
}



void convert_ME_to_RankBoostSample(const ME_Regression_Sample& me,
								   RankBoostSample& rbs)
{
	rbs.clear();
	int i;
	for (i=0; i<me.f_vals.size(); i++)
		rbs.add_real_feature(me.f_vals[i].f_idx,me.f_vals[i].val);
}



void PMCSQS_Scorer::trainPmcRankModels(const Config* config, const SpectraAggregator& sa, 
									   int specificCharge, bool overwrite)
{
	const bool indSampleDiagnostic = false;
	const vector<int>& spectraCounts = sa.getSpectraCountsPerCharge();
	
	
	const int maxSpectraPerFile = 40000;
	
	vector<string> realNames;
	init_PMC_feature_names(realNames);


	// try and read existing pmc model, otherwise init a new one
	string pmcPath = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCR.txt";
	ifstream modelStream(pmcPath.c_str());	

	if (! modelStream.is_open() && modelStream.good())
	{

		modelStream.close();
		string pmcrName = config->get_model_name() + "_PMCR.txt";
		const char* path = pmcPath.c_str();
		read_pmc_rank_models(config,(char *)pmcrName.c_str());
	}
	else
	{
		set_pmc_mass_thresholds();
		set_frag_pair_sum_offset(MASS_PROTON); // b+y - PM+19
		set_bin_increment(0.1);
		pmc_rank_models.resize(pmcMassThresholds_.size());
		pmcMzBiases_.resize(pmcMassThresholds_.size());
	}
	

	maximalChargeWithModels_= pmcMassThresholds_.size()-1;
	const double prop_train = 0.5;


	// It is assumed that the mass thresholds were set according to the training data
	// (this is done manually with values encoded in the set_mass_threhsolds function)
	int charge;
	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		if (specificCharge>0 && charge != specificCharge)
			continue;

		const int num_sizes = pmcMassThresholds_[charge].size();
		pmc_rank_models[charge].resize(num_sizes+1,NULL);
		pmcMzBiases_[charge].resize(num_sizes+1,0);

		
		int sizeIndex;
		for (sizeIndex=0; sizeIndex<=num_sizes; sizeIndex++)
		{
			mass_t minMz =0;
			mass_t maxMz = POS_INF;

			if (sizeIndex>0)
				minMz = pmcMassThresholds_[charge][sizeIndex-1];

			if (sizeIndex<num_sizes)
				maxMz = pmcMassThresholds_[charge][sizeIndex];

			minMz /= charge;
			maxMz /= charge;

			cout << "\nTraining PMC rank model for charge " << charge << " size " << sizeIndex << endl;
			cout << "min m/z: " << minMz << endl;
			cout << "max m/z: " << maxMz << endl;
			if (pmc_rank_models[charge][sizeIndex] && ! overwrite)
			{
				cout << endl << "Already trained..." << endl;
				continue;
			}

			vector<const SingleSpectrumHeader*> testHeaders;
			RankBoostDataset train_ds, test_ds, pos_ds, neg_ds;

		
			// these ranges are given according to pm_with_19
			// so files should be selected through select_files and not
			// select_file_in_mz_range
			SpectraList sl(sa);
			sl.selectHeaders(minMz, maxMz, charge, charge);
		

			if (sl.getNumHeaders()<500)
			{
				cout << "Not enough spectra to train this model, found only " << sl.getNumHeaders() << ", skipping..." << endl;
				continue;
			}
			sl.randomlyReduceListToSize(maxSpectraPerFile);

			cout << "PM m/z range " << minMz << "-" << maxMz << "  using " << sl.getNumHeaders() << " spectra." << endl;

			int numTrainingGroups=0;
			int numTestingGroups=0;
			
			const int numSamples = sl.getNumHeaders();
			
			// first find the bias in number of bins between the true m/z bin and
			// the optimal m/z bin
			vector<bool> skippedIndexes;
			skippedIndexes.resize(numSamples,false);
			int skippedBadMz=0;
			mass_t totalBias=0;
			int i;
			for (i=0; i<numSamples; i++)
			{
				const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
				PeakList pl;

				pl.readPeaksToLocalAllocation(sa, header);
			
				if (header->getPeptideStr().length()<3)
					continue;

				Peptide peptide;
				peptide.parseFromString(config, header->getPeptideStr());
				const mass_t trueMz = (peptide.get_mass()+MASS_H2O+ charge*MASS_PROTON)/(mass_t)charge;

				if (fabs(trueMz - header->getMOverZ())>2.5)
				{
					//cout << setprecision(2) << trueMz << " <---> " << bs.ssf->m_over_z << " skipping" << endl;
					skippedBadMz++;
					skippedIndexes[i]=true;
					continue;
				} 

				initializeForCurrentSpectrum(config, pl);
				calculateCurrentSpectrumPmcValues(pl, bin_increment);
		

				// find the true_mz_bin_idx
				
				const vector<PMCRankStats>& pmc_stats = currentSpectrumPmcTables_[charge];
				int true_mz_bin_idx=0;
				while (true_mz_bin_idx<pmc_stats.size() && pmc_stats[true_mz_bin_idx].m_over_z<trueMz)
					true_mz_bin_idx++;

				if (true_mz_bin_idx == pmc_stats.size())
					true_mz_bin_idx--;

				if (true_mz_bin_idx>0 && pmc_stats[true_mz_bin_idx].m_over_z-trueMz>trueMz-pmc_stats[true_mz_bin_idx-1].m_over_z)
					true_mz_bin_idx--;

				int opt_bin_idx = get_optimal_bin(true_mz_bin_idx, charge);

				if (opt_bin_idx <=0 || opt_bin_idx == pmc_stats.size()-1)
				{
					skippedBadMz++;
					skippedIndexes[i]=true;
					continue;
				}

				totalBias += (pmc_stats[opt_bin_idx].m_over_z - pmc_stats[true_mz_bin_idx].m_over_z);

				if (fabs(pmc_stats[opt_bin_idx].m_over_z - pmc_stats[true_mz_bin_idx].m_over_z)>4.0)
				{
					cout << "opt bin: " << opt_bin_idx << " (" << pmc_stats[opt_bin_idx].m_over_z << ")  ";
					cout << "tru bin: " << true_mz_bin_idx << " ("<< pmc_stats[true_mz_bin_idx].m_over_z << ")" << endl;
				}
			} 

			mass_t mz_bias = totalBias / (mass_t)(numSamples-skippedBadMz);
			pmcMzBiases_[charge][sizeIndex]=mz_bias;

			cout << "m/z bias: " << setprecision(4) << mz_bias << endl;
			cout << "skipped " << skippedBadMz << "/" << numSamples <<
				"  because of m/z more than 2.5 away from observed..." << endl; 

		//	pmcMzBiases_[charge][sizeIndex] = 0;

			for (i=0; i<numSamples; i++)
			{
				if (skippedIndexes[i])
					continue;

				const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
				PeakList pl;

				pl.readPeaksToLocalAllocation(sa, header);
				
				if (header->getPeptideStr().length()<3)
					continue;

				Peptide peptide;
				peptide.parseFromString(config, header->getPeptideStr());
				
				const mass_t trueMz = (peptide.get_mass()+MASS_H2O+ charge*MASS_PROTON)/(mass_t)charge;

				initializeForCurrentSpectrum(config, pl);
				calculateCurrentSpectrumPmcValues(pl, bin_increment);

				// find the true_mz_bin_idx
				
				const vector<PMCRankStats>& pmc_stats = currentSpectrumPmcTables_[charge];
				int true_mz_bin_idx=0;
				while (true_mz_bin_idx<pmc_stats.size() && pmc_stats[true_mz_bin_idx].m_over_z<trueMz)
					true_mz_bin_idx++;

				if (true_mz_bin_idx == pmc_stats.size())
					true_mz_bin_idx--;

				if (true_mz_bin_idx>0 && pmc_stats[true_mz_bin_idx].m_over_z-trueMz>trueMz-pmc_stats[true_mz_bin_idx-1].m_over_z)
					true_mz_bin_idx--;

				int opt_bin_idx = get_optimal_bin(true_mz_bin_idx, charge);

				
				static vector<RankBoostSample> spectrumRankSamples;
				fillRankboostPmcSamples(pl, charge, spectrumRankSamples);

				// select samples and add them to pmc_ds
				int goodIndex=-1;
				vector<int> badIndexes;
				selectTrainingSampleIndexes(charge, spectrumRankSamples, pl, goodIndex, badIndexes);

				const bool indAddToTraining = (myRandom()<prop_train);

				int groupIndex;
				if (indAddToTraining)
				{
					groupIndex= numTrainingGroups++;	
				}
				else
				{
					groupIndex= numTestingGroups++;
					testHeaders.push_back(header);
				}
				
				
				RankBoostDataset& ds = (indAddToTraining ? train_ds : test_ds);

				const int pos_index  = ds.get_num_samples();
				spectrumRankSamples[goodIndex].groupIndex = groupIndex;
				spectrumRankSamples[goodIndex].rank_in_group=0;

				ds.add_sample(spectrumRankSamples[goodIndex]);
				if (indSampleDiagnostic)
					pos_ds.add_sample(spectrumRankSamples[goodIndex]);

				int j;
				for (j=0; j<badIndexes.size(); j++)
				{
					const int bad_idx = badIndexes[j];
					if (bad_idx < 0 || bad_idx>= spectrumRankSamples.size())
						continue;
		
					spectrumRankSamples[bad_idx].groupIndex=groupIndex;
					spectrumRankSamples[bad_idx].rank_in_group=1;

					ds.add_to_phi_vector(ds.get_num_samples(),pos_index);
					ds.add_sample(spectrumRankSamples[bad_idx]);

					if (indSampleDiagnostic)
						neg_ds.add_sample(spectrumRankSamples[bad_idx]);
				}						   
			}

			train_ds.set_num_groups(numTrainingGroups);
			test_ds.set_num_groups(numTestingGroups);
			
			train_ds.compute_total_phi_weight();
			train_ds.initialize_potenital_lists();
			train_ds.initialzie_real_feature_table(realNames.size());

			test_ds.compute_total_phi_weight();

			if (pmc_rank_models[charge][sizeIndex])
				delete pmc_rank_models[charge][sizeIndex];
			
			pmc_rank_models[charge][sizeIndex] = new RankBoostModel;
		

			RankBoostModel* boost = pmc_rank_models[charge][sizeIndex];

			vector<string> empty;
			empty.clear();
			boost->init_rankboost_model_feature_names(empty,realNames);
			boost->init_rankboost_model_for_training(train_ds,100,25);

			train_ds.initialize_real_vote_lists(*boost);

	/*		if (indSampleDiagnostic)
			{
				boost->summarize_features_pos_neg(pos_ds.get_samples(),neg_ds.get_samples());
			}
			else
				boost->summarize_features(train_ds.get_samples());*/

			boost->train_rankboost_model(train_ds,3000,NULL,&test_ds);
			
			boost->ouput_ranked_feature_list();

		//	output_pmc_rank_results(fm,charge,testHeaders);

		//	exit(0);

			ind_initialized_pmcr = true;
		//	string path;
		//	path = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCRtt.txt";
		//	this->write_pmc_rank_models(path.c_str());
			
		}
	}

	string path;
	path = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCR.txt";
	this->write_pmc_rank_models(path.c_str());
	ind_initialized_pmcr = true;
}




/******************************************************************************
Train PMC models from positive example files
*******************************************************************************/
void PMCSQS_Scorer::train_pmc_rank_models(Config *config, const FileManager& fm, 
										  int specificCharge, bool overwrite)
{	
	const bool indSampleDiagnostic = false;
	const vector<int>& spectraCounts = fm.get_spectra_counts();
	
	maximalChargeWithModels_=0;

	int charge;
	for (charge=1; charge<spectraCounts.size(); charge++)
	{
		if (spectraCounts[charge]>=MIN_SPECTRA_FOR_PMCSQS_MODEL)
			maximalChargeWithModels_=charge;
	}

	const int maxSpectraPerFile = 40000;
	
	vector<string> realNames;
	init_PMC_feature_names(realNames);


	// try and read existing pmc model, otherwise init a new one
	string pmcPath = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCR.txt";
	ifstream modelStream(pmcPath.c_str());
	if (modelStream.is_open() && modelStream.good())
	{
		modelStream.close();
		string pmcrName = config->get_model_name() + "_PMCR.txt";
		read_pmc_rank_models(config,pmcrName.c_str());
	}
	else
	{
		set_pmc_mass_thresholds();
	
		set_frag_pair_sum_offset(MASS_PROTON); // b+y - PM+19
		set_bin_increment(0.1);
		pmc_rank_models.resize(pmcMassThresholds_.size());
		pmcMzBiases_.resize(pmcMassThresholds_.size());
	}
	
	const double prop_train = 0.5;


	// It is assumed that the mass thresholds were set according to the training data
	// (this is done manually with values encoded in the set_mass_threhsolds function)
	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		if (specificCharge>0 && charge != specificCharge)
			continue;

		const int num_sizes = pmcMassThresholds_[charge].size();
		pmc_rank_models[charge].resize(num_sizes+1,NULL);
		pmcMzBiases_[charge].resize(num_sizes+1,0);

		
		int sizeIndex;
		for (sizeIndex=0; sizeIndex<=num_sizes; sizeIndex++)
		{
			cout << "\nTraining PMC rank model for charge " << charge << " size " << sizeIndex << endl;
			if (pmc_rank_models[charge][sizeIndex] && ! overwrite)
			{
				cout << endl << "Already trained..." << endl;
				continue;
			}

			vector<SingleSpectrumFile *> testHeaders;
			BasicSpecReader bsr;
			static QCPeak peaks[5000];
			RankBoostDataset train_ds, test_ds, pos_ds, neg_ds;

			mass_t minMass =0;
			mass_t maxMass = POS_INF;

			if (sizeIndex>0)
				minMass = pmcMassThresholds_[charge][sizeIndex-1];

			if (sizeIndex<num_sizes)
				maxMass = pmcMassThresholds_[charge][sizeIndex];

			// these ranges are given according to pm_with_19
			// so files should be selected through select_files and not
			// select_file_in_mz_range
			FileSet fs;		
			fs.select_files(fm,minMass,maxMass,-1,-1,charge);

			if (fs.get_total_spectra()<500)
			{
				cout << "Not enough spectra to train this model, skipping..." << endl;
				continue;
			}

			
			int numTrainingGroups=0;
			int numTestingGroups=0;

			cout << "PM range " << minMass << "-" << maxMass << endl;

			fs.randomly_reduce_ssfs(maxSpectraPerFile);
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
			const int numSamples = all_ssf.size();
			
			// first find the bias in number of bins between the true m/z bin and
			// the optimal m/z bin
			vector<bool> skippedIndexes;
			skippedIndexes.resize(numSamples,false);
			int skippedBadMz=0;
			mass_t totalBias=0;
			int i;
			for (i=0; i<numSamples; i++)
			{
				SingleSpectrumFile* ssf = all_ssf[i];
				BasicSpectrum bs;
			
				bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
				bs.peaks = peaks;
				bs.ssf = ssf;

				ssf->peptide.calc_mass(config);
				
				const mass_t trueMz = (ssf->peptide.get_mass()+MASS_H2O+(mass_t)charge)/(mass_t)charge;

				if (fabs(trueMz - bs.ssf->m_over_z)>2.5)
				{
					//cout << setprecision(2) << trueMz << " <---> " << bs.ssf->m_over_z << " skipping" << endl;
					skippedBadMz++;
					skippedIndexes[i]=true;
					continue;
				} 

				init_for_current_spec(config,bs);
				calculate_curr_spec_pmc_values(bs, bin_increment);

				// find the true_mz_bin_idx
				
				const vector<PMCRankStats>& pmc_stats = currentSpectrumPmcTables_[charge];
				int true_mz_bin_idx=0;
				while (true_mz_bin_idx<pmc_stats.size() && pmc_stats[true_mz_bin_idx].m_over_z<trueMz)
					true_mz_bin_idx++;

				if (true_mz_bin_idx == pmc_stats.size())
					true_mz_bin_idx--;

				if (true_mz_bin_idx>0 && pmc_stats[true_mz_bin_idx].m_over_z-trueMz>trueMz-pmc_stats[true_mz_bin_idx-1].m_over_z)
					true_mz_bin_idx--;

				int opt_bin_idx = get_optimal_bin(true_mz_bin_idx, charge);

				if (opt_bin_idx <=0 || opt_bin_idx == pmc_stats.size()-1)
				{
					skippedBadMz++;
					skippedIndexes[i]=true;
					continue;
				}

				totalBias += (pmc_stats[opt_bin_idx].m_over_z - pmc_stats[true_mz_bin_idx].m_over_z);

				if (fabs(pmc_stats[opt_bin_idx].m_over_z - pmc_stats[true_mz_bin_idx].m_over_z)>4.0)
				{
					cout << "opt bin: " << opt_bin_idx << " (" << pmc_stats[opt_bin_idx].m_over_z << ")  ";
					cout << "tru bin: " << true_mz_bin_idx << " ("<< pmc_stats[true_mz_bin_idx].m_over_z << ")" << endl;
				}
			} 

			mass_t mz_bias = totalBias / (mass_t)(numSamples-skippedBadMz);
			pmcMzBiases_[charge][sizeIndex]=mz_bias;

			cout << "m/z bias: " << setprecision(4) << mz_bias << endl;
			cout << "skipped " << skippedBadMz << "/" << numSamples <<
				"  because of m/z more than 2.5 away from observed..." << endl; 

		//	pmcMzBiases_[charge][sizeIndex] = 0;

			for (i=0; i<numSamples; i++)
			{
				if (skippedIndexes[i])
					continue;

				SingleSpectrumFile* ssf = all_ssf[i];
				BasicSpectrum bs;
			
				bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
				bs.peaks = peaks;
				bs.ssf = ssf;
				const mass_t trueMz = (ssf->peptide.get_mass()+MASS_H2O+(mass_t)charge)/(mass_t)charge;

				init_for_current_spec(config,bs);
				calculate_curr_spec_pmc_values(bs, bin_increment);

				// find the true_mz_bin_idx
				
				const vector<PMCRankStats>& pmc_stats = currentSpectrumPmcTables_[charge];
				int true_mz_bin_idx=0;
				while (true_mz_bin_idx<pmc_stats.size() && pmc_stats[true_mz_bin_idx].m_over_z<trueMz)
					true_mz_bin_idx++;

				if (true_mz_bin_idx == pmc_stats.size())
					true_mz_bin_idx--;

				if (true_mz_bin_idx>0 && pmc_stats[true_mz_bin_idx].m_over_z-trueMz>trueMz-pmc_stats[true_mz_bin_idx-1].m_over_z)
					true_mz_bin_idx--;

				int opt_bin_idx = get_optimal_bin(true_mz_bin_idx, charge);

				
				vector<RankBoostSample> spectrumRankSamples;
				fill_RankBoost_smaples_with_PMC(bs, charge, spectrumRankSamples);

				// select samples and add them to pmc_ds
				int goodIndex;
				vector<int> badIndexes;
				select_training_sample_idxs(charge, spectrumRankSamples, bs,goodIndex,badIndexes);

				const bool indAddToTraining = (myRandom()<prop_train);
				int groupIndex;
				
				if (indAddToTraining)
				{
					groupIndex= numTrainingGroups++;	
				}
				else
				{
					groupIndex= numTestingGroups++;
					testHeaders.push_back(ssf);
				}
				
				
				RankBoostDataset& ds = (indAddToTraining ? train_ds : test_ds);

				const int pos_index  = ds.get_num_samples();
				spectrumRankSamples[goodIndex].groupIndex = groupIndex;
				spectrumRankSamples[goodIndex].rank_in_group=0;

				ds.add_sample(spectrumRankSamples[goodIndex]);
				if (indSampleDiagnostic)
					pos_ds.add_sample(spectrumRankSamples[goodIndex]);

				int j;
				for (j=0; j<badIndexes.size(); j++)
				{
					const int bad_idx = badIndexes[j];
					if (bad_idx < 0 || bad_idx>= spectrumRankSamples.size())
						continue;
		
					spectrumRankSamples[bad_idx].groupIndex=groupIndex;
					spectrumRankSamples[bad_idx].rank_in_group=1;

					ds.add_to_phi_vector(ds.get_num_samples(),pos_index);
					ds.add_sample(spectrumRankSamples[bad_idx]);

					if (indSampleDiagnostic)
						neg_ds.add_sample(spectrumRankSamples[bad_idx]);
				}						   
			}

			train_ds.set_num_groups(numTrainingGroups);
			test_ds.set_num_groups(numTestingGroups);
			
			train_ds.compute_total_phi_weight();
			train_ds.initialize_potenital_lists();
			train_ds.initialzie_real_feature_table(realNames.size());

			test_ds.compute_total_phi_weight();

			if (pmc_rank_models[charge][sizeIndex])
				delete pmc_rank_models[charge][sizeIndex];
			
			pmc_rank_models[charge][sizeIndex] = new RankBoostModel;
		

			RankBoostModel* boost = pmc_rank_models[charge][sizeIndex];

			vector<string> empty;
			empty.clear();
			boost->init_rankboost_model_feature_names(empty,realNames);
			boost->init_rankboost_model_for_training(train_ds,100,25);

			train_ds.initialize_real_vote_lists(*boost);

	/*		if (indSampleDiagnostic)
			{
				boost->summarize_features_pos_neg(pos_ds.get_samples(),neg_ds.get_samples());
			}
			else
				boost->summarize_features(train_ds.get_samples());*/

			boost->train_rankboost_model(train_ds,4000,NULL,&test_ds);
			
			boost->ouput_ranked_feature_list();

		//	output_pmc_rank_results(fm,charge,testHeaders);

		//	exit(0);

			ind_initialized_pmcr = true;
		//	string path;
		//	path = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCRtt.txt";
		//	this->write_pmc_rank_models(path.c_str());
			
		}
	}

	string path;
	path = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCR.txt";
	this->write_pmc_rank_models(path.c_str());
	ind_initialized_pmcr = true;
}


struct offset_pair {
	offset_pair() : offset(POS_INF), inten_sum(0) {};
	offset_pair(mass_t off,float inten) : offset(off), inten_sum(inten) {};
	mass_t offset;
	float inten_sum;
};

bool cmp_offset_pair_offset (const offset_pair& a, const offset_pair& b)
{
	return (a.offset<b.offset);
}

bool cmp_offset_pair_inten (const offset_pair& a, const offset_pair& b)
{
	return (a.inten_sum>b.inten_sum);
}


float calc_mean_abs_offset(const vector<float>& offsets_by_inten)
{
	const float missing_pair_offset = 0.5;
	const int   num_offsets         = 3;

	if (offsets_by_inten.size()==0)
		return 1000;

	float abs_off=0;
	int i;
	for (i=0; i<num_offsets && i<offsets_by_inten.size(); i++)
		abs_off+=fabs(offsets_by_inten[i]);

	abs_off += (3-i)*missing_pair_offset;
	
	return (abs_off/num_offsets);
}


void calculatePmcRankStatisticsForMass(const PeakList& pl, 
									   mass_t SingleChargeFragmentPairSum,
									   mass_t tolerance, 
									   const vector<float>& isotopicLevels,
									   const vector<bool>& strictIsotopicIndicators,
									   PMCRankStats& statistics)
{
	const mass_t minSingleSum = SingleChargeFragmentPairSum - tolerance;
	const mass_t maxSingleSum = SingleChargeFragmentPairSum + tolerance;

	const mass_t minDoubleSum = minSingleSum + MASS_PROTON;
	const mass_t maxDoubleSum = maxSingleSum + MASS_PROTON;
	const mass_t doubleChargeFragmentPairSum = SingleChargeFragmentPairSum + MASS_PROTON;

	const mass_t minSingleSumWithH2OLoss = minSingleSum - MASS_H2O;
	const mass_t maxSingleSumWithH2OLoss = maxSingleSum - MASS_H2O;
	const mass_t single_charge_pair_h2o_sum = SingleChargeFragmentPairSum - MASS_H2O;

	const mass_t min_double_h2o_sum = minDoubleSum - MASS_H2O;
	const mass_t max_double_h2o_sum = maxDoubleSum - MASS_H2O;
	const mass_t double_charge_pair_h2o_sum = doubleChargeFragmentPairSum - MASS_H2O;

	static vector<offset_pair> by_pairs,  strong_pairs;
	static vector<offset_pair> c2_pairs,  strong_c2_pairs;
	static vector<offset_pair> h2o_pairs, c2_h2o_pairs;

	by_pairs.clear();
	strong_pairs.clear();
	c2_pairs.clear();
	strong_c2_pairs.clear();
	h2o_pairs.clear();
	c2_h2o_pairs.clear();

	statistics.clear();

	const Peak* const peaks = pl.getPeaks();
	const int numPeaks = pl.getNumPeaks();
	int forwardIndex = -1;
	int backIndex = numPeaks-1;


	// find pairs of b/y
	while (forwardIndex<backIndex)
	{
		++forwardIndex;
		if (isotopicLevels[forwardIndex]>0)
			continue;

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSum)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSum)
		{
			if (isotopicLevels[backIndex]>0)
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - SingleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			by_pairs.push_back(offset_pair(offset,inten_sum));
			statistics.inten_frag_pairs += inten_sum;

			if (isotopicLevels[forwardIndex] || isotopicLevels[backIndex])
			{
				strong_pairs.push_back(offset_pair(offset,inten_sum));
				statistics.inten_strong_pairs += inten_sum;
			}
		}
	}

	// find pairs b/y2
	forwardIndex = -1;
	backIndex = numPeaks-1;

	const int last_idx = backIndex;
	while (forwardIndex<last_idx && backIndex>=0)
	{
		forwardIndex++;
		if (isotopicLevels[forwardIndex]>0)
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (sum>maxDoubleSum)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minDoubleSum)
		{
			if (isotopicLevels[backIndex]>0)
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs.push_back(offset_pair(offset,inten_sum));
			statistics.inten_c2_pairs += inten_sum;

			if (isotopicLevels[forwardIndex] || isotopicLevels[backIndex])
			{
				strong_c2_pairs.push_back(offset_pair(offset,inten_sum));
				statistics.inten_c2_strong_pairs = inten_sum;
			}
		}
	}

	// find pairs of b/y-H2O
	forwardIndex = -1;
	backIndex = numPeaks-1;

	while (forwardIndex<backIndex)
	{
		forwardIndex++;
		if (isotopicLevels[forwardIndex]>0)
			continue;

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSumWithH2OLoss)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSumWithH2OLoss)
		{
			if (isotopicLevels[backIndex]>0)
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - single_charge_pair_h2o_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			h2o_pairs.push_back(offset_pair(offset,inten_sum));
			statistics.inten_h2o_loss_frag_pairs += inten_sum;
		}
	}

	// find pairs b/y2 - H2O
	forwardIndex = -1;
	backIndex = numPeaks-1;

	while (forwardIndex<last_idx && backIndex>=0)
	{
		forwardIndex++;
		if (isotopicLevels[forwardIndex]>0)
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (sum>max_double_h2o_sum)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > min_double_h2o_sum)
		{
			if (isotopicLevels[backIndex]>0)
				continue;

			const mass_t offset = fabs(sum - double_charge_pair_h2o_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_h2o_pairs.push_back(offset_pair(offset,inten_sum));
			statistics.itnen_h2o_loss_c2_frag_pairs += inten_sum;
		}
	}

	statistics.numFragmentPairs = by_pairs.size();
	statistics.numStrongFragmentPairs = strong_pairs.size();
	statistics.numCharge2FragmentPairs = c2_pairs.size();
	statistics.numStrongCharge2FragmentPairs = strong_c2_pairs.size();
	statistics.num_h2o_loss_frag_pairs = h2o_pairs.size();
	statistics.num_h2o_loss_c2_frag_pairs = c2_h2o_pairs.size();

	int i;

	vector<float>& offset_pairs_ordered_by_inten = statistics.offset_pairs_ordered_by_inten;
	sort(by_pairs.begin(),by_pairs.end(),cmp_offset_pair_inten);
	offset_pairs_ordered_by_inten.resize(by_pairs.size());
	for (i=0; i<by_pairs.size(); i++)
		offset_pairs_ordered_by_inten[i]=by_pairs[i].offset;
	statistics.mean_offset_pairs=calc_mean_abs_offset(offset_pairs_ordered_by_inten);

	vector<float>& strong_offset_pairs_ordered_by_inten = statistics.strong_offset_pairs_ordered_by_inten;
	sort(strong_pairs.begin(),strong_pairs.end(),cmp_offset_pair_inten);
	strong_offset_pairs_ordered_by_inten.resize(strong_pairs.size());
	for (i=0; i<strong_pairs.size(); i++)
		strong_offset_pairs_ordered_by_inten[i]=strong_pairs[i].offset;
	statistics.mean_offset_strong_pairs=calc_mean_abs_offset(strong_offset_pairs_ordered_by_inten);

	vector<float>& c2_offset_pairs_ordered_by_inten = statistics.c2_offset_pairs_ordered_by_inten;
	sort(c2_pairs.begin(),c2_pairs.end(),cmp_offset_pair_inten);
	c2_offset_pairs_ordered_by_inten.resize(c2_pairs.size());
	for (i=0; i<c2_pairs.size(); i++)
		c2_offset_pairs_ordered_by_inten[i]=c2_pairs[i].offset;
	statistics.mean_offset_c2_pairs=calc_mean_abs_offset(c2_offset_pairs_ordered_by_inten);



	// fill in additional iso sum features (look at pairs that sum to expected, expected+1 expected+2)
	
	// find pairs of b/y

	static vector<offset_pair> pairs0,  pairs1, pairs2;
	static vector<offset_pair> c2_pairs0,  c2_pairs1, c2_pairs2;
	
	pairs0.clear();
	forwardIndex = -1;
	backIndex = numPeaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
		if (strictIsotopicIndicators[forwardIndex])
			continue;

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSum)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSum)
		{
			if (strictIsotopicIndicators[backIndex])
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - SingleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs0.push_back(offset_pair(offset,inten_sum));
		}
	}

	pairs1.clear();
	forwardIndex = -1;
	backIndex = numPeaks-1;
	const mass_t max1 = maxSingleSum+1.0;
	const mass_t min1 = minSingleSum+1.0;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>max1)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > min1)
		{
			if (! (strictIsotopicIndicators[backIndex] || strictIsotopicIndicators[forwardIndex]))
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - SingleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs1.push_back(offset_pair(offset,inten_sum));
		}
	}

	pairs2.clear();
	forwardIndex = -1;
	backIndex = numPeaks-1;
	const mass_t max2 = maxSingleSum+2.0;
	const mass_t min2 = minSingleSum+2.0;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>max2)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > min2)
		{
			if (! (strictIsotopicIndicators[backIndex] || strictIsotopicIndicators[forwardIndex]))
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - SingleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs2.push_back(offset_pair(offset,inten_sum));
		}
	}


	c2_pairs0.clear();
	forwardIndex = -1;
	backIndex = numPeaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
		if (strictIsotopicIndicators[forwardIndex])
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxDoubleSum)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minDoubleSum)
		{
			if (strictIsotopicIndicators[backIndex])
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs0.push_back(offset_pair(offset,inten_sum));
		}
	}

	c2_pairs1.clear();
	const mass_t maxc21 = maxDoubleSum + 1.0;
	const mass_t minc21 = minDoubleSum + 1.0;
	forwardIndex = -1;
	backIndex = numPeaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxc21)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minc21)
		{
			if (! (strictIsotopicIndicators[backIndex] || strictIsotopicIndicators[forwardIndex]) )
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs1.push_back(offset_pair(offset,inten_sum));
		}
	}


	c2_pairs2.clear();
	const mass_t maxc22 = maxDoubleSum + 2.0;
	const mass_t minc22 = minDoubleSum + 2.0;
	forwardIndex = -1;
	backIndex = numPeaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxc22)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minc22)
		{
			if (! (strictIsotopicIndicators[backIndex] || strictIsotopicIndicators[forwardIndex]) )
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs2.push_back(offset_pair(offset,inten_sum));
		}
	}

	// use the first 4 peaks
	statistics.inten_strict_pairs0=0;
	statistics.num_strict_pairs0 = pairs0.size();
	sort(pairs0.begin(),pairs0.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs0.size(); i++)
	{
//		statistics.offset_strict_pairs0.push_back(pairs0[i].offset);
		statistics.inten_strict_pairs0+=pairs0[i].inten_sum;
	}


	statistics.inten_strict_pairs1=0;
	statistics.num_strict_pairs1 = pairs1.size();
	sort(pairs1.begin(),pairs1.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs1.size(); i++)
	{
//		statistics.offset_strict_pairs1.push_back(pairs1[i].offset);
		statistics.inten_strict_pairs1+=pairs1[i].inten_sum;
	}

	statistics.inten_strict_pairs2=0;
	statistics.num_strict_pairs2 = pairs2.size();
	sort(pairs2.begin(),pairs2.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs2.size(); i++)
	{
//		statistics.offset_strict_pairs2.push_back(pairs2[i].offset);
		statistics.inten_strict_pairs2+=pairs2[i].inten_sum;
	}
	
	statistics.c2_inten_strict_pairs0=0;
	statistics.c2_num_strict_pairs0 = c2_pairs0.size();
	sort(c2_pairs0.begin(),c2_pairs0.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs0.size(); i++)
	{
//		statistics.c2_offset_strict_pairs0.push_back(c2_pairs0[i].offset);
		statistics.c2_inten_strict_pairs0+=c2_pairs0[i].inten_sum;
	}
	
	statistics.c2_inten_strict_pairs1=0;
	statistics.c2_num_strict_pairs1 = c2_pairs1.size();
	sort(c2_pairs1.begin(),c2_pairs1.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs1.size(); i++)
	{
//		statistics.c2_offset_strict_pairs1.push_back(c2_pairs1[i].offset);
		statistics.c2_inten_strict_pairs1+=c2_pairs1[i].inten_sum;
	}
	
	statistics.c2_inten_strict_pairs2=0;
	statistics.c2_num_strict_pairs2 = c2_pairs2.size();
	sort(c2_pairs2.begin(), c2_pairs2.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs2.size(); i++)
	{
//		statistics.c2_offset_strict_pairs2.push_back(c2_pairs2[i].offset);
		statistics.c2_inten_strict_pairs2+=c2_pairs2[i].inten_sum;
	}

}



void calc_pmc_rank_stats_for_mass(const QCPeak *peaks, 
										  int num_peaks, 
										  mass_t single_charge_pair_sum,
										  mass_t tolerance, 
										  const vector<float>& iso_levels,
										  const vector<bool>& strong_inds,
										  const vector<bool>& strict_iso_inds,
										  PMCRankStats& stats)
{
	const mass_t minSingleSum = single_charge_pair_sum - tolerance;
	const mass_t maxSingleSum = single_charge_pair_sum + tolerance;

	const mass_t minDoubleSum = minSingleSum + 1.0;
	const mass_t maxDoubleSum = maxSingleSum + 1.0;
	const mass_t doubleChargeFragmentPairSum = single_charge_pair_sum +1.0;

	const mass_t minSingleSumWithH2OLoss = minSingleSum - MASS_H2O;
	const mass_t maxSingleSumWithH2OLoss = maxSingleSum - MASS_H2O;
	const mass_t single_charge_pair_h2o_sum = single_charge_pair_sum - MASS_H2O;

	const mass_t min_double_h2o_sum = minDoubleSum - MASS_H2O;
	const mass_t max_double_h2o_sum = maxDoubleSum - MASS_H2O;
	const mass_t double_charge_pair_h2o_sum = doubleChargeFragmentPairSum - MASS_H2O;

	static vector<offset_pair> by_pairs,  strong_pairs;
	static vector<offset_pair> c2_pairs,  strong_c2_pairs;
	static vector<offset_pair> h2o_pairs, c2_h2o_pairs;

	by_pairs.clear();
	strong_pairs.clear();
	c2_pairs.clear();
	strong_c2_pairs.clear();
	h2o_pairs.clear();
	c2_h2o_pairs.clear();

	stats.clear();

	int forwardIndex = -1;
	int backIndex = num_peaks-1;

	// find pairs of b/y
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
		if (iso_levels[forwardIndex]>0)
		{
			continue;
		}

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSum)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSum)
		{
			if (iso_levels[backIndex]>0)
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - single_charge_pair_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			by_pairs.push_back(offset_pair(offset,inten_sum));
			stats.inten_frag_pairs += inten_sum;

			if (strong_inds[forwardIndex] || strong_inds[backIndex])
			{
				strong_pairs.push_back(offset_pair(offset,inten_sum));
				stats.inten_strong_pairs += inten_sum;
			}
		}
	}

	// find pairs b/y2
	forwardIndex = -1;
	backIndex = num_peaks-1;

	const int last_idx =num_peaks-1;
	while (forwardIndex<last_idx && backIndex>=0)
	{
		forwardIndex++;
		if (iso_levels[forwardIndex]>0)
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (sum>maxDoubleSum)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minDoubleSum)
		{
			if (iso_levels[backIndex]>0)
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs.push_back(offset_pair(offset,inten_sum));
			stats.inten_c2_pairs += inten_sum;

			if (strong_inds[forwardIndex] || strong_inds[backIndex])
			{
				strong_c2_pairs.push_back(offset_pair(offset,inten_sum));
				stats.inten_c2_strong_pairs = inten_sum;
			}
		}
	}

	// find pairs of b/y-H2O
	forwardIndex = -1;
	backIndex = num_peaks-1;

	while (forwardIndex<backIndex)
	{
		forwardIndex++;
		if (iso_levels[forwardIndex]>0)
			continue;

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSumWithH2OLoss)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSumWithH2OLoss)
		{
			if (iso_levels[backIndex]>0)
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - single_charge_pair_h2o_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			h2o_pairs.push_back(offset_pair(offset,inten_sum));
			stats.inten_h2o_loss_frag_pairs += inten_sum;
		}
	}

	// find pairs b/y2 - H2O
	forwardIndex = -1;
	backIndex = num_peaks-1;

	while (forwardIndex<last_idx && backIndex>=0)
	{
		forwardIndex++;
		if (iso_levels[forwardIndex]>0)
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (sum>max_double_h2o_sum)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > min_double_h2o_sum)
		{
			if (iso_levels[backIndex]>0)
				continue;

			const mass_t offset = fabs(sum - double_charge_pair_h2o_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_h2o_pairs.push_back(offset_pair(offset,inten_sum));
			stats.itnen_h2o_loss_c2_frag_pairs += inten_sum;
		}
	}

	stats.numFragmentPairs = by_pairs.size();
	stats.numStrongFragmentPairs = strong_pairs.size();
	stats.numCharge2FragmentPairs = c2_pairs.size();
	stats.numStrongCharge2FragmentPairs = strong_c2_pairs.size();
	stats.num_h2o_loss_frag_pairs = h2o_pairs.size();
	stats.num_h2o_loss_c2_frag_pairs = c2_h2o_pairs.size();

	int i;

	vector<float>& offset_pairs_ordered_by_inten = stats.offset_pairs_ordered_by_inten;
	sort(by_pairs.begin(),by_pairs.end(),cmp_offset_pair_inten);
	offset_pairs_ordered_by_inten.resize(by_pairs.size());
	for (i=0; i<by_pairs.size(); i++)
		offset_pairs_ordered_by_inten[i]=by_pairs[i].offset;
	stats.mean_offset_pairs=calc_mean_abs_offset(offset_pairs_ordered_by_inten);

	vector<float>& strong_offset_pairs_ordered_by_inten = stats.strong_offset_pairs_ordered_by_inten;
	sort(strong_pairs.begin(),strong_pairs.end(),cmp_offset_pair_inten);
	strong_offset_pairs_ordered_by_inten.resize(strong_pairs.size());
	for (i=0; i<strong_pairs.size(); i++)
		strong_offset_pairs_ordered_by_inten[i]=strong_pairs[i].offset;
	stats.mean_offset_strong_pairs=calc_mean_abs_offset(strong_offset_pairs_ordered_by_inten);

	vector<float>& c2_offset_pairs_ordered_by_inten = stats.c2_offset_pairs_ordered_by_inten;
	sort(c2_pairs.begin(),c2_pairs.end(),cmp_offset_pair_inten);
	c2_offset_pairs_ordered_by_inten.resize(c2_pairs.size());
	for (i=0; i<c2_pairs.size(); i++)
		c2_offset_pairs_ordered_by_inten[i]=c2_pairs[i].offset;
	stats.mean_offset_c2_pairs=calc_mean_abs_offset(c2_offset_pairs_ordered_by_inten);



	// fill in additional iso sum features (look at pairs that sum to expected, expected+1 expected+2)
	
	// find pairs of b/y

	static vector<offset_pair> pairs0,  pairs1, pairs2;
	static vector<offset_pair> c2_pairs0,  c2_pairs1, c2_pairs2;
	
	pairs0.clear();
	forwardIndex = -1;
	backIndex = num_peaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
		if (strict_iso_inds[forwardIndex])
			continue;

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSum)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSum)
		{
			if (strict_iso_inds[backIndex])
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - single_charge_pair_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs0.push_back(offset_pair(offset,inten_sum));
		}
	}

	pairs1.clear();
	forwardIndex = -1;
	backIndex = num_peaks-1;
	const mass_t max1 = maxSingleSum+1.0;
	const mass_t min1 = minSingleSum+1.0;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>max1)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > min1)
		{
			if (! (strict_iso_inds[backIndex] || strict_iso_inds[forwardIndex]))
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - single_charge_pair_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs1.push_back(offset_pair(offset,inten_sum));
		}
	}

	pairs2.clear();
	forwardIndex = -1;
	backIndex = num_peaks-1;
	const mass_t max2 = maxSingleSum+2.0;
	const mass_t min2 = minSingleSum+2.0;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>max2)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > min2)
		{
			if (! (strict_iso_inds[backIndex] || strict_iso_inds[forwardIndex]))
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - single_charge_pair_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs2.push_back(offset_pair(offset,inten_sum));
		}
	}


	c2_pairs0.clear();
	forwardIndex = -1;
	backIndex = num_peaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
		if (strict_iso_inds[forwardIndex])
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxDoubleSum)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minDoubleSum)
		{
			if (strict_iso_inds[backIndex])
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs0.push_back(offset_pair(offset,inten_sum));
		}
	}

	c2_pairs1.clear();
	const mass_t maxc21 = maxDoubleSum + 1.0;
	const mass_t minc21 = minDoubleSum + 1.0;
	forwardIndex = -1;
	backIndex = num_peaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxc21)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minc21)
		{
			if (! (strict_iso_inds[backIndex] || strict_iso_inds[forwardIndex]) )
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs1.push_back(offset_pair(offset,inten_sum));
		}
	}


	c2_pairs2.clear();
	const mass_t maxc22 = maxDoubleSum + 2.0;
	const mass_t minc22 = minDoubleSum + 2.0;
	forwardIndex = -1;
	backIndex = num_peaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxc22)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minc22)
		{
			if (! (strict_iso_inds[backIndex] || strict_iso_inds[forwardIndex]) )
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs2.push_back(offset_pair(offset,inten_sum));
		}
	}

	// use the first 4 peaks
	stats.inten_strict_pairs0=0;
	stats.num_strict_pairs0 = pairs0.size();
	sort(pairs0.begin(),pairs0.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs0.size(); i++)
	{
//		stats.offset_strict_pairs0.push_back(pairs0[i].offset);
		stats.inten_strict_pairs0+=pairs0[i].inten_sum;
	}


	stats.inten_strict_pairs1=0;
	stats.num_strict_pairs1 = pairs1.size();
	sort(pairs1.begin(),pairs1.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs1.size(); i++)
	{
//		stats.offset_strict_pairs1.push_back(pairs1[i].offset);
		stats.inten_strict_pairs1+=pairs1[i].inten_sum;
	}

	stats.inten_strict_pairs2=0;
	stats.num_strict_pairs2 = pairs2.size();
	sort(pairs2.begin(),pairs2.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs2.size(); i++)
	{
//		stats.offset_strict_pairs2.push_back(pairs2[i].offset);
		stats.inten_strict_pairs2+=pairs2[i].inten_sum;
	}
	
	stats.c2_inten_strict_pairs0=0;
	stats.c2_num_strict_pairs0 = c2_pairs0.size();
	sort(c2_pairs0.begin(),c2_pairs0.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs0.size(); i++)
	{
//		stats.c2_offset_strict_pairs0.push_back(c2_pairs0[i].offset);
		stats.c2_inten_strict_pairs0+=c2_pairs0[i].inten_sum;
	}
	
	stats.c2_inten_strict_pairs1=0;
	stats.c2_num_strict_pairs1 = c2_pairs1.size();
	sort(c2_pairs1.begin(),c2_pairs1.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs1.size(); i++)
	{
//		stats.c2_offset_strict_pairs1.push_back(c2_pairs1[i].offset);
		stats.c2_inten_strict_pairs1+=c2_pairs1[i].inten_sum;
	}
	
	stats.c2_inten_strict_pairs2=0;
	stats.c2_num_strict_pairs2 = c2_pairs2.size();
	sort(c2_pairs2.begin(), c2_pairs2.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs2.size(); i++)
	{
//		stats.c2_offset_strict_pairs2.push_back(c2_pairs2[i].offset);
		stats.c2_inten_strict_pairs2+=c2_pairs2[i].inten_sum;
	}

}


void fillPmcRankStatistics(int charge,
						  const mass_t singleChargePairOffset, // the sum of b+y or c+z
						  mass_t minusRange, 
						  mass_t plusRange,
						  mass_t increment,
						  const Config* config,
						  const PeakList& pl,
						  const vector<bool>&   strongPeakIndicators,
						  const vector<float>&  isotopicLevels,
						  const vector<bool>&   strictIsotopicIndicators,
						  vector<PMCRankStats>& pmcStatisticsVector)
{
	const mass_t tolerance   = config->getTolerance()*0.55;
	const int numBinsPerDa   = static_cast<int>(1.0/increment);
	const int numMinusBins   = static_cast<int>((-minusRange) * numBinsPerDa) + 1;
	const int numPlusBins    = static_cast<int>(plusRange * numBinsPerDa) + 1;
	const mass_t oneOverCharge = 1.0/static_cast<mass_t>(charge);

	const int totalNumBins = numMinusBins + numPlusBins;
	if (pmcStatisticsVector.size() != totalNumBins)
		pmcStatisticsVector.resize(totalNumBins);

	int i;
	for (i=0; i<numMinusBins; i++)
	{	
		const mass_t delta = static_cast<mass_t>(i - numMinusBins)*increment;

		calculatePmcRankStatisticsForMass(pl, singleChargePairOffset + delta, tolerance, 
			isotopicLevels, strictIsotopicIndicators, pmcStatisticsVector[i]);

		pmcStatisticsVector[i].m_over_z = (singleChargePairOffset + delta + (charge-2)*MASS_PROTON)/
										  static_cast<mass_t>(charge);
	}

	const int startBinIndex = i;
	for (i=0; i<numPlusBins; i++)
	{
		const mass_t delta = static_cast<mass_t>(i)*increment;

		calculatePmcRankStatisticsForMass(pl, singleChargePairOffset + delta, tolerance, 
			isotopicLevels, //strongPeakIndicators, 
			strictIsotopicIndicators, pmcStatisticsVector[startBinIndex+i]);

		pmcStatisticsVector[startBinIndex+i].m_over_z = (singleChargePairOffset+delta+(charge-2)*MASS_PROTON)/
														  static_cast<mass_t>(charge);
	}
}


void fill_rank_PMC_stats(int charge,
						  const mass_t single_charge_pair_sum, // the sum of b+y or c+z
						  mass_t minus_range, 
						  mass_t plus_range,
						  mass_t increment,
						  Config *config,
						  const BasicSpectrum& bs,
						  const vector<bool>& strong_inds,
						  const vector<float>& iso_levels,
						  const vector<bool>& iso_inds,
						  vector<PMCRankStats>& pmc_stats_vec)
{

	const mass_t tolerance = config->getTolerance()*0.55;
	const int num_bins_per_Da = (int)(1.0/increment);
	const int num_minus_bins = (int)((-minus_range)*num_bins_per_Da)+1;
	const int num_plus_bins = (int)(plus_range*num_bins_per_Da)+1;
	const mass_t one_over_charge = 1.0/(mass_t)charge;

	const int total_num_bins = num_minus_bins + num_plus_bins;
	if (pmc_stats_vec.size() != total_num_bins)
		pmc_stats_vec.resize(total_num_bins);

	int i;
	for (i=0; i<num_minus_bins; i++)
	{	
		const mass_t delta = (i - num_minus_bins)*increment;

		calc_pmc_rank_stats_for_mass(bs.peaks,bs.num_peaks, single_charge_pair_sum+delta,
			tolerance, iso_levels, strong_inds, iso_inds, pmc_stats_vec[i]);

		pmc_stats_vec[i].m_over_z = (single_charge_pair_sum+delta+(charge-2)*MASS_PROTON)/charge;
	}

	const int start_bin_idx = i;

	for (i=0; i<num_plus_bins; i++)
	{
		const mass_t delta = i*increment;

		calc_pmc_rank_stats_for_mass(bs.peaks,bs.num_peaks, single_charge_pair_sum+delta,
			tolerance, iso_levels, strong_inds, iso_inds, pmc_stats_vec[start_bin_idx+i]);

		pmc_stats_vec[start_bin_idx+i].m_over_z = (single_charge_pair_sum+delta+(charge-2)*MASS_PROTON)/charge;
	}
}




void PMCRankStats::clear()
{
	m_over_z=0;

	rank_score = NEG_INF;

	numFragmentPairs=0;
	numStrongFragmentPairs=0;
	numCharge2FragmentPairs=0;
	numStrongCharge2FragmentPairs=0;
	num_h2o_loss_frag_pairs=0;

	inten_frag_pairs=0;
	inten_strong_pairs=0;
	inten_c2_pairs=0;
	inten_c2_strong_pairs=0;
	inten_h2o_loss_frag_pairs=0;
	itnen_h2o_loss_c2_frag_pairs=0;

	mean_offset_pairs=0;
	mean_offset_strong_pairs=0;
	mean_offset_c2_pairs=0;
	mean_offset_c2_strong_pairs=0;
	mean_offset_h2o_pairs=0;
	mean_offset_c2_h2o_pairs=0;

	ind_pairs_with_min_tol=false;			 
	ind_strong_pairs_with_min_tol=false;
	ind_c2_pairs_with_min_tol=false;
	ind_c2_strong_pairs_with_min_tol=false;
	log_dis_from_pairs_min_tol=0;			 
	log_dis_from_strong_pairs_min_tol=0;
	log_dis_from_c2_pairs_min_tol=0;		 
	log_dis_from_c2_strong_pairs_min_tol=0;

	offset_pairs_ordered_by_inten.clear();
	strong_offset_pairs_ordered_by_inten.clear();
	c2_offset_pairs_ordered_by_inten.clear();


	num_strict_pairs0=0; inten_strict_pairs0=0;
	num_strict_pairs1=0; inten_strict_pairs1=0;
	num_strict_pairs2=0; inten_strict_pairs2=0;
}


void PMCSQS_Scorer::fillRankboostPmcSamples(const PeakList& pl,
								 int charge,
								 vector<RankBoostSample>& samples) const
{
	const int numSamples = currentSpectrumPmcTables_[charge].size();
	const int idx_skip = int((1.0/bin_increment)+0.00001);
	vector<int> idx_offsets;
	int i;

	idx_offsets.clear();
	idx_offsets.push_back(-2*idx_skip);
	idx_offsets.push_back(-1*idx_skip);
	idx_offsets.push_back(idx_skip);
	idx_offsets.push_back(2*idx_skip);

	if (samples.size() != numSamples)
		samples.resize(numSamples);

	for (i=0; i<numSamples; i++)
	{
		const PMCRankStats& stats = currentSpectrumPmcTables_[charge][i];
		RankBoostSample& sam = samples[i];

		const float inten_norm = 1.0/(curr_spec_total_intensity+1.0);
		int r_idx=0;
		const mass_t mz_offset = (stats.m_over_z - pl.getHeader()->getMOverZ());

		sam.clear();
		sam.add_real_feature(r_idx++,mz_offset);

		if (stats.numFragmentPairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.numFragmentPairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.numStrongFragmentPairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

		if (stats.numCharge2FragmentPairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.numCharge2FragmentPairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.numStrongCharge2FragmentPairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

			
	/*	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=2");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS >5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS <4");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS >4");

		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=2");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS >5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS <4");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS >4");*/

		sam.add_real_feature(r_idx++,stats.numFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numStrongFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numCharge2FragmentPairs);
		sam.add_real_feature(r_idx++,stats.numStrongCharge2FragmentPairs);
		sam.add_real_feature(r_idx++,stats.num_h2o_loss_frag_pairs);
		sam.add_real_feature(r_idx++,stats.num_h2o_loss_c2_frag_pairs);

		sam.add_real_feature(r_idx++,stats.inten_frag_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_strong_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_c2_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_c2_strong_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_h2o_loss_frag_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.itnen_h2o_loss_c2_frag_pairs * inten_norm);

		// averages of top k offsets

		float avg=0;
		int j;
		for (j =0; j<7 && j<stats.offset_pairs_ordered_by_inten.size(); j++)
		{
			avg += fabs(stats.offset_pairs_ordered_by_inten[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;

		avg=0;
		for (j =0; j<7 && j<stats.c2_offset_pairs_ordered_by_inten.size(); j++)
		{
			avg += fabs(stats.c2_offset_pairs_ordered_by_inten[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;


		// offset data
	
		if (stats.mean_offset_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_pairs/(1.0+stats.numFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_strong_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_strong_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_strong_pairs/(1.0+stats.numStrongFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_pairs/(1.0+stats.numCharge2FragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_strong_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_strong_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_strong_pairs/(1.0+stats.numStrongCharge2FragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_h2o_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_h2o_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_h2o_pairs/(1.0+stats.num_h2o_loss_frag_pairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_h2o_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_h2o_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_h2o_pairs/(1.0+stats.num_h2o_loss_c2_frag_pairs));
		}
		else
			r_idx+=2;

		// individual offsets
		for (j=0; j<5 && j<stats.offset_pairs_ordered_by_inten.size(); j++)
			sam.add_real_feature(r_idx+j,stats.offset_pairs_ordered_by_inten[j]);
		r_idx+=5;

		for (j=0; j<5 && j<stats.c2_offset_pairs_ordered_by_inten.size(); j++)
			sam.add_real_feature(r_idx+j,stats.c2_offset_pairs_ordered_by_inten[j]);
		r_idx+=5;

	
		// add the +0 +1 +2 strict counts
		sam.add_real_feature(r_idx++,stats.num_strict_pairs0);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.num_strict_pairs1);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.num_strict_pairs2);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs2 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs0);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs1);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs2);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs2 * inten_norm);

		// add comparative features to -2 -1 +1 +2 Da away
		for (j=0; j<idx_offsets.size(); j++)
		{
			const int other_idx = i + idx_offsets[j];
			if (other_idx<0 || other_idx>= samples.size())
			{
				r_idx+=12;
				continue;
			}

			const PMCRankStats& other = currentSpectrumPmcTables_[charge][other_idx];

			sam.add_real_feature(r_idx++,stats.numFragmentPairs - other.numFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numStrongFragmentPairs - other.numStrongFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numCharge2FragmentPairs - other.numCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,stats.numStrongCharge2FragmentPairs - other.numStrongCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,stats.num_h2o_loss_frag_pairs - other.num_h2o_loss_frag_pairs);
			sam.add_real_feature(r_idx++,stats.num_h2o_loss_c2_frag_pairs - other.num_h2o_loss_c2_frag_pairs);

			sam.add_real_feature(r_idx++,(stats.inten_frag_pairs - other.inten_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_strong_pairs - other.inten_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_c2_pairs - other.inten_c2_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_c2_strong_pairs - other.inten_c2_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_h2o_loss_frag_pairs - other.inten_h2o_loss_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.itnen_h2o_loss_c2_frag_pairs - other.itnen_h2o_loss_c2_frag_pairs) * inten_norm);
		}

		const int plus_idx = i + idx_skip;
		const int minus_idx = i-idx_skip;

		if (plus_idx<samples.size() && minus_idx>0)
		{
			const PMCRankStats& plus = currentSpectrumPmcTables_[charge][plus_idx];
			const PMCRankStats& minus = currentSpectrumPmcTables_[charge][minus_idx];

			sam.add_real_feature(r_idx++,plus.numFragmentPairs - minus.numFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numStrongFragmentPairs - minus.numStrongFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numCharge2FragmentPairs - minus.numCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,plus.numStrongCharge2FragmentPairs - minus.numStrongCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,plus.num_h2o_loss_frag_pairs - minus.num_h2o_loss_frag_pairs);
			sam.add_real_feature(r_idx++,plus.num_h2o_loss_c2_frag_pairs - minus.num_h2o_loss_c2_frag_pairs);

			sam.add_real_feature(r_idx++,(plus.inten_frag_pairs - minus.inten_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_strong_pairs - minus.inten_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_c2_pairs - minus.inten_c2_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_c2_strong_pairs - minus.inten_c2_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_h2o_loss_frag_pairs - minus.inten_h2o_loss_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.itnen_h2o_loss_c2_frag_pairs - minus.itnen_h2o_loss_c2_frag_pairs) * inten_norm);
		}
	}
}


/**************************************************************************

  Fills in the RankBoost feature data
***************************************************************************/
void PMCSQS_Scorer::fill_RankBoost_smaples_with_PMC(
									const BasicSpectrum& bs,
									int charge,
									vector<RankBoostSample>& samples) const
{

	const int numSamples = currentSpectrumPmcTables_[charge].size();
	const int idx_skip = int((1.0/bin_increment)+0.00001);
	vector<int> idx_offsets;
	int i;

	idx_offsets.clear();
	idx_offsets.push_back(-2*idx_skip);
	idx_offsets.push_back(-1*idx_skip);
	idx_offsets.push_back(idx_skip);
	idx_offsets.push_back(2*idx_skip);

	if (samples.size() != numSamples)
		samples.resize(numSamples);

	for (i=0; i<numSamples; i++)
	{
		const PMCRankStats& stats = currentSpectrumPmcTables_[charge][i];
		RankBoostSample& sam = samples[i];

		const float inten_norm = 1.0/(curr_spec_total_intensity+1.0);
		int r_idx=0;
		const mass_t mz_offset = (stats.m_over_z - bs.ssf->m_over_z);

		sam.clear();
		sam.add_real_feature(r_idx++,mz_offset);

		if (stats.numFragmentPairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.numFragmentPairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.numStrongFragmentPairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

		if (stats.numCharge2FragmentPairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.numCharge2FragmentPairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.numStrongCharge2FragmentPairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

			
	/*	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=2");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS >5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS <4");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS >4");

		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=2");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS >5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS <4");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS >4");*/

		sam.add_real_feature(r_idx++,stats.numFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numStrongFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numCharge2FragmentPairs);
		sam.add_real_feature(r_idx++,stats.numStrongCharge2FragmentPairs);
		sam.add_real_feature(r_idx++,stats.num_h2o_loss_frag_pairs);
		sam.add_real_feature(r_idx++,stats.num_h2o_loss_c2_frag_pairs);

		sam.add_real_feature(r_idx++,stats.inten_frag_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_strong_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_c2_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_c2_strong_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_h2o_loss_frag_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.itnen_h2o_loss_c2_frag_pairs * inten_norm);

		// averages of top k offsets

		float avg=0;
		int j;
		for (j =0; j<7 && j<stats.offset_pairs_ordered_by_inten.size(); j++)
		{
			avg += fabs(stats.offset_pairs_ordered_by_inten[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;

		avg=0;
		for (j =0; j<7 && j<stats.c2_offset_pairs_ordered_by_inten.size(); j++)
		{
			avg += fabs(stats.c2_offset_pairs_ordered_by_inten[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;


		// offset data
	
		if (stats.mean_offset_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_pairs/(1.0+stats.numFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_strong_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_strong_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_strong_pairs/(1.0+stats.numStrongFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_pairs/(1.0+stats.numCharge2FragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_strong_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_strong_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_strong_pairs/(1.0+stats.numStrongCharge2FragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_h2o_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_h2o_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_h2o_pairs/(1.0+stats.num_h2o_loss_frag_pairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_h2o_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_h2o_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_h2o_pairs/(1.0+stats.num_h2o_loss_c2_frag_pairs));
		}
		else
			r_idx+=2;

		// individual offsets
		for (j=0; j<5 && j<stats.offset_pairs_ordered_by_inten.size(); j++)
			sam.add_real_feature(r_idx+j,stats.offset_pairs_ordered_by_inten[j]);
		r_idx+=5;

		for (j=0; j<5 && j<stats.c2_offset_pairs_ordered_by_inten.size(); j++)
			sam.add_real_feature(r_idx+j,stats.c2_offset_pairs_ordered_by_inten[j]);
		r_idx+=5;

	
		// add the +0 +1 +2 strict counts
		sam.add_real_feature(r_idx++,stats.num_strict_pairs0);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.num_strict_pairs1);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.num_strict_pairs2);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs2 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs0);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs1);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs2);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs2 * inten_norm);

		// add comparative features to -2 -1 +1 +2 Da away
		for (j=0; j<idx_offsets.size(); j++)
		{
			const int other_idx = i + idx_offsets[j];
			if (other_idx<0 || other_idx>= samples.size())
			{
				r_idx+=12;
				continue;
			}

			const PMCRankStats& other = currentSpectrumPmcTables_[charge][other_idx];

			sam.add_real_feature(r_idx++,stats.numFragmentPairs - other.numFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numStrongFragmentPairs - other.numStrongFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numCharge2FragmentPairs - other.numCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,stats.numStrongCharge2FragmentPairs - other.numStrongCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,stats.num_h2o_loss_frag_pairs - other.num_h2o_loss_frag_pairs);
			sam.add_real_feature(r_idx++,stats.num_h2o_loss_c2_frag_pairs - other.num_h2o_loss_c2_frag_pairs);

			sam.add_real_feature(r_idx++,(stats.inten_frag_pairs - other.inten_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_strong_pairs - other.inten_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_c2_pairs - other.inten_c2_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_c2_strong_pairs - other.inten_c2_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_h2o_loss_frag_pairs - other.inten_h2o_loss_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.itnen_h2o_loss_c2_frag_pairs - other.itnen_h2o_loss_c2_frag_pairs) * inten_norm);
		}

		const int plus_idx = i + idx_skip;
		const int minus_idx = i-idx_skip;

		if (plus_idx<samples.size() && minus_idx>0)
		{
			const PMCRankStats& plus = currentSpectrumPmcTables_[charge][plus_idx];
			const PMCRankStats& minus = currentSpectrumPmcTables_[charge][minus_idx];

			sam.add_real_feature(r_idx++,plus.numFragmentPairs - minus.numFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numStrongFragmentPairs - minus.numStrongFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numCharge2FragmentPairs - minus.numCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,plus.numStrongCharge2FragmentPairs - minus.numStrongCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,plus.num_h2o_loss_frag_pairs - minus.num_h2o_loss_frag_pairs);
			sam.add_real_feature(r_idx++,plus.num_h2o_loss_c2_frag_pairs - minus.num_h2o_loss_c2_frag_pairs);

			sam.add_real_feature(r_idx++,(plus.inten_frag_pairs - minus.inten_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_strong_pairs - minus.inten_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_c2_pairs - minus.inten_c2_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_c2_strong_pairs - minus.inten_c2_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_h2o_loss_frag_pairs - minus.inten_h2o_loss_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.itnen_h2o_loss_c2_frag_pairs - minus.itnen_h2o_loss_c2_frag_pairs) * inten_norm);
		}
	}
}


void init_PMC_feature_names(vector<string>& names)
{
	names.clear();
	int i;

	names.push_back("OFFSET FROM MEASURED M/Z");

	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=2");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS >4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS <3");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS >=3");

	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=2");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS >4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS <3");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS >=3");
	
	names.push_back("# PAIRS");
	names.push_back("# STRONG PAIRS");
	names.push_back("# C2 PAIRS");
	names.push_back("# STRONG C2 PAIRS");
	names.push_back("# H2O PAIRS");
	names.push_back("# C2 H2O PAIRS");

	names.push_back("INTEN PAIRS");
	names.push_back("INTEN STRONG PAIRS");
	names.push_back("INTEN C2 PAIRS");
	names.push_back("INTEN STRONG C2 PAIRS");
	names.push_back("INTEN H2O PAIRS");
	names.push_back("INTEN C2 H2O PAIRS");

	for (i=2; i<7; i++)
	{
		char name[64];
		sprintf(name,"AVG OFFSET TOP (STRONG %d)",i);
		names.push_back(name);
	}

	for (i=2; i<7; i++)
	{
		char name[64];
		sprintf(name,"AVG OFFSET TOP C2 (STRONG %d)",i);
		names.push_back(name);
	}

	names.push_back("MEAN OFFSET PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET PAIRS");

	names.push_back("MEAN OFFSET STRONG PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET STRONG PAIRS");

	names.push_back("MEAN OFFSET C2 PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET C2 PAIRS");

	names.push_back("MEAN OFFSET STRONG C2 PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET STRONG C2 PAIRS");

	names.push_back("MEAN OFFSET H2O PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET H2O PAIRS");

	names.push_back("MEAN OFFSET C2 H2O PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET C2 H2O PAIRS");

	for (i=0; i<5; i++)
	{
		char name[64];
		sprintf(name,"PAIR OFFSET (STRONG %d)",i+1);
		names.push_back(name);
	}

	for (i=0; i<5; i++)
	{
		char name[64];
		sprintf(name,"C2 PAIR OFFSET (STRONG %d)",i+1);
		names.push_back(name);
	}



	names.push_back("NUM STRICT 0");
	names.push_back("INTEN STRICT 0");

	names.push_back("NUM STRICT 1");
	names.push_back("INTEN STRICT 1");

	names.push_back("NUM STRICT 2");
	names.push_back("INTEN STRICT 2");

	names.push_back("NUM C2 STRICT 0");
	names.push_back("INTEN C2 STRICT 0");

	names.push_back("NUM C2 STRICT 1");
	names.push_back("INTEN C2 STRICT 1");

	names.push_back("NUM C2 STRICT 2");
	names.push_back("INTEN C2 STRICT 2");

	// diff features with -2 -1 +1 +2
	const string dis_labels[]={"-2","-1","+1","+2"};
	for (i=0; i<4; i++)
	{
		const string prefix = "DIFF WITH "+dis_labels[i]+" ";

		names.push_back(prefix+"# PAIRS");
		names.push_back(prefix+"# STRONG PAIRS");
		names.push_back(prefix+"# C2 PAIRS");
		names.push_back(prefix+"# STRONG C2 PAIRS");
		names.push_back(prefix+"# H2O PAIRS");
		names.push_back(prefix+"# C2 H2O PAIRS");

		names.push_back(prefix+"INTEN PAIRS");
		names.push_back(prefix+"INTEN STRONG PAIRS");
		names.push_back(prefix+"INTEN C2 PAIRS");
		names.push_back(prefix+"INTEN STRONG C2 PAIRS");
		names.push_back(prefix+"INTEN H2O PAIRS");
		names.push_back(prefix+"INTEN C2 H2O PAIRS");
	}

	names.push_back("DIFF +1/-1 # PAIRS");
	names.push_back("DIFF +1/-1 # STRONG PAIRS");
	names.push_back("DIFF +1/-1 # C2 PAIRS");
	names.push_back("DIFF +1/-1 # STRONG C2 PAIRS");
	names.push_back("DIFF +1/-1 # H2O PAIRS");
	names.push_back("DIFF +1/-1 # C2 H2O PAIRS");

	names.push_back("DIFF +1/-1 INTEN PAIRS");
	names.push_back("DIFF +1/-1 INTEN STRONG PAIRS");
	names.push_back("DIFF +1/-1 INTEN C2 PAIRS");
	names.push_back("DIFF +1/-1 INTEN STRONG C2 PAIRS");
	names.push_back("DIFF +1/-1 INTEN H2O PAIRS");
	names.push_back("DIFF +1/-1 INTEN C2 H2O PAIRS");
	cout << "Initialized: " << names.size() << " real feature names..." << endl;
}




void PMCSQS_Scorer::output_pmc_rank_results(const FileManager& fm,
											int charge,
											const vector<SingleSpectrumFile *>& testHeaders) 
{
	BasicSpecReader bsr;
	static QCPeak peaks[5000];

	vector<int> org_offset_counts, new_offset_counts;
	org_offset_counts.resize(201,0);
	new_offset_counts.resize(201,0);

	vector<mass_t> org_offsets;
	vector<mass_t> corr_offsets;

	org_offsets.clear();
	corr_offsets.clear();

	int i;
	for (i=0; i<testHeaders.size(); i++)
	{
		SingleSpectrumFile* ssf = testHeaders[i];
		BasicSpectrum bs;
	
	//	bs.num_peaks = bsr.read_basic_spec(config_,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.ssf = ssf;

	//	init_for_current_spec(config_,bs);
		calculate_curr_spec_pmc_values(bs, bin_increment);

		PmcSqsChargeRes res;
		find_best_mz_values_from_rank_model(bs, charge, config_->get_pm_tolerance(),res);

		ssf->peptide.calc_mass(config_);
		mass_t trueMz = (ssf->peptide.get_mass() + 18.01 + charge)/charge;

		org_offsets.push_back(trueMz - ssf->m_over_z);
		corr_offsets.push_back(trueMz - res.mz1);
	}

	mass_t m_org,sd_org,m_corr,sd_corr;
	calc_mean_sd(org_offsets,&m_org, &sd_org);
	calc_mean_sd(corr_offsets,&m_corr,&sd_corr);

	cout << "CHARGE: " << charge << endl;
	cout << "ORG:  mean " << m_org << " " << sd_org << endl;
	cout << "CORR: mean " << m_corr << " " << sd_corr << endl;

	for (i=0; i<org_offsets.size(); i++)
	{
		int org_idx = 100 + int(org_offsets[i] * 20);
		if (org_idx<0)
			org_idx = 0;
		if (org_idx>200)
			org_idx=200;
		org_offset_counts[org_idx]++;

		int new_idx = 100 + int(corr_offsets[i] * 20);
		if (new_idx<0)
			new_idx = 0;
		if (new_idx>200)
			new_idx=200;
		new_offset_counts[new_idx]++;
	}

	int cum_org=0;
	int cum_new=0;
	for (i=0; i<=200; i++)
	{

		if (org_offset_counts[i]==0 && new_offset_counts[i]==0)
			continue;
		
		cum_org+=org_offset_counts[i];
		cum_new+=new_offset_counts[i];
		cout << fixed << setprecision(3) << i*0.05 - 5.0 << "\t" <<
			org_offset_counts[i]/(float)org_offsets.size() << "\t" <<
			new_offset_counts[i]/(float)corr_offsets.size() << "\t" <<
			cum_org/(float)org_offsets.size() << "\t"<<
			cum_new/(float)corr_offsets.size() << endl;

	}


}


/**********************************************************************************
Finds the best m/z values for a given charge
If the pm_tolerance is low (less than 0.5) then the m/z vlaues need to reflect +-X Da
from the recorded m/z
***********************************************************************************/
void PMCSQS_Scorer::computeBestMzValuesForCharge(const PeakList& pl, 
	     						  int charge,
								  mass_t precursorMassTolerance,
								  PmcSqsChargeRes& result)
{
	if (charge<1)
	{
		cout << "Error: trying to find m/z values which carge < 1 !" << endl;
		exit(1);
	}

	static vector<RankBoostSample> spectrumRankSamples;
	static vector<float> rankingScores;

	const mass_t spectrumMOverZ  = pl.getHeader()->getMOverZ();
	const mass_t maxMzDiffAllowed = precursorMassTolerance / static_cast<mass_t>(charge);

	fillRankboostPmcSamples(pl, charge, spectrumRankSamples);
//	fill_RankBoost_smaples_with_PMC(bs, charge, spectrumRankSamples);

	if (rankingScores.size() != spectrumRankSamples.size())
		rankingScores.resize(spectrumRankSamples.size(),NEG_INF);

	const mass_t pm_with_19 = spectrumMOverZ * charge - (charge + 1);
	const int sizeIndex = get_rank_model_size_idx(charge, pm_with_19);
	
	int best_idx=-1;
	float best_score=NEG_INF;

	if (charge>= pmc_rank_models.size() ||
		sizeIndex>= pmc_rank_models[charge].size() ||
		! pmc_rank_models[charge][sizeIndex])
	{
		//
	}
	else
	{
	//	cout << "spec: " << spectrumMOverZ << "  charge: " << charge << endl;
		int i;
		for (i=0; i<spectrumRankSamples.size(); i++)
		{
			// make sure the m/z is in an allowed range (i.e., it has a mass that can possibly be a +-X Da
			// shift from the true pm). Assume shift can be at most
			if (precursorMassTolerance<0.5)
			{
				const mass_t table_mz = currentSpectrumPmcTables_[charge][i].m_over_z;
				const mass_t one_over_charge = 1.0033 / (mass_t)charge; // ~mass of isotopic peak difference
				int d;
				for (d=-4; d<=4; d++)
				{
					const mass_t mz_diff = fabs(table_mz + d*one_over_charge - spectrumMOverZ);
					if (mz_diff<maxMzDiffAllowed)
						break;
				}

				if (d>4)
				{
					rankingScores[i] = NEG_INF;
					continue;
				}
			//	cout << "ok : " << table_mz << endl;
			}

			rankingScores[i]=pmc_rank_models[charge][sizeIndex]->calc_rank_score(spectrumRankSamples[i]);
			currentSpectrumPmcTables_[charge][i].rank_score = rankingScores[i];
			if (rankingScores[i]>best_score)
			{
				best_score=rankingScores[i];
				best_idx = i;
			}
		}
	}

	// no suitable models were found for this spectrum
	if (best_idx<0)
	{
		result.mz1 = pl.getHeader()->getMOverZ();
		result.score1 = 10.0;
		result.mz2 = NEG_INF;
		result.score2 = NEG_INF;
		return;
	}

	
	result.mz1 = currentSpectrumPmcTables_[charge][best_idx].m_over_z;
	result.score1 = best_score;

	// look for additional m/z
	int second_best_idx=-1;
	float second_best_score=NEG_INF;

	const mass_t mz_diff = currentSpectrumPmcTables_[charge][1].m_over_z - 
						   currentSpectrumPmcTables_[charge][0].m_over_z;

	const int idx_diff = (int)(0.45/(charge * mz_diff));

	int i;
	for (i=0; i<spectrumRankSamples.size(); i++)
	{
		if (rankingScores[i]>NEG_INF && abs(i-best_idx)<idx_diff)
			continue;

		if (rankingScores[i]>second_best_score)
		{
			second_best_score=rankingScores[i];
			second_best_idx = i;
		}

	}
 
	if (second_best_idx>=0)
	{
		result.mz2 = currentSpectrumPmcTables_[charge][second_best_idx].m_over_z;
		result.score2 = second_best_score;
	} 
	else
	{
		result.mz2 = NEG_INF;
		result.score2 = NEG_INF;
	}
}

/**********************************************************************************
Finds the best m/z values for a given charge
If the pm_tolerance is low (less than 0.5) then the m/z vlaues need to reflect +-X Da
from the recorded m/z
***********************************************************************************/
void PMCSQS_Scorer::find_best_mz_values_from_rank_model(
										const BasicSpectrum& bs, 
										int charge,
										mass_t pm_tolerance,
										PmcSqsChargeRes& res)
{
	static vector<RankBoostSample> spectrumRankSamples;
	static vector<float> rankingScores;

	const mass_t spectrumMOverZ = bs.ssf->m_over_z;
	const mass_t maxMzDiffAllowed = pm_tolerance / charge;

	fill_RankBoost_smaples_with_PMC(bs, charge, spectrumRankSamples);

	if (rankingScores.size() != spectrumRankSamples.size())
		rankingScores.resize(spectrumRankSamples.size(),NEG_INF);

	const mass_t pm_with_19 = bs.ssf->m_over_z * charge - (charge + 1);
	const int sizeIndex = get_rank_model_size_idx(charge, pm_with_19);
	
	int best_idx=-1;
	float best_score=NEG_INF;

	if (charge>= pmc_rank_models.size() ||
		sizeIndex>= pmc_rank_models[charge].size() ||
		! pmc_rank_models[charge][sizeIndex])
	{
		//
	}
	else
	{
	//	cout << "spec: " << spectrumMOverZ << "  charge: " << charge << endl;
		int i;
		for (i=0; i<spectrumRankSamples.size(); i++)
		{
			// make sure the m/z is in an allowed range (i.e., it has a mass that can possibly be a +-X Da
			// shift from the true pm). Assume shift can be at most
			if (pm_tolerance<0.5)
			{
				const mass_t table_mz = currentSpectrumPmcTables_[charge][i].m_over_z;
				const mass_t one_over_charge = 1.0033 / (mass_t)charge; // ~mass of isotopic peak difference
				int d;
				for (d=-4; d<=4; d++)
				{
					const mass_t mz_diff = fabs(table_mz + d*one_over_charge - spectrumMOverZ);
					if (mz_diff<maxMzDiffAllowed)
						break;
				}

				if (d>4)
				{
					rankingScores[i] = NEG_INF;
					continue;
				}
			//	cout << "ok : " << table_mz << endl;
			}

			rankingScores[i]=pmc_rank_models[charge][sizeIndex]->calc_rank_score(spectrumRankSamples[i]);
			currentSpectrumPmcTables_[charge][i].rank_score = rankingScores[i];
			if (rankingScores[i]>best_score)
			{
				best_score=rankingScores[i];
				best_idx = i;
			}
		}
	}

	// no suitable models were found for this spectrum
	if (best_idx<0)
	{
		res.mz1 = bs.ssf->m_over_z;
		res.score1 = 10.0;
		res.mz2 = NEG_INF;
		res.score2 = NEG_INF;
		return;
	}

	
	res.mz1 = currentSpectrumPmcTables_[charge][best_idx].m_over_z;
	res.score1 = best_score;

	// look for additional m/z
	int second_best_idx=-1;
	float second_best_score=NEG_INF;

	const mass_t mz_diff = currentSpectrumPmcTables_[charge][1].m_over_z - 
						   currentSpectrumPmcTables_[charge][0].m_over_z;

	const int idx_diff = (int)(0.45/(charge * mz_diff));

	int i;
	for (i=0; i<spectrumRankSamples.size(); i++)
	{
		if (rankingScores[i]>NEG_INF && abs(i-best_idx)<idx_diff)
			continue;

		if (rankingScores[i]>second_best_score)
		{
			second_best_score=rankingScores[i];
			second_best_idx = i;
		}

	}
 
	if (second_best_idx>=0)
	{
		res.mz2 = currentSpectrumPmcTables_[charge][second_best_idx].m_over_z;
		res.score2 = second_best_score;
	} 
	else
	{
		res.mz2 = NEG_INF;
		res.score2 = NEG_INF;
	}

//	const int sizeIndex = this->get_rank_model_sizeIndex(charge,res.mz1*charge-charge+1);
//	res.mz1 -= pmcMzBiases_[charge][sizeIndex];
//	res.mz2 -= pmcMzBiases_[charge][sizeIndex];

//	cout << charge << " ]\t" << res.mz1 << "\t" << res.prob1 << "\t" << res.mz2 << "\t" << res.prob2 << endl;


}




