#include "AllScoreModels.h"

AllScoreModels::AllScoreModels()
{
	ind_was_initialized = false;
	
	setTrainingStageNames();

	int i;
	for (i=0; i<3; i++)
		peptideRankModels_[i]=NULL;

	for (i=0; i<10; i++)
	{
		tagRankModels_[i]=NULL;
		cumulativeSequenceModels_[i]=NULL;
	}

	const int numPeakModels = sizeof(this->peak_prediction_models)/sizeof(PeakRankModel*);
	for (i=0; i<numPeakModels; i++)
		peak_prediction_models[i]=NULL;
}

void TrainingStage::writeNameHeader() const
{
	ostringstream oss;
	oss << "STAGE " << index << " : " << name;
	cout << endl << oss.str() << endl;
	int i;
	for (i=0; i<oss.str().length(); i++)
		cout << "=";
	cout << endl << endl;

}

void AllScoreModels::setTrainingStageNames()
{
	trainingStages_.clear();
	trainingStages_.push_back(TrainingStage("Partitioning according to size and charge",0));
	trainingStages_.push_back(TrainingStage("Choosing set of fragment ion types",1));
	trainingStages_.push_back(TrainingStage("Precusror ion and fragment ion mass tolerances",2));
	trainingStages_.push_back(TrainingStage("Sequence Quality Score models (SQS)",3));
	trainingStages_.push_back(TrainingStage("Precursor Mass Correction models (PMCR)",4));
	trainingStages_.push_back(TrainingStage("Breakage score models (PRM node scores)",5));
	trainingStages_.push_back(TrainingStage("PRM Score normalizers",6));
	trainingStages_.push_back(TrainingStage("Edge score models",7));
	trainingStages_.push_back(TrainingStage("Peak rank predictions (complete sequences)",8));
	trainingStages_.push_back(TrainingStage("Peak rank predictions (de novo sequences)",9));
	trainingStages_.push_back(TrainingStage("Rerank models (database score reranker)",10));
	trainingStages_.push_back(TrainingStage("Rerank models (denovo results reranker)",11));
	trainingStages_.push_back(TrainingStage("Rerank models (complete denovo sequence reranker)",12));
	trainingStages_.push_back(TrainingStage("Rerank models (peptide sequence tag reranker)",13));	
	trainingStages_.push_back(TrainingStage("Amino acid probability models for de novo sequences and tags",14));
	trainingStages_.push_back(TrainingStage("Cumulative probability models for de novo sequences and tags",15));
}

// checks that all model components were initialized
bool AllScoreModels::checkModelsInitializationStatus(bool report)
{
	int i;
	for (i=0; i<trainingStages_.size(); i++)
		trainingStages_[i].indWasInitialized = false;

	if (config_.getIndWasInitialized())
	{
		if (config_.get_size_thresholds().size()>0)
			trainingStages_[0].indWasInitialized = true; // size and charge thresholds

		if (config_.get_all_fragments().size()>0)
			trainingStages_[1].indWasInitialized = true; // fragments

		if (config_.getTolerance()>0 && config_.get_pm_tolerance()>0) // tolerances
			trainingStages_[2].indWasInitialized = true;
	}

	if (pmcsqs.getIndInitializedSqs())
		trainingStages_[3].indWasInitialized = true;

	if (pmcsqs.getIndInitializedPmcr())
		trainingStages_[4].indWasInitialized = true;

	if (this->prmNodeScoreModel_.get_ind_was_initialized())
		trainingStages_[5].indWasInitialized = true;

	if (prmNodeScoreModel_.getIndNormalizersInitialized())
		trainingStages_[6].indWasInitialized = true;

	if (this->edgeModel_.get_ind_was_initialized())
		trainingStages_[7].indWasInitialized = true;

	if (this->get_peak_prediction_model_ptr(3))
		trainingStages_[8].indWasInitialized = true;

	if (this->get_peak_prediction_model_ptr(4))
		trainingStages_[9].indWasInitialized = true;

	if (this->get_rank_model_ptr(0))
		trainingStages_[10].indWasInitialized = true;

	if (this->get_rank_model_ptr(1))
		trainingStages_[11].indWasInitialized = true;

	if (this->get_rank_model_ptr(2))
		trainingStages_[12].indWasInitialized = true;

	if (this->get_rank_tag_model_ptr(3))
		trainingStages_[13].indWasInitialized = true;

	if (this->amino_acid_probs.get_ind_initialized())
		trainingStages_[14].indWasInitialized = true;

	if (cumulativeSequenceModels_ && cumulativeSequenceModels_[0] &&
		cumulativeSequenceModels_[0]->get_ind_initialized())
		trainingStages_[15].indWasInitialized = true;


	if (report)
	{
		cout << "Model initialzation report:" << endl;
		cout << "---------------------------" << endl;
		int i;
		for (i=0; i<trainingStages_.size(); i++)
			cout << i << "\t" << (trainingStages_[i].indWasInitialized ? "  +" : "  -") << "\t" <<
				trainingStages_[i].name << endl;
		cout << endl;
	}

	for (i=0; i<trainingStages_.size(); i++)
		if (! trainingStages_[i].indWasInitialized)
			break;

	return (i<trainingStages_.size());
}


/*************************************************************************
This function performs the entire training process of the model
Allows for training in stages, gives better output and checks that
previous stages are intialized
**************************************************************************/
void AllScoreModels::trainModelsInStages(
					const char *newModelName, 
					const FileManager& fm,
					const SpectraAggregator& sa,
					mass_t initialToleranceEstimate, 
					int startTrainingStage,
					int endTrainingStage,
					int specificCharge, 
					int specificSize, 
					int specificRegion,
					const char *pathNegativeSpectra)
{
	setTrainingStageNames();


	if (endTrainingStage >= trainingStages_.size())
		endTrainingStage = trainingStages_.size() -1;
	
	// check what stages were initialized already
	checkModelsInitializationStatus(true);
	
	model_name = newModelName;
	config_.set_model_name(string(newModelName));

	// check if starting stage dependencies
	if (startTrainingStage>0)
	{
		int i;
		for (i=0; i<startTrainingStage && i<7; i++)
			if (! trainingStages_[i].indWasInitialized)
				break;

		if (i<startTrainingStage && i<7)
		{
			cout << "Stage " << i << " (" << trainingStages_[i].name << ") was not initialized!" << endl;
			cout << "Starting training at this stage!" << endl;
			startTrainingStage=i;
		}
	}


	// start training according to stages

	// partition according to sizes
	if (startTrainingStage<=0 && endTrainingStage>=0)
	{
		trainingStages_[0].writeNameHeader();
		config_.computeSizeThresholds(sa);
	}

	// selection of fragments
	if (startTrainingStage<=1 && endTrainingStage>=1)
	{
		trainingStages_[1].writeNameHeader();
		config_.setTolerances(initialToleranceEstimate);
		config_.selectFragmentIonTypes(sa, 16, 0.05);
	}

	// fragment tolerance and precursor tolerance
	if (startTrainingStage<=2 && endTrainingStage>=2)
	{
		trainingStages_[2].writeNameHeader();
		
		config_.learnTolerancesFromData(sa, initialToleranceEstimate);
		write_model();
	}

	// SQS - spectra quality score
	if (startTrainingStage<=3 && endTrainingStage>=3)
	{
		trainingStages_[3].writeNameHeader();

		if (! pathNegativeSpectra)
		{
			cout << "Error: to train SQS models you must supply a file with negative spectra samples (with the -neg_spec_list flag)." << endl;
			write_model();
			exit(1);
		}

		pmcsqs.trainSqsModels(&config_, sa, pathNegativeSpectra, specificCharge);
		write_model();
	}


	// PMCR - precursor mass correction
	if (startTrainingStage<=4 && endTrainingStage>=4)
	{
		trainingStages_[4].writeNameHeader();

		pmcsqs.trainPmcRankModels(&config_, sa, specificCharge);
		write_model();
	}


	// Node scores
	if (startTrainingStage<=5 && endTrainingStage>=5)
	{
		trainingStages_[5].writeNameHeader();
		prmNodeScoreModel_.trainNodeScoreModels(static_cast<void*>(this), newModelName, sa, 
												specificCharge, specificSize, specificRegion);	
		write_model();
	}

	// Node score normalizer
	if (startTrainingStage<=6 && endTrainingStage>=6)
	{
		trainingStages_[6].writeNameHeader();
		prmNodeScoreModel_.learnPrmNormalizerValue(static_cast<void*>(this), sa);
		prmNodeScoreModel_.write_prm_normalizer_values();
	}

	// Edge scores
	if (startTrainingStage<=7 && endTrainingStage>=7)
	{
		trainingStages_[7].writeNameHeader();
		edgeModel_.train_all_edge_models(sa,static_cast<void*>(this),specificCharge);
		write_model();
	}

	// rerank model (database scores)
	if (startTrainingStage<=8 && endTrainingStage>=7)
	{
		trainingStages_[8].writeNameHeader();
	}



	
	exit(0);

	cout << endl << "STAGE 7: Train Amino Acid models" << endl;
	cout <<         "********************************" << endl << endl;
	if (startTrainingStage>7)
	{
		cout << endl << "Already done." << endl;
	}
	else
	{
		if (specificCharge>0)
			cout << "+++ Only specified charge " <<  specificCharge << endl << endl;

		amino_acid_probs.train_amino_acid_prob_models(fm,this,specificCharge, specificSize);
	}


/*	cout << endl << "STAGE 8: Train Cumulative de novo probability models" << endl;
	cout <<         "****************************************************" << endl << endl;
	if (startTrainingStage>8)
	{
		cout << endl << "Already done." << endl;
	}
	else
	{
		if (specific_charge>0)
			cout << "+++ Only specified charge " <<  specific_charge << endl << endl;

		.train_seq_prob_models(fm,this,specific_charge,specific_size);


	} */


	if (endTrainingStage<=8)
	{
		write_model();
		exit(0);
	}


	exit(0);
}


struct PeakTuple {
	bool operator< (const PeakTuple& rhs) const
	{
		return (score>rhs.score);
	}
	size_t frag_idx;
	size_t pos;
	int	   rank;
	float  score;
};

void AllScoreModels::predict_fragmentation( const char* input_file, size_t num_peaks)
{
	FILE* stream = fopen(input_file,"r");
	if (! stream)
	{
		cout << "Error: couldn't open file for reading: " << input_file << endl;
		exit(1);
	}

	PeptideRankScorer *dnv_rank = (PeptideRankScorer *)get_rank_model_ptr(1);
	PeakRankModel *prm = this->get_peak_prediction_model_ptr(3);

	char buffer[128];
	char pep_str[128];
	while (fgets(buffer,128,stream))
	{
		int charge;

		if (sscanf(buffer,"%s %d",pep_str,&charge) != 2)
			continue;

		cout << ">> " << pep_str << "\t" << charge << endl;
		if (charge<1 || charge>=prm->get_size_thresholds().size())
		{
			cout << "Invalid charge!" << endl;
			continue;
		}

		Peptide pep;
		pep.parseFromString(&config_,static_cast<string>(pep_str));


		PeptideSolution sol;
		sol.pep = pep;
		sol.reaches_n_terminal=true;
		sol.reaches_c_terminal=true;
		sol.charge = charge;
		sol.pm_with_19 = pep.get_mass_with_19();

		PeptidePeakPrediction ppp;
		prm->calc_peptide_predicted_scores(sol, ppp);

		const size_t num_frags = ppp.frag_idxs.size();
		vector< vector<int> > predicted_ranks;
		calc_combined_peak_ranks(ppp.rank_scores, predicted_ranks);

		vector<PeakTuple> tuples;
		for (size_t f=0; f<num_frags; f++)
			for (size_t i=0; i<ppp.rank_scores[f].size(); i++)
				if (predicted_ranks[f][i]<999)
				{
					PeakTuple pt;
					pt.frag_idx = f;
					pt.pos =i;
					pt.rank = predicted_ranks[f][i];
					pt.score = ppp.rank_scores[f][i];
					tuples.push_back(pt);
				}
		
		sort(tuples.begin(),tuples.end());

		if (tuples.size()<1)
			continue;

		const size_t num_aas = pep.get_num_aas();
		vector<mass_t> breakage_masses;
		pep.calc_expected_breakage_masses(&config_, breakage_masses);

		cout << fixed << "Rank\tIon\tm/z\tScore" << endl;
		for (size_t i=0; i<num_peaks && i<tuples.size(); i++)
		{
			PeakTuple pt = tuples[i];
			cout << i+1 << "\t";
			const FragmentType& ft = config_.get_fragment(ppp.frag_idxs[pt.frag_idx]);
			cout << ft.label << ":" << (ft.orientation == PREFIX ? pt.pos : num_aas - pt.pos) << "\t";

			mass_t mz =  ft.calc_expected_mass(breakage_masses[pt.pos],pep.get_mass_with_19());
			cout << setprecision(2);
			if (mz<100)
				cout << " ";
			if (mz<1000)
				cout << " ";

			cout << mz << "\t";
			cout << setprecision(3) << pt.score << endl;
		}
		cout << endl;
	}

	fclose(stream);
}

