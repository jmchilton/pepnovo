#include "PMCSQS.h"
#include "SpectraList.h"
#include "auxfun.h"

const char * SQS_var_names[]={
"SQS_CONST",	"SQS_PEAK_DENSITY",	"SQS_PROP_UPTO2G",	"SQS_PROP_UPTO5G",	"SQS_PROP_UPTO10G",	"SQS_PROP_MORE10G",	"SQS_PROP_INTEN_UPTO2G",	"SQS_PROP_INTEN_UPTO5G",	"SQS_PROP_INTEN_MORE5G",	"SQS_PROP_ISO_PEAKS",	"SQS_PROP_STRONG_WITH_ISO_PEAKS",	"SQS_PROP_ALL_WITH_H2O_LOSS",	"SQS_PROP_ALL_WITH_NH3_LOSS",	"SQS_PROP_ALL_WITH_CO_LOSS",	"SQS_PROP_STRONG_WITH_H2O_LOSS",	"SQS_PROP_STRONG_WITH_NH3_LOSS",	
"SQS_PROP_STRONG_WITH_CO_LOSS",	"SQS_C2_PROP_ALL_WITH_H2O_LOSS",	"SQS_C2_PROP_ALL_WITH_NH3_LOSS",	"SQS_C2_PROP_ALL_WITH_CO_LOSS",	"SQS_C2_PROP_STRONG_WITH_H2O_LOSS",	"SQS_C2_PROP_STRONG_WITH_NH3_LOSS",	"SQS_C2_PROP_STRONG_WITH_CO_LOSS",	"SQS_DIFF_ALL_WITH_H2O_LOSS",	"SQS_DIFF_ALL_WITH_NH3_LOSS",	"SQS_DIFF_ALL_WITH_CO_LOSS",	"SQS_DIFF_STRONG_WITH_H2O_LOSS",	"SQS_DIFF_STRONG_WITH_NH3_LOSS",	
"SQS_DIFF_STRONG_WITH_CO_LOSS",	"SQS_PROP_PEAKS_WITH_C1C2",	"SQS_PROP_STRONG_PEAKS_WITH_C1C2",	"SQS_PROP_INTEN_WITH_C1C2",	"SQS_IND_MAX_TAG_LENGTH_ABOVE_4",	"SQS_IND_MAX_TAG_LENGTH_BELOW_4",	"SQS_MAX_TAG_LENGTH_ABOVE_4",	"SQS_MAX_TAG_LENGTH_BELOW_4",	"SQS_PROP_INTEN_IN_TAGS",	"SQS_PROP_TAGS1",	"SQS_PROP_STRONG_PEAKS_IN_TAG1",	"SQS_PROP_INTEN_TAG1",	"SQS_IND_PROP_STRONG_BELOW30_TAG1",	
"SQS_PROP_TAGS2",	"SQS_PROP_STRONG_PEAKS_IN_TAG2",	"SQS_PROP_INTEN_TAG2",	"SQS_IND_PROP_STRONG_BELOW20_TAG2",	"SQS_PROP_TAGS3",	"SQS_PROP_STRONG_PEAKS_IN_TAG3",	"SQS_PROP_INTEN_TAG3",	"SQS_IND_PROP_STRONG_BELOW10_TAG3",	"SQS_C2_IND_MAX_TAG_LENGTH_ABOVE_4",	"SQS_C2IND_MAX_TAG_LENGTH_BELOW_4",	"SQS_C2MAX_TAG_LENGTH_ABOVE_4",	"SQS_C2MAX_TAG_LENGTH_BELOW_4",	"SQS_C2PROP_INTEN_IN_TAGS",	
"SQS_C2PROP_TAGS1",	"SQS_C2PROP_STRONG_PEAKS_IN_TAG1",	"SQS_C2PROP_INTEN_TAG1",	"SQS_IND_C2PROP_STRONG_BELOW30_TAG1",	"SQS_C2PROP_TAGS2",	"SQS_C2PROP_STRONG_PEAKS_IN_TAG2",	"SQS_C2PROP_INTEN_TAG2",	"SQS_IND_C2PROP_STRONG_BELOW20_TAG2",	"SQS_C2PROP_TAGS3",	"SQS_C2PROP_STRONG_PEAKS_IN_TAG3",	"SQS_C2PROP_INTEN_TAG3",	"SQS_IND_C2PROP_STRONG_BELOW10_TAG3",	"SQS_DIFF_MAX_TAG_LENGTH",	
"SQS_DIFF_PROP_INTEN_IN_TAGS",	"SQS_DIFF_PROP_TAGS1",	"SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG1",	"SQS_DIFF_PROP_INTEN_TAG1",	"SQS_DIFF_PROP_TAGS2",	"SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG2",	"SQS_DIFF_PROP_INTEN_TAG2",	"SQS_DIFF_PROP_TAGS3",	"SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG3",	"SQS_DIFF_PROP_INTEN_TAG3",	"SQS_PEAK_DENSE_T1",	"SQS_PEAK_DENSE_T2",	"SQS_PEAK_DENSE_T3",	"SQS_INTEN_DENSE_T1",	
"SQS_INTEN_DENSE_T2",	"SQS_INTEN_DENSE_T3",	"SQS_PEAK_DENSE_H1",	"SQS_PEAK_DENSE_H2",	"SQS_INTEN_DENSE_H1",	"SQS_INTEN_DENSE_H2",	"SQS_PROP_MZ_RANGE_WITH_33_INTEN",	"SQS_PROP_MZ_RANGE_WITH_50_INTEN",	"SQS_PROP_MZ_RANGE_WITH_75_INTEN",	"SQS_PROP_MZ_RANGE_WITH_90_INTEN",	"SQS_NUM_FRAG_PAIRS_1",	"SQS_NUM_STRONG_FRAG_PAIRS_1",	"SQS_NUM_C2_FRAG_PAIRS_1",	"SQS_NUM_STRONG_C2_FRAG_PAIRS_1",	
"SQS_NUM_FRAG_PAIRS_2",	"SQS_NUM_STRONG_FRAG_PAIRS_2",	"SQS_NUM_C2_FRAG_PAIRS_2",	"SQS_NUM_STRONG_C2_FRAG_PAIRS_2",	"SQS_NUM_FRAG_PAIRS_3",	"SQS_NUM_STRONG_FRAG_PAIRS_3",	"SQS_NUM_C2_FRAG_PAIRS_3",	"SQS_NUM_STRONG_C2_FRAG_PAIRS_3",	"SQS_PROP_OF_MAX_FRAG_PAIRS_1",	"SQS_PROP_OF_MAX_STRONG_FRAG_PAIRS_1",	"SQS_PROP_OF_MAX_C2_FRAG_PAIRS_1",	"SQS_PROP_OF_MAX_STRONG_C2_FRAG_PAIRS_1",	
"SQS_PROP_OF_MAX_FRAG_PAIRS_2",	"SQS_PROP_OF_MAX_STRONG_FRAG_PAIRS_2",	"SQS_PROP_OF_MAX_C2_FRAG_PAIRS_2",	"SQS_PROP_OF_MAX_STRONG_C2_FRAG_PAIRS_2",	"SQS_PROP_OF_MAX_FRAG_PAIRS_3",	"SQS_PROP_OF_MAX_STRONG_FRAG_PAIRS_3",	"SQS_PROP_OF_MAX_C2_FRAG_PAIRS_3",	"SQS_PROP_OF_MAX_STRONG_C2_FRAG_PAIRS_3",	"SQS_PROP_FRAG_PAIRS_1",	"SQS_PROP_STRONG_FRAG_PAIRS_1",	"SQS_PROP_C2_FRAG_PAIRS_1",	
"SQS_PROP_STRONG_C2_FRAG_PAIRS_1",	"SQS_PROP_FRAG_PAIRS_2",	"SQS_PROP_STRONG_FRAG_PAIRS_2",	"SQS_PROP_C2_FRAG_PAIRS_2",	"SQS_PROP_STRONG_C2_FRAG_PAIRS_2",	"SQS_PROP_FRAG_PAIRS_3",	"SQS_PROP_STRONG_FRAG_PAIRS_3",	"SQS_PROP_C2_FRAG_PAIRS_3",	"SQS_PROP_STRONG_C2_FRAG_PAIRS_3",	"SQS_DIFF_NUM_FRAG_PAIRS_23",	"SQS_DIFF_NUM_STRONG_FRAG_PAIRS_23",	"SQS_DIFF_NUM_C2_FRAG_PAIRS_23",	"SQS_DIFF_NUM_STRONG_C2_FRAG_PAIRS_23",	
"SQS_DIFF_PROP_OF_MAX_FRAG_PAIRS_23",	"SQS_DIFF_PROP_OF_MAX_STRONG_FRAG_PAIRS_23",	"SQS_DIFF_PROP_OF_MAX_C2_FRAG_PAIRS_23",	"SQS_DIFF_PROP_OF_MAX_STRONG_C2_FRAG_PAIRS_23",	"SQS_DIFF_PROP_FRAG_PAIRS_23",	"SQS_DIFF_PROP_STRONG_FRAG_PAIRS_23",	"SQS_DIFF_PROP_C2_FRAG_PAIRS_23",	"SQS_DIFF_PROP_STRONG_C2_FRAG_PAIRS_23",	"SQS_NUM_FIELDS",	"SQS_Fields"};



void PMCSQS_Scorer::setClassWeightsAccordingToData(const SpectraAggregator& sa, 
												   vector< vector<float> >& classWeights) const
{
	vector< vector<int> > counters;
	const int numSizes = getNumSizes();

	classWeights.resize(maximalChargeWithModels_ + 1);
	counters.resize(maximalChargeWithModels_ + 1);

	int c;
	for (c=1; c<=maximalChargeWithModels_; c++)
	{
		classWeights[c].resize(numSizes,1.0);
		counters[c].resize(numSizes,0);

		SpectraList sl(sa);
		sl.selectHeaders(0, POS_INF, c, c);

		int i;
		for (i=0; i<sl.getNumHeaders(); i++)
		{
			const mass_t mz = sl.getSpectrumHeader(i)->getMOverZ();
			const int sizeIndex = this->getSqsSizeIndex(mz);
			counters[c][sizeIndex]++;
		}
	}

	// assign weights so the probability for each size (across the different charges equals 1)
	int i;
	for (i=0; i<numSizes; i++)
	{
		int total=0;
		int c;
		for (c=1; c<=maximalChargeWithModels_; c++) 
			total+=counters[c][i];

		if (total<1000)
			continue;

		for (c=1; c<=maximalChargeWithModels_; c++)
			classWeights[c][i]=counters[c][i]/static_cast<float>(total);

		// make sure that minimal class weight is 0.2/ max_charge
		const float minWeight = 0.2 / static_cast<float>(maximalChargeWithModels_);
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (classWeights[c][i] > 0.0005 && classWeights[c][i] < minWeight)
				classWeights[c][i] = minWeight;
			
	}



	cout << endl << "SQS class weights as determined from the data:" << endl;
	cout << "Charge";
	for (i=0; i<numSizes; i++)
		cout << "\tSize " << i;
	cout << endl;
	for (c=1; c<=maximalChargeWithModels_; c++)
	{
		cout << c;
		for (i=0; i<numSizes; i++)
			cout << setprecision(3) << fixed << "\t" << classWeights[c][i];
		cout << endl;
	}
	cout << endl;
}



void PMCSQS_Scorer::trainSqsModels(const Config* config, 
								   const SpectraAggregator& positiveSpectra,
								   const char* pathNegativeSpectraList,
								   int specificCharge,
								   vector< vector<float> >* inputWeights)
{
	// TODO add weight file that can be read from outside to set the weights... ?

	vector< vector< vector<ME_Regression_Sample> > > samples; //  first dim: neg, +1, +2, +3
															  // second dim: sizeIndex


	
	maximalChargeWithModels_ = (inputWeights ? inputWeights->size()-1 : 3);

	set_frag_pair_sum_offset(MASS_PROTON); // b+y - PM+19
	set_bin_increment(0.1);

	set_sqs_mass_thresholds();
	if (pmcMassThresholds_.size() == 0)
	{
		pmcMassThresholds_=config->get_size_thresholds();
	}

	vector<vector<float> > classWeights;
	if (inputWeights)
	{
		classWeights = *inputWeights;
	}
	else
		setClassWeightsAccordingToData(positiveSpectra, classWeights);


	const int numSizes = sqsMassThresholds_.size();
	cout << "number of sizes for SQS models " << numSizes+1 << endl;

	samples.resize(maximalChargeWithModels_+1);
	
	SpectraAggregator negativeSpectra;
	negativeSpectra.initializeFromTextFile(pathNegativeSpectraList, config);
	const int maxHeadersPerModel = 8000;
	
	// read all samples
	size_t charge;
	for (charge=0; charge<=maximalChargeWithModels_; charge++)
	{
		if (charge>0 && specificCharge>0 && charge != specificCharge)
			continue; 

		samples[charge].resize(numSizes+1);

		size_t sizeIndex;
		for (sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
		{	
			const mass_t minMass = (sizeIndex == 0 ? 0 : sqsMassThresholds_[sizeIndex-1]);
			const mass_t maxMass = (sizeIndex == numSizes ? POS_INF : sqsMassThresholds_[sizeIndex]);

			const SpectraAggregator& sa = (charge == 0 ? negativeSpectra : positiveSpectra);
			SpectraList sl(sa);

			if (charge == 0)
			{
				sl.selectHeaders(minMass, maxMass);
			}
			else
				sl.selectHeaders(minMass, maxMass, charge, charge);

			cout << "Found " << sl.getNumHeaders() << " for charge " << charge << " ranges:" <<
				minMass << " - " << maxMass << endl;

			sl.randomlyReduceListToSize(maxHeadersPerModel);

			
			const int label = (charge == 0 ? 1 : 0);	
			samples[charge][sizeIndex].resize(sl.getNumHeaders());
			int i;
			for (i=0; i<sl.getNumHeaders(); i++)
			{
				const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
				PeakList pl;

				pl.readPeaksToLocalAllocation(sa,header);
				pl.initializePeakList(config, true);
			
				initializeForCurrentSpectrum(config, pl);

				calculateCurrentSpectrumPmcValues(pl, bin_increment);
			
				fillSqsMeSample(pl, samples[charge][sizeIndex][i]);
				samples[charge][sizeIndex][i].label = label;
			}
		}
	}

	// cout sample composition
	cout << "Sample composition:" << endl;
	for (charge=0; charge<=maximalChargeWithModels_; charge++)
	{
		cout << charge;
		size_t i;
		for (i=0; i<samples[charge].size(); i++)
			cout << "\t" << samples[charge][i].size();
		cout << endl;
	}

	// create SQS models
	sqs_models.resize(maximalChargeWithModels_+1);
	for (charge =0; charge<=maximalChargeWithModels_; charge++)
	{
		sqs_models[charge].resize(maximalChargeWithModels_+1);
		int j;
		for (j=0; j<sqs_models[charge].size(); j++)
			sqs_models[charge][j].resize(numSizes+1,NULL);
	}



	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		int sizeIndex;
		for (sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
		{
			cout << endl << "CHARGE " << charge << " SIZE " << sizeIndex << endl;

			
			ME_Regression_DataSet ds;
			ds.num_classes=2;
			ds.num_features=SQS_NUM_FIELDS;
			ds.add_samples(samples[0][sizeIndex]);
			ds.add_samples(samples[charge][sizeIndex]);
			ds.tally_samples();

			if (ds.class_weights[0]<0.0001 || ds.class_weights[1]<0.0001)
			{
				cout << "Warning: insufficient number of samples, not trianing model for this charge " << charge <<
					" size " << sizeIndex << endl;
				continue;
			}

			const double pos_weight = 0.2 + classWeights[charge][sizeIndex]*0.3;

			ds.randomly_remove_samples_with_activated_feature(1,SQS_IND_MAX_TAG_LENGTH_ABOVE_4,0.5);

			ds.calibrate_class_weights(pos_weight); // charge vs bad spectra
			ds.print_feature_summary(cout,SQS_var_names);

			sqs_models[charge][0][sizeIndex]=new ME_Regression_Model;
			sqs_models[charge][0][sizeIndex]->train_cg(ds,250);
			sqs_models[charge][0][sizeIndex]->print_ds_probs(ds);
		
		}
	}

		
	////////////////////////////////////////////
	// train model vs. model if charge1>charge2
	if (1)
	{
		int charge1,charge2;
		for (charge1=2; charge1<=maximalChargeWithModels_; charge1++)
		{
			for (charge2=1; charge2<charge1; charge2++)
			{
				int sizeIndex;
				for (sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
				{
					ME_Regression_DataSet ds;

					ds.num_classes=2;
					ds.num_features=SQS_NUM_FIELDS;

					ds.add_samples(samples[charge1][sizeIndex]);

					int i;
					for (i=0; i<samples[charge2][sizeIndex].size(); i++)
					{
						samples[charge2][sizeIndex][i].label=1;
						ds.add_sample(samples[charge2][sizeIndex][i]);
						samples[charge2][sizeIndex][i].label=0;
					}

					float relative_weight = classWeights[charge1][sizeIndex]/
						(classWeights[charge1][sizeIndex]+classWeights[charge2][sizeIndex]);

					ds.tally_samples();

					if (ds.class_weights[0]<0.0001 || ds.class_weights[1]<0.0001)
					{
						cout << "Warning: insufficient number of samples, not trianing model for charge " << charge1 <<
							" vs charge " << charge2<< " (size " << sizeIndex << ")" << endl;
						continue;
					}

					ds.calibrate_class_weights(relative_weight);

					sqs_models[charge1][charge2][sizeIndex] = new ME_Regression_Model;

					cout << endl << "CHARGE " << charge1 << " vs " << charge2 << "  size " << sizeIndex << endl;
					cout << "Relative weights: " << charge1 << "/(" << charge1 << "+" <<
						charge2 << "): " << relative_weight << endl;
				
					ds.print_feature_summary(cout,SQS_var_names);

					sqs_models[charge1][charge2][sizeIndex]->train_cg(ds,300);
					sqs_models[charge1][charge2][sizeIndex]->print_ds_probs(ds);
				}
			}
		}
	}

	init_sqs_correct_factors(maximalChargeWithModels_, sqsMassThresholds_.size());

	////////////////////////////////////////////
	// final report on datasets
	cout << endl;

	int sizeIndex;
	for (sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
	{
		cout << endl << "SIZE: " << sizeIndex << endl;
		cout << "--------" << endl;
		float p_thresh = 0.05;
		int d;
		for (d=0; d<=maximalChargeWithModels_; d++)
		{
			vector<int> counts;
			vector<int> max_counts;
			counts.resize(maximalChargeWithModels_+1,0);
			max_counts.resize(maximalChargeWithModels_+1,0);

			int i;
			for (i=0; i<samples[d][sizeIndex].size(); i++)
			{
				bool above_thresh=false;
				float max_prob=0;
				int   max_class=0;
				int c;
				for (c=1; c<=maximalChargeWithModels_; c++)
				{
					if (! sqs_models[c][0][sizeIndex])
						continue;

					float prob = sqs_models[c][0][sizeIndex]->p_y_given_x(0,samples[d][sizeIndex][i]);
					if (prob>p_thresh)
					{
						counts[c]++;
						above_thresh=true;
						if (prob>max_prob)
						{
							max_prob=prob;
							max_class=c;
						}
					}
				}
				max_counts[max_class]++;

				if (! above_thresh)
					counts[0]++;
			}

			cout << d << "\t";
			for (i=0; i<=maximalChargeWithModels_; i++)
				cout << fixed << setprecision(4) << max_counts[i]/(float)samples[d][sizeIndex].size() << "\t";
			cout << endl;
		}
	}



	ind_initialized_sqs = true;

	string path;
	path = config->get_resource_dir() + "/" + config->get_model_name() + "_SQS.txt";
	write_sqs_models(path.c_str());
}


/**********************************************************************************
***********************************************************************************/
void PMCSQS_Scorer::train_sqs_models(Config *config, 
									 const FileManager& fm_pos, 
									 const char *neg_list,
									 int specificCharge, 
									 vector<vector<float> > *inputWeights)
{
	vector< vector< vector<ME_Regression_Sample> > > samples; //  neg, p1, p2, p3 / sizeIndex
	FileManager fm_neg;

	const vector<int>& spectra_counts = fm_pos.get_spectra_counts();
	maximalChargeWithModels_ = (inputWeights ? inputWeights->size()-1 : 3);
	int charge;

	set_frag_pair_sum_offset(MASS_PROTON); // b+y - PM+19
	set_bin_increment(0.1);
	this->set_sqs_mass_thresholds();

	if (this->pmcMassThresholds_.size() == 0)
	{
		pmcMassThresholds_=config->get_size_thresholds();
	}

	vector<vector<float> > classWeights;
	if (inputWeights)
	{
		classWeights = *inputWeights;
	}
	else
	{
		classWeights.resize(maximalChargeWithModels_+1);
		int i;
		for (i=0; i<classWeights.size(); i++)
			classWeights[i].resize(maximalChargeWithModels_+1,1.0);
	}

	const int numSizes = this->sqsMassThresholds_.size();
	cout << "NUM SIZE MODELS: " << numSizes+1 << endl;

	samples.resize(maximalChargeWithModels_+1);

	fm_neg.init_from_list_file(config, neg_list);
	const int max_to_read_per_file = 8000;

	for (charge=0; charge<=maximalChargeWithModels_; charge++)
	{
		if (charge>0 && specificCharge>0 && charge != specificCharge)
			continue; 

		int sizeIndex;
		for (sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
		{	
			const mass_t minMass = (sizeIndex == 0 ? 0 : sqsMassThresholds_[sizeIndex-1]);
			const mass_t maxMass = (sizeIndex == numSizes ? POS_INF : sqsMassThresholds_[sizeIndex]);

			samples[charge].resize(numSizes+1);

			BasicSpecReader bsr;
			QCPeak peaks[5000]; 

			FileSet fs;
			if (charge == 0)
			{
				fs.select_files_in_mz_range(fm_neg,minMass, maxMass,0);	
			}
			else
			{
				fs.select_files_in_mz_range(fm_pos, minMass, maxMass, charge);
			}

			cout << "Found " << fs.get_total_spectra() << " for charge " << charge << " ranges:" <<
				minMass << " - " << maxMass << endl;

			fs.randomly_reduce_ssfs(max_to_read_per_file);
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
			const int label = (charge == 0 ? 1 : 0);
			const int num_samples =  all_ssf.size();
						
			samples[charge][sizeIndex].resize(num_samples);

			
			int i;
			for (i=0; i<num_samples; i++)
			{
				SingleSpectrumFile* ssf = all_ssf[i];
				BasicSpectrum bs;

				bs.peaks = peaks;
				bs.ssf = ssf;
			
				if (charge==0)
				{
					bs.num_peaks = bsr.read_basic_spec(config,fm_neg,ssf,peaks);
					bs.ssf->charge=0;
				}
				else
					bs.num_peaks = bsr.read_basic_spec(config,fm_pos,ssf,peaks);

				init_for_current_spec(config,bs);
				calculate_curr_spec_pmc_values(bs, bin_increment);
			
				fill_fval_vector_with_SQS(bs, samples[charge][sizeIndex][i]);
				
				samples[charge][sizeIndex][i].label = label;
			}
		}
	}

	// cout sample composition
	cout << "Sample composition:" << endl;
	for (charge=0; charge<=maximalChargeWithModels_; charge++)
	{
		cout << charge;
		int i;
		for (i=0; i<samples[charge].size(); i++)
			cout << "\t" << samples[charge][i].size();
		cout << endl;
	}

	// create SQS models
	this->sqs_models.resize(maximalChargeWithModels_+1);
	for (charge =0; charge<=maximalChargeWithModels_; charge++)
	{
		sqs_models[charge].resize(maximalChargeWithModels_+1);
		int j;
		for (j=0; j<sqs_models[charge].size(); j++)
			sqs_models[charge][j].resize(numSizes+1,NULL);
	}



	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		int sizeIndex;
		for (sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
		{
			ME_Regression_DataSet ds;

			cout << endl << "CHARGE " << charge << " SIZE " << sizeIndex << endl;
			ds.num_classes=2;
			ds.num_features=SQS_NUM_FIELDS;
			ds.add_samples(samples[0][sizeIndex]);
			ds.add_samples(samples[charge][sizeIndex]);
			ds.tally_samples();

			if (ds.class_weights[0]<0.0001 || ds.class_weights[1]<0.0001)
			{
				cout << "Warning: insufficient number of samples, not trianing model for this charge " << charge <<
					" size " << sizeIndex << endl;
				continue;
			}

			const double pos_weight = 0.2 + classWeights[charge][sizeIndex]*0.3;

			ds.randomly_remove_samples_with_activated_feature(1,SQS_IND_MAX_TAG_LENGTH_ABOVE_4,0.5);

			ds.calibrate_class_weights(pos_weight); // charge vs bad spectra
			ds.print_feature_summary(cout,SQS_var_names);

			sqs_models[charge][0][sizeIndex]=new ME_Regression_Model;

			sqs_models[charge][0][sizeIndex]->train_cg(ds,250);

			sqs_models[charge][0][sizeIndex]->print_ds_probs(ds);

		}
	}

		
	////////////////////////////////////////////
	// train model vs. model if charge1>charge2
	if (1)
	{
		int charge1,charge2;
		for (charge1=2; charge1<=maximalChargeWithModels_; charge1++)
		{
			for (charge2=1; charge2<charge1; charge2++)
			{
				int sizeIndex;
				for (sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
				{
					ME_Regression_DataSet ds;

					ds.num_classes=2;
					ds.num_features=SQS_NUM_FIELDS;

					ds.add_samples(samples[charge1][sizeIndex]);

					int i;
					for (i=0; i<samples[charge2][sizeIndex].size(); i++)
					{
						samples[charge2][sizeIndex][i].label=1;
						ds.add_sample(samples[charge2][sizeIndex][i]);
						samples[charge2][sizeIndex][i].label=0;
					}

					float relative_weight = classWeights[charge1][sizeIndex]/
						(classWeights[charge1][sizeIndex]+classWeights[charge2][sizeIndex]);

					ds.tally_samples();

					if (ds.class_weights[0]<0.0001 || ds.class_weights[1]<0.0001)
					{
						cout << "Warning: insufficient number of samples, not trianing model for charge " << charge1 <<
							" vs charge " << charge2<< " (size " << sizeIndex << ")" << endl;
						continue;
					}

					ds.calibrate_class_weights(relative_weight);

					sqs_models[charge1][charge2][sizeIndex] = new ME_Regression_Model;

					cout << endl << "CHARGE " << charge1 << " vs " << charge2 << "  size " << sizeIndex << endl;
					cout << "Relative weights: " << charge1 << "/(" << charge1 << "+" <<
						charge2 << "): " << relative_weight << endl;
				
					ds.print_feature_summary(cout,SQS_var_names);

					sqs_models[charge1][charge2][sizeIndex]->train_cg(ds,300);
					sqs_models[charge1][charge2][sizeIndex]->print_ds_probs(ds);
				}
			}
		}
	}

	init_sqs_correct_factors(maximalChargeWithModels_, sqsMassThresholds_.size());

	////////////////////////////////////////////
	// final report on datasets
	cout << endl;

	int sizeIndex;
	for (sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
	{
		cout << endl << "SIZE: " << sizeIndex << endl;
		cout << "--------" << endl;
		float p_thresh = 0.05;
		int d;
		for (d=0; d<=maximalChargeWithModels_; d++)
		{
			vector<int> counts;
			vector<int> max_counts;
			counts.resize(maximalChargeWithModels_+1,0);
			max_counts.resize(maximalChargeWithModels_+1,0);

			int i;
			for (i=0; i<samples[d][sizeIndex].size(); i++)
			{
				bool above_thresh=false;
				float max_prob=0;
				int   max_class=0;
				int c;
				for (c=1; c<=maximalChargeWithModels_; c++)
				{
					if (! sqs_models[c][0][sizeIndex])
						continue;

					float prob = sqs_models[c][0][sizeIndex]->p_y_given_x(0,samples[d][sizeIndex][i]);
					if (prob>p_thresh)
					{
						counts[c]++;
						above_thresh=true;
						if (prob>max_prob)
						{
							max_prob=prob;
							max_class=c;
						}
					}
				}
				max_counts[max_class]++;

				if (! above_thresh)
					counts[0]++;
			}

			cout << d << "\t";
			for (i=0; i<=maximalChargeWithModels_; i++)
				cout << fixed << setprecision(4) << max_counts[i]/(float)samples[d][sizeIndex].size() << "\t";
			cout << endl;
		}
	}



	ind_initialized_sqs = true;

	string path;
	path = config->get_resource_dir() + "/" + config->get_model_name() + "_SQS.txt";
	write_sqs_models(path.c_str());
}


void PMCSQS_Scorer::write_sqs_models(const char *path) const
{
	ofstream out_stream(path,ios::out);
	if (! out_stream.good())
	{
		cout << "Error: couldn't open pmc model for writing: " << path << endl;
		exit(1);
	}
	int i;

	out_stream << sqs_models.size() << endl;
	out_stream << this->sqsMassThresholds_.size() << setprecision(2) << fixed;
	for (i=0; i<sqsMassThresholds_.size(); i++)
		out_stream << " " << sqsMassThresholds_[i];
	out_stream << endl;
		

	const int numSizes = sqsMassThresholds_.size();

	for (i=0; i<sqs_models.size(); i++)
	{
		out_stream << this->sqs_correction_factors[i].size() << setprecision(4);
		int j;
		for (j=0; j<sqs_correction_factors[i].size(); j++)
			out_stream << " " << sqs_correction_factors[i][j] << " " << sqs_mult_factors[i][j];
		out_stream << endl;
	}
	

	
	// write ME models
	
	for (i=0; i<sqs_models.size(); i++)
	{
		int j;
		for (j=0; j<sqs_models[i].size(); j++)
		{
			int k;
			for (k=0; k<sqs_models[i][j].size(); k++)
			{
				if (sqs_models[i][j][k])
				{
					out_stream << i << " " << j << " " << k << endl;
					sqs_models[i][j][k]->write_regression_model(out_stream);
				}
			}
		}
	}	
	out_stream.close();
}


bool PMCSQS_Scorer::read_sqs_models(Config *_config, char *file)
{
	config_ = _config;

	string path;
	path = config_->get_resource_dir() + "/" + string(file);


	ifstream in_stream(path.c_str(),ios::in);
	if (! in_stream.good())
	{
		cout << "Warning: couldn't open sqs model for reading: " << path << endl;
		return false;
	}

	int i;
	char buff[512];
	int num_sqs_charges=-1;

	in_stream.getline(buff,256);
	istringstream iss(buff);

	iss >> num_sqs_charges;
	
	in_stream.getline(buff,256);
	istringstream iss1(buff);
	int numSizes=0;
	iss1 >> numSizes;

	this->sqsMassThresholds_.resize(numSizes,POS_INF);
	for (i=0; i<numSizes; i++)
		iss1 >> sqsMassThresholds_[i];

	this->sqs_correction_factors.resize(num_sqs_charges);
	this->sqs_mult_factors.resize(num_sqs_charges);
	for (i=0; i<num_sqs_charges; i++)
	{
		in_stream.getline(buff,512);
		istringstream iss(buff);
		int num_threshes = 0;
		iss >> num_threshes;

		if (num_threshes>0)
		{
			sqs_correction_factors[i].resize(num_threshes+1,0);
			sqs_mult_factors[i].resize(num_threshes+1,1.0);
			int j;
			for (j=0; j<=num_threshes; j++)
				iss >> sqs_correction_factors[i][j] >> sqs_mult_factors[i][j];
		}
	}

	sqs_models.resize(num_sqs_charges);
	for (i=0; i<num_sqs_charges; i++)
	{
		sqs_models[i].resize(num_sqs_charges);
		int j;
		for (j=0; j<num_sqs_charges; j++)
			sqs_models[i][j].resize(numSizes+1,NULL);
	}

	
	// read ME models
	
	while (in_stream.getline(buff,128))
	{
		int charge1=-1,charge2=-1, sizeIndex=-1;
		sscanf(buff,"%d %d %d",&charge1,&charge2,&sizeIndex); 

		if (charge1<1 || charge2<0 || charge1>=num_sqs_charges || charge2>=charge1 || sizeIndex<0)
		{
			cout << "\nError: reading SQS, bad charge numbers in line: " << endl << buff << endl;
			exit(1);
		}

		sqs_models[charge1][charge2][sizeIndex] = new ME_Regression_Model;
		sqs_models[charge1][charge2][sizeIndex]->read_regression_model(in_stream);
		continue;
	}

	if (maximalChargeWithModels_>0 && maximalChargeWithModels_ >=num_sqs_charges)
	{
		maximalChargeWithModels_ = num_sqs_charges - 1;
	}
	else
		maximalChargeWithModels_ = num_sqs_charges - 1;

	in_stream.close();
	this->ind_initialized_sqs = true;
	return true;
}






/****************************************************************************
Finds the bin which has the optimal values (look for the maximal number of pairs).
Performs search near the peptide's true m/z value to compensate for systematic bias
in the precursor mass.
*****************************************************************************/
int PMCSQS_Scorer::get_optimal_bin(int true_mz_bin, int charge) const
{
	const int max_bin_offset = 6-charge; // look in the range +- of this value
	const vector<PMCRankStats>& pmcStatistics = currentSpectrumPmcTables_[charge];
	const int min_bin_idx = (true_mz_bin - max_bin_offset>=0 ? true_mz_bin - max_bin_offset : 0);
	const int max_bin_idx = (true_mz_bin + max_bin_offset>= pmcStatistics.size() ? pmcStatistics.size()-1 :
								true_mz_bin + max_bin_offset);

	if (pmcStatistics[true_mz_bin].numFragmentPairs==0 &&
		pmcStatistics[true_mz_bin].numCharge2FragmentPairs==0)
		return true_mz_bin;
	
	int   optimal_bin_idx=NEG_INF;
	float max_num_pairs=NEG_INF;
	float best_offset=POS_INF;

	if (pmcStatistics[true_mz_bin].numFragmentPairs>=pmcStatistics[true_mz_bin].numCharge2FragmentPairs)
	{
		float max_num_pairs=0;
		int bin_idx;
		for (bin_idx = min_bin_idx; bin_idx<=max_bin_idx; bin_idx++)
			if (pmcStatistics[bin_idx].numFragmentPairs > max_num_pairs)
				max_num_pairs = pmcStatistics[bin_idx].numFragmentPairs;

		// find minimal offset
		for (bin_idx = min_bin_idx; bin_idx<=max_bin_idx; bin_idx++)
			if (pmcStatistics[bin_idx].numFragmentPairs == max_num_pairs &&
				pmcStatistics[bin_idx].mean_offset_pairs < best_offset)
			{
				optimal_bin_idx = bin_idx;
				best_offset = pmcStatistics[bin_idx].mean_offset_pairs;
			}

		return optimal_bin_idx;
		
	}
	else
	// use the charge 2 fragment pairs
	{
		float max_num_pairs=0; 
		int bin_idx;
		for (bin_idx = min_bin_idx; bin_idx<=max_bin_idx; bin_idx++)
			if (pmcStatistics[bin_idx].numCharge2FragmentPairs > max_num_pairs)
				max_num_pairs = pmcStatistics[bin_idx].numCharge2FragmentPairs;

		// find minimal offset
		for (bin_idx = min_bin_idx; bin_idx<=max_bin_idx; bin_idx++)
			if (pmcStatistics[bin_idx].numCharge2FragmentPairs == max_num_pairs &&
				pmcStatistics[bin_idx].mean_offset_c2_pairs < best_offset)
			{
				optimal_bin_idx = bin_idx;
				best_offset = pmcStatistics[bin_idx].mean_offset_c2_pairs;
			}

		return optimal_bin_idx;	
	}


	return -1;
}



/*********************************************************************************
Takes a set of samples around the correct mass ([-3+5] every 0.1 Da.)
Selects the bin of the correct mass as positive and a set from offseted m/z
as negative samples. 
**********************************************************************************/
void PMCSQS_Scorer::selectTrainingSampleIndexes(
					int charge,
					const vector<RankBoostSample>& samples,
					const PeakList& pl,
					int& correctIndex,
					vector<int>& badPmcIndexes) const
{
	const vector<PMCRankStats>& pmcStatistics = currentSpectrumPmcTables_[charge];

	
	if (pl.getHeader()->getPeptideStr().length()<3)
	{
		cout << "Error: supplied training spectrum without peptide!" << endl;
		exit(1);
	}

	Peptide peptide;
	peptide.parseFromString(config_, pl.getHeader()->getPeptideStr());
	const mass_t peptideMass = peptide.get_mass()+MASS_H2O;
	const int sizeIndex = this->get_rank_model_size_idx(charge, peptideMass);
	const mass_t true_mz = (peptideMass + charge)/(float)charge + this->pmcMzBiases_[charge][sizeIndex];
	const mass_t observed_mz = pl.getHeader()->getMOverZ();

	// check that the training sample has an ok offset
	if (fabs(true_mz-observed_mz)>10.0)
	{
		
		cout << "Erorr in m/z offsets (remove this spectrum from training set): " << endl;
		cout << fixed << setprecision(2) << "file m/z: " << observed_mz << "\t" << 
			"\"true\" m/z: " << true_mz << "\t peptide: " << peptide.as_string(config_) << endl;
		cout << "spectrum: " << pl.getHeader()->getTitle() << endl;
		
		cout << "Mass Cys = " << this->config_->get_aa2mass()[Cys] << endl;

		exit(1);
	}

	// find the entry with the correct m/z

	int idx=0;
	while (idx<pmcStatistics.size() && pmcStatistics[idx].m_over_z<true_mz)
		idx++;

	if (idx>= pmcStatistics.size())
		idx--;

	if (idx>0 && pmcStatistics[idx].m_over_z-true_mz>true_mz-pmcStatistics[idx-1].m_over_z)
		idx--;

//	correctIndex=get_optimal_bin(idx,charge);
	correctIndex = idx;

	vector<int> idxs;
	idxs.clear();
	badPmcIndexes.clear();


	idxs.push_back(correctIndex+4);
	idxs.push_back(correctIndex+5);
	idxs.push_back(correctIndex+7);
	idxs.push_back(correctIndex+9);
	idxs.push_back(correctIndex+10);
	idxs.push_back(correctIndex+15);

	idxs.push_back(correctIndex+19);
	idxs.push_back(correctIndex+20);
	
	idxs.push_back(correctIndex-4);
	idxs.push_back(correctIndex-5);
	idxs.push_back(correctIndex-7);
	idxs.push_back(correctIndex-9);
	idxs.push_back(correctIndex-10);;
	idxs.push_back(correctIndex-15);
	idxs.push_back(correctIndex-19);
	idxs.push_back(correctIndex-20);

	




	// select upto 5 random samples (make sure they are not close to the correct one)
	int i;
	for (i=0; i<5; i++)
	{
		int idx = (int)(myRandom()*pmcStatistics.size());
		if (abs(correctIndex-idx)<6)
			continue;

		idxs.push_back(idx);
	}

	sort(idxs.begin(),idxs.end());
	for (i=0; i<idxs.size(); i++)
		if (idxs[i]>=0 && idxs[i]<pmcStatistics.size())
			badPmcIndexes.push_back(idxs[i]);

}

/*********************************************************************************
Takes a set of samples around the correct mass ([-3+5] every 0.1 Da.)
Selects the bin of the correct mass as positive and a set from offseted m/z
as negative samples. 
**********************************************************************************/
void PMCSQS_Scorer::select_training_sample_idxs(
		int charge,
		const vector<RankBoostSample>& spec_samples,
		const BasicSpectrum& bs,
		int& correctIndex,
		vector<int>& badPmcIndexes) const
{
	const vector<PMCRankStats>& pmcStatistics = currentSpectrumPmcTables_[charge];

	bs.ssf->peptide.calc_mass(config_);
	const mass_t peptideMass = bs.ssf->peptide.get_mass()+MASS_H2O;
	const int sizeIndex = this->get_rank_model_size_idx(charge,peptideMass);
	const mass_t true_mz = (peptideMass + charge)/(float)charge + this->pmcMzBiases_[charge][sizeIndex];
	const mass_t observed_mz = bs.ssf->m_over_z;

	// check that the training sample has an ok offset
	if (fabs(true_mz-observed_mz)>10.0)
	{
		
		cout << "Erorr in m/z offsets (remove this spectrum from training set): " << endl;
		cout << fixed << setprecision(2) << "file m/z: " << observed_mz << "\t" << 
			"\"true\" m/z: " << true_mz << "\t peptide: " << bs.ssf->peptide.as_string(config_) << endl;
		cout << "spectrum: " << bs.ssf->single_name << endl;
		
		cout << "Mass Cys = " << this->config_->get_aa2mass()[Cys] << endl;

		exit(1);
	}

	// find the entry with the correct m/z

	int idx=0;
	while (idx<pmcStatistics.size() && pmcStatistics[idx].m_over_z<true_mz)
		idx++;

	if (idx>= pmcStatistics.size())
		idx--;

	if (idx>0 && pmcStatistics[idx].m_over_z-true_mz>true_mz-pmcStatistics[idx-1].m_over_z)
		idx--;

//	correctIndex=get_optimal_bin(idx,charge);
	correctIndex = idx;

	vector<int> idxs;
	idxs.clear();
	badPmcIndexes.clear();


	idxs.push_back(correctIndex+4);
	idxs.push_back(correctIndex+5);
	idxs.push_back(correctIndex+7);
	idxs.push_back(correctIndex+9);
	idxs.push_back(correctIndex+10);
	idxs.push_back(correctIndex+15);

	idxs.push_back(correctIndex+19);
	idxs.push_back(correctIndex+20);
	
	idxs.push_back(correctIndex-4);
	idxs.push_back(correctIndex-5);
	idxs.push_back(correctIndex-7);
	idxs.push_back(correctIndex-9);
	idxs.push_back(correctIndex-10);;
	idxs.push_back(correctIndex-15);
	idxs.push_back(correctIndex-19);
	idxs.push_back(correctIndex-20);

	




	// select upto 5 random samples (make sure they are not close to the correct one)
	int i;
	for (i=0; i<5; i++)
	{
		int idx = (int)(myRandom()*pmcStatistics.size());
		if (abs(correctIndex-idx)<6)
			continue;

		idxs.push_back(idx);
	}

	sort(idxs.begin(),idxs.end());
	for (i=0; i<idxs.size(); i++)
		if (idxs[i]>=0 && idxs[i]<pmcStatistics.size())
			badPmcIndexes.push_back(idxs[i]);

}



/*************************************************************************
Tests the performance of precursor mass correction
**************************************************************************/
void PMCSQS_Scorer::test_pmc(Config *config, char *specs_file, int charge, 
							 mass_t minMass, mass_t maxMass)
{
	BasicSpecReader bsr;
	static QCPeak peaks[5000];

	FileManager fm;
	FileSet fs;
		
	fm.init_from_file(config,specs_file);
	fs.select_files_in_mz_range(fm,minMass,maxMass,charge);

	const int max_to_read_per_file = 5000;

	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	const int num_samples = (all_ssf.size()<max_to_read_per_file ? all_ssf.size() :
									max_to_read_per_file);
	
	vector<mass_t> org_offsets;
	vector<mass_t> corr_offsets;

	vector<int> ssf_idxs;
	if (num_samples<all_ssf.size())
	{
		chooseKFromN(num_samples,all_ssf.size(),ssf_idxs);
	}
	else
	{
		int i;
		ssf_idxs.resize(all_ssf.size());
		for (i=0; i<all_ssf.size(); i++)
			ssf_idxs[i]=i;
	}

	vector<SingleSpectrumFile *> ssfs;
	int i;
	for (i=0; i<num_samples; i++)
		ssfs.push_back( all_ssf[ssf_idxs[i]]);
	
	output_pmc_rank_results(fm,charge,ssfs);
}


/***********************************************************************************

Functions for training set.


************************************************************************************/

struct ScanPair {
	ScanPair(int f,int sc, string& se) : file_idx(f), scan(sc), seq(se) {};
	ScanPair(int f,int s) : file_idx(f), scan(s) {};
	ScanPair() : file_idx(-1), scan(-1) {};

	bool operator< (const ScanPair& other) const
	{
		return (file_idx<other.file_idx || 
			    (file_idx == other.file_idx && scan<other.scan));
	}

	bool operator == (const ScanPair& other) const
	{
		return (file_idx == other.file_idx && scan == other.scan);
	}


	int file_idx;
	int scan;
	string seq;
};

void read_idxs_from_file(char *file, vector<ScanPair>& final_pairs, int max_size)
{
	ifstream inp(file,ios::in);
	
	if (! inp.good())
	{	
		cout << "Error opening: " << file << endl;
		exit(1);
	}

	vector<ScanPair> pairs;
	pairs.clear();

	char buff[256];
	while (inp.getline(buff,256))
	{
		istringstream iss(buff);
		int f,s;
		string seq;

		iss >> f >> s >> seq;

		
		if (f>=0 && s>=0)
		{
			if (seq.length()>2)
			{
				pairs.push_back(ScanPair(f,s,seq));
			}
			else
				pairs.push_back(ScanPair(f,s));
		}


	}
	inp.close();

	if (pairs.size() > max_size)
	{
		vector<int> idxs;
		chooseKFromN(max_size,pairs.size(),idxs);
		final_pairs.resize(max_size);
		int i;
		for (i=0; i<max_size; i++)
			final_pairs[i]=pairs[idxs[i]];
	}
	else
	{
		final_pairs=pairs;
	}

	sort(final_pairs.begin(),final_pairs.end());
}


void create_training_files(Config *config)
{
	char mzxml_list[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\HEK293_mzxml_list.txt"};
	char idxs_neg_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul_neg_samples.txt"};
//	char idxs1_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul_pos_samples.1.txt"};
//	char idxs2_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul_pos_samples.2.txt"};
//	char idxs2_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\Len10_pos_samples.2.txt"};
	char idxs1_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\sqs_train_pos_samples.1.txt"};
	char idxs2_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\sqs_train_pos_samples.2.txt"};
	char idxs3_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul_pos_samples.3.txt"};

	char out_base[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\sqs_train"};
	string out_neg (out_base); 
	string out1=out_neg;
	string out2=out_neg;
	string out3=out_neg;

	out_neg += "_neg.mgf";
	out1 += "_1.mgf";
	out2 += "_2.mgf";
	out3 += "_3.mgf";

	ofstream stream_neg (out_neg.c_str(),ios::out);
	ofstream stream1(out1.c_str(),ios::out);
	ofstream stream2(out2.c_str(),ios::out);
	ofstream stream3(out3.c_str(),ios::out);

	vector<ScanPair> neg_pairs, pairs1,pairs2,pairs3;


	read_idxs_from_file(idxs_neg_file,neg_pairs,12000);
	read_idxs_from_file(idxs1_file,pairs1,12000);
	read_idxs_from_file(idxs2_file,pairs2,12000);
	read_idxs_from_file(idxs3_file,pairs3,8000);

	cout << "Read " << neg_pairs.size() << " neg idxs\n";
	cout << "Read " << pairs1.size() << " pos1 idxs\n";
	cout << "Read " << pairs2.size() << " pos2 idxs\n";
	cout << "Read " << pairs3.size() << " pos3 idxs\n";

	vector<bool> file_inds;
	file_inds.resize(10000,false);
	int i;

	for (i=0; i<neg_pairs.size(); i++)
		file_inds[neg_pairs[i].file_idx]=true;

	for (i=0; i<pairs1.size(); i++)
		file_inds[pairs1[i].file_idx]=true;

	for (i=0; i<pairs2.size(); i++)
		file_inds[pairs2[i].file_idx]=true;

	for (i=0; i<pairs3.size(); i++)
		file_inds[pairs3[i].file_idx]=true;

	
	FileManager fm;
	FileSet fs;

	fm.init_from_list_file(config,mzxml_list,file_inds);
	fs.select_all_files(fm);
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

	

	// read spectra
	BasicSpecReader bsr;
	QCPeak peaks[5000];

	int num_out_neg=0, num_out1=0, num_out2=0, num_out3=0;
	int neg_idx=0,c1_idx=0,c2_idx=0,c3_idx=0;


	for (i=0; i<all_ssf.size(); i++)
	{
		MZXML_single *ssf = (MZXML_single *)all_ssf[i];
		ScanPair ssf_pair(ssf->file_idx,ssf->scan_number);
		string seq="";

		int out_dest=-1;

		while (neg_idx<neg_pairs.size() && neg_pairs[neg_idx]<ssf_pair)
			neg_idx++;

		if (neg_idx<neg_pairs.size() && neg_pairs[neg_idx]==ssf_pair)
			out_dest=0;


		while (c1_idx<pairs1.size() && pairs1[c1_idx]<ssf_pair)
			c1_idx++;
		if (c1_idx<pairs1.size() && pairs1[c1_idx]==ssf_pair)
		{
			seq = pairs1[c1_idx].seq;
			out_dest=1;
		}


		while (c2_idx<pairs2.size() && pairs2[c2_idx]<ssf_pair)
			c2_idx++;
		if (c2_idx<pairs2.size() && pairs2[c2_idx]==ssf_pair)
		{
			seq = pairs2[c2_idx].seq;
			out_dest=2;
		}


		while (c3_idx<pairs3.size() && pairs3[c3_idx]<ssf_pair)
			c3_idx++;
		if (c3_idx<pairs3.size() && pairs3[c3_idx]==ssf_pair)
		{
			seq = pairs3[c3_idx].seq;
			out_dest=3;
		}

		if (out_dest<0)
			continue;

		BasicSpectrum bs;
		bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.ssf = ssf;

	//	if (out_dest>0)
	//		bs.ssf->peptide.parse_from_string(config,seq);
	
		char name_buff[64];
		if (out_dest==0)
		{
			sprintf(name_buff,"train_neg_%d_%d_%d",num_out_neg,ssf->file_idx,ssf->scan_number);
			bs.ssf->single_name = string(name_buff);
			bs.output_to_mgf(stream_neg,config);
			num_out_neg++;
			continue;
		}

		if (out_dest==1)
		{
			sprintf(name_buff,"train_pos1_%d_%d_%d",num_out1,ssf->file_idx,ssf->scan_number);
			bs.ssf->single_name = string(name_buff);
			bs.output_to_mgf(stream1,config,seq.c_str());
			num_out1++;
			continue;
		}

		if (out_dest==2)
		{
			sprintf(name_buff,"train_pos2_%d_%d_%d",num_out2,ssf->file_idx,ssf->scan_number);
			bs.ssf->single_name = string(name_buff);
			bs.output_to_mgf(stream2,config,seq.c_str());
			num_out2++;
			continue;
		}

		if (out_dest==3)
		{
			sprintf(name_buff,"train_pos3_%d_%d_%d",num_out3,ssf->file_idx,ssf->scan_number);
			bs.ssf->single_name = string(name_buff);
			bs.output_to_mgf(stream3,config,seq.c_str());
			num_out3++;
			continue;
		}
	}

	cout << "Wrote: " << endl;
	cout << "Neg " << num_out_neg << " / " << neg_pairs.size() << endl;
	cout << "Pos1 " << num_out1 << " / " << pairs1.size() << endl;
	cout << "Pos2 " << num_out2 << " / " << pairs2.size() << endl;
	cout << "Pos3 " << num_out3 << " / " << pairs3.size() << endl;

	stream_neg.close();
	stream1.close();
	stream2.close();
	stream3.close();
	
}





void PMCSQS_Scorer::print_spec(const BasicSpectrum& bs) const
{
	cout << bs.ssf->single_name << endl;
	int i;
	for (i=0; i<bs.num_peaks; i++)
	{
		cout << setprecision(2) << fixed << bs.peaks[i].mass << "\t" << bs.peaks[i].intensity << "\t";
		if (currentSpectrumIsotopeLevels_[i]>0)
			cout << " ISO ";
		if (currentSpectrumStrongPeakFlags_[i])
			cout << " STRONG ";
		cout << endl;
	}
}





