#include "PMCSQS.h"
#include "auxfun.h"
#include "SpectraList.h"
#include "AnnotatedSpectrum.h"

extern const char * SQS_var_names[];

/********************************************************************
Select a set of mzs nd charges that have high probabilities.
If unless the use_spectrum_mz is set to 1 in the config, it adds -1,+1  
(if there are too many mzs, not all these offfsets will be added)
Uses emprically set probability threhsolds to chose corrected pms.
*********************************************************************/
void PMCSQS_Scorer::selectPrecursorMassesAndCharges(
								const Config *config, 
								const PeakList& pl,
								vector<mass_t>& precursorMassesWith19,
								vector<int>&    charges,
								vector<PmcSqsChargeRes>* allResults)
{
	static const mass_t offsets[]={-MASS_ISO,MASS_ISO};
	static const int numAllOffsets = sizeof(offsets)/sizeof(mass_t);
	static const float minimalCompareProbForAdding  = 0.04;
	static const float minimalSqsProbForAdding		= 0.05;
	const mass_t halfTolerance = config->getTolerance() * 0.5;
	
	const SingleSpectrumHeader* const header = pl.getHeader();
	int		spectrumCharge	= header->getCharge();
	mass_t	spectrumMOverZ	= header->getMOverZ();
	const int originalSpectrumCharge = spectrumCharge;
	const mass_t originalSpectrumMOverZ = spectrumMOverZ;
	
	precursorMassesWith19.clear();
	charges.clear();

	int    specificCharge=0;
	mass_t specificMOverZ=0;
	
	if (config->get_use_spectrum_charge())
		specificCharge=spectrumCharge;

	static vector<PmcSqsChargeRes> results;
	int maximalprobabilityCharge=0;
	float maximalProbability = 0;

	if (spectrumCharge>0 && 
		config->get_use_spectrum_charge() && 
		config->get_use_spectrum_mz())
	{
		precursorMassesWith19.push_back( spectrumMOverZ * spectrumCharge - 
										(spectrumCharge-1)*MASS_PROTON );
		charges.push_back(spectrumCharge);
		maximalprobabilityCharge=spectrumCharge;
		if (allResults)
			allResults->clear();
		return;
	}
	else
	{
		if (! this->ind_initialized_pmcr)
		{
			cout << "Error: no PMC model was read!" << endl;
			if (header->getCharge() <=0)
				cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
			exit(1);
		}

		if (config->get_use_spectrum_charge() && spectrumCharge >= this->pmc_rank_models.size())
		{
			// adjust charge and m/z to lowest charge modeled
			spectrumCharge = pmc_rank_models.size()-1;
			mass_t mass_with_19 = maximalprobabilityCharge * spectrumMOverZ - (maximalprobabilityCharge-1)*MASS_PROTON;
			spectrumMOverZ = (mass_with_19 + (spectrumCharge -1 ) * MASS_PROTON) / spectrumCharge;

			// TODO there should be a better way to do this correction
			SingleSpectrumHeader *nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
			nonConstHeader->setCharge(spectrumCharge);
			nonConstHeader->setMOverZ(spectrumMOverZ);
		}
		

	//	get_pmcsqs_results_for_spectrum(config,bs,res);
		if (computePmcsqsResultsForSpectrum(config, pl, results) <= 0.0)
		{
			precursorMassesWith19.push_back(-1.0);
			charges.push_back(0);

			if (allResults)
				allResults->clear();
			

			return;
		}

		if (specificCharge>0)
		{
			maximalprobabilityCharge = spectrumCharge;
			maximalProbability = results[spectrumCharge].min_comp_prob;
		}
		else
		{
			int c;
			for (c=1; c<results.size(); c++)
			{
				if (results[c].min_comp_prob>maximalProbability)
				{
					maximalprobabilityCharge = c;
					maximalProbability = results[c].min_comp_prob;
				}
			}
		}

		precursorMassesWith19.push_back(results[maximalprobabilityCharge].mz1 * 
										maximalprobabilityCharge - (maximalprobabilityCharge-1)*MASS_PROTON);
		charges.push_back(maximalprobabilityCharge);
		if (results[maximalprobabilityCharge].mz2>0)
		{
			precursorMassesWith19.push_back(results[maximalprobabilityCharge].mz2 * 
										maximalprobabilityCharge - (maximalprobabilityCharge-1)*MASS_PROTON);
			charges.push_back(maximalprobabilityCharge);
		}

	}

	// only added to the first pm_with_19, check that it doesn't overlap with others
	if (! config->get_use_spectrum_mz())
	{
		int numOffsets = numAllOffsets;
		int i;
		for (i=0; i<numOffsets; i++)
		{
			const mass_t pm_with_19 = precursorMassesWith19[0] + offsets[i];
			int j;
			for (j=0; j<precursorMassesWith19.size(); j++)
				if (charges[j]==maximalprobabilityCharge && 
					fabs(pm_with_19 - precursorMassesWith19[j]) < halfTolerance)
					break;
	
			if (j==precursorMassesWith19.size())
			{
				precursorMassesWith19.push_back(pm_with_19);
				charges.push_back(maximalprobabilityCharge);
			}
		}

		// add tol to -1 of 2nd
		if (charges[0]>=2)
		{
			const mass_t pm_with_19 = precursorMassesWith19[1] - MASS_PROTON;
			int j;
			for (j=0; j<precursorMassesWith19.size(); j++)
				if (charges[j]==maximalprobabilityCharge &&
					fabs(pm_with_19 - precursorMassesWith19[j]) < halfTolerance)
						break;
			
			if (j==precursorMassesWith19.size())
			{
				precursorMassesWith19.push_back(pm_with_19);
				charges.push_back(maximalprobabilityCharge);
			}
		}
		

		if (charges[0]>=3)
		{
			const mass_t pm_with_19 = precursorMassesWith19[1] + MASS_PROTON;
			int j;
			for (j=0; j<precursorMassesWith19.size(); j++)
				if (charges[j]==maximalprobabilityCharge &&
					fabs(pm_with_19 - precursorMassesWith19[j]) < halfTolerance)
						break;

			if (j==precursorMassesWith19.size())
			{
				precursorMassesWith19.push_back(pm_with_19);
				charges.push_back(maximalprobabilityCharge);
			}
		}
	}

	// add other charges if their comp probability and sqs probs are high enough
	if (specificCharge==0)
	{
		// find best charge
		int c; 
		float max_prob=-1.0;

		for (c=1; c<results.size(); c++)
		{
			if (c==maximalprobabilityCharge)
				continue;

			if (results[c].min_comp_prob > minimalCompareProbForAdding &&
				results[c].sqs_prob > minimalSqsProbForAdding)
			{
				precursorMassesWith19.push_back(results[c].mz1 * c - (c-1)*MASS_PROTON);
				charges.push_back(c);
			}
		}
	}

	if (allResults)
		*allResults = results;

	if (spectrumCharge != maximalprobabilityCharge)
	{
		// TODO there should be a better way to do this correction
		SingleSpectrumHeader *nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
		nonConstHeader->setCharge(maximalprobabilityCharge);
		nonConstHeader->setMOverZ(originalSpectrumMOverZ);
	}
}


/***************************************************************************
Gives detailed m/z and prob values for the different charges.
Returns the highest sqs probability found for a spectrum.
****************************************************************************/
float PMCSQS_Scorer::computePmcsqsResultsForSpectrum(const Config* config,
													 const PeakList& pl,
													 vector<PmcSqsChargeRes>& results)
{
	ME_Regression_Sample sqs_sam;

	if (! initializeForCurrentSpectrum(config, pl))
		return 0; // corrupt file

//	calculate_curr_spec_pmc_values(bs,bin_increment);
	calculateCurrentSpectrumPmcValues(pl, bin_increment);

	const int sqs_size_idx = getSqsSizeIndex(pl.getHeader()->getMOverZ());

	if (results.size() <= maximalChargeWithModels_)
		results.resize(maximalChargeWithModels_+1);

	
	vector<float> minComparisonProbs;	// the minimal probability when used in a comparison model
	minComparisonProbs.resize(maximalChargeWithModels_+1,2.0);
	minComparisonProbs[0]=0;

	if (ind_initialized_sqs)
	{
		const int maxSqsCharge = ( sqs_models.size()<maximalChargeWithModels_ ?
			sqs_models.size() : maximalChargeWithModels_);

	//	fill_fval_vector_with_SQS(bs, sqs_sam);
		fillSqsMeSample(pl, sqs_sam);
	
		int charge;
		for (charge=1; charge<=maxSqsCharge; charge++)
		{
			if (sqs_models[charge].size()<1 ||
				sqs_models[charge][0].size() <= sqs_size_idx)
				continue;
			
			float prob = sqs_models[charge][0][sqs_size_idx]->p_y_given_x(0,sqs_sam);

		//	cout << "SQS: " << charge << "\t" << prob;
			// correct prob
			prob += sqs_correction_factors[charge][sqs_size_idx];
			prob *= sqs_mult_factors[charge][sqs_size_idx];
			if (prob<0.0)
				prob=0.0;

			results[charge].sqs_prob=prob;

		//	 cout << "\t" << prob << endl;
		}
		
		int c;
		for (c=1; c<maxSqsCharge; c++)
		{
			if (sqs_models[c].size()<=0)
				continue;

			float comp_prob = 2.0;
			int d;
			for (d=c+1; d<=maxSqsCharge; d++)
			{
				if (sqs_models[d][c].size() <= sqs_size_idx ||
					! sqs_models[d][c][sqs_size_idx])
					continue;

			//	if (c==2 && d == 3)
			//		sqs_sam.print(SQS_var_names);

				comp_prob=sqs_models[d][c][sqs_size_idx]->p_y_given_x(0,sqs_sam);

				if (comp_prob<minComparisonProbs[d])
					minComparisonProbs[d]=comp_prob;

				const float one_minus_prob = 1.0-comp_prob;
				if (one_minus_prob<minComparisonProbs[c])
					minComparisonProbs[c]=one_minus_prob;

			//	cout << "COMP:  " << c << ": " << minComparisonProbs[c] << "\t" << d << ": " << minComparisonProbs[d] << endl;
			}
		}
		
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (minComparisonProbs[c]>1)
				minComparisonProbs[c]=0;
	}
	else // give the max prob charge to the input charge, the rest get 0
	{
		const int charge = pl.getHeader()->getCharge();
		if (charge<=0)
		{
			cout << "Error: no SQS model was read, so charge must be supplied in the spectrum!" << endl;
			exit(1);
		}

		int c;
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (c != charge)
				minComparisonProbs[c]=0;
	}

	
	float max_prob=0;
	int charge;

	const int max_pmc_charge = (pmc_rank_models.size() < maximalChargeWithModels_ ? 
		pmc_rank_models.size() : maximalChargeWithModels_);

	for (charge=1; charge<=max_pmc_charge; charge++)
	{
		const float prob = minComparisonProbs[charge];
		if (prob>max_prob)
			max_prob=prob;

		results[charge].min_comp_prob = prob;

		if (prob>0.02)
		{
		//	find_best_mz_values_from_rank_model(bs,charge,config->get_pm_tolerance(),results[charge]);
			computeBestMzValuesForCharge(pl, charge, config->get_pm_tolerance(), results[charge]);

		//	cout << "RES: " << charge << "\t" << results[charge].score1 << "\t" << results[charge].mz1 << endl;
		//	cout << "RES: " << charge << "\t" << results[charge].score2 << "\t" << results[charge].mz2 << endl;
		}
		else // only give one guess, the second one will come from another charge
		{
			results[charge].mz1= pl.getHeader()->getMOverZ();
			results[charge].score1=0;
			results[charge].mz2=NEG_INF;
			results[charge].score2=NEG_INF;
		}
	}
	

	return max_prob;
}


bool PMCSQS_Scorer::initializeForCurrentSpectrum(const Config* config, const PeakList& pl)
{
	config_ = config;

	pl.calculateIsotopicLevels(config->getTolerance(), currentSpectrumIsotopeLevels_);
	if (! pl.selectStrongPeakIndexesForPmcsqs(currentSpectrumIsotopeLevels_, currentSpectrumStrongPeakFlags_))
		return false;

	pl.markAllPossibleIsotopicPeaks(config->getTolerance(),currentSpectrumStrictIsotopeFlags_);

/*	bs.calc_peak_isotope_levels(config->getTolerance(),this->currentSpectrumIsotopeLevels_);
	if (! bs.select_strong_peak_idxs(this->currentSpectrumIsotopeLevels_,this->currentSpectrumStrongPeakFlags_))
		return false;
	bs.mark_all_possible_isotope_peaks(config->getTolerance(),currentSpectrumStrictIsotopeFlags_);*/

	curr_spec_total_intensity=0;
	curr_spec_strong_intensity=0;
	curr_spec_num_strong=0;
	int i;
	for (i=0; i<pl.getNumPeaks(); i++)
	{
		if (currentSpectrumIsotopeLevels_[i]>0)
			continue;

		curr_spec_total_intensity += pl.getPeakIntensity(i);
		if (currentSpectrumStrongPeakFlags_[i])
		{
			curr_spec_strong_intensity += pl.getPeakIntensity(i);
			curr_spec_num_strong++;
		}
	}


	if (currentSpectrumPmcTables_.size() <= maximalChargeWithModels_+1)
	{
		currentSpectrumPmcTables_.clear();
		currentSpectrumPmcTables_.resize(maximalChargeWithModels_+1);

		currentSpectrumBackgroundStats_.clear();
		currentSpectrumBackgroundStats_.resize(maximalChargeWithModels_+1);

		currentSpectrumMaximalValues_.clear();
		currentSpectrumMaximalValues_.resize(maximalChargeWithModels_+1);
	}

	return true;
}

void PMCSQS_Scorer::calculateCurrentSpectrumPmcValues(const PeakList& pl, mass_t binIncrement)
{
	if (frag_pair_sum_offset<-999)
	{
		set_frag_pair_sum_offset(MASS_PROTON); // b+y - PM+19
		set_bin_increment(0.1);
	}

	const mass_t mOverZ = pl.getHeader()->getMOverZ();
	int charge;
	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		const mass_t thisChargePmWith19 = mOverZ* charge - MASS_PROTON*(charge - 1);
		const mass_t frag_pair_sum = thisChargePmWith19 + frag_pair_sum_offset;
		const int size_idx = this->get_rank_model_size_idx(charge,thisChargePmWith19);

		float bias = 0;
		if (pmcMzBiases_.size()>0 && pmcMzBiases_[charge].size()>size_idx)
			bias = this->pmcMzBiases_[charge][size_idx];

		if (config_->get_pm_tolerance()<0.075) // don't use bias if dealing with accurate pm
			bias=0;

		int bin_range =3*charge;
		if (charge>=2)
			bin_range = int(2*charge);
		if (charge>=3)
			bin_range = int(1.5*charge);


		fillPmcRankStatistics(charge,
							  frag_pair_sum + bias,
							  -bin_range-1, 
							  bin_range-1,
							  binIncrement,
							  config_,
							  pl,
							  currentSpectrumStrongPeakFlags_,
							  currentSpectrumIsotopeLevels_,
							  currentSpectrumStrictIsotopeFlags_,
							  currentSpectrumPmcTables_[charge]);


		calculateBackgroundStatistics(frag_pair_sum, config_, pl, 
			currentSpectrumStrongPeakFlags_, currentSpectrumIsotopeLevels_, 
			this->currentSpectrumStrictIsotopeFlags_, currentSpectrumBackgroundStats_[charge]);


		// find maximal values
		int i;
		PMCRankStats& maximal = currentSpectrumMaximalValues_[charge];
		maximal = currentSpectrumPmcTables_[charge][0];

		for (i=1; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			const PMCRankStats& curr_stats = currentSpectrumPmcTables_[charge][i];

			// mathced pairs look for maximum number

			if (curr_stats.numCharge2FragmentPairs>maximal.numCharge2FragmentPairs)
				maximal.numCharge2FragmentPairs = curr_stats.numCharge2FragmentPairs;

			if (curr_stats.numFragmentPairs>maximal.numFragmentPairs)
				maximal.numFragmentPairs = curr_stats.numFragmentPairs;

			if (curr_stats.numStrongFragmentPairs>maximal.numStrongFragmentPairs)
				maximal.numStrongFragmentPairs = curr_stats.numStrongFragmentPairs;

			if (curr_stats.numStrongCharge2FragmentPairs > maximal.numStrongCharge2FragmentPairs)
				maximal.numStrongCharge2FragmentPairs = curr_stats.numStrongCharge2FragmentPairs;

			// intensity look for maximal numbers

			if (curr_stats.inten_frag_pairs > maximal.inten_frag_pairs)
				maximal.inten_frag_pairs = curr_stats.inten_frag_pairs;
			
			if (curr_stats.inten_strong_pairs> maximal.inten_strong_pairs)
				maximal.inten_strong_pairs = curr_stats.inten_strong_pairs;

			if (curr_stats.inten_c2_pairs> maximal.inten_c2_pairs)
				maximal.inten_c2_pairs = curr_stats.inten_c2_pairs;

			if (curr_stats.inten_c2_strong_pairs > maximal.inten_c2_strong_pairs)
				maximal.inten_c2_strong_pairs = curr_stats.inten_c2_strong_pairs;


	
		}

		// find indicators for min tolerance and max pairs
		float tol_pairs =           POS_INF;
		float tol_strong_pairs =    POS_INF;
		float tol_c2_pairs =        POS_INF;
		float tol_c2_strong_pairs = POS_INF;

		int   idx_pairs=0;
		int	  idx_strong_pairs=0;
		int   idx_c2_pairs=0;
		int   idx_c2_strong_pairs=0;

		for (i=0; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			const PMCRankStats& curr_stats = currentSpectrumPmcTables_[charge][i];

			if (curr_stats.numFragmentPairs == maximal.numFragmentPairs &&
				curr_stats.mean_offset_pairs < tol_pairs)
			{
				idx_pairs=i;
				tol_pairs=curr_stats.mean_offset_pairs;
			}

			if (curr_stats.numStrongFragmentPairs == maximal.numStrongFragmentPairs &&
				curr_stats.mean_offset_c2_strong_pairs < tol_strong_pairs)
			{
				idx_strong_pairs=i;
				tol_strong_pairs=curr_stats.mean_offset_c2_strong_pairs;
			}

			if (curr_stats.numCharge2FragmentPairs == maximal.numCharge2FragmentPairs &&
				curr_stats.mean_offset_c2_pairs < tol_c2_pairs)
			{
				idx_c2_pairs=i;
				tol_c2_pairs=curr_stats.mean_offset_c2_pairs;
			}

			if (curr_stats.numStrongCharge2FragmentPairs == maximal.numStrongCharge2FragmentPairs &&
				curr_stats.mean_offset_c2_strong_pairs < tol_c2_strong_pairs)
			{
				idx_c2_strong_pairs=i;
				tol_c2_strong_pairs=curr_stats.mean_offset_c2_strong_pairs;
			}

		}

		currentSpectrumPmcTables_[charge][idx_pairs].ind_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_strong_pairs].ind_strong_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_c2_pairs].ind_c2_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_c2_strong_pairs].ind_c2_strong_pairs_with_min_tol=true;


		static vector<float> log_distances;
		if (log_distances.size()<currentSpectrumPmcTables_[charge].size())
		{
			log_distances.resize(currentSpectrumPmcTables_[charge].size(),0);
			int i;
			for (i=1; i<log_distances.size(); i++)
				log_distances[i]=log(1.0+(float)i);
		}

		for (i=0; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			PMCRankStats& curr_stats = currentSpectrumPmcTables_[charge][i];
			curr_stats.log_dis_from_pairs_min_tol = log_distances[abs(i-idx_pairs)];
			curr_stats.log_dis_from_strong_pairs_min_tol = log_distances[abs(i-idx_strong_pairs)];
			curr_stats.log_dis_from_c2_pairs_min_tol = log_distances[abs(i-idx_c2_pairs)];
			curr_stats.log_dis_from_c2_strong_pairs_min_tol = log_distances[abs(i-idx_c2_strong_pairs)];
		}

	}
}



PMCRankStats::PMCRankStats(const PMCRankStats & Source)
{
	*this = Source;
//	offset_pairs_ordered_by_inten		 = Source.offset_pairs_ordered_by_inten;
//	strong_offset_pairs_ordered_by_inten = Source.strong_offset_pairs_ordered_by_inten;
//  c2_offset_pairs_ordered_by_inten	 =  Source.c2_offset_pairs_ordered_by_inten;
}


void PMCSQS_Scorer::set_sqs_mass_thresholds()
{
	sqsMassThresholds_.clear();
	sqsMassThresholds_.push_back(800);
	sqsMassThresholds_.push_back(1200);
}


void PMCSQS_Scorer::set_default_sqs_correct_factors()
{
	sqs_correction_factors.resize(4);
	sqs_mult_factors.resize(4);

	sqs_correction_factors[1].push_back(0);
	sqs_correction_factors[1].push_back(-0.1);
	sqs_correction_factors[1].push_back(-0.15);

	sqs_correction_factors[2].push_back(0);
	sqs_correction_factors[2].push_back(0.05);
	sqs_correction_factors[2].push_back(0.05);

	sqs_correction_factors[2].push_back(0);
	sqs_correction_factors[2].push_back(0.05);
	sqs_correction_factors[2].push_back(0.05);

	int i;
	for (i=0; i<4; i++)
	{
		int j;
		for (j=0; j<sqs_correction_factors.size(); j++)
			sqs_mult_factors[i].push_back(1.0/(1.0+sqs_correction_factors[i][j]));
	}
}


void PMCSQS_Scorer::init_sqs_correct_factors(int max_charge, int num_sizes)
{
	sqs_correction_factors.resize(max_charge+max_charge);
	sqs_mult_factors.resize(max_charge+1);

	int i;
	for (i=1; i<=max_charge; i++)
	{
		sqs_correction_factors[i].clear();
		sqs_correction_factors[i].resize(num_sizes,0);
		sqs_mult_factors[i].clear();
		sqs_mult_factors[i].resize(num_sizes,1);
	}
}





/***************************************************************************
Gives detailed m/z and prob values for the different charges.
Returns the highest sqs probability found for a spectrum.
****************************************************************************/
float PMCSQS_Scorer::get_pmcsqs_results_for_spectrum(const Config *config, 
													 const BasicSpectrum& bs,
													 vector<PmcSqsChargeRes>& res)
{
	static ME_Regression_Sample sqs_sam;

	if (! init_for_current_spec(config,bs))
		return 0; // corrupt file

	calculate_curr_spec_pmc_values(bs,bin_increment);

	const int sqs_size_idx = this->getSqsSizeIndex(bs.ssf->m_over_z);

	if (res.size()<=maximalChargeWithModels_)
		res.resize(maximalChargeWithModels_+1);

	
	vector<float> minComparisonProbs;		// the minimal probability when used in a comparison model
	minComparisonProbs.resize(maximalChargeWithModels_+1,2.0);
	minComparisonProbs[0]=0;

	if (ind_initialized_sqs)
	{
		const int maxSqsCharge = ( sqs_models.size()<maximalChargeWithModels_ ?
			sqs_models.size() : maximalChargeWithModels_);

		fill_fval_vector_with_SQS(bs, sqs_sam);
	
		int charge;
		for (charge=1; charge<=maxSqsCharge; charge++)
		{
			if (sqs_models[charge].size()<1 ||
				sqs_models[charge][0].size() <= sqs_size_idx)
				continue;
			
			float prob = sqs_models[charge][0][sqs_size_idx]->p_y_given_x(0,sqs_sam);

			// correct prob
			prob += sqs_correction_factors[charge][sqs_size_idx];
			prob *= sqs_mult_factors[charge][sqs_size_idx];
			if (prob<0.0)
				prob=0.0;

			res[charge].sqs_prob=prob;
		}
		
		int c;
		for (c=1; c<maxSqsCharge; c++)
		{
			if (sqs_models[c].size()<=0)
				continue;

			float comp_prob = 2.0;
			int d;
			for (d=c+1; d<=maxSqsCharge; d++)
			{
				if (sqs_models[d][c].size() <= sqs_size_idx ||
					! sqs_models[d][c][sqs_size_idx])
					continue;

				comp_prob=sqs_models[d][c][sqs_size_idx]->p_y_given_x(0,sqs_sam);

				if (comp_prob<minComparisonProbs[d])
					minComparisonProbs[d]=comp_prob;

				const float one_minus_prob = 1.0-comp_prob;
				if (one_minus_prob<minComparisonProbs[c])
					minComparisonProbs[c]=one_minus_prob;
			}
		}
		
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (minComparisonProbs[c]>1)
				minComparisonProbs[c]=0;
	}
	else // give the max prob charge to the input charge, the rest get 0
	{
		if (bs.ssf->charge<=0)
		{
			cout << "Error: no SQS model was read, so charge must be supplied in the spectrum!" << endl;
			exit(1);
		}

		int c;
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (c != bs.ssf->charge)
				minComparisonProbs[c]=0;
	}

	
	float max_prob=0;
	int charge;

	const int max_pmc_charge = (pmc_rank_models.size() < maximalChargeWithModels_ ? 
		pmc_rank_models.size() : maximalChargeWithModels_);

	for (charge=1; charge<=max_pmc_charge; charge++)
	{
		const float prob = minComparisonProbs[charge];
		if (prob>max_prob)
			max_prob=prob;

		res[charge].min_comp_prob = prob;

		if (prob>0.02)
		{
			find_best_mz_values_from_rank_model(bs,charge,config->get_pm_tolerance(),res[charge]);
		}
		else // only give one guess, the second one will come from another charge
		{
			res[charge].mz1=bs.ssf->m_over_z;
			res[charge].score1=0;
			res[charge].mz2=NEG_INF;
			res[charge].score2=NEG_INF;
		}
	}
	

	return max_prob;
}







float PMCSQS_Scorer::computeBestMzAndChargeForSpectrum(const Config *config, const PeakList& pl,
							mass_t* mz1, int* charge1, float *prob1,
							mass_t* mz2, int* charge2, float *prob2,
							vector<PmcSqsChargeRes>* all_res)
{
	static vector<PmcSqsChargeRes> res;

	if (! this->ind_initialized_pmcr)
	{
		cout << "Error: no PMC model was read!" << endl;
		if (pl.getHeader()->getCharge() <=0)
			cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
		exit(1);
	}


	if (! computePmcsqsResultsForSpectrum(config, pl, res))
	{
		*mz1=-1.0;
		*charge1=0;
		*mz2=-1.0;
		*charge2=0;
		return 0.0;
	}

	float best_prob=-1;
	int best_charge=0;

	int charge;
	for (charge =1; charge<=maximalChargeWithModels_; charge++)
		if (res[charge].min_comp_prob>best_prob)
		{
			best_charge = charge;
			best_prob = res[charge].min_comp_prob;
		}

	const PmcSqsChargeRes& cr = res[best_charge];
	float second_best_prob=1.0 - (cr.score1-cr.score2)/(fabs(cr.score1)+fabs(cr.score2)+2.0);
	int   second_best_charge = best_charge;
	
	for (charge =1; charge<=maximalChargeWithModels_; charge++)
		if (charge != best_charge && res[charge].min_comp_prob>second_best_prob)
		{
			second_best_charge = charge;
			second_best_prob = res[charge].min_comp_prob;
		}

	*mz1 = (mass_t)res[best_charge].mz1;
	*charge1 = best_charge;
	*prob1 = best_prob;

	if (mz2 && charge2)
	{
		*charge2 = second_best_charge;
		*prob2   = second_best_prob;
		if (second_best_charge == best_charge)
		{
			*mz2 = (mass_t)res[best_charge].mz2;
		}
		else
			*mz2 = (mass_t)res[second_best_charge].mz1;
	}

	if (all_res)
		*all_res = res;

	return best_prob;
}

/******************************************************************************
Selects the the two best values of charges 1,2,3
returns the max prob found
*******************************************************************************/
float PMCSQS_Scorer::get_best_mz_charge(const Config *config, const BasicSpectrum& bs, 
						   mass_t* mz1, int* charge1, float *prob1,
						   mass_t* mz2, int* charge2, float *prob2,
						   vector<PmcSqsChargeRes>* all_res)
{
	static vector<PmcSqsChargeRes> res;

	if (! this->ind_initialized_pmcr)
	{
		cout << "Error: no PMC model was read!" << endl;
		if (bs.ssf->charge <=0)
			cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
		exit(1);
	}

	get_pmcsqs_results_for_spectrum(config,bs,res);

	float best_prob=-1;
	int best_charge=0;

	int charge;
	for (charge =1; charge<=maximalChargeWithModels_; charge++)
		if (res[charge].min_comp_prob>best_prob)
		{
			best_charge = charge;
			best_prob = res[charge].min_comp_prob;
		}

	const PmcSqsChargeRes& cr = res[best_charge];
	float second_best_prob=1.0 - (cr.score1-cr.score2)/(fabs(cr.score1)+fabs(cr.score2)+2.0);
	int   second_best_charge = best_charge;
	
	for (charge =1; charge<=maximalChargeWithModels_; charge++)
		if (charge != best_charge && res[charge].min_comp_prob>second_best_prob)
		{
			second_best_charge = charge;
			second_best_prob = res[charge].min_comp_prob;
		}

	*mz1 = (mass_t)res[best_charge].mz1;
	*charge1 = best_charge;
	*prob1 = best_prob;

	if (mz2 && charge2)
	{
		*charge2 = second_best_charge;
		*prob2   = second_best_prob;
		if (second_best_charge == best_charge)
		{
			*mz2 = (mass_t)res[best_charge].mz2;
		}
		else
			*mz2 = (mass_t)res[second_best_charge].mz1;
	}

	if (all_res)
		*all_res = res;

	return best_prob;
}



/********************************************************************
Select a set of mzs nd charges that have high probabilities.
If unless the use_spectrum_mz is set to 1 in the config, it adds -1,+1  
(if there are too many mzs, not all these offfsets will be added)
Uses emprically set probability threhsolds to chose corrected pms.
*********************************************************************/
void PMCSQS_Scorer::select_pms_and_charges(const Config *config, 
								const BasicSpectrum& bs,
								vector<mass_t>& pms_with_19,
								vector<int>&    charges,
								vector<PmcSqsChargeRes>* all_res)
{
	static const mass_t offsets[]={-MASS_ISO,MASS_ISO};
	static const int num_all_offsets = sizeof(offsets)/sizeof(mass_t);
	static const float min_comp_prob_for_adding_second = 0.04;
	static const float min_sqs_prob_for_adding_second = 0.05;
	const mass_t smallerTolerance = config->getTolerance() * 0.5;

	int spec_charge = bs.ssf->charge;
	mass_t spec_mz = bs.ssf->m_over_z;
	int org_spec_charge = spec_charge;
	mass_t org_spec_mz = spec_mz;
	
	pms_with_19.clear();
	charges.clear();
	int    specific_charge=0;
	mass_t specific_mz=0;
	

	if (config->get_use_spectrum_charge())
		specific_charge=spec_charge;

	static vector<PmcSqsChargeRes> res;
	int max_charge=0;
	float max_prob = 0;
	if (spec_charge>0 && 
		config->get_use_spectrum_charge() && 
		config->get_use_spectrum_mz())
	{
		pms_with_19.push_back(spec_mz * spec_charge - (spec_charge-1)*MASS_PROTON);
		charges.push_back(spec_charge);
		max_charge=spec_charge;
	}
	else
	{
		if (! this->ind_initialized_pmcr)
		{
			cout << "Error: no PMC model was read!" << endl;
			if (bs.ssf->charge <=0)
				cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
			exit(1);
		}

		if (config->get_use_spectrum_charge() && spec_charge>= this->pmc_rank_models.size())
		{
			// adjust charge and m/z to lowest charge modeled
			spec_charge = pmc_rank_models.size()-1;
			mass_t mass_with_19 = org_spec_charge * spec_mz - (org_spec_charge-1)*MASS_PROTON;
			spec_mz = (mass_with_19 + (spec_charge -1 ) * MASS_PROTON) / spec_charge;

			bs.ssf->charge = spec_charge;
			bs.ssf->m_over_z = spec_mz;
		}
		

		get_pmcsqs_results_for_spectrum(config,bs,res);

		if (specific_charge>0)
		{
			max_charge = spec_charge;
			max_prob = res[spec_charge].min_comp_prob;
		}
		else
		{
			int c;
			for (c=1; c<res.size(); c++)
			{
				if (res[c].min_comp_prob>max_prob)
				{
					max_prob = res[c].min_comp_prob;
					max_charge = c;
				}
			}
		}

		pms_with_19.push_back(res[max_charge].mz1 * max_charge - (max_charge-1)*MASS_PROTON);
		charges.push_back(max_charge);
		if (res[max_charge].mz2>0)
		{
			pms_with_19.push_back(res[max_charge].mz2 * max_charge - (max_charge-1)*MASS_PROTON);
			charges.push_back(max_charge);
		}

	}

	// only added to the first pm_with_19, check that it doesn't overlap with others
	if (! config->get_use_spectrum_mz())
	{
		int numOffsets = num_all_offsets;
		int i;
		for (i=0; i<numOffsets; i++)
		{
			const mass_t pm_with_19 = pms_with_19[0] + offsets[i];
			int j;
			for (j=0; j<pms_with_19.size(); j++)
				if (charges[j]==max_charge && fabs(pm_with_19 - pms_with_19[j])<smallerTolerance)
					break;
	
			if (j==pms_with_19.size())
			{
				pms_with_19.push_back(pm_with_19);
				charges.push_back(max_charge);
			}
		}

		// add tol to -1 of 2nd
		if (charges[0]>=2)
		{
			const mass_t pm_with_19 = pms_with_19[1] - MASS_PROTON;
			int j;
			for (j=0; j<pms_with_19.size(); j++)
				if (charges[j]==max_charge &&
					fabs(pm_with_19 - pms_with_19[j])<smallerTolerance)
						break;
			
			if (j==pms_with_19.size())
			{
				pms_with_19.push_back(pm_with_19);
				charges.push_back(max_charge);
			}
		}
		

		if (charges[0]>=3)
		{
			const mass_t pm_with_19 = pms_with_19[1] + MASS_PROTON;
			int j;
			for (j=0; j<pms_with_19.size(); j++)
				if (charges[j]==max_charge &&
					fabs(pm_with_19 - pms_with_19[j])<smallerTolerance)
						break;

			if (j==pms_with_19.size())
			{
				pms_with_19.push_back(pm_with_19);
				charges.push_back(max_charge);
			}
		}
	}

	// add other charges if their comp probability and sqs probs are high enough
	if (specific_charge==0)
	{
		// find best charge
		int c; 
		float max_prob=-1.0;

		for (c=1; c<res.size(); c++)
		{
			if (c==max_charge)
				continue;

			if (res[c].min_comp_prob > min_comp_prob_for_adding_second &&
				res[c].sqs_prob > min_sqs_prob_for_adding_second)
			{
				pms_with_19.push_back(res[c].mz1 * c - (c-1)*MASS_PROTON);
				charges.push_back(c);
			}
		}
	}

	if (all_res)
		*all_res = res;

	if (spec_charge != org_spec_charge)
	{
		bs.ssf->charge = org_spec_charge;
		bs.ssf->m_over_z = org_spec_mz;
	}
}



void PMCSQS_Scorer::benchmark_pm_selection(Config *config, FileManager& fm, mass_t pm_val_tol)
{
	const vector< vector< mass_t > >& threshes = config->get_size_thresholds();
	int c=1;

	for (c=1; c<threshes.size(); c++)
	{
		int s;
		for (s=0; s<threshes[c].size(); s++)
		{
			mass_t min_mz = (s>0 ? threshes[c][s-1]/c : 0);
			mass_t max_mz = threshes[c][s]/c;

			FileSet fs;
			fs.select_files_in_mz_range(fm,min_mz,max_mz,c);

			if (fs.get_total_spectra()<200)
				continue;

			cout << "CHARGE " << c <<" size " << s << "  (" << fs.get_total_spectra() << " spectra)" << endl;

			const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();
			BasicSpecReader bsr;
			vector<QCPeak> peaks;
			peaks.resize(10000);

			vector<int> correct_counts;
			correct_counts.resize(8,0);

			int num_correct=0;
			int num_wrong_charge=0;
			int num_diff_charge=0;

			int i;
			for (i=0; i<all_ssfs.size(); i++)
			{
				SingleSpectrumFile *ssf=all_ssfs[i];
				BasicSpectrum bs;

				bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
				bs.peaks = &peaks[0];
				bs.ssf = ssf;

				const mass_t true_mass = ssf->peptide.get_mass_with_19();

				vector<mass_t> pms_with_19;
				vector<int>    charges;
				select_pms_and_charges(config,bs,pms_with_19,charges);

				if (charges[0] != c)
					num_wrong_charge++;

				bool got_diff_charge=false;
				int j;
				for (j=1; j<charges.size(); j++)
					if (charges[j] != charges[0])
						got_diff_charge=true;

				if (got_diff_charge)
					num_diff_charge++;

				for (j=0; j<pms_with_19.size(); j++)
					if (fabs(pms_with_19[j]-true_mass)<pm_val_tol)
						break;
				
				if (j==pms_with_19.size())
					continue;
				

				num_correct++;

				if (got_diff_charge && j == pms_with_19.size()-1)
				{
					correct_counts[7]++;
				}
				else
					correct_counts[j]++;
			}

			double num_total = (double)all_ssfs.size();
			cout << "Had correct       " << fixed << setprecision(4) << num_correct/num_total << endl;
			cout << "First correct     " << correct_counts[0]/num_total << endl;
			cout << "Second correct    " << correct_counts[1]/num_total << endl;
			cout << "Off-1  correct    " << correct_counts[2]/num_total << endl;
			cout << "Off+1  correct    " << correct_counts[3]/num_total << endl;
			cout << "Off-1  2nd        " << correct_counts[4]/num_total << endl;
			cout << "Off+1  2nd        " << correct_counts[5]/num_total << endl;
			cout << "Diff Ch correct   " << correct_counts[7]/num_total << endl;
			cout << "With wrong charge " << num_wrong_charge/num_total << endl;
			cout << "With diff charge  " << num_diff_charge/num_total << endl << endl;
		}
	}
}


const int DPColumnBytes = sizeof(int)*(Val+1);


// for each peak,aa holds the idx of the previous aa if they have a mass
// diff of that aa (within tolerance)
// entry 0 in each column holds an indicator if peak is in aa diff
// entry 1 in each column holds an indicator if peak has a 
void PMCSQS_Scorer::fillSqsDynamicProgrammingTable(const PeakList& pl, 
												   vector<DPColumn>& dp, 
												   int fargmentCharge ) const
{
	const Peak* const peaks = pl.getPeaks();
	const size_t numPeaks   = pl.getNumPeaks();
	const vector<mass_t>& aa2mass = config_->get_aa2mass();
	const mass_t tagTolerance    = (config_->getTolerance()<0.15 ? 
									config_->getTolerance() : config_->getTolerance()*0.33);

	const mass_t multVal = (1.0 / static_cast<mass_t>(fargmentCharge));
	dp.resize(numPeaks);
	int i;
	for (i=0; i<numPeaks; i++)
	{
		dp[i].pointers[0]=0;
		dp[i].fPointers[0]=0;
		int j;
		for (j=1; j<=Val; j++)
		{
			dp[i].pointers[j]=-1;
			dp[i].fPointers[j]=-1;
		}
	}

	int aa;
	for (aa=Ala; aa<=Val; aa++)
	{
		if (aa==Ile || aa==Xle)
			continue;

		const mass_t aamass = multVal * aa2mass[aa];
		const mass_t minimalOffset = aamass-tagTolerance;
		const mass_t maximalOffset = aamass+tagTolerance;

		size_t trailIndex=0;
		size_t leadIndex=1;

		while (leadIndex<numPeaks)
		{
			if (currentSpectrumIsotopeLevels_[leadIndex]>0)
			{
				leadIndex++;
				continue;
			}

			while (peaks[leadIndex].mass-peaks[trailIndex].mass>maximalOffset)
				trailIndex++;

			if (currentSpectrumIsotopeLevels_[trailIndex]==0 && 
				peaks[leadIndex].mass-peaks[trailIndex].mass>minimalOffset)
			{
				dp[leadIndex].pointers[aa]=trailIndex;
				dp[leadIndex].pointers[0]=1;
				dp[trailIndex].pointers[0]=1;
				dp[trailIndex].fPointers[aa]=leadIndex;
			}
			else 
				dp[leadIndex].pointers[aa]=-1;

			leadIndex++;
		}
	}
}



// for each peak,aa holds the idx of the previous aa if they have a mass
// diff of that aa (within tolerance)
// entry 0 in each column holds an indicator if peak is in aa diff
// entry 1 in each column holds an indicator if peak has a 
void PMCSQS_Scorer::fill_SQS_DP(const BasicSpectrum& bs, vector<DPColumn>& dp, int charge ) const
{
	const QCPeak *peaks = bs.peaks;
	const int num_peaks = bs.num_peaks;
	const vector<mass_t>& aa2mass = config_->get_aa2mass();
	const mass_t tagTolerance = (config_->getTolerance()<0.15 ? 
									config_->getTolerance() : config_->getTolerance()*0.33);

	
	const mass_t multVal = (1.0 / charge);
	dp.resize(num_peaks);
	int i;
	for (i=0; i<num_peaks; i++)
	{
		dp[i].pointers[0]=0;
		int j;
		for (j=1; j<=Val; j++)
			dp[i].pointers[j]=-1;
	}

	int aa;
	for (aa=Ala; aa<=Val; aa++)
	{
		if (aa==Ile || aa==Xle)
			continue;

		const mass_t aamass = multVal * aa2mass[aa];
		const mass_t minimalOffset = aamass-tagTolerance;
		const mass_t maximalOffset = aamass+tagTolerance;

		int trailIndex=0;
		int leadIndex=1;

		while (leadIndex<num_peaks)
		{
			if (currentSpectrumIsotopeLevels_[leadIndex]>0)
			{
				leadIndex++;
				continue;
			}

			while (peaks[leadIndex].mass-peaks[trailIndex].mass>maximalOffset)
				trailIndex++;

			if (currentSpectrumIsotopeLevels_[trailIndex]==0 && 
				peaks[leadIndex].mass-peaks[trailIndex].mass>minimalOffset)
			{
				dp[leadIndex].pointers[aa]=trailIndex;
				dp[leadIndex].pointers[0]=1;
				dp[trailIndex].pointers[0]=1;

			//	cout << "Off: " << peaks[leadIndex].mass-peaks[trailIndex].mass - minimalOffset<< endl;
			}
			else 
				dp[leadIndex].pointers[aa]=-1;

			leadIndex++;
		}

	}
}


/*****************************************************************************
Does the raw calculations required for processing the spectrum.
******************************************************************************/
bool PMCSQS_Scorer::init_for_current_spec(const Config *_config, 
										  const BasicSpectrum& bs)
{
	config_ = _config;

	bs.calc_peak_isotope_levels(config_->getTolerance(),this->currentSpectrumIsotopeLevels_);
	if (! bs.select_strong_peak_idxs(this->currentSpectrumIsotopeLevels_,this->currentSpectrumStrongPeakFlags_))
		return false;

	bs.mark_all_possible_isotope_peaks(config_->getTolerance(),currentSpectrumStrictIsotopeFlags_);

	curr_spec_total_intensity=0;
	curr_spec_strong_intensity=0;
	curr_spec_num_strong=0;
	int i;
	for (i=0; i<bs.num_peaks; i++)
	{
		if (currentSpectrumIsotopeLevels_[i]>0)
			continue;

		curr_spec_total_intensity+=bs.peaks[i].intensity;
		if (currentSpectrumStrongPeakFlags_[i])
		{
			curr_spec_strong_intensity+=bs.peaks[i].intensity;
			curr_spec_num_strong++;
		}
	}


	if (currentSpectrumPmcTables_.size()<=maximalChargeWithModels_+1)
	{
		currentSpectrumPmcTables_.clear();
		currentSpectrumPmcTables_.resize(maximalChargeWithModels_+1);

		currentSpectrumBackgroundStats_.clear();
		currentSpectrumBackgroundStats_.resize(maximalChargeWithModels_+1);

		currentSpectrumMaximalValues_.clear();
		currentSpectrumMaximalValues_.resize(maximalChargeWithModels_+1);
	}

	return true;
}


void PMCSQS_Scorer::fillSqsMeSample(const PeakList& pl, ME_Regression_Sample& sample) const
{
	const mass_t tolerance = config_->getTolerance();
	const mass_t fragmentTolerance = (tolerance<0.15 ? tolerance : tolerance * 0.5);
	const Peak* const peaks  = pl.getPeaks();
	const int  numPeaks		 = pl.getNumPeaks();
	const mass_t mOverZ		 = pl.getHeader()->getMOverZ();
	const int numStrongPeaks=curr_spec_num_strong;
	const float totalIntensity=curr_spec_total_intensity;
	const float oneOverNumberOfPeaks = 10.0 / (numPeaks+1);
	const float oneOverNumberOfStrongPeaks = 5.0/ (numStrongPeaks + 1);
	const float oneOverTotalIntensity = 1.0 / (totalIntensity + 1.0);
	const mass_t maximalPeakMass = (numPeaks <5 ? POS_INF : peaks[numPeaks-1].mass);

	sample.f_vals.clear();
	sample.f_vals.push_back(fval(SQS_CONST,1.0));
	sample.f_vals.push_back(fval(SQS_PEAK_DENSITY, (float)numPeaks/maximalPeakMass));

	float grassLevelIntensity=NEG_INF;
	
	// calculate grass level peaks
	if (1)
	{
		int i;
		vector<float> peakIntensities;
		peakIntensities.resize(numPeaks);
		for (i=0; i<numPeaks; i++)
			peakIntensities[i]=peaks[i].intensity;

		sort(peakIntensities.begin(),peakIntensities.end());

		int indexForGrass = numPeaks/3;
		grassLevelIntensity = peakIntensities[indexForGrass];

		float intensityUpto2G=0;
		float intensityOf2G=2.0 * grassLevelIntensity;
		int index2G = indexForGrass;
		while (index2G<numPeaks && peakIntensities[index2G]<intensityOf2G)
		{
			intensityUpto2G+=peakIntensities[index2G];
			index2G++;
		}

		float intensityUpto5G=0;
		float intensityOf5G = 5.0 * grassLevelIntensity;
		int index5G = index2G;
		while (index5G<numPeaks && peakIntensities[index5G]<intensityOf5G)
		{
			intensityUpto5G+=peakIntensities[index5G];
			index5G++;
		}

		float intensityUpto10G = 10.0 * grassLevelIntensity;
		int index10G = index5G;
		while (index10G<numPeaks && peakIntensities[index10G]<intensityUpto10G)
			index10G++;

		sample.f_vals.push_back(fval(SQS_PROP_UPTO2G,(float)index2G*oneOverNumberOfPeaks));
		sample.f_vals.push_back(fval(SQS_PROP_UPTO5G,(float)(index5G-index2G)*oneOverNumberOfPeaks));
		sample.f_vals.push_back(fval(SQS_PROP_UPTO10G,(float)(index10G-index5G)*oneOverNumberOfPeaks)); 
		sample.f_vals.push_back(fval(SQS_PROP_MORE10G,(float)(numPeaks-index10G)*oneOverNumberOfPeaks));
		sample.f_vals.push_back(fval(SQS_PROP_INTEN_UPTO2G,intensityUpto2G*oneOverTotalIntensity));
		sample.f_vals.push_back(fval(SQS_PROP_INTEN_UPTO5G,intensityUpto5G*oneOverTotalIntensity));
		sample.f_vals.push_back(fval(SQS_PROP_INTEN_MORE5G, (totalIntensity-intensityUpto2G-intensityUpto5G)*oneOverTotalIntensity));
	}

	// isotope features
	if (1)
	{
		int i;
		int numPeaksWithIso=0;
		int numStrongPeaksWithIso=0;
		for (i=1; i<numPeaks; i++)
			if (currentSpectrumIsotopeLevels_[i]>0)
			{
				numPeaksWithIso++;
				if (currentSpectrumIsotopeLevels_[i-1]==0 && currentSpectrumStrongPeakFlags_[i-1])
					numStrongPeaksWithIso++;
			}

		sample.f_vals.push_back(fval(SQS_PROP_ISO_PEAKS,numPeaksWithIso*oneOverNumberOfPeaks)); 
		sample.f_vals.push_back(fval(SQS_PROP_STRONG_WITH_ISO_PEAKS,numStrongPeaksWithIso*oneOverNumberOfStrongPeaks)); 

	//	cout << "PROP WITH ISO  : " << numPeaksWithIso/(float)numPeaks << endl;
	//	cout << "STRONG WITH IOS: " << numStrongPeaksWithIso/(float)numStrongPeaks << endl;
	}


	// neutral loss features
	if (1)
	{
		vector<float> tmpVals;
		tmpVals.clear();
		int charge;
		for (charge=1; charge<=2; charge++)
		{
			const mass_t offsets[3]={MASS_H2O/charge, MASS_NH3/charge, MASS_CO/charge};
			const int numOffsets = 3;
			
			const int indexDelta = (charge-1)*numOffsets * 2;
			int i;
			for (i=0; i<numOffsets; i++)
			{
				const mass_t minimalOffset = offsets[i]-fragmentTolerance;
				const mass_t maximalOffset = offsets[i]+fragmentTolerance;

				int numPairs=0;
				int numStrongPairs=0;

				int trailIndex=0;
				int leadIndex=1;

				while (leadIndex<numPeaks)
				{
					while (peaks[leadIndex].mass-peaks[trailIndex].mass>maximalOffset)
						trailIndex++;

					if (peaks[leadIndex].mass-peaks[trailIndex].mass>minimalOffset)
					{
						numPairs++;
						if (currentSpectrumStrongPeakFlags_[leadIndex])
							numStrongPairs++;
					}
					leadIndex++;
				}

				const float propAllPeaksWithLoss = numPairs*oneOverNumberOfPeaks;
				const float propStrongPeaksWithLoss = numStrongPairs*oneOverNumberOfStrongPeaks;
				sample.f_vals.push_back(fval(SQS_PROP_ALL_WITH_H2O_LOSS+indexDelta+i,propAllPeaksWithLoss));
				sample.f_vals.push_back(fval(SQS_PROP_STRONG_WITH_H2O_LOSS+indexDelta+i,propStrongPeaksWithLoss));

				tmpVals.push_back(propAllPeaksWithLoss);
				tmpVals.push_back(propStrongPeaksWithLoss);

			//	cout << "OFF REG " << i << " " << numPairs/(float)numPeaks << endl;
			//	cout << "OFF STR " << i << " " << numStrongPairs/(float)numStrongPeaks << endl;
			}
		}

		int i;
		const int half_size = tmpVals.size()/2; 
		for (i=0; i<half_size; i++)
			sample.f_vals.push_back(fval(SQS_DIFF_ALL_WITH_H2O_LOSS+i,tmpVals[i]-tmpVals[half_size+i]));

	//	currentSpectrumStrongPeakFlags_
		
		

	//	
		const mass_t halfMaxMass = maximalPeakMass*0.5;
		const mass_t smallerTolerance = tolerance * 0.66;
		
		int numPairs=0;
		int numStrongPairs=0;
		float pairIntensity=0;
		int doublePeakIndex=0;
		for (i=0; i<numPeaks; i++)
		{
			const mass_t peakMass = peaks[i].mass;
			if (peakMass>halfMaxMass)
				break;

			const mass_t doublePeakMass = peakMass * 2 - MASS_PROTON;
			const mass_t minimalMass = doublePeakMass - smallerTolerance;
			const mass_t maximalMass = doublePeakMass + smallerTolerance;
			while (doublePeakIndex<numPeaks && peaks[doublePeakIndex].mass<minimalMass)
				doublePeakIndex++;

			if (doublePeakIndex==numPeaks)
				break;
			
			if (peaks[doublePeakIndex].mass<maximalMass)
			{
				numPairs++;
				if (currentSpectrumStrongPeakFlags_[i] || currentSpectrumStrongPeakFlags_[doublePeakIndex])
					numStrongPairs++;
				pairIntensity+= peaks[i].intensity + peaks[doublePeakIndex].intensity;
			}
		}

		sample.f_vals.push_back(fval(SQS_PROP_PEAKS_WITH_C1C2,numPairs*oneOverNumberOfPeaks));
		sample.f_vals.push_back(fval(SQS_PROP_STRONG_PEAKS_WITH_C1C2,numStrongPairs*oneOverNumberOfStrongPeaks));
		sample.f_vals.push_back(fval(SQS_PROP_INTEN_WITH_C1C2,pairIntensity*oneOverTotalIntensity));
	}

	// tag features
	if (1)
	{
		vector<DPColumn> dp;

		vector<float> tmpVals;
		int charge;
		for (charge=1; charge<=2; charge++)
		{

			fillSqsDynamicProgrammingTable(pl, dp, charge);
			float intensityInTags[3]={0,0,0};
			int numPeaksInTags[3]={0,0,0};
			int numStrongPeaksInTags[3]={0,0,0};

			vector<int> tagLengths;
			tagLengths.resize(numPeaks,0);

			int longest=0;
			int i;
			for (i=1; i<numPeaks; i++)
			{
				int maxTagLength=-1;
				int aa;
				for (aa=Ala; aa<=Val; aa++)
				{
					const int prev = dp[i].pointers[aa];
					if (prev>=0 && tagLengths[prev]>maxTagLength)
						maxTagLength = tagLengths[prev];
				}

				tagLengths[i]=maxTagLength+1;
				if (tagLengths[i]>longest)
					longest = tagLengths[i];

				int k;
				for (k=0; k<3; k++)
				{
					if (tagLengths[i]==k)
						break;
					
					intensityInTags[k]+=peaks[i].intensity;
					numPeaksInTags[k]++;
					if (currentSpectrumStrongPeakFlags_[i])
						numStrongPeaksInTags[k]++;
				}		
			}

			
			const int idx_off = (SQS_C2_IND_MAX_TAG_LENGTH_ABOVE_4 - SQS_IND_MAX_TAG_LENGTH_ABOVE_4)*(charge-1);
			if (longest>=4)
			{
				sample.f_vals.push_back(fval(SQS_IND_MAX_TAG_LENGTH_ABOVE_4+idx_off,1.0));
				sample.f_vals.push_back(fval(SQS_MAX_TAG_LENGTH_ABOVE_4+idx_off,(float)longest-4.0));
			}
			else
			{
				sample.f_vals.push_back(fval(SQS_IND_MAX_TAG_LENGTH_BELOW_4+idx_off,1.0));
				sample.f_vals.push_back(fval(SQS_MAX_TAG_LENGTH_BELOW_4+idx_off,(float)longest));
			}

	

			float inten_in_tags_both_sides=0;
			for (i=0; i<numPeaks; i++)
				if (dp[i].pointers[0])
					inten_in_tags_both_sides+=peaks[i].intensity;

			sample.f_vals.push_back(fval(SQS_PROP_INTEN_IN_TAGS+idx_off,inten_in_tags_both_sides*oneOverTotalIntensity));

			tmpVals.push_back(longest);
			tmpVals.push_back(inten_in_tags_both_sides*oneOverTotalIntensity);
		//	cout << "MAX TAG: " << longest << endl;
		//	cout << "Inten in TAG: " << inten_in_tags_both_sides/totalIntensity << endl;

			const float strong_threshes[]={0.3,0.2,0.1};
			int k;
			for (k=0; k<3; k++)
			{
				const int pos_off = 4*k;
				const float prop_tags = (float)numPeaksInTags[k]*oneOverNumberOfPeaks;
				const float propStrongPeaksWithLoss = (float)numStrongPeaksInTags[k]*oneOverNumberOfStrongPeaks;
				const float prop_inten  = (float)intensityInTags[k]*oneOverTotalIntensity;
				sample.f_vals.push_back(fval(SQS_PROP_TAGS1+idx_off+pos_off,prop_tags));
				sample.f_vals.push_back(fval(SQS_PROP_STRONG_PEAKS_IN_TAG1+idx_off+pos_off,propStrongPeaksWithLoss));
				sample.f_vals.push_back(fval(SQS_PROP_INTEN_TAG1+idx_off+pos_off,prop_inten));

				if (propStrongPeaksWithLoss<strong_threshes[k])
					sample.f_vals.push_back(fval(SQS_IND_PROP_STRONG_BELOW30_TAG1+idx_off+pos_off,1.0));

				// save vals for diff features
				tmpVals.push_back(prop_tags);
				tmpVals.push_back(propStrongPeaksWithLoss);
				tmpVals.push_back(prop_inten);

			}
		}

		int pos_off = tmpVals.size()/2;
		int i;
		for (i=0; i<pos_off; i++)
			sample.f_vals.push_back(fval(SQS_DIFF_MAX_TAG_LENGTH+i,tmpVals[i]-tmpVals[pos_off+i]));
		/*
		SQS_DIFF_MAX_TAG_LENGTH, SQS_DIFF_PROP_INTEN_IN_TAGS,
		SQS_DIFF_PROP_TAGS1, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG1, SQS_DIFF_PROP_INTEN_TAG1,
		SQS_DIFF_PROP_TAGS2, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG2, SQS_DIFF_PROP_INTEN_TAG2,
		SQS_DIFF_PROP_TAGS3, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG3, SQS_DIFF_PROP_INTEN_TAG3,
		*/
	}


	// density features
	if (1)
	{
		mass_t t1=mOverZ*0.6666;
		mass_t t2=mOverZ*1.3333;
		mass_t h1=mOverZ;

		float intensityT1=0,intensityT2=0, intensityH1=0;
		int	  numPeaksT1=0, numPeaksT2=0, numPeaksH1=0;

		int i;
		for (i=0; i<numPeaks && peaks[i].mass<t1; i++)
			intensityT1+=peaks[i].intensity;
		
		numPeaksT1=i;

		for ( ; i<numPeaks && peaks[i].mass<h1; i++)
			intensityH1+=peaks[i].intensity;

		numPeaksH1=i;

		intensityT2=intensityH1;

		for ( ; i<numPeaks && peaks[i].mass<t2; i++)
			intensityT2+=peaks[i].intensity;

		numPeaksT2=i-numPeaksT1;

		int numPeaksT3= numPeaks - numPeaksT2 - numPeaksT1;
		float intensityT3 = totalIntensity - intensityT2 - intensityT1;

	
		sample.f_vals.push_back(fval(SQS_PEAK_DENSE_T1,(float)(numPeaksT1)*oneOverNumberOfPeaks));
		sample.f_vals.push_back(fval(SQS_PEAK_DENSE_T2,(float)(numPeaksT2)*oneOverNumberOfPeaks));
		sample.f_vals.push_back(fval(SQS_PEAK_DENSE_T3,(float)(numPeaksT3)*oneOverNumberOfPeaks));

	//	cout << "T1 PEAKS: " << (float)(numPeaksT1)/(float)num_peaks << endl;
	//	cout << "T2 PEAKS: " << (float)(numPeaksT2)/(float)num_peaks << endl;
	//	cout << "T3 PEAKS: " << (float)(numPeaksT3)/(float)num_peaks << endl;


		sample.f_vals.push_back(fval(SQS_INTEN_DENSE_T1,intensityT1*oneOverTotalIntensity));
		sample.f_vals.push_back(fval(SQS_INTEN_DENSE_T2,intensityT2*oneOverTotalIntensity));
		sample.f_vals.push_back(fval(SQS_INTEN_DENSE_T3,intensityT3*oneOverTotalIntensity));

	//	cout << "T1 INTEN: " << intensityT1/totalIntensity << endl;
	//	cout << "T2 INTEN: " << intensityT2/totalIntensity << endl;
	//	cout << "T3 INTEN: " << intensityT3/totalIntensity << endl;


		int numPeaksH2 = numPeaks - numPeaksH1;
		float intensityH2 = totalIntensity - intensityH1;

		sample.f_vals.push_back(fval(SQS_PEAK_DENSE_H1,(numPeaksH1)*oneOverNumberOfPeaks));
		sample.f_vals.push_back(fval(SQS_PEAK_DENSE_H2,(numPeaksH2)*oneOverNumberOfPeaks));

		sample.f_vals.push_back(fval(SQS_INTEN_DENSE_H1,intensityH1*oneOverTotalIntensity));
		sample.f_vals.push_back(fval(SQS_INTEN_DENSE_H2,intensityH2*oneOverTotalIntensity));
	
		const float inten33=0.333* totalIntensity;
		const float inten50=0.5  * totalIntensity;
		const float inten75=0.75 * totalIntensity;
		const float inten90=0.90 * totalIntensity;

		float cumulativeIntensity=0;
		const mass_t oneOverMz = 1.0 / (mOverZ + 0.1);

		for (i=0; i<numPeaks && cumulativeIntensity<inten33; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_33_INTEN,(peaks[i].mass*oneOverMz)-0.333));

		for ( ; i<numPeaks && cumulativeIntensity<inten50; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_50_INTEN,(peaks[i].mass*oneOverMz)-0.5));

		for ( ; i<numPeaks && cumulativeIntensity<inten75; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_75_INTEN,(peaks[i].mass*oneOverMz)-0.75));

		for ( ; i<numPeaks && cumulativeIntensity<inten90; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_90_INTEN,(peaks[i].mass*oneOverMz)-0.9));		
	}

	// pmc features
	if (1)
	{
		static vector< vector<float> > sqsFeatures;
		getSqsFeaturesFromPmcTabels(pl, sqsFeatures);
	//	get_sqs_features_from_pmc_tables(bs,sqsFeatures);

		vector<float> maximalValues;
		maximalValues.resize(4,0.000001);

		vector< vector<float> > tmpVals;
		tmpVals.resize(4);

		// first give absolute counts
		int charge;
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqsFeatures[charge];
			const int idx_off = 4 *( charge-1);

			int i;
			for (i=0; i<4; i++)
			{
				sample.f_vals.push_back(fval(SQS_NUM_FRAG_PAIRS_1+idx_off+i,vals[i]));	
				tmpVals[charge].push_back(vals[i]);
				if (vals[i]>maximalValues[i])
					maximalValues[i]=vals[i];
			}	
		}

		// add features for prop of max
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqsFeatures[charge];
			const int idx_off = 4 *( charge-1);

			int i;
			for (i=0; i<4; i++)
			{
				const float ratio = (maximalValues[i]>0 ? vals[i]/maximalValues[i] : 0);
				sample.f_vals.push_back(fval(SQS_PROP_OF_MAX_FRAG_PAIRS_1+idx_off+i,ratio));			
				tmpVals[charge].push_back(ratio);
			}
		}

		// conver to proportions by dividing by the number of peaks/strong_peaks
		// and subtract the background levels (first)
		const float one_over_peaks = 1.0/((float)numPeaks+0.1);
		const float one_over_strong = 1.0/((float)numStrongPeaks+0.1);

		
	
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqsFeatures[charge];
			const int idx_off = 4 *( charge-1);

			// normalize to get proportions of total number of peaks/strong
			vals[0]-= this->currentSpectrumBackgroundStats_[charge].numFragmentPairs;
			vals[1]-= this->currentSpectrumBackgroundStats_[charge].numStrongFragmentPairs;
			vals[2]-= this->currentSpectrumBackgroundStats_[charge].numCharge2FragmentPairs;
			vals[3]-= this->currentSpectrumBackgroundStats_[charge].numStrongCharge2FragmentPairs;

			vals[0]*=one_over_peaks;
			vals[1]*=one_over_strong;
			vals[2]*=one_over_peaks;
			vals[3]*=one_over_strong;

			int i;
			for (i=0; i<4; i++)
			{
				sample.f_vals.push_back(fval(SQS_PROP_FRAG_PAIRS_1+idx_off+i,vals[i]));
				tmpVals[charge].push_back(vals[i]);
			}
		}

		// diff is between values for c2 and c3
		int i;
		for (i=0; i<tmpVals[2].size(); i++)
			sample.f_vals.push_back(fval(SQS_DIFF_NUM_FRAG_PAIRS_23+i,tmpVals[2][i]-tmpVals[3][i]));
	}

	
	sort(sample.f_vals.begin(),sample.f_vals.end());
}



/**************************************************************************
Fills a sample with the SQS features (listed in sqsFeatures.txt)
***************************************************************************
void PMCSQS_Scorer::fillSqsSample(const PeakList& pl, Sample& sample) const
{
	const mass_t tolerance = config_->getTolerance();
	const mass_t fragmentTolerance = (tolerance<0.15 ? tolerance : tolerance * 0.5);
	const Peak* const peaks  = pl.getPeaks();
	const int  numPeaks		 = pl.getNumPeaks();
	const mass_t mOverZ		 = pl.getHeader()->getMOverZ();
	const int numStrongPeaks = curr_spec_num_strong;
	const float totalIntensity = curr_spec_total_intensity;
	const float oneOverNumberOfPeaks = 1.0 / (numPeaks+1);
	const float oneOverNumberOfStrongPeaks = 1.0/ (numStrongPeaks + 1);
	const float oneOverTotalIntensity = 1.0 / (totalIntensity + 1.0);
	const mass_t maximalPeakMass = (numPeaks <5 ? POS_INF : peaks[numPeaks-1].mass);

	sample.clearPairs();
	sample.addPair(0,1.0); // SQS_CONST
	sample.addPair(1, static_cast<value_t>(numPeaks)/maximalPeakMass); // SQS_PEAK_DENSITY


	float grassLevelIntensity=0;
	
	// calculate grass level peaks and add intensity features
	if (1)
	{
		int i;
		vector<float> peakIntensities;
		peakIntensities.resize(numPeaks);
		for (i=0; i<numPeaks; i++)
			peakIntensities[i]=peaks[i].intensity;

		sort(peakIntensities.begin(),peakIntensities.end());

		size_t indexForGrass = numPeaks/3;
		grassLevelIntensity = peakIntensities[indexForGrass];

		float intensityUpto1G=0.0;
		size_t index=0;
		while (index<numPeaks && peakIntensities[index]<grassLevelIntensity)
			intensityUpto1G+=peakIntensities[index++];

		float intensityUpto2G=0.0;
		float intensityOf2G = 2.0 * grassLevelIntensity;
		size_t index2G = indexForGrass;
		while (index2G<numPeaks && peakIntensities[index2G]<intensityOf2G)
			intensityUpto2G+=peakIntensities[index2G++];
		
		float intensityUpto5G=0.0;
		float intensityOf5G = 5.0 * grassLevelIntensity;
		size_t index5G = index2G;
		while (index5G<numPeaks && peakIntensities[index5G]<intensityOf5G)
			intensityUpto5G+=peakIntensities[index5G++];
		
		float intensityUpto10G=0.0;
		float intensityOf10G = 10.0 * grassLevelIntensity;
		size_t index10G = index5G;
		while (index10G<numPeaks && peakIntensities[index10G]<intensityOf10G)
			intensityUpto10G+=peakIntensities[index10G++];

		sample.addPair(2, static_cast<value_t>(indexForGrass)* oneOverNumberOfPeaks);  // SQS_PROP_GRASS
		sample.addPair(3, static_cast<value_t>(index2G-indexForGrass)* oneOverNumberOfPeaks); // SQS_PROP_PEAKS_1GTO2G
		sample.addPair(4, static_cast<value_t>(index5G-index2G)* oneOverNumberOfPeaks); // SQS_PROP_PEAKS_2GTO5G
		sample.addPair(5, static_cast<value_t>(index10G-index5G)* oneOverNumberOfPeaks); // SQS_PROP_PEAKS_5GTO10G
		sample.addPair(6, static_cast<value_t>(numPeaks-index10G)* oneOverNumberOfPeaks); // SQS_PROP_PEAKS_MORE_10G

		sample.addPair(7, intensityUpto1G * oneOverTotalIntensity); // SQS_PROP_INTENSITY_GRASS
		sample.addPair(8, intensityUpto2G * oneOverTotalIntensity); // SQS_PROP_INTENSITY_2G
		sample.addPair(9, intensityUpto5G * oneOverTotalIntensity); // SQS_PROP_INTENSITY_5G
		sample.addPair(10, intensityUpto10G * oneOverTotalIntensity); // SQS_PROP_INTENSITY_10G
		sample.addPair(11, (totalIntensity-intensityUpto1G-intensityUpto2G-intensityUpto5G-intensityUpto10G) 
					  * oneOverTotalIntensity); // SQS_PROP_INTENSITY_MORE_10G
	}

	// isotope features
	if (1)
	{
		size_t i;
		size_t numPeaksWithIso=0;
		size_t numStrongPeaksWithIso=0;
		for (i=1; i<numPeaks; i++)
			if (currentSpectrumIsotopeLevels_[i]>0)
			{
				numPeaksWithIso++;
				if (currentSpectrumIsotopeLevels_[i-1]==0.0 && currentSpectrumStrongPeakFlags_[i-1])
					numStrongPeaksWithIso++;
			}

		sample.addPair(12, numPeaksWithIso*oneOverNumberOfPeaks); // SQS_PROP_PEAKS_WITH_ISO
		sample.addPair(13, numStrongPeaksWithIso*oneOverNumberOfStrongPeaks); // PROP_STRONG_PEAKS_WITH_ISO
	}


	// neutral loss features
	if (1)
	{
		float propPeaksWithAnyLoss=0;
		float propStrongPeaksWithAnyLoss=0;

		size_t charge;
		for (charge=1; charge<=2; charge++)
		{
			const mass_t offsets[3]={MASS_H2O/charge, MASS_NH3/charge, MASS_CO/charge};
			const size_t numOffsets = sizeof(offsets)/sizeof(mass_t);

			const size_t indexStart = 14+(charge-1)*(numOffsets * 2 + 2);

			size_t i;
			for (i=0; i<numOffsets; i++)
			{
				const mass_t minimalOffset = offsets[i] - fragmentTolerance;
				const mass_t maximalOffset = offsets[i] + fragmentTolerance;

				size_t numPairs=0;
				size_t numStrongPairs=0;

				size_t trailIndex=0;
				size_t leadIndex=1;

				while (leadIndex<numPeaks)
				{
					while (peaks[leadIndex].mass-peaks[trailIndex].mass>maximalOffset)
						trailIndex++;

					if (peaks[leadIndex].mass-peaks[trailIndex].mass>minimalOffset)
					{
						numPairs++;
						if (currentSpectrumStrongPeakFlags_[leadIndex])
							numStrongPairs++;
					}
					leadIndex++;
				}

				const float propAllPeaksWithLoss = numPairs*oneOverNumberOfPeaks;
				const float propStrongPeaksWithLoss = numStrongPairs*oneOverNumberOfStrongPeaks;

				propPeaksWithAnyLoss += propAllPeaksWithLoss; // don't care if same peak gets recounted
				propStrongPeaksWithAnyLoss += propStrongPeaksWithLoss;

				sample.addPair(indexStart + 2*i,  propAllPeaksWithLoss);  // LOSS
				sample.addPair(indexStart + 2*i+1,propStrongPeaksWithLoss);	// STRONG LOSS
			}
			sample.addPair(indexStart + 2*numOffsets,   propPeaksWithAnyLoss); // ANY_LOSS
			sample.addPair(indexStart + 2*numOffsets+1, propStrongPeaksWithAnyLoss); // ANY STRONG LOSS 
		}
	}
		
	//	double/single complement features
	if (1)
	{
		const mass_t halfMaxMass = maximalPeakMass*0.5;
		const mass_t smallerTolerance = tolerance * 0.66;
		
		size_t numPairs=0;
		size_t numStrongPairs=0;
		float pairIntensity=0;
		size_t doublePeakIndex=0;
		size_t i;
		for (i=0; i<numPeaks; i++)
		{
			const mass_t peakMass = peaks[i].mass;
			if (peakMass>halfMaxMass)
				break;

			const mass_t doublePeakMass = peakMass * 2 - MASS_PROTON;
			const mass_t minimalMass = doublePeakMass - smallerTolerance;
			const mass_t maximalMass = doublePeakMass + smallerTolerance;
			while (doublePeakIndex<numPeaks && peaks[doublePeakIndex].mass<minimalMass)
				doublePeakIndex++;

			if (doublePeakIndex==numPeaks)
				break;
			
			if (peaks[doublePeakIndex].mass<maximalMass)
			{
				numPairs++;
				if (currentSpectrumStrongPeakFlags_[i] || currentSpectrumStrongPeakFlags_[doublePeakIndex])
					numStrongPairs++;
				pairIntensity+= peaks[i].intensity + peaks[doublePeakIndex].intensity;
			}
		}

		sample.addPair(30, numPairs*oneOverNumberOfPeaks);	//	SQS_PROP_PEAKS_WITH_C1C2_PAIRS
		sample.addPair(31, numStrongPairs*oneOverNumberOfStrongPeaks); // SQS_PROP_STRONG_WITH_C1C2_PAIRS
		sample.addPair(32, pairIntensity*oneOverTotalIntensity); // SQS_PROP_INTEN_IN_C1C2_PAIRS
	}

	// tag features
	if (1)
	{
		vector<DPColumn> dp;
		vector<float> tmpVals;
		int charge;
		for (charge=1; charge<=2; charge++)
		{
			fillSqsDynamicProgrammingTable(pl, dp, charge);

			vector<size_t> forwardTagLengths;
			vector<size_t> backwardsTagLengths;

			forwardTagLengths.resize(numPeaks,0);
			backwardsTagLengths.resize(numPeaks,0);

			size_t i;
			for (i=1; i<numPeaks; i++)
			{
				int maxTagLength = -1;
				int aa;
				for (aa=Ala; aa<=Val; aa++)
				{
					const int prev = dp[i].pointers[aa];
					if (prev>=0 && backwardsTagLengths[prev]>maxTagLength)
						maxTagLength = backwardsTagLengths[prev];
				}
				backwardsTagLengths[i]=static_cast<size_t>(maxTagLength+1);
			}

			for (i=numPeaks-1;  ; i--)
			{
				int maxTagLength = -1;
				int aa;
				for (aa=Ala; aa<=Val; aa++)
				{
					const int prev = dp[i].fPointers[aa];
					if (prev>=0 && forwardTagLengths[prev]>maxTagLength)
						maxTagLength = forwardTagLengths[prev];
				}
				forwardTagLengths[i]=static_cast<size_t>(maxTagLength+1);
				if (i == 0)
					break;
			}

			vector<size_t> countsPerLength(6,0);
			vector<size_t> strongCountsPerLength(6,0);
			size_t maxTagLength=0;

			for (i=0; i<numPeaks; i++)
			{
				const size_t len = backwardsTagLengths[i]+forwardTagLengths[i];
				if (len>maxTagLength)
					maxTagLength=len;

				if (len>=5)
				{
					countsPerLength[5]++;
					if (currentSpectrumStrongPeakFlags_[i])
						strongCountsPerLength[5]++;
				}
				else
				{
					countsPerLength[len]++;
					if (currentSpectrumStrongPeakFlags_[i])
						strongCountsPerLength[len]++;
				}
			}

			const size_t startIndex = (charge == 1 ? 33 : 42);
			sample.addPair(startIndex, maxTagLength);
			sample.addPair(startIndex+1, countsPerLength[0]*oneOverNumberOfPeaks);
			sample.addPair(startIndex+2, strongCountsPerLength[0]*oneOverNumberOfStrongPeaks);
			sample.addPair(startIndex+3, (countsPerLength[1]+countsPerLength[2])*oneOverNumberOfPeaks);
			sample.addPair(startIndex+4, (strongCountsPerLength[1]+strongCountsPerLength[2])*oneOverNumberOfStrongPeaks);
			sample.addPair(startIndex+5, (countsPerLength[3]+countsPerLength[4])*oneOverNumberOfPeaks);
			sample.addPair(startIndex+6, (strongCountsPerLength[3]+strongCountsPerLength[4])*oneOverNumberOfStrongPeaks);
			sample.addPair(startIndex+7, countsPerLength[5]*oneOverNumberOfPeaks);
			sample.addPair(startIndex+8, strongCountsPerLength[5]*oneOverNumberOfStrongPeaks);
		}
	}


	// regional density features
	if (1)
	{
		mass_t t1=mOverZ*0.6666;
		mass_t t2=mOverZ*1.3333;
		mass_t h1=mOverZ;

		float intensityT1=0,intensityT2=0, intensityH1=0;
		size_t numPeaksT1=0, numPeaksT2=0, numPeaksH1=0;

		size_t i;
		for (i=0; i<numPeaks && peaks[i].mass<t1; i++)
			intensityT1+=peaks[i].intensity;
		
		numPeaksT1=i;

		for ( ; i<numPeaks && peaks[i].mass<h1; i++)
			intensityH1+=peaks[i].intensity;

		numPeaksH1=i;

		intensityT2=intensityH1;

		for ( ; i<numPeaks && peaks[i].mass<t2; i++)
			intensityT2+=peaks[i].intensity;

		numPeaksT2=i-numPeaksT1;

		const size_t numPeaksT3 = numPeaks - numPeaksT2 - numPeaksT1;
		const float intensityT3 = totalIntensity - intensityT2 - intensityT1;

		sample.addPair(51, static_cast<value_t>(numPeaksT1) * oneOverNumberOfPeaks);
		sample.addPair(52, static_cast<value_t>(numPeaksT2) * oneOverNumberOfPeaks);
		sample.addPair(53, static_cast<value_t>(numPeaksT3) * oneOverNumberOfPeaks);

		sample.addPair(54,	intensityT1 * oneOverTotalIntensity);
		sample.addPair(55,	intensityT2 * oneOverTotalIntensity);
		sample.addPair(56,	intensityT3 * oneOverTotalIntensity);

		const size_t numPeaksH2 = numPeaks - numPeaksH1;
		const float intensityH2 = totalIntensity - intensityH1;

		sample.addPair(57, static_cast<value_t>(numPeaksH1) * oneOverNumberOfPeaks);
		sample.addPair(58, static_cast<value_t>(numPeaksH2) * oneOverNumberOfPeaks);
		sample.addPair(59,	intensityH1 * oneOverTotalIntensity);
		sample.addPair(60,	intensityH2 * oneOverTotalIntensity);
	
		const float inten33=0.333* totalIntensity;
		const float inten50=0.5  * totalIntensity;
		const float inten75=0.75 * totalIntensity;
		const float inten90=0.90 * totalIntensity;

		float cumulativeIntensity=0;
		const mass_t oneOverMz = 1.0 / (mOverZ + 0.1);

		for (i=0; i<numPeaks && cumulativeIntensity<inten33; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.addPair(61, (peaks[i].mass*oneOverMz)-0.333 ); // SQS_PROP_MZ_RANGE_WITH_33_INTEN

		for ( ; i<numPeaks && cumulativeIntensity<inten50; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.addPair(62, (peaks[i].mass*oneOverMz)-0.5 ); //  SQS_PROP_MZ_RANGE_WITH_50_INTEN


		for ( ; i<numPeaks && cumulativeIntensity<inten75; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.addPair(63, (peaks[i].mass*oneOverMz)-0.75 ); //  SQS_PROP_MZ_RANGE_WITH_75_INTEN

		for ( ; i<numPeaks && cumulativeIntensity<inten90; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.addPair(64, (peaks[i].mass*oneOverMz)-0.9 ); //  SQS_PROP_MZ_RANGE_WITH_90_INTEN	
	}

	// pmc features
	if (1 && numPeaks>0 && numStrongPeaks>0)
	{
		vector< vector<value_t> > sqsFeatures;
		getSqsFeaturesFromPmcTabels(pl, sqsFeatures);
	
		vector<size_t> maxValueCharges(4,1);
		
		// find what charges give the best fragment pair counts
		size_t charge;
		for (charge=2; charge<=maximalChargeWithModels_; charge++)
		{
			size_t i;
			for (i=0; i<4; i++)
				if (sqsFeatures[charge][i]>sqsFeatures[maxValueCharges[i]][i])
					maxValueCharges[i]=charge;
		}

		// find charges that give the second best counts
		vector<size_t> secondBestCharges(4,0);
		for (charge=1; charge<=maximalChargeWithModels_; charge++)
		{
			size_t i;
			for (i=0; i<4; i++)
			{
				if (charge == maxValueCharges[i])
					continue;

				if (secondBestCharges[i]==0 || 
					sqsFeatures[charge][i] > sqsFeatures[secondBestCharges[i]][i])
					secondBestCharges[i]=charge;
			}
		}

		vector<value_t> maxVals(4), secondVals(4);
		size_t i;
		for (i=0; i<4; i++)
		{
			maxVals[i]=sqsFeatures[maxValueCharges[i]][i];
			secondVals[i]=sqsFeatures[secondBestCharges[i]][i];
		}

		sample.addPair(65, maxVals[0]); // SQS_CHARGE_MAX_PAIRS	
		sample.addPair(66, maxVals[1]); // SQS_CHARGE_MAX_STRONG_PAIRS
		sample.addPair(67, maxVals[2]); // SQS_CHARGE_MAX_C2_PAIRS
		sample.addPair(68, maxVals[3]); // SQS_CHARGE_MAX_STRONG_C2_PAIRS

		sample.addPair(69, secondVals[0]); // SQS_2ND_CHARGE_MAX_PAIRS
		sample.addPair(70, secondVals[1]); // SQS_2ND_CHARGE_MAX_STRONG_PAIRS
		sample.addPair(71, secondVals[2]); // SQS_2ND_CHARGE_MAX_C2_PAIRS
		sample.addPair(72, secondVals[3]); // SQS_2ND_CHARGE_MAX_STRONG_C2_PAIRS


		// the following features measure the difference in pair counts
		// the vlaues are first normalized to account for peak frequency issues
		const value_t oneOverNumberOfPeaks       = 1.0/static_cast<value_t>(numPeaks);
		const value_t oneOverNumberOfStrongPeaks = 1.0/static_cast<value_t>(numStrongPeaks);

		sample.addPair(73, (maxVals[0]-secondVals[0]) * oneOverNumberOfPeaks); // SQS_DIFF_MAX_PAIRS
		sample.addPair(74, (maxVals[1]-secondVals[1]) * oneOverNumberOfStrongPeaks); // SQS_DIFF_MAX_STRONG_PAIRS
		sample.addPair(75, (maxVals[2]-secondVals[2]) * oneOverNumberOfPeaks); // SQS_DIFF_MAX_C2_PAIRS
		sample.addPair(76, (maxVals[3]-secondVals[3]) * oneOverNumberOfStrongPeaks); // SQS_DIFF_MAX_STRONG_C2_PAIRS

		// add features that compare the counts to the background count of pairs (i.e., what we'd
		// epect to find at random). These values are also normalized

		// SQS_DIFF_BK_PROP_PAIRS
		sample.addPair(77, oneOverNumberOfPeaks * 
			(maxVals[0] - currentSpectrumBackgroundStats_[maxValueCharges[0]].numFragmentPairs));
		
		// SQS_DIFF_BK_PROP_STRONG_PAIRS
		sample.addPair(78, oneOverNumberOfStrongPeaks *
			(maxVals[1] - currentSpectrumBackgroundStats_[maxValueCharges[1]].numStrongFragmentPairs));

		// SQS_DIFF_BK_PROP_C2_PAIRS
		sample.addPair(79, oneOverNumberOfPeaks * 
			(maxVals[2] - currentSpectrumBackgroundStats_[maxValueCharges[2]].numCharge2FragmentPairs));

		// SQS_DIFF_BK_PROP_C2_STRONG_PAIRS
		sample.addPair(80, oneOverNumberOfStrongPeaks *
			(maxVals[3] - currentSpectrumBackgroundStats_[maxValueCharges[3]].numStrongCharge2FragmentPairs));
	}

	#ifdef MY_DEBUG_CHECKS
	/*	if (! sample.checkConsistency())
		{
			sample.print();
			error("Index mismatch in sample!");
		}
	#endif
}*/


/****************************************************************************
*****************************************************************************/
void PMCSQS_Scorer::fill_fval_vector_with_SQS(const BasicSpectrum& bs,
								   ME_Regression_Sample& sam) const
{
	const mass_t tolerance = config_->getTolerance();
	const mass_t fragmentTolerance = (tolerance<0.15 ? tolerance : tolerance * 0.5);
	const QCPeak *peaks  = bs.peaks;
	const int  num_peaks = bs.num_peaks;
	const mass_t m_over_z = bs.ssf->m_over_z;
	const int numStrongPeaks=curr_spec_num_strong;
	const float totalIntensity=curr_spec_total_intensity;
	const float oneOverNumberOfPeaks = 10.0 / (num_peaks+1);
	const float oneOverNumberOfStrongPeaks = 5.0/ (numStrongPeaks + 1);
	const float oneOverTotalIntensity = 1.0 / (totalIntensity + 1.0);
	const mass_t maximalPeakMass = (num_peaks <5 ? POS_INF : peaks[num_peaks-1].mass);

	sam.f_vals.clear();
	sam.f_vals.push_back(fval(SQS_CONST,1.0));
	sam.f_vals.push_back(fval(SQS_PEAK_DENSITY, (float)num_peaks/maximalPeakMass));

	float grassLevelIntensity=NEG_INF;
	
	// calculate grass level peaks
	if (1)
	{
		int i;
		vector<float> peakIntensities;
		peakIntensities.resize(num_peaks);
		for (i=0; i<num_peaks; i++)
			peakIntensities[i]=peaks[i].intensity;

		sort(peakIntensities.begin(),peakIntensities.end());

		int indexForGrass = num_peaks/3;
		grassLevelIntensity = peakIntensities[indexForGrass];

		float intensityUpto2G=0;
		float intensityOf2G=2.0 * grassLevelIntensity;
		int index2G = indexForGrass;
		while (index2G<num_peaks && peakIntensities[index2G]<intensityOf2G)
		{
			intensityUpto2G+=peakIntensities[index2G];
			index2G++;
		}

		float intensityUpto5G=0;
		float intensityOf5G = 5.0 * grassLevelIntensity;
		int index5G = index2G;
		while (index5G<num_peaks && peakIntensities[index5G]<intensityOf5G)
		{
			intensityUpto5G+=peakIntensities[index5G];
			index5G++;
		}

		float intensityUpto10G = 10.0 * grassLevelIntensity;
		int index10G = index5G;
		while (index10G<num_peaks && peakIntensities[index10G]<intensityUpto10G)
			index10G++;

		sam.f_vals.push_back(fval(SQS_PROP_UPTO2G,(float)index2G*oneOverNumberOfPeaks));
		sam.f_vals.push_back(fval(SQS_PROP_UPTO5G,(float)(index5G-index2G)*oneOverNumberOfPeaks));
		sam.f_vals.push_back(fval(SQS_PROP_UPTO10G,(float)(index10G-index5G)*oneOverNumberOfPeaks)); 
		sam.f_vals.push_back(fval(SQS_PROP_MORE10G,(float)(num_peaks-index10G)*oneOverNumberOfPeaks));
		sam.f_vals.push_back(fval(SQS_PROP_INTEN_UPTO2G,intensityUpto2G*oneOverTotalIntensity));
		sam.f_vals.push_back(fval(SQS_PROP_INTEN_UPTO5G,intensityUpto5G*oneOverTotalIntensity));
		sam.f_vals.push_back(fval(SQS_PROP_INTEN_MORE5G, (totalIntensity-intensityUpto2G-intensityUpto5G)*oneOverTotalIntensity));
	}

	// isotope features
	if (1)
	{
		int i;
		int numPeaksWithIso=0;
		int numStrongPeaksWithIso=0;
		for (i=1; i<num_peaks; i++)
			if (currentSpectrumIsotopeLevels_[i]>0)
			{
				numPeaksWithIso++;
				if (currentSpectrumIsotopeLevels_[i-1]==0 && currentSpectrumStrongPeakFlags_[i-1])
					numStrongPeaksWithIso++;
			}

		sam.f_vals.push_back(fval(SQS_PROP_ISO_PEAKS,numPeaksWithIso*oneOverNumberOfPeaks)); 
		sam.f_vals.push_back(fval(SQS_PROP_STRONG_WITH_ISO_PEAKS,numStrongPeaksWithIso*oneOverNumberOfStrongPeaks)); 

	//	cout << "PROP WITH ISO  : " << numPeaksWithIso/(float)num_peaks << endl;
	//	cout << "STRONG WITH IOS: " << numStrongPeaksWithIso/(float)numStrongPeaks << endl;
	}


	// neutral loss features
	if (1)
	{
		vector<float> tmpVals;
		tmpVals.clear();
		int charge;
		for (charge=1; charge<=2; charge++)
		{
			const mass_t offsets[3]={MASS_H2O/charge, MASS_NH3/charge, MASS_CO/charge};
			const int numOffsets = 3;
			
			const int indexDelta = (charge-1)*numOffsets * 2;
			int i;
			for (i=0; i<numOffsets; i++)
			{
				const mass_t minimalOffset = offsets[i]-fragmentTolerance;
				const mass_t maximalOffset = offsets[i]+fragmentTolerance;

				int numPairs=0;
				int numStrongPairs=0;

				int trailIndex=0;
				int leadIndex=1;

				while (leadIndex<num_peaks)
				{
					while (peaks[leadIndex].mass-peaks[trailIndex].mass>maximalOffset)
						trailIndex++;

					if (peaks[leadIndex].mass-peaks[trailIndex].mass>minimalOffset)
					{
						numPairs++;
						if (currentSpectrumStrongPeakFlags_[leadIndex])
							numStrongPairs++;
					}
					leadIndex++;
				}

				const float propAllPeaksWithLoss = numPairs*oneOverNumberOfPeaks;
				const float propStrongPeaksWithLoss = numStrongPairs*oneOverNumberOfStrongPeaks;
				sam.f_vals.push_back(fval(SQS_PROP_ALL_WITH_H2O_LOSS+indexDelta+i,propAllPeaksWithLoss));
				sam.f_vals.push_back(fval(SQS_PROP_STRONG_WITH_H2O_LOSS+indexDelta+i,propStrongPeaksWithLoss));

				tmpVals.push_back(propAllPeaksWithLoss);
				tmpVals.push_back(propStrongPeaksWithLoss);

			//	cout << "OFF REG " << i << " " << numPairs/(float)num_peaks << endl;
			//	cout << "OFF STR " << i << " " << numStrongPairs/(float)numStrongPeaks << endl;
			}
		}

		int i;
		const int half_size = tmpVals.size()/2; 
		for (i=0; i<half_size; i++)
			sam.f_vals.push_back(fval(SQS_DIFF_ALL_WITH_H2O_LOSS+i,tmpVals[i]-tmpVals[half_size+i]));

	//	currentSpectrumStrongPeakFlags_
		
		

	//	
		const mass_t halfMaxMass = maximalPeakMass*0.5;
		const mass_t smallerTolerance = tolerance * 0.66;
		
		int numPairs=0;
		int numStrongPairs=0;
		float pairIntensity=0;
		int doublePeakIndex=0;
		for (i=0; i<num_peaks; i++)
		{
			const mass_t peakMass = peaks[i].mass;
			if (peakMass>halfMaxMass)
				break;

			const mass_t doublePeakMass = peakMass * 2 - MASS_PROTON;
			const mass_t minimalMass = doublePeakMass - smallerTolerance;
			const mass_t maximalMass = doublePeakMass + smallerTolerance;
			while (doublePeakIndex<num_peaks && peaks[doublePeakIndex].mass<minimalMass)
				doublePeakIndex++;

			if (doublePeakIndex==num_peaks)
				break;
			
			if (peaks[doublePeakIndex].mass<maximalMass)
			{
				numPairs++;
				if (currentSpectrumStrongPeakFlags_[i] || currentSpectrumStrongPeakFlags_[doublePeakIndex])
					numStrongPairs++;
				pairIntensity+= peaks[i].intensity + peaks[doublePeakIndex].intensity;
			}
		}

		sam.f_vals.push_back(fval(SQS_PROP_PEAKS_WITH_C1C2,numPairs*oneOverNumberOfPeaks));
		sam.f_vals.push_back(fval(SQS_PROP_STRONG_PEAKS_WITH_C1C2,numStrongPairs*oneOverNumberOfStrongPeaks));
		sam.f_vals.push_back(fval(SQS_PROP_INTEN_WITH_C1C2,pairIntensity*oneOverTotalIntensity));
	}

	// tag features
	if (1)
	{
		vector<DPColumn> dp;

		vector<float> tmpVals;
		int charge;
		for (charge=1; charge<=2; charge++)
		{
			fill_SQS_DP(bs,dp,charge);
			float intensityInTags[3]={0,0,0};
			int numPeaksInTags[3]={0,0,0};
			int numStrongPeaksInTags[3]={0,0,0};

			vector<int> tagLengths;
			tagLengths.resize(num_peaks,0);

			int longest=0;
			int i;
			for (i=1; i<num_peaks; i++)
			{
				int maxTagLength=-1;
				int aa;
				for (aa=Ala; aa<=Val; aa++)
				{
					const int prev = dp[i].pointers[aa];
					if (prev>=0 && tagLengths[prev]>maxTagLength)
						maxTagLength = tagLengths[prev];
				}

				tagLengths[i]=maxTagLength+1;
				if (tagLengths[i]>longest)
					longest = tagLengths[i];

				int k;
				for (k=0; k<3; k++)
				{
					if (tagLengths[i]==k)
						break;
					
					intensityInTags[k]+=peaks[i].intensity;
					numPeaksInTags[k]++;
					if (currentSpectrumStrongPeakFlags_[i])
						numStrongPeaksInTags[k]++;
				}		
			}

			
			const int idx_off = (SQS_C2_IND_MAX_TAG_LENGTH_ABOVE_4 - SQS_IND_MAX_TAG_LENGTH_ABOVE_4)*(charge-1);
			if (longest>=4)
			{
				sam.f_vals.push_back(fval(SQS_IND_MAX_TAG_LENGTH_ABOVE_4+idx_off,1.0));
				sam.f_vals.push_back(fval(SQS_MAX_TAG_LENGTH_ABOVE_4+idx_off,(float)longest-4.0));
			}
			else
			{
				sam.f_vals.push_back(fval(SQS_IND_MAX_TAG_LENGTH_BELOW_4+idx_off,1.0));
				sam.f_vals.push_back(fval(SQS_MAX_TAG_LENGTH_BELOW_4+idx_off,(float)longest));
			}

	

			float inten_in_tags_both_sides=0;
			for (i=0; i<num_peaks; i++)
				if (dp[i].pointers[0])
					inten_in_tags_both_sides+=peaks[i].intensity;

			sam.f_vals.push_back(fval(SQS_PROP_INTEN_IN_TAGS+idx_off,inten_in_tags_both_sides*oneOverTotalIntensity));

			tmpVals.push_back(longest);
			tmpVals.push_back(inten_in_tags_both_sides*oneOverTotalIntensity);
		//	cout << "MAX TAG: " << longest << endl;
		//	cout << "Inten in TAG: " << inten_in_tags_both_sides/totalIntensity << endl;

			const float strong_threshes[]={0.3,0.2,0.1};
			int k;
			for (k=0; k<3; k++)
			{
				const int pos_off = 4*k;
				const float prop_tags = (float)numPeaksInTags[k]*oneOverNumberOfPeaks;
				const float propStrongPeaksWithLoss = (float)numStrongPeaksInTags[k]*oneOverNumberOfStrongPeaks;
				const float prop_inten  = (float)intensityInTags[k]*oneOverTotalIntensity;
				sam.f_vals.push_back(fval(SQS_PROP_TAGS1+idx_off+pos_off,prop_tags));
				sam.f_vals.push_back(fval(SQS_PROP_STRONG_PEAKS_IN_TAG1+idx_off+pos_off,propStrongPeaksWithLoss));
				sam.f_vals.push_back(fval(SQS_PROP_INTEN_TAG1+idx_off+pos_off,prop_inten));

				if (propStrongPeaksWithLoss<strong_threshes[k])
					sam.f_vals.push_back(fval(SQS_IND_PROP_STRONG_BELOW30_TAG1+idx_off+pos_off,1.0));

				// save vals for diff features
				tmpVals.push_back(prop_tags);
				tmpVals.push_back(propStrongPeaksWithLoss);
				tmpVals.push_back(prop_inten);

			}
		}

		int pos_off = tmpVals.size()/2;
		int i;
		for (i=0; i<pos_off; i++)
			sam.f_vals.push_back(fval(SQS_DIFF_MAX_TAG_LENGTH+i,tmpVals[i]-tmpVals[pos_off+i]));
		/*
		SQS_DIFF_MAX_TAG_LENGTH, SQS_DIFF_PROP_INTEN_IN_TAGS,
		SQS_DIFF_PROP_TAGS1, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG1, SQS_DIFF_PROP_INTEN_TAG1,
		SQS_DIFF_PROP_TAGS2, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG2, SQS_DIFF_PROP_INTEN_TAG2,
		SQS_DIFF_PROP_TAGS3, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG3, SQS_DIFF_PROP_INTEN_TAG3,
		*/
	}


	// density features
	if (1)
	{
		mass_t t1=m_over_z*0.6666;
		mass_t t2=m_over_z*1.3333;
		mass_t h1=m_over_z;

		float intensityT1=0,intensityT2=0, intensityH1=0;
		int	  numPeaksT1=0, numPeaksT2=0, numPeaksH1=0;

		int i;
		for (i=0; i<num_peaks && peaks[i].mass<t1; i++)
			intensityT1+=peaks[i].intensity;
		
		numPeaksT1=i;

		for ( ; i<num_peaks && peaks[i].mass<h1; i++)
			intensityH1+=peaks[i].intensity;

		numPeaksH1=i;

		intensityT2=intensityH1;

		for ( ; i<num_peaks && peaks[i].mass<t2; i++)
			intensityT2+=peaks[i].intensity;

		numPeaksT2=i-numPeaksT1;

		int numPeaksT3= num_peaks - numPeaksT2 - numPeaksT1;
		float intensityT3 = totalIntensity - intensityT2 - intensityT1;

	
		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_T1,(float)(numPeaksT1)*oneOverNumberOfPeaks));
		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_T2,(float)(numPeaksT2)*oneOverNumberOfPeaks));
		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_T3,(float)(numPeaksT3)*oneOverNumberOfPeaks));

	//	cout << "T1 PEAKS: " << (float)(numPeaksT1)/(float)num_peaks << endl;
	//	cout << "T2 PEAKS: " << (float)(numPeaksT2)/(float)num_peaks << endl;
	//	cout << "T3 PEAKS: " << (float)(numPeaksT3)/(float)num_peaks << endl;


		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_T1,intensityT1*oneOverTotalIntensity));
		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_T2,intensityT2*oneOverTotalIntensity));
		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_T3,intensityT3*oneOverTotalIntensity));

	//	cout << "T1 INTEN: " << intensityT1/totalIntensity << endl;
	//	cout << "T2 INTEN: " << intensityT2/totalIntensity << endl;
	//	cout << "T3 INTEN: " << intensityT3/totalIntensity << endl;


		int numPeaksH2 = num_peaks - numPeaksH1;
		float intensityH2 = totalIntensity - intensityH1;

		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_H1,(numPeaksH1)*oneOverNumberOfPeaks));
		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_H2,(numPeaksH2)*oneOverNumberOfPeaks));

		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_H1,intensityH1*oneOverTotalIntensity));
		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_H2,intensityH2*oneOverTotalIntensity));
	
		const float inten33=0.333* totalIntensity;
		const float inten50=0.5  * totalIntensity;
		const float inten75=0.75 * totalIntensity;
		const float inten90=0.90 * totalIntensity;

		float cumulativeIntensity=0;
		const mass_t oneOverMz = 1.0 / (m_over_z + 0.1);

		for (i=0; i<num_peaks && cumulativeIntensity<inten33; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==num_peaks)
			i--;

		sam.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_33_INTEN,(peaks[i].mass*oneOverMz)-0.333));

		for ( ; i<num_peaks && cumulativeIntensity<inten50; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==num_peaks)
			i--;

		sam.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_50_INTEN,(peaks[i].mass*oneOverMz)-0.5));

		for ( ; i<num_peaks && cumulativeIntensity<inten75; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==num_peaks)
			i--;

		sam.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_75_INTEN,(peaks[i].mass*oneOverMz)-0.75));

		for ( ; i<num_peaks && cumulativeIntensity<inten90; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==num_peaks)
			i--;

		sam.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_90_INTEN,(peaks[i].mass*oneOverMz)-0.9));		
	}

	// pmc features
	if (1)
	{
		static vector< vector<float> > sqsFeatures;
		get_sqs_features_from_pmc_tables(bs,sqsFeatures);

		vector<float> maximalValues;
		maximalValues.resize(4,0.000001);

		vector< vector<float> > tmpVals;
		tmpVals.resize(4);

		// first give absolute counts
		int charge;
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqsFeatures[charge];
			const int idx_off = 4 *( charge-1);

			int i;
			for (i=0; i<4; i++)
			{
				sam.f_vals.push_back(fval(SQS_NUM_FRAG_PAIRS_1+idx_off+i,vals[i]));	
				tmpVals[charge].push_back(vals[i]);
				if (vals[i]>maximalValues[i])
					maximalValues[i]=vals[i];
			}	
		}

		// add features for prop of max
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqsFeatures[charge];
			const int idx_off = 4 *( charge-1);

			int i;
			for (i=0; i<4; i++)
			{
				const float ratio = (maximalValues[i]>0 ? vals[i]/maximalValues[i] : 0);
				sam.f_vals.push_back(fval(SQS_PROP_OF_MAX_FRAG_PAIRS_1+idx_off+i,ratio));			
				tmpVals[charge].push_back(ratio);
			}
		}

		// conver to proportions by dividing by the number of peaks/strong_peaks
		// and subtract the background levels (first)
		const float one_over_peaks = 1.0/((float)num_peaks+0.1);
		const float one_over_strong = 1.0/((float)numStrongPeaks+0.1);

		
	
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqsFeatures[charge];
			const int idx_off = 4 *( charge-1);

			// normalize to get proportions of total number of peaks/strong
			vals[0]-= this->currentSpectrumBackgroundStats_[charge].numFragmentPairs;
			vals[1]-= this->currentSpectrumBackgroundStats_[charge].numStrongFragmentPairs;
			vals[2]-= this->currentSpectrumBackgroundStats_[charge].numCharge2FragmentPairs;
			vals[3]-= this->currentSpectrumBackgroundStats_[charge].numStrongCharge2FragmentPairs;

			vals[0]*=one_over_peaks;
			vals[1]*=one_over_strong;
			vals[2]*=one_over_peaks;
			vals[3]*=one_over_strong;

			int i;
			for (i=0; i<4; i++)
			{
				sam.f_vals.push_back(fval(SQS_PROP_FRAG_PAIRS_1+idx_off+i,vals[i]));
				tmpVals[charge].push_back(vals[i]);
			}
		}

		// diff is between values for c2 and c3
		int i;
		for (i=0; i<tmpVals[2].size(); i++)
			sam.f_vals.push_back(fval(SQS_DIFF_NUM_FRAG_PAIRS_23+i,tmpVals[2][i]-tmpVals[3][i]));
	}

	
	sort(sam.f_vals.begin(),sam.f_vals.end());
}




void calculateBackgroundStatistics(const mass_t single_charge_pair_sum, // the sum of b+y or c+z
						  const Config *config,
						  const PeakList& pl,
						  const vector<bool>& strong_inds,
						  const vector<float>& iso_levels,
						  const vector<bool>& strict_iso_inds,
						  PMCRankStats& pmc_stats_total)
{
	const mass_t tolerance = config->getTolerance();
	const mass_t offsets[]={-22.0,-10.0,-8.5,8.5,12.0,22.5};
	const int numOffsets = sizeof(offsets)/sizeof(mass_t);
	const float one_over_num_offsets = 1.0 / (float)numOffsets;

	int i;
	for (i=0; i<numOffsets; i++)
	{
		static PMCRankStats pmc_stats; 
		
		calculatePmcRankStatisticsForMass(pl, single_charge_pair_sum+offsets[i], 
			1.5*tolerance,iso_levels, strict_iso_inds,pmc_stats);

		if (i==0)
		{
			pmc_stats_total = pmc_stats;
		}
		else
		{
			pmc_stats_total.numCharge2FragmentPairs        += pmc_stats.numCharge2FragmentPairs;
			pmc_stats_total.numFragmentPairs           += pmc_stats.numFragmentPairs;
			pmc_stats_total.numStrongCharge2FragmentPairs += pmc_stats.numStrongCharge2FragmentPairs;
			pmc_stats_total.numStrongFragmentPairs	 += pmc_stats.numStrongFragmentPairs;

			pmc_stats_total.inten_frag_pairs		 += pmc_stats.inten_frag_pairs;
			pmc_stats_total.inten_strong_pairs		 += pmc_stats.inten_strong_pairs;
			pmc_stats_total.inten_c2_pairs			 += pmc_stats.inten_c2_pairs;
			pmc_stats_total.inten_c2_strong_pairs	 += pmc_stats.inten_c2_strong_pairs;

			
		}
	}

	pmc_stats_total.numCharge2FragmentPairs		 *= one_over_num_offsets;
	pmc_stats_total.numFragmentPairs			 *= one_over_num_offsets;
	pmc_stats_total.numStrongCharge2FragmentPairs *= one_over_num_offsets;
	pmc_stats_total.numStrongFragmentPairs	 *= one_over_num_offsets;

	pmc_stats_total.inten_frag_pairs		 *= one_over_num_offsets;
	pmc_stats_total.inten_strong_pairs		 *= one_over_num_offsets;
	pmc_stats_total.inten_c2_pairs			 *= one_over_num_offsets;
	pmc_stats_total.inten_c2_strong_pairs	 *= one_over_num_offsets;

}

/*******************************************************************************

  Calculates the average statistics observed for the background using different
  mass offsets

********************************************************************************/
void calc_background_stats(const mass_t single_charge_pair_sum, // the sum of b+y or c+z
						  Config *config,
						  const QCPeak *peaks, 
						  const int num_peaks,
						  const vector<bool>& strong_inds,
						  const vector<float>& iso_levels,
						  const vector<bool>& strict_iso_inds,
						  PMCRankStats& pmc_stats_total)
{
	const mass_t tolerance = config->getTolerance();
	const mass_t offsets[]={-22.0,-10.0,-8.5,8.5,12.0,22.5};
	const int numOffsets = sizeof(offsets)/sizeof(mass_t);
	const float one_over_num_offsets = 1.0 / (float)numOffsets;

	int i;
	for (i=0; i<numOffsets; i++)
	{
		static PMCRankStats pmc_stats; 
		calc_pmc_rank_stats_for_mass(peaks,num_peaks, single_charge_pair_sum+offsets[i], 
			1.5*tolerance,iso_levels,strong_inds, strict_iso_inds,pmc_stats);

		if (i==0)
		{
			pmc_stats_total = pmc_stats;
		}
		else
		{
			pmc_stats_total.numCharge2FragmentPairs        += pmc_stats.numCharge2FragmentPairs;
			pmc_stats_total.numFragmentPairs           += pmc_stats.numFragmentPairs;
			pmc_stats_total.numStrongCharge2FragmentPairs += pmc_stats.numStrongCharge2FragmentPairs;
			pmc_stats_total.numStrongFragmentPairs	 += pmc_stats.numStrongFragmentPairs;

			pmc_stats_total.inten_frag_pairs		 += pmc_stats.inten_frag_pairs;
			pmc_stats_total.inten_strong_pairs		 += pmc_stats.inten_strong_pairs;
			pmc_stats_total.inten_c2_pairs			 += pmc_stats.inten_c2_pairs;
			pmc_stats_total.inten_c2_strong_pairs	 += pmc_stats.inten_c2_strong_pairs;

			
		}
	}

	pmc_stats_total.numCharge2FragmentPairs		 *= one_over_num_offsets;
	pmc_stats_total.numFragmentPairs			 *= one_over_num_offsets;
	pmc_stats_total.numStrongCharge2FragmentPairs *= one_over_num_offsets;
	pmc_stats_total.numStrongFragmentPairs	 *= one_over_num_offsets;

	pmc_stats_total.inten_frag_pairs		 *= one_over_num_offsets;
	pmc_stats_total.inten_strong_pairs		 *= one_over_num_offsets;
	pmc_stats_total.inten_c2_pairs			 *= one_over_num_offsets;
	pmc_stats_total.inten_c2_strong_pairs	 *= one_over_num_offsets;


}





/*********************************************************************************
Fills the PMC sample features
Assumes the frag_pair_sum_offset is set.
Sets the PMC features and also creates vetors of features that are to be added
to the SQS sampels.
**********************************************************************************/
void PMCSQS_Scorer::calculate_curr_spec_pmc_values(const BasicSpectrum& bs,
												   mass_t increment)
{
	const mass_t m_over_z = bs.ssf->m_over_z;


	if (frag_pair_sum_offset<-999)
	{
		set_frag_pair_sum_offset(MASS_PROTON); // b+y - PM+19
		set_bin_increment(0.1);
	}

	int charge;
	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		const mass_t org_pm_with_19 = bs.ssf->m_over_z * charge - MASS_PROTON*(charge - 1);
		const mass_t frag_pair_sum = org_pm_with_19 + frag_pair_sum_offset;
		const int size_idx = this->get_rank_model_size_idx(charge,org_pm_with_19);

		float bias = 0;
		if (pmcMzBiases_.size()>0 && pmcMzBiases_[charge].size()>size_idx)
			bias = this->pmcMzBiases_[charge][size_idx];

		if (config_->get_pm_tolerance()<0.075) // don't use bias if dealing with accurate pm
			bias=0;

		int bin_range =3*charge;
		if (charge>=2)
			bin_range = int(2*charge);
		if (charge>=3)
			bin_range = int(1.5*charge);

	/*	fill_rank_PMC_stats(  charge,
							  frag_pair_sum + bias,
							  -bin_range-1, 
							  bin_range-1,
							  increment,
							  config_,
							  bs,
							  currentSpectrumStrongPeakFlags_,
							  currentSpectrumIsotopeLevels_,
							  currentSpectrumStrictIsotopeFlags_,
							  currentSpectrumPmcTables_[charge]);*/


/*		calc_background_stats(frag_pair_sum,config_,bs.peaks,bs.num_peaks, 
			currentSpectrumStrongPeakFlags_, currentSpectrumIsotopeLevels_, 
			this->currentSpectrumStrictIsotopeFlags_, currentSpectrumBackgroundStats_[charge]);*/


		// find maximal values
		int i;
		PMCRankStats& maximal = currentSpectrumMaximalValues_[charge];
		maximal = currentSpectrumPmcTables_[charge][0];

		for (i=1; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			const PMCRankStats& curr_stats = currentSpectrumPmcTables_[charge][i];

			// mathced pairs look for maximum number

			if (curr_stats.numCharge2FragmentPairs>maximal.numCharge2FragmentPairs)
				maximal.numCharge2FragmentPairs = curr_stats.numCharge2FragmentPairs;

			if (curr_stats.numFragmentPairs>maximal.numFragmentPairs)
				maximal.numFragmentPairs = curr_stats.numFragmentPairs;

			if (curr_stats.numStrongFragmentPairs>maximal.numStrongFragmentPairs)
				maximal.numStrongFragmentPairs = curr_stats.numStrongFragmentPairs;

			if (curr_stats.numStrongCharge2FragmentPairs > maximal.numStrongCharge2FragmentPairs)
				maximal.numStrongCharge2FragmentPairs = curr_stats.numStrongCharge2FragmentPairs;

			// intensity look for maximal numbers

			if (curr_stats.inten_frag_pairs > maximal.inten_frag_pairs)
				maximal.inten_frag_pairs = curr_stats.inten_frag_pairs;
			
			if (curr_stats.inten_strong_pairs> maximal.inten_strong_pairs)
				maximal.inten_strong_pairs = curr_stats.inten_strong_pairs;

			if (curr_stats.inten_c2_pairs> maximal.inten_c2_pairs)
				maximal.inten_c2_pairs = curr_stats.inten_c2_pairs;

			if (curr_stats.inten_c2_strong_pairs > maximal.inten_c2_strong_pairs)
				maximal.inten_c2_strong_pairs = curr_stats.inten_c2_strong_pairs;


	
		}

		// find indicators for min tolerance and max pairs
		float tol_pairs =           POS_INF;
		float tol_strong_pairs =    POS_INF;
		float tol_c2_pairs =        POS_INF;
		float tol_c2_strong_pairs = POS_INF;

		int   idx_pairs=0;
		int	  idx_strong_pairs=0;
		int   idx_c2_pairs=0;
		int   idx_c2_strong_pairs=0;

		for (i=0; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			const PMCRankStats& curr_stats = currentSpectrumPmcTables_[charge][i];

			if (curr_stats.numFragmentPairs == maximal.numFragmentPairs &&
				curr_stats.mean_offset_pairs < tol_pairs)
			{
				idx_pairs=i;
				tol_pairs=curr_stats.mean_offset_pairs;
			}

			if (curr_stats.numStrongFragmentPairs == maximal.numStrongFragmentPairs &&
				curr_stats.mean_offset_c2_strong_pairs < tol_strong_pairs)
			{
				idx_strong_pairs=i;
				tol_strong_pairs=curr_stats.mean_offset_c2_strong_pairs;
			}

			if (curr_stats.numCharge2FragmentPairs == maximal.numCharge2FragmentPairs &&
				curr_stats.mean_offset_c2_pairs < tol_c2_pairs)
			{
				idx_c2_pairs=i;
				tol_c2_pairs=curr_stats.mean_offset_c2_pairs;
			}

			if (curr_stats.numStrongCharge2FragmentPairs == maximal.numStrongCharge2FragmentPairs &&
				curr_stats.mean_offset_c2_strong_pairs < tol_c2_strong_pairs)
			{
				idx_c2_strong_pairs=i;
				tol_c2_strong_pairs=curr_stats.mean_offset_c2_strong_pairs;
			}

		}

		currentSpectrumPmcTables_[charge][idx_pairs].ind_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_strong_pairs].ind_strong_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_c2_pairs].ind_c2_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_c2_strong_pairs].ind_c2_strong_pairs_with_min_tol=true;


		static vector<float> log_distances;
		if (log_distances.size()<currentSpectrumPmcTables_[charge].size())
		{
			log_distances.resize(currentSpectrumPmcTables_[charge].size(),0);
			int i;
			for (i=1; i<log_distances.size(); i++)
				log_distances[i]=log(1.0+(float)i);
		}

		for (i=0; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			PMCRankStats& curr_stats = currentSpectrumPmcTables_[charge][i];
			curr_stats.log_dis_from_pairs_min_tol = log_distances[abs(i-idx_pairs)];
			curr_stats.log_dis_from_strong_pairs_min_tol = log_distances[abs(i-idx_strong_pairs)];
			curr_stats.log_dis_from_c2_pairs_min_tol = log_distances[abs(i-idx_c2_pairs)];
			curr_stats.log_dis_from_c2_strong_pairs_min_tol = log_distances[abs(i-idx_c2_strong_pairs)];
		}

	}
}
	


/************************************************************************
Looks at the best m/z value for each charge and returns
the number of fragment pairs, strong fragment pairs, 
charge 2 fragment pairs, and strong charge 2 fragment pairs
*************************************************************************/
void PMCSQS_Scorer::getSqsFeaturesFromPmcTabels(const PeakList& pl, 
									            vector< vector<float> >& sqsFeatures) const
{
	const mass_t originalMOverZ = pl.getHeader()->getMOverZ();
	float maxNumStrongPairs=0;
	int   bestTableIndex=0;
	mass_t mzDiff = POS_INF;

	if (sqsFeatures.size() != maximalChargeWithModels_ +1)
		sqsFeatures.resize(maximalChargeWithModels_+1);

	size_t charge;
	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		// find entry which has the maximal number of strong pairs, while maining the minimial
		// m/z distance shift from the original
		size_t i;
		for (i=1; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			const float currentNumStrong = currentSpectrumPmcTables_[charge][i].numStrongFragmentPairs;
			
			if (currentNumStrong >= maxNumStrongPairs)
			{
				const mass_t distance = fabs(currentSpectrumPmcTables_[charge][i].m_over_z - originalMOverZ);
				if (currentNumStrong == maxNumStrongPairs && distance>= mzDiff)
					continue;
				
				maxNumStrongPairs = currentNumStrong;
				mzDiff			  = distance;
				bestTableIndex    = i;
			}
		}

		const PMCRankStats& bestStatistics = currentSpectrumPmcTables_[charge][bestTableIndex];

		sqsFeatures[charge].resize(4);
		sqsFeatures[charge][0]=bestStatistics.numFragmentPairs;
		sqsFeatures[charge][1]=bestStatistics.numStrongFragmentPairs;
		sqsFeatures[charge][2]=bestStatistics.numCharge2FragmentPairs;
		sqsFeatures[charge][3]=bestStatistics.numStrongCharge2FragmentPairs;
	}
}

/*****************************************************************************************

******************************************************************************************/
void PMCSQS_Scorer::get_sqs_features_from_pmc_tables(const BasicSpectrum& bs,
						vector< vector<float> >& sqsFeatures) const
{
	const mass_t org_m_over_z = bs.ssf->m_over_z;
	float max_num_strong_pairs=0;
	int   best_table_idx=0;
	mass_t mz_diff = 99999;

	if (sqsFeatures.size() != maximalChargeWithModels_ +1)
		sqsFeatures.resize(maximalChargeWithModels_+1);

	int charge;
	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		// find entry which has the maximal number of strong pairs, while maining the minimial
		// m/z distance shift from the original
		int i;
		for (i=1; i<this->currentSpectrumPmcTables_[charge].size(); i++)
		{
			const float curr_num_strong = currentSpectrumPmcTables_[charge][i].numStrongFragmentPairs;
			
			if (curr_num_strong>=max_num_strong_pairs)
			{
				const mass_t distance = fabs(currentSpectrumPmcTables_[charge][i].m_over_z - org_m_over_z);
				if (curr_num_strong == max_num_strong_pairs && distance>= mz_diff)
					continue;
				
				max_num_strong_pairs = curr_num_strong;
				mz_diff = distance;
				best_table_idx = i;
			}
		}

		const PMCRankStats& best_stats = currentSpectrumPmcTables_[charge][best_table_idx];

		sqsFeatures[charge].clear();

		sqsFeatures[charge].push_back(best_stats.numFragmentPairs);

		sqsFeatures[charge].push_back(best_stats.numStrongFragmentPairs);

		sqsFeatures[charge].push_back(best_stats.numCharge2FragmentPairs);

		sqsFeatures[charge].push_back(best_stats.numStrongCharge2FragmentPairs);
	}
}








void PMCSQS_Scorer::create_filtered_peak_list_for_sqs(
									  QCPeak *org_peaks, int num_org_peaks,
									  QCPeak *new_peaks, int& num_new_peaks) const
{
	static vector<MassInten> peak_list;
	static int peak_list_size =0;
	int i;

	if (num_org_peaks>peak_list_size)
	{
		peak_list_size = (int)(num_org_peaks * 1.5);
		if (peak_list_size<2000)
			peak_list_size = 2000;

		peak_list.resize(peak_list_size);
	}

	// copy org_peaks to the temporary peak_list
	int f_idx=0;
	for (i=0; i<num_org_peaks; i++)
	{
		peak_list[i].mass=org_peaks[i].mass;
		peak_list[i].intensity=org_peaks[i].intensity;
	}

	const mass_t tolerance = config_->getTolerance();
	const mass_t join_tolerance = (tolerance < 0.05 ? tolerance : 0.5 * tolerance);
	int p_idx=0;
	i=1;
	while (i<num_org_peaks)
	{
		if (peak_list[i].mass - peak_list[p_idx].mass<=join_tolerance)
		{
			intensity_t inten_sum = peak_list[i].intensity + peak_list[p_idx].intensity;
			mass_t new_mass = (peak_list[i].intensity * peak_list[i].mass + 
							   peak_list[p_idx].intensity * peak_list[p_idx].mass ) / inten_sum;

			peak_list[p_idx].mass = new_mass;
			peak_list[p_idx].intensity = inten_sum;	
		}
		else
		{
			peak_list[++p_idx]=peak_list[i];
		}
		i++;
	}
	int num_peaks = p_idx+1;


	// filter low intensity noise
	const mass_t half_window_size = 0.5 * config_->get_local_window_size();
	const int num_peaks_in_window = config_->get_max_number_peaks_per_local_window();
	const int max_peak_idx = num_peaks -1;
	int min_idx=1;
	int max_idx=1;
	p_idx =1;

	
	new_peaks[0].mass      = peak_list[0].mass;
	new_peaks[0].intensity = peak_list[0].intensity;
	f_idx=1;

	// check the rest of the peaks
	for (i=1; i<max_peak_idx; i++)
	{
		const mass_t& peakMass=peak_list[i].mass;
		mass_t minimalMass = peak_list[min_idx].mass;
		mass_t maximalMass = peak_list[max_idx].mass;

		// advance min/max pointers
		while (peakMass-minimalMass > half_window_size)
			minimalMass=peak_list[++min_idx].mass;

		while (max_idx < max_peak_idx && maximalMass - peakMass <= half_window_size)
			maximalMass=peak_list[++max_idx].mass;

		if (maximalMass - peakMass > half_window_size)
			max_idx--;

		// if there are less than the maximum number of peaks in the window, keep it.
		if (max_idx-min_idx < num_peaks_in_window)
		{
			new_peaks[f_idx].mass = peak_list[i].mass;
			new_peaks[f_idx].intensity = peak_list[i].intensity;
			f_idx++;
			continue;
		}

		// check if this is one of the top peaks in the window
		int higher_count=0;
		for (int j=min_idx; j<=max_idx; j++)
			if (peak_list[j].intensity > peak_list[i].intensity)
				higher_count++;

		if (higher_count < num_peaks_in_window)
		{
			new_peaks[f_idx].mass = peak_list[i].mass;
			new_peaks[f_idx].intensity = peak_list[i].intensity;
			f_idx++;
		}
	}
	new_peaks[f_idx].mass = peak_list[i].mass;
	new_peaks[f_idx].intensity = peak_list[i].intensity;
	f_idx++;


	num_new_peaks = f_idx;

	// normalize intensities

	if (1)
	{
		intensity_t total_inten=0;

		for (i=1; i<num_new_peaks; i++)
			total_inten+=new_peaks[i].intensity;

		const mass_t one_over_total_inten = (1000.0 / total_inten);

		for (i=1; i<num_new_peaks; i++)
			new_peaks[i].intensity *= one_over_total_inten; 
	}
}



/*****************************************************************
Creates for each input file an mgf file that holds the spectra
that passed quality filtering does not correct PM and charge in the 
mgf files.
******************************************************************/
void PMCSQS_Scorer::output_filtered_spectra_to_mgfs(
									 Config *config,
									 const vector<string>& files,
									 char *out_dir,
									 float filter_prob, 
									 int& total_num_written, 
									 int& total_num_read)
{
	total_num_read = 0;
	total_num_read = 0;
	int f;
	for (f=0; f<files.size(); f++)
	{
		const char *spectra_file = files[f].c_str();
		
		string fname, mgf_name, map_name;
		getFileNameWithoutExtension(files[f].c_str(),fname);

		mgf_name = string(out_dir) + "/" + fname + "_fil.mgf";
		map_name = string(out_dir) + "/" + fname + "_map.txt";

		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(spectra_file, config);

		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();

		int  num_spec_written=0;
		bool first=true;
		ofstream out_stream, map_stream;

		int sc;
		for (sc=0; sc<sl.getNumHeaders(); sc++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(sc);
		
			AnnotatedSpectrum as;
			as.readSpectrum(sa, header, false);

		
			vector<mass_t> pms_with_19;
			vector<int>    charges;
			vector<PmcSqsChargeRes> pmcsqs_res;
			selectPrecursorMassesAndCharges(config, as, pms_with_19, charges, &pmcsqs_res);

			float prob = 0.0;
			if (charges.size()>0 && pmcsqs_res.size()>charges[0] && pmcsqs_res[charges[0]].sqs_prob<filter_prob)
			{
				continue;
			}

			if (charges.size()>0 && pmcsqs_res.size()>charges[0])
				prob = pmcsqs_res[charges[0]].sqs_prob;

			if (first)
			{
				out_stream.open(mgf_name.c_str(),ios::out);
				map_stream.open(map_name.c_str(),ios::out);

				cout << "Filtering spectra to minumum quality score: " << filter_prob << endl;
				cout << "Writing spectra info to:" << endl;
				cout << mgf_name << endl << map_name << endl;

				if (! out_stream.is_open() || ! out_stream.good())
				{
					cout << "Error: couldn\'t open for out mgf stream for writing: " <<
						endl << mgf_name << endl;
					exit(1);
				}
				first = false;
			}

			char single_name[64];
			sprintf(single_name,"%d:%d",f,as.getHeader()->getScanNumber());
			SingleSpectrumHeader* ssh = const_cast<SingleSpectrumHeader*>(as.getHeader());
			
			ssh->setTitle(single_name);
			as.outputToMgfStream(out_stream);
			if (prob>1.0)
				prob=1.0;
			map_stream << num_spec_written++ << "\t" << as.getHeader()->getScanNumber() << "\t" << fixed << prob << endl;
		}
		out_stream.close();
		map_stream.close();

		total_num_read+= sl.getNumHeaders();
		total_num_written += num_spec_written;
	}
}



	
