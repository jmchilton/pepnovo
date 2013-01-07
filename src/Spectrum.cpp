#include "Spectrum.h"
#include "Isotopes.h"
#include "QuickClustering.h"
#include "FileManagement.h"
#include "auxfun.h"


bool Spectrum::readSpectrum(const SpectraAggregator& sa,	
							const SingleSpectrumHeader* header, 
							bool indFilterSpectrum)
{

	config_ = sa.getConfig();

	const int numPeaksRead = readPeaksToLocalAllocation(sa, header);
	if (numPeaksRead<5)
		return false;

	copyHeaderInformation();

	initializePeakList(config_, indFilterSpectrum);

	if (! sanityCheck())
		return false;

	initializeSpectrum();

	return (numPeaks_>0);
}

void Spectrum::copyHeaderInformation()
{
	if (! header_)
	{
		cout << "Error: must set header field in spectrum before copying!" << endl;
		exit(1);
	}

	numPeaks_		  = header_->getOriginalNumPeaks();
	mOverZ_			  = header_->getMOverZ();
	originalPmWith19_ = header_->getOriginalPmWith19();
	charge_			  = header_->getCharge();
	scanNumber_		  = header_->getScanNumber();
	retentionTime_	  = header_->getRetentionTime();
	clusterSize_	  = header_->getClusterSize();

	if (header_->getPeptideStr().length()>2)
		peptide_.parseFromString(config_,header_->getPeptideStr());
}



void Spectrum::initializeSpectrum(void *ssf, bool indFilterSpectrum)
{
	SingleSpectrumFile* ssfPtr = static_cast<SingleSpectrumFile*>(ssf);
	SingleSpectrumHeader* h = new SingleSpectrumHeader;
	h->intializeFromSingleSpectrumFile(ssfPtr);
	h->setPeptideStr(ssfPtr->peptide.as_string(config_));
	header_ = h;
	

//	header_ = new SingleSpectrumHeader;
//	header_->intializeFromSingleSpectrumFile(static_cast<SingleSpectrumFile*>(ssf));
	//join_adjacent_peaks();

	this->joinAdjacentPeaks(config_->getTolerance());

	if (indFilterSpectrum)
	{
//		filter_peaks();
		this->filterWeakPeaks(config_, (correctedPmWith19_>0 ? correctedPmWith19_ : originalPmWith19_));
	}

//	init_index_array();
//	calc_ranks();
//	calc_log_local_ranks();

	this->normalizePeakIntensities(); 

	this->createIndexArray(indexArray_, maximalPeakMassToConsider_);
	this->calculatePeakRanks(ranks_);
	this->calculateLogLocalRanks(config_->get_local_window_size(), logLocalRanks_);
	


//	normalize_intensities();


	this->computeLogIntensities(logIntensities_);

//	calc_isotope_levels();
	this->calculateIsotopicLevels(config_->getTolerance(), isotopicLevels_);
//	this->calculateMonoisotpicRanks(isotopicLevels_, 

//	select_strong_peaks();
	this->selectStrongPeakIndexes(config_->get_number_of_strong_peaks_per_local_window(),
								  logLocalRanks_, isotopicLevels_, strongPeakIndexes_);
//	set_log_random_probs();
	this->calculateLogRandomProbabilities(logIntensities_, logRandomProbabilities_);
}

void Spectrum::initializeSpectrum()
{
	

	this->computeLogIntensities(logIntensities_); 
	this->createIndexArray(indexArray_, maximalPeakMassToConsider_);
	this->calculatePeakRanks(ranks_);
	this->calculateLogLocalRanks(config_->get_local_window_size(), logLocalRanks_);
	this->calculateIsotopicLevels(config_->getTolerance(), isotopicLevels_);
	this->selectStrongPeakIndexes(config_->get_number_of_strong_peaks_per_local_window(),
								  logLocalRanks_, isotopicLevels_, strongPeakIndexes_);
	this->calculateLogRandomProbabilities(logIntensities_, logRandomProbabilities_);

	if (numPeaks_ >0)
	{
		minimalPeakMass_ = peaks_[0].mass -1; // margin 1 Dalton
		maximalPeakMass_ = peaks_[numPeaks_-1].mass +1;
	
		if (charge_>=2)
		{
			maximalPeakMassToConsider_ =  originalPmWith19_ > maximalPeakMass_ ? 
						originalPmWith19_ : maximalPeakMass_;
		}
		else
			maximalPeakMassToConsider_ = maximalPeakMass_;
	}

	maximalPeakMassToConsider_ += 10.0; // margin of error

	sizeIndex_ = config_->calc_size_idx(charge_, originalPmWith19_);

}




void Spectrum::init_spectrum(bool perform_filtering)
{


	join_adjacent_peaks();
	if (perform_filtering)
		filter_peaks();

	init_index_array();
	calc_ranks();
	calc_log_local_ranks();
	normalize_intensities();
	calc_isotope_levels();
	select_strong_peaks();
	set_log_random_probs();	
}





/*********************************************************************
Initializes the index array.
For each rounded off Dalton m, it gives the index of the closest peak i
with mass m_i>= m.
**********************************************************************/
void Spectrum::init_index_array()
{
	int i,c,size=(int)(maximalPeakMassToConsider_+27.0);
	const int max_peak_idx = numPeaks_-1;
	
	indexArray_.clear();
	indexArray_.resize(size+1,max_peak_idx);
	
	i=0;
	int m=(int)peaks_[0].mass;
	while (i<m)
		indexArray_[i++]=0;

	c=0;
	while (c< max_peak_idx)
	{
		int curr_m=(int)peaks_[c].mass;
		int next_m = curr_m;
		int next_c = c;

		while (next_m == curr_m && next_c<max_peak_idx)
			next_m=(int)peaks_[++next_c].mass;

		if (next_m>=size)
			next_m = size-1;

		while (i<next_m)
			indexArray_[i++]=c;
		
		c=next_c;
	}
}




/**********************************************************************
It then joins all pairs that are less than max_proximity away from each
other.
**********************************************************************/
void Spectrum::join_adjacent_peaks()
{
	vector<Peak> new_peaks;
	mass_t max_proximity = config_->getTolerance() * 0.5;
	int i;

	if (numPeaks_ == 0)
		return;

	new_peaks.clear();
	new_peaks.push_back(peaks_[0]);
	int prev_idx=0;
	for (i=1; i<numPeaks_; i++)
	{
		if 	(peaks_[i].mass - new_peaks[prev_idx].mass < max_proximity)
		{
			// join peaks with proportion to their intensities

			int tt = new_peaks.size();

			intensity_t inten_sum=(new_peaks[prev_idx].intensity + peaks_[i].intensity);
			mass_t ratio = new_peaks[prev_idx].intensity/inten_sum;
			mass_t new_mass = ratio *new_peaks[prev_idx].mass + (1-ratio)*peaks_[i].mass;
			
			new_peaks[prev_idx].intensity = inten_sum;
			new_peaks[prev_idx].mass = new_mass;
		}
		else
		{
			new_peaks.push_back(peaks_[i]);
			prev_idx++;
		}
	}


	if (new_peaks.size() < numPeaks_)
	{
		memcpy(peaks_, &new_peaks[0], (new_peaks.size() * sizeof(Peak)));
		numPeaks_ = new_peaks.size();
	}
}



/*******************************************************************
Filters the number of peaks in the spectra. Keeps only the highest
intensity peaks in the window. Uses values for number of peaks per
windows of 100 Da. from the Config.
Also keeps pairs of peaks that add up to mass+20 (b,y pairs)
Also keeps neutral losses for kept peaks (-H2O and -NH3)
********************************************************************/
void Spectrum::filter_peaks()
{
	const mass_t max_allowed_peak_mass = maximalPeakMassToConsider_;
	const mass_t window_size = 0.5 * config_->get_local_window_size();
	const int num_peaks_in_window = config_->get_max_number_peaks_per_local_window();
	const mass_t tolerance = config_->getTolerance();
	int i,j,min_idx,max_idx;
	vector<bool> keep_peaks;	
	int new_num_peaks=0;
	vector<Peak> new_peaks;

	if (numPeaks_<5)
		return;

	int max_peak_idx = numPeaks_ -1;

	keep_peaks.resize(numPeaks_,false);

	// keep first peak and last peak
	keep_peaks[0]=true;
	if (peaks_[max_peak_idx].mass<max_allowed_peak_mass)
		keep_peaks[max_peak_idx]=true;
	
	min_idx=1;
	max_idx=1;

	// check the rest of the peaks
	for (i=1; i<max_peak_idx; i++)
	{
		mass_t peak_mass=peaks_[i].mass;
		mass_t min_mass=peaks_[min_idx].mass;
		mass_t max_mass=peaks_[max_idx].mass;

		if (peaks_[i].mass > max_allowed_peak_mass)
			break;

		// advance min/max pointers
		while (peak_mass-min_mass > window_size)
			min_mass=peaks_[++min_idx].mass;

		while (max_idx < max_peak_idx && max_mass - peak_mass <= window_size)
			max_mass=peaks_[++max_idx].mass;

		if (max_mass - peak_mass > window_size)
			max_idx--;

		// this peak might already be marked for keeping (isotpoic peak)
		if (keep_peaks[i])
			continue;

		// if there are less than the maximum number of peaks in the window, keep it.
		if (max_idx-min_idx < num_peaks_in_window)
		{
			keep_peaks[i]=true;
			continue;
		}

		// check if this is one of the top peaks in the window
		int higher_count=0;
		for (int j=min_idx; j<=max_idx; j++)
			if (peaks_[j].intensity > peaks_[i].intensity)
				higher_count++;

		if (higher_count < num_peaks_in_window)
		{
			keep_peaks[i]=true;
		}
	}



	// look for b/y pairs
	mass_t pm_with_20 = (originalPmWith19_>0 ? correctedPmWith19_ : originalPmWith19_) + 
						 MASS_PROTON ;

	mass_t pm_with_20_upper  = pm_with_20 + tolerance;
	mass_t pm_with_20_lower = pm_with_20 - tolerance;

	int f_idx =0;
	int b_idx = numPeaks_-1;
	while (f_idx<numPeaks_ && b_idx>=0)
	{
		if (! keep_peaks[f_idx])
		{
			f_idx++;
			continue;
		}

		while (b_idx>=0 && peaks_[f_idx].mass + peaks_[b_idx].mass > pm_with_20_upper )
			b_idx--;

		if (b_idx<0)
			break;

		mass_t mass_sum = peaks_[f_idx].mass + peaks_[b_idx].mass;
		if (mass_sum > pm_with_20_lower && mass_sum < pm_with_20_upper)
		{
			keep_peaks[f_idx]=true;
			keep_peaks[b_idx]=true;
		}
		f_idx++;
	}


	// look for -H2O -NH3 peaks
	// 17.0265, 18.0105
	const mass_t frag_tolerance = tolerance * 0.6;
	int p_idx = 0;
	while (p_idx<numPeaks_)
	{
		if (! keep_peaks[p_idx])
		{
			p_idx++;
			continue;
		}

		const mass_t upper_H2O = peaks_[p_idx].mass + MASS_H2O + frag_tolerance;
		const mass_t lower_H2O = peaks_[p_idx].mass + MASS_H2O - frag_tolerance;
		const mass_t upper_NH3 = peaks_[p_idx].mass + MASS_NH3 + frag_tolerance;
		const mass_t lower_NH3 = peaks_[p_idx].mass + MASS_NH3 - frag_tolerance;

		int f_idx;
		for (f_idx=p_idx+1; f_idx<numPeaks_; f_idx++)
		{
			if (peaks_[f_idx].mass>upper_H2O)
				break;

			if (peaks_[f_idx].mass>lower_H2O)
			{
				keep_peaks[f_idx]=true;
				break;
			}

			if (peaks_[f_idx].mass>=lower_NH3 && peaks_[f_idx].mass<=upper_NH3)
				keep_peaks[f_idx]=true;
		}
		p_idx++;
	}




	new_num_peaks=0;
	for (i=0; i<numPeaks_; i++)
		if (keep_peaks[i])
			new_num_peaks++;

	new_peaks.resize(new_num_peaks);
	
	j=0;
	for (i=0; i<numPeaks_; i++)
		if (keep_peaks[i])
			new_peaks[j++]=peaks_[i];
		
	numPeaks_ = j;

	if (numPeaks_ >0)
	{
		minimalPeakMass_ = peaks_[0].mass -1; 
		maximalPeakMass_ = peaks_[numPeaks_-1].mass +1;
	}
}



struct peak_pair {
	bool operator< (const peak_pair& other) const
	{
		return inten>other.inten;
	}
	int idx;
	intensity_t inten;
};

void Spectrum::calc_ranks()
{
	int i;
	vector<peak_pair> pairs;
	pairs.resize(numPeaks_);

	for (i=0; i<numPeaks_; i++)
	{
		pairs[i].idx=i;
		pairs[i].inten=peaks_[i].intensity;
	}

	sort(pairs.begin(), pairs.end());
	ranks_.resize(pairs.size());

	for (i=0; i<pairs.size(); i++)
		ranks_[pairs[i].idx]=i+1;
}


/*********************************************************************
// gives each peak it log local rank
// good be done more effciently...
**********************************************************************/
void Spectrum::calc_log_local_ranks()
{
	const mass_t half_window_size = config_->get_local_window_size() * 0.5;
	
	logLocalRanks_.resize(numPeaks_);	
	
	int i;
	for (i=0; i<numPeaks_; i++)
	{
		
		const mass_t peak_mass = peaks_[i].mass;

		const PeakRange pr= findPeaksInRange(peak_mass - half_window_size,
											 peak_mass + half_window_size);
		int above=0;
		int j;
		for (j=pr.low_idx; j<=pr.high_idx && j<numPeaks_; j++)
			if (peaks_[j].intensity>peaks_[i].intensity)
				above++;

		logLocalRanks_[i] = log(1.0 + static_cast<float>(above));
	}	
}



/**************************************************************
Calculates the normalized intensity for each peaks. 
The normalization is done so that the new sum of all intensities 
equals 1000
***************************************************************/
void Spectrum::normalize_intensities()
{
	int i;
	
	logIntensities_.clear();
	logIntensities_.resize(numPeaks_,0);

	if (! config_->get_need_to_normalize())
		return;

	// remove the intensity of pm+20 / 2 and pm+2/2 from the total intensity
	// if this is an issue, normalization might have to be done twice,
	// before and after parent mass correction

	if (totalPeakIntensity_<=0)
	{
		totalPeakIntensity_ =0;
		for (i=0; i<numPeaks_; i++)
			totalPeakIntensity_ += peaks_[i].intensity;
	}
	
	const float norm_val = 1000.0 / totalPeakIntensity_;
	totalPeakIntensity_ = 1000.0;
	for (i=0; i<numPeaks_; i++)
	{
		peaks_[i].intensity *= norm_val;
		logIntensities_[i] = log(1.0 + peaks_[i].intensity);
	}
}


struct rank_pair {
	bool operator< (const rank_pair& other) const
	{
		return inten > other.inten;
	}

	int idx;
	intensity_t inten;
};


/**************************************************************
Sets the iso_level for each peak. Iso level 0 means that there 
is no evidence that this peaks is an isotopic peak. The higher
the level, the more this looks like an isotopic peak
***************************************************************/
void Spectrum::calc_isotope_levels()
{
	int i;
	const mass_t iso_tolerance = (config_->getTolerance()<0.2) ?
								  config_->getTolerance() : 0.2;

	const int last_peak_idx = numPeaks_-1;

	isotopicLevels_.clear();
	isotopicLevels_.resize(numPeaks_,0);

	for (i=0; i<last_peak_idx; i++)
	{	
		if (isotopicLevels_[i]>0)
			continue;

		// look for +1 peak
		int idx1 = findPeakWithMaxIntensity(peaks_[i].mass + MASS_PROTON, iso_tolerance);
		if (idx1<0)  
			continue;

		float one_over_intensity = 1.0 / peaks_[i].intensity;
		float ratio1 = peaks_[idx1].intensity * one_over_intensity;

		// ignore strong +1
		if ( ratio1 > 3)
			continue;

		// examine ratios
		vector<float> expected_ratios, observed_ratios, relative_ratios;
		vector<int> iso_idxs;
		observed_ratios.resize(6);
		observed_ratios[0]=1.0;
		observed_ratios[1]= ratio1;

		iso_idxs.resize(6);
		iso_idxs[0]=i;
		iso_idxs[1]=idx1;

		// find additional peaks
		int j;
		for (j=2; j<=5; j++)
		{
			int idx = findPeakWithMaxIntensity(peaks_[i].mass + j, iso_tolerance);
			if (idx<0)
				break;
			observed_ratios[j] = peaks_[idx].mass * one_over_intensity;
			iso_idxs[j]=idx;
		}
		int last_iso = j-1;

		// get expected iso ratios
		calc_expected_iso_ratios(peaks_[i].mass,expected_ratios,j);

		// calc ratios between observed and expected		
		relative_ratios.resize(j);
		relative_ratios[0]=1;
		for (j=1; j<=last_iso; j++)
			relative_ratios[j]=observed_ratios[j] / expected_ratios[j];

		float level_mul=1.0;
		for (j=1; j<= last_iso; j++)
		{
			float iso_level;

			if (relative_ratios[j]>= 0.75 && relative_ratios[j]<=1.333)
			{
				iso_level=2.0;
			}
			else if (relative_ratios[j] >= 0.5 && relative_ratios[j] <=2)
			{
				iso_level=1.3333;
			}
			else if (relative_ratios[j] >= 0.3333 && relative_ratios[j] <=3)
			{
				iso_level=0.6666;
			}
			else if (relative_ratios[j] >= 0.25 && relative_ratios[j] <= 4)
			{
				iso_level=0.3333;
			}
			else
				break;

		//	if (relative_ratios[j] / relative_ratios[j-1] > 3)
		//		break;
			
			isotopicLevels_[iso_idxs[j]] = isotopicLevels_[iso_idxs[j-1]] + level_mul * iso_level;

		//	peaks[iso_idxs[j]].IsotopicLevel = peaks[iso_idxs[j-1]].IsotopicLevel + 
		//									   level_mul * iso_level;
			level_mul *= 0.5;
		}
	}

	vector<rank_pair> mono_pairs, iso_pairs;
	for (i=0; i<numPeaks_; i++)
	{
		rank_pair p;
		p.idx=i;
		p.inten = peaks_[i].intensity;
		if (isotopicLevels_[i] <= 0)
		{
			mono_pairs.push_back(p);
		}
		else
			iso_pairs.push_back(p);

	}

	sort(mono_pairs.begin(),mono_pairs.end());
	sort(iso_pairs.begin(),iso_pairs.end());

/*	ranked_mono_peaks.resize(numPeaks_);
	
	for (i=0; i<mono_pairs.size(); i++)
		ranked_mono_peaks[i]=mono_pairs[i].idx;

	for (i=0; i<iso_pairs.size(); i++)
		ranked_mono_peaks[mono_pairs.size()+i]=iso_pairs[i].idx;
*/
}



/*****************************************************************
Find the strongest peaks by compairng their log local rank
to the number from the config file
******************************************************************/
void Spectrum::select_strong_peaks()
{
	int i;
	float thresh_log_level = log((float)(config_->get_number_of_strong_peaks_per_local_window()));

	strongPeakIndexes_.clear();
	for (i=0 ;i<numPeaks_; i++)
		if (logLocalRanks_[i]<=thresh_log_level && isotopicLevels_[i]<=0.0)
			strongPeakIndexes_.push_back(i);
}





/***********************************************************************
calcs for each peak the probability of observing it at random (based on
 the neighbor's distribution. Assumes the log_intens are distributed
 according to a normal disribution.
************************************************************************/
void Spectrum::set_log_random_probs()
{
	if (numPeaks_<2)
		return;
	const float log_add = log(1.2);
	const float one_over_sqr_2pi = 1.0 / sqrt(2*3.1415927);
	const mass_t peak_window_size = 0.6; // this is fixed and independent of the tolerance!
	const mass_t margin      = 25.0;
	const mass_t window_size = 100.0;
	const mass_t min_mass = peaks_[0].mass;
	const mass_t max_mass = peaks_[numPeaks_-1].mass;
	const mass_t viz_range = (max_mass - min_mass);

	logRandomProbabilities_.resize(numPeaks_);	
	int i;
	for (i=0; i<numPeaks_; i++)
	{
		const mass_t peak_mass    = peaks_[i].mass;
		const mass_t rel_position = (peak_mass-min_mass)/viz_range;
		const mass_t left_window  = margin + rel_position * window_size;
		const mass_t right_window = margin + window_size - left_window;
		const PeakRange pr = findPeaksInRange(peak_mass-left_window,peak_mass+right_window);
		const float peak_window_prob = peak_window_size /(left_window + right_window);

		// some freak cases have 0 peak counts (only in unix)
		const int num_peaks_in_range = (pr.num_peaks>0 ? pr.num_peaks : 1);
		const float zero_prob = pow((1.0 - peak_window_prob),num_peaks_in_range); 

		if (pr.num_peaks<5)
		{
			logRandomProbabilities_[i] = log(1.0-zero_prob) + log_add;
		}
		else
		{
			vector<float> log_intens;
			int j;
			for (j=pr.low_idx; j<=pr.high_idx; j++)
				log_intens.push_back(logIntensities_[j]);
		
			float mean=0,sd=1;
			calc_mean_sd(log_intens,&mean,&sd);
			const float e = (logIntensities_[i] - mean)/sd;
			if (e<0)
			{
				logRandomProbabilities_[i] = log(1 - zero_prob) + log_add;
			}
			else
			{
				const float norm = (one_over_sqr_2pi/ sd) * exp(-0.5*e*e);
				const float norm_const = (1 - zero_prob) / (one_over_sqr_2pi/ sd); //
				logRandomProbabilities_[i] = log(norm*norm_const);
			}
		} 
	}
}
								 






/****************************************************************
Calculates different pm_with_19 values, for different charges
(calculates assumed pm_with_19 from the m/z value)
*****************************************************************/
void Spectrum::calc_original_pm_with_19_for_different_charges(
								   vector<mass_t>& pms_with_19) const
{
	int i;
	mass_t mz = (originalPmWith19_ - MASS_PROTON) / charge_;

	pms_with_19.resize(5);
	pms_with_19[0]=-1;
	for (i=1; i<=4; i++)
		pms_with_19[i] = mz * i + MASS_PROTON;
}


/*********************************************************
Reads the spectrum info from a dta file.
returns false if there was a problem
**********************************************************/
bool Spectrum::read_from_dta( const Config* _config, const char *dta_name, 
							 int _charge , bool verbose)
{
	ifstream fs(dta_name,ios::in);
	if (! fs)
		return false;

	return false;

/*	peaks.clear();
	config_ = _config;
	file_name = dta_name;

	if (verbose)
		cout << "Reading: " << file_name << endl;

	char buff[1024];

	fs.getline(buff,1024);
	if (fs.bad())
		return false;

	while (buff[0] =='#') 
	{
		istringstream is(buff);
		string arg1,arg2;

		is >> arg1 >> arg2;

		if (! arg1.compare(0,5,"#FILE"))
		{
			file_name = arg2;
			if (verbose)
				cout << "#FILE " << file_name << endl;
		}
		else 
		if (! arg1.compare(0,4,"#SEQ"))
		{
			peptide_.parse_from_string(config_,arg2);
			if (verbose)
				cout << "#SEQ " << peptide_.as_string(config_) << endl;

		}
		fs.getline(buff,1024);
	}
	
	istringstream is(buff);
	is >> originalPmWith19_ >> charge_;

	if (_charge>=0)
		charge_ = _charge;

//	charge =0;

	if (charge_ == 0)
	{
		mOverZ_ = originalPmWith19_;
	}
	else
		mOverZ_ = (originalPmWith19_ + (charge_ - 1) * 1.0078)  / charge_;

	if (verbose)
		cout << originalPmWith19_ << " " << charge_ << endl;


	while (fs.getline(buff,1024))
	{
		istringstream is(buff);
		Peak p;
		is >> p.mass >> p.intensity;
	
		if (p.mass <0 || p.intensity==0)   // the peak probably got rounded down
			continue;

		peaks.push_back(p);
		total_original_intensity+=p.intensity;

		if (verbose)
			cout << p.mass << " " << p.intensity << endl;
	}

	if (peaks.size() >0)
	{
		min_peak_mass = peaks[0].mass -1; // margin 1 Dalton
		max_peak_mass = peaks[peaks.size()-1].mass +1;
	}
	else
		return false;


	if (charge_>=2)
	{
		max_mass =  originalPmWith19_ > max_peak_mass ? 
					originalPmWith19_ : max_peak_mass;
	}
	else
		max_mass = max_peak_mass;

	max_mass += 11.0; // margin of error

	size_idx = config_->calc_size_idx(charge_, originalPmWith19_);

	return true;*/
}





bool Spectrum::read_from_MGF_stream(const Config *_config, FILE *stream, int _charge, bool verbose)
{
	config_ = _config;
	totalPeakIntensity_ = 0;
	char buff[1024];

	peaks_ = new Peak[5000];
	localAllocationSize_ = 5000;
	
/*	while (1)
	{
		if( ! fgets(buff, 1024, stream) )
			return false;

		if (strlen(buff)<3)
			continue;

		if (strncmp(buff,"BEGIN IONS",10) )
			continue;

		break;
	}*/
	scanNumber_		=-1;
	retentionTime_	=-1;
	clusterSize_	=-1;
	charge_			=-1;
	mOverZ_			=-1;
	originalPmWith19_=-1;




	// read header info and first peak
	while (1)
	{
		if( ! fgets(buff, 1024, stream))
			return false;

		if (strlen(buff)<3)
			continue;

		if (strncmp(buff,"BEGIN IONS",10) )
			continue;

		if (! strncmp(buff,"END IONS",7))
			return false;

		if (! strncmp(buff,"TITLE=",6) )
		{
			buff[strlen(buff)-1]='\0';
			title_ = static_cast<string>(buff+6);
			continue;
		}
		else
		if (! strncmp(buff,"SEQ=",4) )
		{
			string seq_string(buff+4);
			peptide_.parseFromString(config_,seq_string);
			peptide_.calc_mass(config_);
			continue;		
		}
		else
		if (! strncmp(buff,"SCAN=",5) )
		{
			if (sscanf(buff+5,"%d",&scanNumber_) != 1)
			{
				cout << "Error: couldn't read scan number!" << endl;
				exit(1);
			}
			continue;
		}
		if (! strncmp(buff,"SCANS=",6) )
		{
			if (sscanf(buff+6,"%d",&scanNumber_) != 1)
			{
				cout << "Error: couldn't read scan number!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buff,"RT=",3) )
		{
			if (sscanf(buff+3,"%f",&retentionTime_) != 1)
			{
				cout << "Error: couldn't read retention_time!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buff,"CLUSTER_SIZE=",13) )
		{
			if (sscanf(buff+13,"%d",&clusterSize_) != 1)
			{
				cout << "Error: couldn't read cluster size!" << endl;
				exit(1);
			}
			continue;
		}
		else	
		if ( ! strncmp(buff,"CHARGE=",6))
		{
			if (sscanf(buff,"CHARGE=%d",&charge_) != 1)
			{
				cout <<  "Error: couldn't read charge!" << endl;
				return false;
			}
		}
		else
		if (! strncmp(buff,"PEPMASS=",8))
		{
			istringstream is(buff+8);
			is >> mOverZ_;
			
			if (mOverZ_ < 0)
			{
				cout << "Error: reading pepmass:" << mOverZ_ << endl;
				return false;
			}		
		}
		else
		{
			istringstream is(buff);
			Peak p;
	
			is >> p.mass >> p.intensity;

			if (p.mass >0 && p.intensity>0)
			{
				peaks_[0]=p;
				break;
			}
		}
	}



	if (_charge > 0)
			charge_ = _charge;

	originalPmWith19_ = mOverZ_ * charge_ + MASS_PROTON * (1 - charge_);
	if (originalPmWith19_<0)
		originalPmWith19_=-1;



	// read rest of peaks

	int peakIndex = 1;
	
	while (fgets(buff, 1024, stream))
	{
		istringstream is(buff);
		Peak p;

		if (! strncmp("END IONS",buff,7))
			break;

		is >> p.mass >> p.intensity;
	
		if (p.mass <0 || p.intensity==0)   // the peak probably got rounded down
		{
			continue;
		}

		peaks_[peakIndex++]=p;
		
		if (verbose)
			cout << p.mass << " " << p.intensity << endl;
	}

	numPeaks_ = peakIndex;
	if (numPeaks_ >0)
	{
		minimalPeakMass_ = peaks_[0].mass -1.0; // margin 1 Dalton
		maximalPeakMass_ = peaks_[numPeaks_-1].mass +1.0;
	}
	else
		return false;


	maximalPeakMassToConsider_ =      originalPmWith19_ > maximalPeakMass_ ? 
					originalPmWith19_ : maximalPeakMass_;

	maximalPeakMassToConsider_ += 11.0; // margin of error

	sizeIndex_ = config_->calc_size_idx(charge_, originalPmWith19_);

	if (peptide_.get_num_aas()>0 && charge_>0)
	{
		mass_t org_diff = originalPmWith19_ - peptide_.get_mass() - 19.0183;

		if (fabs(org_diff)>7.0)
		{
			// try and correct charge!
			for (charge_=1; charge_<=4; charge_++)
			{

				mass_t new_org_pm_with_19 =  mOverZ_ * charge_ + MASS_PROTON * (1 - charge_);
				mass_t diff = fabs(new_org_pm_with_19 - peptide_.get_mass() - MASS_OHHH);
				if (diff<6)
				{
					originalPmWith19_ = new_org_pm_with_19;
					return true;
				}
			}
		
			cout << "Error: sequence mass doesn't add up: " << title_ << " (diff: "
				 << setprecision(3) << fixed << org_diff << ")" << endl;
			cout << "Pepitde: " << peptide_.as_string(config_) << endl;
			cout << "Mass Cys = " << config_->get_session_tables().get_aa2mass(Cys) << endl;
			exit(1);
		}
	}

	return true;
}


bool Spectrum::read_from_peak_arrays(const Config* _config, Peptide& _pep, mass_t _pm_with_19, 
		int _charge, int _num_peaks, mass_t *_masses, intensity_t *_intensities)
{
	
	config_ = _config;
	title_ = "xxx";

	peptide_ =_pep;
	originalPmWith19_ = _pm_with_19;
	charge_ = _charge;
	mOverZ_ = (originalPmWith19_ + (charge_ - 1) * MASS_PROTON ) / charge_;

	localAllocationSize_ = _num_peaks;
	peaks_ = new Peak[localAllocationSize_];

	int i;
	for (i=0; i<_num_peaks; i++)
	{
		peaks_[i].mass = _masses[i];
		peaks_[i].intensity = _intensities[i];
	}

	if (numPeaks_ >0)
	{
		minimalPeakMass_ = peaks_[0].mass -1; // margin 1 Dalton
		maximalPeakMass_ = peaks_[numPeaks_-1].mass +1;
	}
	else
		return false;


	if (charge_>=2)
	{
		maximalPeakMassToConsider_ =  originalPmWith19_ > maximalPeakMass_ ? 
					originalPmWith19_ : maximalPeakMass_;
	}
	else
		maximalPeakMassToConsider_ = maximalPeakMass_;

	maximalPeakMassToConsider_ += 11.0; // margin of error

	sizeIndex_ = config_->calc_size_idx(charge_, originalPmWith19_);

	return true;
}

bool Spectrum::init_from_QCPeaks(const Config* _config, 
								 void *QCPeaks_ptr, 
								 int num_peaks, 
								 void *ssf_ptr, 
								 bool ind_perfrom_filtering)
{
	SingleSpectrumFile *ssf = (SingleSpectrumFile *)ssf_ptr;
	QCPeak     *QCpeaks = (QCPeak *)QCPeaks_ptr;

	config_ = _config;
	title_ = "xxx";

	charge_ = (ssf->charge>0 ? ssf->charge : 2);
	mOverZ_ = ssf->m_over_z;
	originalPmWith19_ = mOverZ_ * charge_ - ( charge_ - 1)* MASS_PROTON;
	peptide_=ssf->peptide;

	if (ssf->single_name.length()>0)
		title_ = ssf->single_name;

	localAllocationSize_ = num_peaks;
	peaks_ = new Peak[localAllocationSize_];
	numPeaks_ = num_peaks;
	int i;
	for (i=0; i<num_peaks; i++)
	{
		peaks_[i].mass = QCpeaks[i].mass;
		peaks_[i].intensity = QCpeaks[i].intensity;
	}

	if (numPeaks_ >0)
	{
		minimalPeakMass_ = peaks_[0].mass -1; // margin 1 Dalton
		maximalPeakMass_ = peaks_[numPeaks_-1].mass +1;
	}
	else
		return false;


//	if (charge>=2)
//	{
		maximalPeakMassToConsider_ =  originalPmWith19_ > maximalPeakMass_ ? 
					originalPmWith19_ : maximalPeakMass_;
//	}
//	else
//		max_mass = max_peak_mass;

	maximalPeakMassToConsider_ += 2.0; // margin of error

	sizeIndex_ = 0;

//	init_spectrum(ind_perfrom_filtering);
	this->initializeSpectrum(ssf_ptr,ind_perfrom_filtering);

	return true;
}


void Spectrum::print_spectrum(ostream& os) const
{
	int i;
	if (title_.length() > 0)
		os << "#TITLE " << title_ << endl;
	if (peptide_.get_num_aas()>0)
		os << "#SEQ " << peptide_.as_string(config_) << endl;

	os << fixed << setprecision(2) << originalPmWith19_ << " " << charge_ << endl;

	for (i=0; i<numPeaks_; i++)
	{
		os << left << setw(5) << i << setw(8) << setprecision(NUM_SIG_DIGITS) 
			<< fixed << right  << peaks_[i].mass;
		os << setw(12) << right  << setprecision(1) << peaks_[i].intensity << " " << setw(4) 
			<< setprecision(1) << isotopicLevels_[i] <<  setw(4) 
			<< setprecision(1) << logLocalRanks_[i] << "\t" << setprecision(3) << exp(logRandomProbabilities_[i]) << endl;
	}
}


void Spectrum::output_as_MGF(ostream& os) const
{
	os << "BEGIN IONS" << endl;
	os << "TITLE=" << title_ << endl;
	
	if (peptide_.get_num_aas()>0)
		os << "SEQ=" << peptide_.as_string(config_) << endl;
	
	if (scanNumber_>=0 && clusterSize_ == 1)
		os << "SCANS=" << scanNumber_ << endl;
	
//	if (retention_time >=0)
//		os << "RT=" << retention_time << endl;

	if (clusterSize_>1)
		os << "CLUSTER_SIZE=" << clusterSize_ << endl;

	os << "CHARGE=+" << charge_ << endl;

	os << "PEPMASS=" << fixed << setprecision(NUM_SIG_DIGITS) << mOverZ_ << endl;

	int i;
	for (i=0; i<numPeaks_; i++)
		os << fixed << setprecision(NUM_SIG_DIGITS) << peaks_[i].mass << " " 
		   << fixed << setprecision(NUM_SIG_DIGITS) << peaks_[i].intensity << endl;
	
	os << "END IONS" << endl << endl;
}



void Spectrum::print_expected_by(ostream& os) const
{
	vector<string> labels;
	labels.push_back("b");
	labels.push_back("y");
	print_expected_fragment_peaks(labels,os);
}

void Spectrum::print_expected_fragment_peaks(vector<string>& frag_labels, ostream& os) const
{
	int i;
	const vector<FragmentType>& all_fragments = config_->get_all_fragments();
	vector<mass_t> break_masses;
	
	const int frag_label_width =20;

	if (peptide_.get_num_aas() == 0)
		return;

	vector<string> pre_frags, suf_frags;
	vector<int> pre_frag_idxs, suf_frag_idxs;
	
	for (i=0; i<frag_labels.size(); i++)
	{
		if (frag_labels[i].length()>0 &&
			(frag_labels[i][0] == 'a' || frag_labels[i][0]== 'b' || frag_labels[i][0] == 'c') )
		{
			pre_frags.push_back(frag_labels[i]);
			pre_frag_idxs.push_back(config_->get_frag_idx_from_label(frag_labels[i]));
		}
		else
		{
			suf_frags.push_back(frag_labels[i]);
			suf_frag_idxs.push_back(config_->get_frag_idx_from_label(frag_labels[i]));
		}
	}


	peptide_.calc_expected_breakage_masses(config_,break_masses);
	mass_t true_mass_with_19 = peptide_.get_mass() + MASS_OHHH;
	os << peptide_.as_string(config_) << " (" << fixed << setprecision(4) << true_mass_with_19 << ")" << endl;

	if ( suf_frags.size() == 0)
	{
		os << setw(frag_label_width*frag_labels.size()+3)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		os << "|  |";
		for (i=0; i<frag_labels.size(); i++)
		{
			int w = (frag_label_width+6- frag_labels[i].length())/2;
			os << setw(w) << " " << setw(frag_label_width-w-1) << left << frag_labels[i] << "|";
		}
		os<<endl;

		os << setw(frag_label_width*frag_labels.size()+3)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		
		for (i=0; i<break_masses.size(); i++)
		{
			int region = config_->calc_region_idx(break_masses[i],
				true_mass_with_19, charge_, minimalPeakMass_, maximalPeakMass_);

			const RegionalFragments& rf = config_->get_regional_fragments(charge_, sizeIndex_, region);

			if (rf.get_frag_type_idxs().size() == 0)
			{
				cout << "Error: no fragments selected for region " << charge_ << " " <<
						sizeIndex_ << " " << region << endl;
				exit(1);
			}

			cout << "|" << setw(2) << left << i <<"|";
			
			int j;
			for (j=0; j<frag_labels.size(); j++)
			{
				int idx=pre_frag_idxs[j];
				if (rf.get_position_of_frag_type_idx(idx)<0)
					idx=-1;
				
				if (idx<0)
				{
					os << setw(frag_label_width-2) <<  left << "      - ";
				}
				else
				{
					mass_t exp_mass = all_fragments[idx].calc_expected_mass(break_masses[i],
										true_mass_with_19);
					int p_idx=findPeakWithMaxIntensity(exp_mass,config_->getTolerance());
					os << " " << setw(6) << right << fixed << setw(9) << setprecision(3) << exp_mass;
					if (p_idx>=0)
					{
						os << " " << setw(6) << fixed << setprecision(3) << peaks_[p_idx].mass - exp_mass;
					}
					else
						os << "       ";
				}
				os << " |";
			}
			os << endl;
		}
		os << setw(frag_label_width*frag_labels.size()+5)<< setfill('-') << right << " " << endl;
		os << setfill(' ');
	}
	else
	if ( pre_frags.size() == 0)
	{
		os << setw(frag_label_width*frag_labels.size()+5)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		os << "|  |";
		for (i=0; i<frag_labels.size(); i++)
		{
			int w = (frag_label_width - 1- frag_labels[i].length())/2;
			os << setw(w) << " " << setw(frag_label_width-1-w) << left << frag_labels[i] << "|";
		}
		os<<endl;

		os << setw(frag_label_width*frag_labels.size()+5)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		for (i=0; i<break_masses.size(); i++)
		{
			int region = config_->calc_region_idx(break_masses[i],
				true_mass_with_19, charge_, minimalPeakMass_, maximalPeakMass_);

			const RegionalFragments& rf = config_->get_regional_fragments(charge_, sizeIndex_,region);

			cout << "|" << setw(2) << left << break_masses.size() - i - 1<<"|";
			
			int j;
			for (j=0; j<frag_labels.size(); j++)
			{
				int idx=suf_frag_idxs[j];
				if (rf.get_position_of_frag_type_idx(idx)<0)
					idx=-1;

				if (idx<0)
				{
					os << setw(13) <<  left << "      - ";
				}
				else
				{
					mass_t exp_mass = all_fragments[idx].calc_expected_mass(break_masses[i],
										true_mass_with_19);
					int p_idx=findPeakWithMaxIntensity(exp_mass,config_->getTolerance());
					os << " " << setw(6) << right << fixed << setw(9) << setprecision(3) << exp_mass;
					if (p_idx>=0)
					{
						os << " " << setw(6) << fixed << setprecision(3) << peaks_[p_idx].mass - exp_mass;
					}
					else
						os << "       ";
				}
				os << " |";
			}
			os << endl;
		}

		os << setw(frag_label_width*frag_labels.size()+5)<< setfill('-') << right << " " << endl;
		os << setfill(' ');
	}
	else // have both prefix and suffix fragments, separate between them
	{
		os << setw(frag_label_width*frag_labels.size()+7)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		os << "|  |";
		for (i=0; i<pre_frags.size(); i++)
		{
			int w = (frag_label_width - 2 - pre_frags[i].length())/2;
			os << setw(w) << " " << setw(frag_label_width-2-w) << left << pre_frags[i] << "|";
		}
		os << "   |";
		for (i=0; i<suf_frags.size(); i++)
		{
			int w = (frag_label_width - 2 - suf_frags[i].length())/2;
			os << setw(w) << " " << setw(frag_label_width-2-w) << left << suf_frags[i] << "|";
		}
		os<<endl;

		os << setw(frag_label_width*frag_labels.size()+7)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		for (i=0; i<break_masses.size(); i++)
		{
			int region = config_->calc_region_idx(break_masses[i],
				true_mass_with_19, charge_, minimalPeakMass_, maximalPeakMass_);

			const RegionalFragments& rf = config_->get_regional_fragments(charge_, sizeIndex_, region);


			cout << "|" << setw(2) << left << i <<"|";
			
			int j;
			for (j=0; j<pre_frags.size(); j++)
			{
				int idx=pre_frag_idxs[j];
				if (rf.get_position_of_frag_type_idx(idx)<0)
					idx=-1;

				if (idx<0)
				{
					os << setw(13) <<  left << "      - ";
				}
				else
				{
					mass_t exp_mass = all_fragments[idx].calc_expected_mass(break_masses[i],
										true_mass_with_19);
					int p_idx=findPeakWithMaxIntensity(exp_mass,config_->getTolerance());
					os << " " << setw(6) << right << fixed << setw(9) << setprecision(3) << exp_mass;
					if (p_idx>=0)
					{
						os << " " << setw(6) << fixed << setprecision(3) << peaks_[p_idx].mass - exp_mass;
					}
					else
						os << "       ";
				}
				os << " |";
			}

			// output suffix fragments

			cout << " " << setw(2) << left << break_masses.size() - i - 1<<"|";
			for (j=0; j<suf_frags.size(); j++)
			{
				int idx=suf_frag_idxs[j];
				if (rf.get_position_of_frag_type_idx(idx)<0)
					idx=-1;

				if (idx<0)
				{
					os << setw(13) <<  left << "      - ";
				}
				else
				{
					mass_t exp_mass = all_fragments[idx].calc_expected_mass(break_masses[i],
										true_mass_with_19);
					int p_idx=findPeakWithMaxIntensity(exp_mass,config_->getTolerance());
					os << " " << setw(6) << right << fixed << setw(9) << setprecision(3) << exp_mass;
					if (p_idx>=0)
					{
						os << " " << setw(6) << fixed << setprecision(3) << peaks_[p_idx].mass - exp_mass;
					}
					else
						os << "       ";
				}
				os << " |";
			}


			os << endl;
		}

		os << setw(frag_label_width*frag_labels.size()+7)<< setfill('-') << right << " " << endl;
		os << setfill(' ');
	}

	// print pm stats

	

	os << "true - org_pm:    " << fixed << setprecision(4) << true_mass_with_19 << " - " <<
		originalPmWith19_ << "  = " << true_mass_with_19 - originalPmWith19_ << endl;
	
	if (this->correctedPmWith19_>0)
	{
		os << "true - cor1_pm:   "  << fixed << setprecision(4) << true_mass_with_19 << " - " <<
		correctedPmWith19_ << "  = " << true_mass_with_19 - correctedPmWith19_ << endl;
	}
	if (secondaryPmWith19_>0)
	{
		os << "true - cor2_pm: "  << fixed << setprecision(4) << true_mass_with_19 << " - " <<
		secondaryPmWith19_ << "  = " << true_mass_with_19 - secondaryPmWith19_ << endl;
	}
	os << endl;
	
}


// checks several charges to see which one has a good m over z match
// writes the expected b/y ions if finds a match
bool Spectrum::check_m_over_z_and_sequence(ostream& os)
{
	int c;
	mass_t pep_mass = peptide_.get_mass() + MASS_OHHH;

	for (c=1; c<=4; c++)
	{
		mass_t pm_19 = mOverZ_ * c - (c-1)*MASS_PROTON;

	//	cout <<"c: " <<c << "   pm19: " << pm_19 << endl;

		if (fabs(pm_19-pep_mass)<7)
		{
			set_org_pm_with_19(pm_19);
			setCharge(c);
		//	print_expected_by(os);
		//	os << endl;
			return true;
		}
	}

	for (c=1; c<=4; c++)
	{
		mass_t pm_19 = mOverZ_ * c - (c-1)*MASS_PROTON;

		cout <<"c: " <<c << "   pm19: " << pm_19 << endl;
	}

	os << "No charge matches: " << peptide_.as_string(config_) << "  (" <<
		peptide_.get_mass() + 19.083 << ") !!!! " << endl << endl;

	return false;
}


ostream& operator << (ostream& os, const Spectrum& spec)
{
	const Peak* const peaks = spec.getPeaks();
	const string& title = spec.getTitle();
	const Peptide& peptide = spec.getPeptide();
	const Config *config = spec.getConfig();

	if (title.length() > 0)
		os << "#TITLE " << title << endl;
	if (peptide.get_num_aas()>0)
		os << "#SEQ " << peptide.as_string(config) << endl;

	os << spec.get_org_pm_with_19() << " " << spec.getCharge() << endl;
	int i;
	for (i=0; i<spec.getNumPeaks(); i++)
		os << setw(8) << left << peaks[i].mass << " \t" << peaks[i].intensity << 
		" \t" << exp(spec.getPeakLogRandomProbability(i)) << endl;

	return os;
}





