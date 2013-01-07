#include "QuickClustering.h"
#include "auxfun.h"

struct ppair {
	bool operator< (const ppair& other) const
	{
		return (inten>other.inten);
	}

	intensity_t inten;
	int idx;
};






 
// sets the adjusted intensity of the peaks
void set_adjusted_inten(QCPeak *peaks, int num_peaks)
{
	int i;
	float total_intensity = 0.0;
	for (i=0; i<num_peaks; i++)
		total_intensity+= peaks[i].intensity;

	if (total_intensity <=0)
		return;

	const float norm_factor = 1000.0/total_intensity;

	for (i=0; i<num_peaks; i++)
		peaks[i].adjusted_inten = log(norm_factor * peaks[i].intensity);

	return; 
} 






double num_similarity_calculations=0;




/************************************************************************
// Calculates the cosine dot product of to sets of peaks (designated by the idxs)
// the values that are used are the ones in the adjusted_inten field
*************************************************************************/
float calc_selected_dot_prod(mass_t tolerance, 
							 const QCPeak *pa, int na, const vector<int>& peak_idxs_a,
	  					     const QCPeak *pb, int nb, const vector<int>& peak_idxs_b,
							 bool verbose)
{	
	float top_sum = 0;
	float sum_a_sqr = 0;
	float sum_b_sqr = 0;

	num_similarity_calculations++;

	const int num_a = peak_idxs_a.size();
	const int num_b = peak_idxs_b.size();

	int a=0,b=0;

	float top=0;
	while (a<num_a && b<num_b)
	{
		const int a_idx = peak_idxs_a[a];
		const int b_idx = peak_idxs_b[b];

		if (a_idx<0 || a_idx>=na)
		{
			int i;
			for (i=0; i<num_a; i++)
				cout << "a " << i << "\t" << peak_idxs_a[i] << endl;
			int qq=1;
		}

		if (b_idx<0 || b_idx>=nb)
		{
			int i;
			for (i=0; i<num_b; i++)
				cout << "b " << i << "\t" << peak_idxs_b[i] << endl;
			int qq=1;
		}
		const QCPeak& peak_a = pa[a_idx];
		const QCPeak& peak_b = pb[b_idx];

		if (fabs(peak_a.mass - peak_b.mass)<=tolerance)
		{
			top += peak_a.adjusted_inten * peak_b.adjusted_inten;
			a++;
			b++;
			continue;

		}

		if (peak_a.mass<peak_b.mass)
		{
			a++;
		}
		else
			b++;
	}

	float sqr_a=0;
	for (a=0; a<num_a; a++)
	{
		const float& inten = pa[peak_idxs_a[a]].adjusted_inten;
		sqr_a += (inten*inten);
	}

	float sqr_b=0;
	for (b=0; b<num_b; b++)
	{
		const float& inten = pb[peak_idxs_b[b]].adjusted_inten;
		sqr_b += (inten*inten);
	}

	if (top<=0)
		return 0;

	float simz = (float)top/sqrt(sqr_a*sqr_b);
	if (simz>1)
		simz=1;

	return simz; 
}



// scans the peaks in the spectrum to select the idxs of the peaks with the
// highest ranks (both absolute and using a local window). Invalidates certian 
// mass such as M M-H2O, etc.
void select_top_peak_idxs(QCPeak *peaks, 
						  int num_peaks, 
						  mass_t m_over_z, 
						  mass_t tolerance, 
						  vector<int>& top_ranked_peak_idxs, 
						  float top_x_masses[NUM_TOP_CLUSTER_PEAKS],
						  int top_peaks_per_1000da,
						  const Config *config)
{
	static vector<bool> inv_peak_indicators;
	static vector<bool> high_inten_indicators;
	static vector<ppair> pairs;
	const mass_t local_window_size = 75.0; // also look at the strongest peak in the local window of this size
	const mass_t iso_mass_margin = 1.003 + tolerance * 0.6;
	const mass_t min_exclude_mass = (config ? config->get_min_exclude_range() : 99999);
	const mass_t max_exclude_mass = (config ? config->get_max_exclude_range() : -1 );
	int i;

	int num_global_high_inten = (int)(top_peaks_per_1000da * peaks[num_peaks-1].mass * 0.001);

	if (NUM_TOP_CLUSTER_PEAKS<1)
		num_global_high_inten = NUM_TOP_CLUSTER_PEAKS;

	// reallocate static memory
	if (num_peaks<1 || num_peaks>1E5)
	{
		cout << "Error: invalid number of peaks : " << num_peaks << endl;
		exit(1);
	}
	
	inv_peak_indicators.clear();
	high_inten_indicators.clear();

	inv_peak_indicators.resize(num_peaks,false);
	high_inten_indicators.resize(num_peaks,false);

	// mark invalid peaks such as isotopes and derivatives of the parent ion
	// set indicators if there is a peak ahead
	for (i=0; i<num_peaks; i++)
	{
		const mass_t forward_mass = peaks[i].mass + iso_mass_margin;
		int j=i+1;
		while (j<num_peaks && peaks[j].mass < forward_mass)
		{
			if (peaks[j].intensity>peaks[i].intensity)
			{
				inv_peak_indicators[i]=true;
				break;
			}
			++j;
		}

		const mass_t backwards_mass = peaks[i].mass - iso_mass_margin;
		j=i-1;
		while (j>=0 && peaks[j].mass > backwards_mass)
		{
			if (peaks[j].intensity>peaks[i].intensity)
			{
				inv_peak_indicators[i]=true;
				break;
			}
			--j;
		}
	}

	// find highest intensity peaks
	pairs.clear();	
	for (i=0; i<num_peaks; i++)
	{
		if (inv_peak_indicators[i])
			continue;

		if (max_exclude_mass>0)
		{
			const mass_t peak_mass = peaks[i].mass;
			if (peak_mass<=max_exclude_mass &&
				peak_mass>=min_exclude_mass &&
				config->check_if_mass_is_in_exclude_range(peak_mass))
			{
//				cout << "Excluded: " << peak_mass << endl;
				continue;
			}
		}

		ppair p;
		p.idx=i;
		p.inten = peaks[i].intensity;
		pairs.push_back(p);
	}
	sort(pairs.begin(),pairs.end());

	
	for (i=0; i<num_global_high_inten && i<pairs.size(); i++)
		high_inten_indicators[pairs[i].idx]=true;
		
	if (top_x_masses)
	{
		vector<float> top_masses;

		for (i=0; i<NUM_TOP_CLUSTER_PEAKS && i<pairs.size(); i++)
			top_masses.push_back(peaks[pairs[i].idx].mass);
	
		sort(top_masses.begin(),top_masses.end());
		for (i=0; i<top_masses.size(); i++)
			top_x_masses[i]=top_masses[i];
	}
	
	// collect the selected peak idxs
	top_ranked_peak_idxs.clear();
	for (i=0; i<num_peaks; i++)
		if (high_inten_indicators[i])
			top_ranked_peak_idxs.push_back(i);
}






// functions for dot product statistics















