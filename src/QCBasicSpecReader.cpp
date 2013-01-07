#include "QuickClustering.h"
#include "Isotopes.h"



int BasicSpecReader::read_basic_spec(const Config* config, 
									 const FileManager& fm,
									 SingleSpectrumFile *ssf, 
									 QCPeak* peaks, 
									 bool override_file_idx,
									 bool no_processing)
{
	const mass_t tolerance = config->getTolerance();
	int i,num_peaks = -1;

	if (ssf->num_peaks > peak_list.size())
	{
		max_peak_list_size = ssf->num_peaks * 2;
		if (max_peak_list_size<2000)
			max_peak_list_size = 2000;
	
		peak_list.resize(max_peak_list_size);
	}

//	cout << ssf->file_idx << " " << ssf->single_name  << " NUM PEAKS: " << ssf->num_peaks << endl;

	if (ssf->type == DTA)
	{
		num_peaks = get_peak_list_from_DTA(ssf->single_name.c_str());
		if (num_peaks != ssf->num_peaks)
		{
			cout << "Heh heh ... there is an error!!!" << endl;
			exit(1);
		}
	}
	else if (ssf->type == MGF)
	{
		if (this->current_mgf_file_idx != ssf->file_idx)
		{
			if (mgf_stream)
			{
				fclose(mgf_stream);
			}

			int file_idx = ssf->file_idx;
			if (override_file_idx)
				file_idx = 0;

			const string& fname = fm.get_mgf_file(file_idx).mgf_name.c_str();
			
			mgf_stream=fopen(fname.c_str(),"r");
			if (! mgf_stream)
			{
				cout << "Error: couldn't open mgf: " << fname << endl;
				exit(1);
				
			}
			this->current_mgf_file_idx = ssf->file_idx;
		}

		if ( fseek(mgf_stream,ssf->file_pos,SEEK_SET) )
		{
			cout << "Error: could not skip in file!" << endl;
			exit(1);
		}

		num_peaks = get_peak_list_from_MGF(mgf_stream);
		MGF_single *mgf_ssf = (MGF_single *)ssf;
		if (num_peaks>0 && peak_list[0].mass != mgf_ssf->first_peak_mass)
		{
			mgf_ssf->print_ssf_stats(config);
			cout << "Error: mismatch in first peak masses: " << setprecision(5) << mgf_ssf->first_peak_mass << " " << peak_list[0].mass << endl;
			cout << "This error could arise because of Unix/Windows issues." << endl;
			cout << "Try running unix2dos on the input files." << endl;
			exit(1);
		}
	}
	else if (ssf->type == MZXML)
	{
		// overide if the mzXML was read entirely into the MZXML_file
		// via peak buff. Should know if this happend because 
		MZXML_single *mzxml_ssf = (MZXML_single *)ssf;
		if (mzxml_ssf->peak_buff_start_idx>=0)
		{
			int i;
			float *ptr;
			fm.copy_mzxml_peak_buff_ptr(&ptr);

			float *mzxml_peak_buff = ptr + mzxml_ssf->peak_buff_start_idx;
			int idx=0;
			for (i=0; i<mzxml_ssf->num_peaks; i++)
			{
				peaks[i].mass     =mzxml_peak_buff[idx++];
				peaks[i].intensity=mzxml_peak_buff[idx++];
			}
			// peaks already filtered
			return mzxml_ssf->num_peaks;
		}
		else
		{
			if (this->current_mzxml_file_idx != ssf->file_idx)
			{
				if (mzxml_stream)
				{
					fclose(mzxml_stream);
				}

				int file_idx = ssf->file_idx;
				if (override_file_idx)
					file_idx = 0;

				const string& fname = fm.get_mzxml_file(file_idx).mzxml_name.c_str();
				
				mzxml_stream=fopen(fname.c_str(),"r");
				if (! mzxml_stream)
				{
					cout << "Error: couldn't open mzxml: " << fname << endl;
					exit(1);
					
				}
				this->current_mzxml_file_idx = ssf->file_idx;
			}

			if ( fseek(mzxml_stream,ssf->file_pos,0) )
			{
				cout << "Error: could not skip in file!" << endl;
				exit(1);
			}
			num_peaks = get_peak_list_from_MZXML(mzxml_stream);
		}
	}
	else if (ssf->type == DAT)
	{
		if (this->current_dat_file_idx != ssf->file_idx)
		{
			if (dat_file.is_open())
			{
				dat_file.close();
				
			}

			int file_idx = ssf->file_idx;
			if (override_file_idx)
				file_idx = 0;

			const string& fname = fm.get_dat_file(file_idx).dat_name.c_str();
			
			dat_file.open(fname.c_str(), ios::in | ios::binary);
	
			if (! dat_file.is_open())
			{
				cout << "Error: couldn't open dat: " << fname << endl;
				exit(1);
			}
			this->current_dat_file_idx = ssf->file_idx;
		}

		dat_file.seekg(ssf->file_pos);


		// If the peaks are read from a DAT file, there is no need to perform a joining of the peaks
	//	num_peaks = get_peak_list_from_DAT(dat_file, &peak_list[0]);
		num_peaks = get_peak_list_from_DAT(dat_file, peaks);
		if (num_peaks<3)
			num_peaks=-1;
		return num_peaks;	
	}
	else if (ssf->type == PKL)
	{
		int file_idx = ssf->file_idx;
		if (override_file_idx)
				file_idx = 0;

		const string& pkl_name = fm.get_pkl_dir(file_idx).dir_path;
		string file_path = pkl_name + "/" + ssf->single_name;

		num_peaks = get_peak_list_from_PKL(file_path);
	}
	else if (ssf->type == MS2)
	{
		if (this->current_ms2_file_idx != ssf->file_idx)
		{
			if (ms2_stream)
			{
				fclose(ms2_stream);
			}

			const string& fname = fm.get_ms2_file(ssf->file_idx).ms2_name.c_str();
			
			ms2_stream=fopen(fname.c_str(),"rb");
			if (! mgf_stream)
			{
				cout << "Error: couldn't open ms2: " << fname << endl;
				exit(1);
				
			}
			this->current_ms2_file_idx = ssf->file_idx;
		}

		if ( fseek(ms2_stream,ssf->file_pos,SEEK_SET) )
		{
			cout << "Error: could not skip in file!" << endl;
			exit(1);
		}

		num_peaks = get_peak_list_from_MS2(ms2_stream);
	
	}
	else
	{
		cout << "Error: invalid file type:" << ssf->type <<  endl;
		exit(1);
	}

	if (num_peaks<3)
		num_peaks=-1;


	// copy peaks as is
	if (no_processing)
	{
		int i;
		for (i=0; i<num_peaks; i++)
			peaks[i]=peak_list[i];

		return num_peaks;
	}


	// join adjacent peaks
	const mass_t join_tolerance = (tolerance < 0.05 ? tolerance : 0.6 * tolerance);
	int p_idx=0;
	i=1;
	while (i<num_peaks)
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

	num_peaks = p_idx+1;
	


	// filter low intensity noise
	// and mark those that are good peaks
	int max_peak_idx = num_peaks -1;
	int min_idx=1;
	int max_idx=1;
	p_idx =1;

	vector<bool> indicators;

	mark_top_peaks_with_sliding_window(&peak_list[0], 
									   num_peaks,
									   config->get_local_window_size(),
									   config->get_max_number_peaks_per_local_window(),
									   indicators);

	if (ssf->precursor_intensity <=0)
	{

	//	cout << "PRecursor was zero!" << endl;
	//	exit(1);

		ssf->precursor_intensity=0;
		for (i=0; i<num_peaks; i++)
			ssf->precursor_intensity+=peak_list[i].intensity;
	}

	p_idx=0;
	for (i=0; i<num_peaks; i++)
		if (indicators[i])
			peaks[p_idx++]=peak_list[i];
	
	num_peaks = p_idx;



	// normalize intensities

	if (config->get_need_to_normalize() && ! config->get_itraq_mode())
	{
		intensity_t total_inten=0;
		for (i=0; i<num_peaks; i++)
			total_inten+=peaks[i].intensity;

		const mass_t one_over_total_inten = (1000.0 / total_inten);

		for (i=0; i<num_peaks; i++)
			peaks[i].intensity *= one_over_total_inten;
	}


	// sanity check
	if (peaks[0].mass<=0 && peaks[0].intensity<0)
	{
		cout << "Invalid mass or intensity at first peak: " << peaks[0].mass << "\t" << peaks[0].intensity << endl;
		exit(1);
	}

	for (i=1; i<num_peaks; i++)
	{
	
		if (peaks[i].intensity<0 || peaks[i].mass <0 || peaks[i].mass < peaks[i-1].mass)
		{
			cout << "Error: peak mismatches in file! (peak " << i << "/" << num_peaks-1 << ")" << endl;

			cout << fixed << setprecision(4);
			int j;
			for (j=0; j<num_peaks; j++)
				cout << peaks[j].mass << " " << peaks[j].intensity << endl;
			exit(1);
		}
	}

//	for (i=0; i<num_peaks; i++)
//		cout << "]] " << i << "\t" << peaks[i].mass << "\t" << peaks[i].intensity << endl;

	return num_peaks;
}


// these functions just extract the peak list from the spectrum file, return the actual
// number of peaks (after joining)

int BasicSpecReader::get_peak_list_from_DTA(const char* dta_name)
{
	ifstream fs(dta_name,ios::in);
	if (! fs)
		return -1;

	char buff[256];

	fs.getline(buff,256);

	if (fs.bad())
		return false;

	while (buff[0] =='#') 
	{
		fs.getline(buff,256);
	}
	
	

	int p_count=0;
	while (fs.getline(buff,256))
	{
		istringstream is(buff);

		mass_t& mass = peak_list[p_count].mass;
		intensity_t& intensity = peak_list[p_count].intensity;
		is >> mass >> intensity;

		if (mass <0 || intensity==0)   // the peak probably got rounded down
			continue;

		p_count++;
	
	}

	fs.close();
	return p_count;
}

int BasicSpecReader::get_peak_list_from_MGF(FILE *mgf_stream)
{
	char buff[1024];
	
	mass_t peak_mass=-1;
	intensity_t peak_intensity=-1;

	while (1)
	{

		if( ! fgets(buff, 1024, mgf_stream) )
			return -1;

		if (! strncmp("END IONS",buff,7))
			return -1;

		istringstream is(buff);
		is >> peak_mass >> peak_intensity;

		if (peak_mass>0 &&  peak_intensity>0)
			break;
	}
	
	peak_list[0].mass = peak_mass;
	peak_list[0].intensity = peak_intensity;

	int p_count = 1;

	while (fgets(buff, 256, mgf_stream))
	{
		istringstream is(buff);
		
		if (! strncmp("END IONS",buff,7))
			break;

		
		mass_t mass=-1;
		intensity_t intensity=-1;

		is >> mass >> intensity;
	
		if (mass <=0 || intensity<=0)   // the peak probably got rounded down
			continue;

		peak_list[p_count].mass		= mass;
		peak_list[p_count].intensity= intensity;

		p_count++;
	}
	return p_count;
}



// 
int BasicSpecReader::get_peak_list_from_DAT(ifstream& dat_file, QCPeak *peaks)
{
	const int num_header_bytes_to_read = sizeof(mass_t) + 2* sizeof(float) + 4 * sizeof(int);
	const int num_peaks_pos = sizeof(mass_t) + 3 *sizeof(int);
	static char      header_buff[num_header_bytes_to_read];
	static vector<float> peak_buff;
	static int    peak_buff_size = 0;

	memset(header_buff,0,num_header_bytes_to_read);
	
	// read header info

	if (! dat_file.read(header_buff,num_header_bytes_to_read))
	{
		cout << "Error: reading header info from DAT file!" << endl;
		exit(1);
	}

	const int num_peaks = *(int *)(header_buff + num_peaks_pos);
	if (num_peaks<3)
		return -1;

	if (num_peaks>peak_buff_size)
	{
		if (num_peaks> 100000)
		{
			cout << "Error: too many peaks: " << num_peaks << endl;
			exit(1);
		}
		peak_buff_size = num_peaks * 2;
		if (peak_buff_size<8000)
			peak_buff_size = 8000;

		peak_buff.resize(peak_buff_size * 2);
	}

	// read peaks
	int num_peak_bytes_to_read = 2*sizeof(float) * num_peaks;
	if (! dat_file.read((char *)&peak_buff[0],num_peak_bytes_to_read) )
	{
		cout << "Error reading peak info from dat file! (" << 
			num_peak_bytes_to_read << " bytes)" <<  endl;
		exit(1);
	}
	
	// copy peaks to peak list
	int i;
	int pos=0;
	for (i=0; i<num_peaks; i++)
	{
		peaks[i].mass = peak_buff[pos++];
		peaks[i].intensity = peak_buff[pos++];
	}

	// sanity check
	for (i=1; i<num_peaks; i++)
	{
	
		if (peaks[i].intensity<0 || peaks[i].mass <0 || peaks[i].mass < peaks[i-1].mass)
		{
			cout << "Error: peak mismatches in DAT file! (peak " << i << "/" << num_peaks-1 << ")" << endl;

			cout << fixed << setprecision(4);
			int j;
			for (j=0; j<num_peaks; j++)
				cout << peaks[j].mass << " " << peaks[j].intensity << endl;
			exit(1);
		}
	}

	return num_peaks;
}




int BasicSpecReader::get_peak_list_from_MS2(FILE *ms2_stream)
{
	char buff[128];
	
	int p_count=0;
	while (fgets(buff, 128, ms2_stream))
	{
		istringstream is(buff);
	
		if (strlen(buff)<3 || buff[0] == ':') // try reading until we get to an empty line
			break;
	
		mass_t& mass = peak_list[p_count].mass;
		intensity_t& intensity = peak_list[p_count].intensity;
		is >> mass >> intensity;
	
		if (mass <0 || intensity==0)   // the peak probably got rounded down
			continue;

		p_count++;
	}
	return p_count;
}


int BasicSpecReader::get_peak_list_from_PKL(const string& pkl_path)
{
	char buff[128];
	
	FILE *pkl_stream = fopen(pkl_path.c_str(),"r");
	if (! pkl_stream)
	{
		cout << "Error: couldn't open pkl sinlge file for reading: " << endl;
		cout << pkl_path << endl;
		exit(1);
	}

	fgets(buff, 128, pkl_stream); // skip first line

	int p_count=0;
	while (fgets(buff, 128, pkl_stream))
	{
		istringstream is(buff);
	
		mass_t mass = -1;
		intensity_t intensity =-1;
		
		is >> mass >> intensity;
	
		if (mass <0 || intensity==0)   // the peak probably got rounded down
			continue;

		peak_list[p_count].mass = mass;
		peak_list[p_count].intensity = intensity;

		p_count++;
	}

	fclose(pkl_stream);

	return p_count;
}


struct ppair {
	bool operator< (const ppair& other) const
	{
		return (inten>other.inten);
	}

	intensity_t inten;
	int idx;
};




int BasicSpectrum::get_number_of_matching_peaks(mass_t tolerance, const vector<mass_t>& masses) const
{
	int count=0;
	int p_idx=0;
	int m_idx=0;

	while (m_idx<masses.size())
	{
		while (p_idx<num_peaks && masses[m_idx]>peaks[p_idx].mass)
			p_idx++;

		if (p_idx == num_peaks)
			return count;
	
		if ((peaks[p_idx].mass - masses[m_idx] < tolerance) ||
			(p_idx>0 && masses[m_idx] - peaks[p_idx-1].mass < tolerance))
			count++;

		m_idx++;
	}
	return count;
}


/*****************************************************************
******************************************************************/
float BasicSpectrum::calc_signal_level()
{
	int i;
	intensity_t total=0;

	vector<intensity_t> intensities;
	intensities.resize(num_peaks);
	for (i=0; i<num_peaks; i++)
	{
		intensities[i]=peaks[i].intensity;
		total+=peaks[i].intensity;
	}

	sort(intensities.begin(),intensities.end());

	int num_signal_peaks = (int)(0.04 * peaks[num_peaks-1].mass);

/*	float max_inten = intensities[num_peaks-1];
	if (num_peaks>num_signal_peaks)
	{
		intensity_t grass_times_20 = intensities[num_peaks-num_signal_peaks] * 20;
		if (max_inten>grass_times_20)
			max_inten=grass_times_20;
	}

	signal_level = total / max_inten;

	return signal_level;

  */


	int half_num_peaks = (int)(num_peaks*0.5);
	
	if (num_signal_peaks>half_num_peaks)
		num_signal_peaks = half_num_peaks;

	intensity_t strong_inten = 0;
	const int strong_start_idx = num_peaks - num_signal_peaks;
	const int strong_end_idx   = num_peaks;
	for (i = strong_start_idx; i<strong_end_idx; i++)
		strong_inten += intensities[i];

	if (strong_inten<=0)
		return 0;

	signal_level = strong_inten /
		(num_signal_peaks*intensities[strong_start_idx]);
	return signal_level;


	intensity_t low_inten = 0;
	int low_start_idx = num_peaks/2 - num_signal_peaks;
	if (low_start_idx<0)
		low_start_idx=0;
	const int low_end_idx = low_start_idx + num_signal_peaks;
	for (i = low_start_idx; i<low_end_idx; i++)
		low_inten+=intensities[i];
	
	signal_level = strong_inten/low_inten;
	if (signal_level<1.0)
		signal_level =1.0;

	return signal_level;
}


typedef float IsoIntens[6];
typedef int	  IsoIdxs[6];

void fill_iso_intens(mass_t iso_tolerance, QCPeak *peaks, int num_peaks, 
					 int idx, IsoIntens intens, IsoIdxs idxs)
{
	int i;
	for (i=0; i<6; i++)
	{
		intens[i]=0;
		idxs[i]=-1;
	}

	intens[0]=peaks[idx].intensity;
	idxs[0]=idx;
	mass_t base = peaks[idx].mass;

	int p_idx=idx+1;
	for (i=1; i<6; i++)
	{
		mass_t min = base + i - iso_tolerance;
		mass_t max = base + i + iso_tolerance;
		
		while (p_idx<num_peaks && peaks[p_idx].mass<min)
			p_idx++;

		if (p_idx==num_peaks || peaks[p_idx].mass>max)
			break;
		
		intens[i]=peaks[p_idx].intensity;
		idxs[i]=p_idx;
	}
}

void  BasicSpectrum::calc_peak_isotope_levels(mass_t tolerance, vector<float>& iso_levels) const
{
	int i;
	const mass_t iso_tolerance = (tolerance< 0.25) ? tolerance : 0.25;

	iso_levels.clear();
	iso_levels.resize(num_peaks,0);

	const int last_peak_idx = num_peaks-1;

	for (i=0; i<last_peak_idx; i++)
	{	
		if (iso_levels[i]>0)
			continue;

		// look for +1 peak
		if (peaks[i].intensity <=0)
		{
			iso_levels[i]=1;
			continue;
		}
		
		if (peaks[i+1].mass - peaks[i].mass>1.2)
			continue;

		IsoIntens iso_intens;
		IsoIdxs   iso_idxs;
		fill_iso_intens(iso_tolerance,peaks,num_peaks,i,iso_intens,iso_idxs);
		if (iso_intens[1]<=0)
			continue;

		float one_over_intensity = 1.0 / peaks[i].intensity;
		float ratio1 = 	one_over_intensity * iso_intens[1];
	
		// ignore strong +1
		if ( ratio1 > 3.5)
			continue;

		// examine ratios
		vector<float> relative_ratios,observed_ratios,expected_ratios;
		
		observed_ratios.resize(6,0);
		observed_ratios[0]=1.0;
		int j;
		for (j=1; j<6 && iso_idxs[j]>=0; j++)
			observed_ratios[j]=iso_intens[j]*one_over_intensity;
		
		int last_iso = (j<6 ? j : 5);
		// get expected iso ratios
		calc_expected_iso_ratios(peaks[i].mass,expected_ratios,last_iso);

		// calc ratios between observed and expected		
		relative_ratios.resize(j+1,0);
		relative_ratios[0]=1;
		for (j=1; j<=last_iso; j++)
			relative_ratios[j]=observed_ratios[j] / expected_ratios[j];

		iso_levels[i]=0;

		float level_mul=1.0;
		for (j=1; j<= last_iso; j++)
		{
			const float rel_ratio = relative_ratios[j];
			float iso_level;

			if (rel_ratio>= 0.75 && rel_ratio<=1.333)
			{
				iso_level=2.0;
			}
			else if (rel_ratio >= 0.5 && rel_ratio <=2)
			{
				iso_level=1.3333;
			}
			else if (rel_ratio >= 0.3333 && rel_ratio <=3)
			{
				iso_level=0.6666;
			}
			else if (rel_ratio >= 0.25 && rel_ratio <= 4)
			{
				iso_level=0.3333;
			}
			else
				break;
			
			iso_levels[iso_idxs[j]] = iso_levels[iso_idxs[j-1]] + level_mul * iso_level;
			level_mul *= 0.5;
		}
	}
}


void  BasicSpectrum::mark_all_possible_isotope_peaks(mass_t tolerance, vector<bool>& iso_inds) const
{
	int i;
	const mass_t iso_tolerance = (tolerance< 0.25) ? tolerance : 0.25;
	const mass_t max_iso_diff = 1.0 + iso_tolerance;
	const mass_t min_iso_diff = 1.0 - iso_tolerance;

	iso_inds.clear();
	iso_inds.resize(num_peaks,false);

	const int last_peak_idx = num_peaks-1;

	for (i=0; i<last_peak_idx; i++)
	{	
		// look for +1 peak
		if (peaks[i].intensity <=0)
		{
			iso_inds[i]=true;
			continue;
		}
		
		if (peaks[i+1].mass - peaks[i].mass>1.25)
			continue;

		int forward_idx=i+1;
		while (forward_idx<num_peaks && peaks[forward_idx].mass - peaks[i].mass<max_iso_diff)
		{
			if (peaks[forward_idx].mass - peaks[i].mass>min_iso_diff)
				iso_inds[i+1]=true;
			forward_idx++;
		}
	}
}



struct PeakPair {
	PeakPair() {};
	PeakPair(int i, float in) : idx(i), intensity(in) {};

	bool operator< (const PeakPair& other) const
	{
		return intensity>other.intensity;
	}

	int idx;
	float intensity;
};

bool  BasicSpectrum::select_strong_peak_idxs(const vector<float>& iso_levels, 
											 vector<bool>& strong_indicators) const
{
	vector<bool> inds;

	if ( ! mark_top_peaks_with_sliding_window(peaks, num_peaks, 120, 3, inds))
		return false;

	// also mark the top 20 peaks (non isotopic)

	vector<PeakPair> pairs;
	pairs.resize(num_peaks);
	int i;
	for (i=0; i<num_peaks; i++)
		pairs[i]=PeakPair(i,peaks[i].intensity);
	
	sort(pairs.begin(),pairs.end());
	const int half_num_peaks = num_peaks/2;
	const int max_to_mark = (half_num_peaks<20 ? half_num_peaks : 20);
	for (i=0; i<max_to_mark; i++)
		inds[pairs[i].idx]=true;

	strong_indicators.resize(num_peaks,false);

	for (i=0; i<num_peaks; i++)
		strong_indicators[i]=(inds[i] && iso_levels[i]==0);

	return true;

}



