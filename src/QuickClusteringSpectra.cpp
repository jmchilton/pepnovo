#include "QuickClustering.h"
#include "auxfun.h"
#include "AnnotatedSpectrum.h"


void read_mzXML_annotations(char *mzXML_list, 
							char *ann_file, 
							vector< vector<int> >& annotation_idxs, 
							vector<mzXML_annotation>& annotations,
							int max_ann_size) 
{
	int i;
	char buff[256];
	FILE *mzxml_stream = fopen(mzXML_list,"r");
	if (! mzxml_stream)
	{
		cout << "Error: couldn't open annotation file for mzXML run!: " << mzXML_list << endl;
		exit(1);
	}

	int count=0;
	while (fgets(buff,256,mzxml_stream))
	{
		count++;
	}
	fclose(mzxml_stream);

	annotation_idxs.resize(count+1);
	for (i=0; i<=count; i++)
		annotation_idxs[i].resize(max_ann_size,-1);

	FILE *ann_stream = fopen(ann_file,"r");
	if (! ann_stream)
	{
		cout << "Error: couldn't open annotation files for run!: " << ann_file << endl;
		exit(1);
	}

/*	b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2-19.mzXML'...
20989 spectra...
255 Parse spectra from 'C:/Work/Data/Briggs\H293b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2\H293
b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2-20.mzXML'...
20712 spectra...
256 Parse spectra from 'C:/Work/Data/Briggs\H293b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2\H293
b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2-21.mzXML'...
20603 spectra...*/

	

	i=0;
	while (fgets(buff,256,ann_stream))
	{
		int file_idx=-1, mzXML_idx=-1; 
		char only_peptide[128];
		int scan=-1,charge=0;

		if (sscanf(buff,"%d %d %d %d %s",&file_idx,&mzXML_idx,&scan,&charge,only_peptide)<5)
			continue;

		mzXML_annotation ann;
		ann.charge =charge;
		ann.mzXML_file_idx = file_idx;
		ann.pep = only_peptide;


		annotation_idxs[file_idx][scan]=annotations.size();
		annotations.push_back(ann);
		
	//	cout << file_idx << " : >> " << scan << "   " << only_peptide << " " << charge << endl;

	//	if (++i>=200)
	//		break;
	
	}
}


struct peak_idx_pair
{
	bool operator< (const peak_idx_pair& other) const
	{
		return (intensity>other.intensity);
	}
	int idx;
	intensity_t intensity;
};

/***********************************************************************
Uses a heuristic approach jumps every half window
************************************************************************/
bool mark_top_peaks_with_sliding_window(const QCPeak *peaks, 
										int num_peaks, 
										mass_t window_size, 
										int num_peaks_per_window, 
										vector<bool>& indicators)
{
	// filter low intensity noise
	// and mark those that are good peaks
	const mass_t half_window_size = 0.5 * window_size;
	const int max_peak_idx = num_peaks -1;

	if (num_peaks<=5)
	{
		indicators.resize(num_peaks,true);
		return false;
	}
	int i;
	for (i=0; i<5; i++)
	{
		if (peaks[i].scaled_intensity<=0)
			break;
	}

	const bool use_scaled_intensity = (i==5);
	int start_window_idx =0;

	indicators.resize(num_peaks,false);
	indicators[0]=true;
	indicators[max_peak_idx]=true;

	while (start_window_idx<max_peak_idx)
	{
		const mass_t max_window_mass = peaks[start_window_idx].mass + window_size;

		int end_window_idx=start_window_idx;
		while (end_window_idx<max_peak_idx && peaks[end_window_idx].mass<max_window_mass)
			end_window_idx++;


		if (end_window_idx - start_window_idx>num_peaks_per_window)
		{
			const int num_peaks_in_window = end_window_idx - start_window_idx+1;
			vector<peak_idx_pair> pairs;
			pairs.resize(num_peaks_in_window);

			if (use_scaled_intensity)
			{
				int i;
				for (i=0; i<num_peaks_in_window; i++)
				{
					const int peak_idx = i+start_window_idx;
					peak_idx_pair& pair = pairs[i];
					pair.idx = peak_idx ;
					pair.intensity = peaks[peak_idx].scaled_intensity;
				}
			}
			else
			{
				int i;
				for (i=0; i<num_peaks_in_window; i++)
				{
					const int peak_idx = i+start_window_idx;
					peak_idx_pair& pair = pairs[i];
					pair.idx = peak_idx ;
					pair.intensity = peaks[peak_idx].intensity;
				}
			}

			sort(pairs.begin(),pairs.end());

			if (pairs[0].intensity<pairs[1].intensity)
			{
				printf("Error: with peak intensity order (possible corruption in the files)!\n");
			//	int i;
			//	for (i=0; i<pairs.size(); i++)
			//		cout << i << " " << pairs[i].intensity << endl;
			//	exit(1);
				return false;
			}

			int i;
			for (i=0; i<num_peaks_per_window; i++)
				indicators[pairs[i].idx]=true;	
		}
		else 
		{
			int i;
			for (i=start_window_idx; i<=end_window_idx; i++)
				indicators[i]=true;
		}

		// advance half a window
		const mass_t mid_mass = peaks[start_window_idx].mass + half_window_size;
		start_window_idx++;
		while (start_window_idx<max_peak_idx && peaks[start_window_idx].mass<mid_mass)
			start_window_idx++;

	}

	return true;
}




void BasicSpectrum::output_to_mgf(ostream& mgf, const Config *config, const char *seq) const
{
	mgf << "BEGIN IONS" << endl;
	mgf << "TITLE=" <<  ssf->single_name << endl;
	
	if (ssf->peptide.get_num_aas()>0)
	{
		mgf << "SEQ=" << ssf->peptide.as_string(config) << endl;
	}
	else if (seq && strlen(seq)>2)
		mgf << "SEQ=" << seq << endl;
	
	if (ssf->type == MZXML)
	{
		MZXML_single *mzxml_single = (MZXML_single *)ssf;
		if (mzxml_single->scan_number>=0)
			mgf << "SCAN=" <<mzxml_single->scan_number << endl;

//		if (mzxml_single->retention_time>=0)
//			mgf << "RT=" << mzxml_single->retention_time << endl;
	}

	mgf << "CHARGE=+" << ssf->charge << endl;
		
	mgf << "PEPMASS=" << ssf->m_over_z << endl;
	
	int i;
	for (i=0; i<this->num_peaks; i++)
		mgf << fixed << setprecision(3) << peaks[i].mass << " " << peaks[i].intensity << endl;

	mgf << "END IONS" << endl << endl;
}






