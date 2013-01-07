#include "QuickClustering.h"







/***************************************************************************
	This function reads an mgf file into a basic spectra structure
	It allocates memory and stores the spectra in the supplied vector.
	returns false if unsuccsedul in read
****************************************************************************/
bool read_mgf_file_into_basic_spectra(Config *config, 
									  char *mgf_file, 
									  QCPeak *basic_peaks,
									  vector<BasicSpectrum>& basic_spectra)
{
	char buff[256];
	int max_num_spectra,max_num_lines;
	QCPeak *peak_storage;
	MGF_single *mgf_singles;

	examine_mgf_file(mgf_file,&max_num_spectra,&max_num_lines);

	mgf_singles = new MGF_single[max_num_spectra];
	basic_peaks = new QCPeak[max_num_lines -4*max_num_spectra];

//	cout << "Num spectra: " << max_num_spectra << endl;
//	cout << "Num lines:   " << max_num_lines << endl;

	peak_storage = basic_peaks;
	basic_spectra.clear();

	FILE *stream = fopen(mgf_file,"r");
	if (! stream)
	{
		cout << "Error: couldn't open file for reading: " << mgf_file << endl;
		exit(1);
	}



	int num_spectra=0;
	while (fgets(buff, 256, stream))
	{
		if (strncmp("BEGIN IONS",buff,10))
			continue;

		BasicSpectrum bs;
		MGF_single *ssf = &mgf_singles[num_spectra];
		bs.ssf = ssf;

		while (fgets(buff, 256, stream))
		{
			if (! strncmp(buff,"TITLE=",6) )
			{
				mass_t exp_pm;

				buff[strlen(buff)-1]='\0';
				ssf->single_name = buff + 6;
				
				if( ! fgets(buff, 1024, stream) )
					return false;

				if (! strncmp(buff,"SEQ=",4) )
				{
					string seq_string(buff+4);
					ssf->peptide.parseFromString(config,seq_string);

					if( ! fgets(buff, 1024, stream) )
						return false;
				}

				if (! strncmp(buff,"SCAN=",5) )
				{
					if (sscanf(buff+5,"%d",&ssf->scan_number) != 1 )
					{
						cout << "Error: couldn't read scan number!" << endl;
						exit(1);
					}
					if( ! fgets(buff, 1024, stream) )
						return false;
				}

				if (! strncmp(buff,"RT=",3) )
				{
					if (sscanf(buff+3,"%f",&ssf->retention_time) != 1)
					{
						cout << "Error: couldn't read retention_time!" << endl;
						exit(1);
					}
					if( ! fgets(buff, 1024, stream) )
						return false;
				}


				if (! strncmp(buff,"CLUSTER_SIZE=",13) )
				{
					if (sscanf(buff+13,"%d",&ssf->cluster_size) != 1 )
					{
						cout << "Error: couldn't read num single spectra!" << endl;
						exit(1);
					}
					if( ! fgets(buff, 1024, stream) )
						return false;
				}


				
				if ( sscanf(buff,"CHARGE=%d",&ssf->charge) != 1)
				{
					cout << "Error: couldn't read charge!" << endl;
					exit(1);
				}

				if( ! fgets(buff, 1024, stream) )
					return false;

				if (! strncmp(buff,"PEPMASS=",8) )
				{
					istringstream is(buff+8);
					is >> exp_pm;

					if (exp_pm<0)
					{
						cout << "Error: couldn't read mass!" << endl;
						exit(1);
					}
				}
				else
				{
					cout << "Error: couldn't read mass!" << endl;
					exit(1);
				}

				ssf->org_pm_with_19 = exp_pm * (ssf->charge) - (ssf->charge - 1)*MASS_PROTON;
				if (ssf->org_pm_with_19<0)
					ssf->org_pm_with_19 = -1;

				ssf->m_over_z = exp_pm;

				if (ssf->org_pm_with_19<5)
					ssf->org_pm_with_19 = ssf->m_over_z;

				int p_count = 0;
				while (fgets(buff, 256, stream))
				{
					istringstream is(buff);
					
					if (! strncmp("END IONS",buff,7))
						break;

					mass_t& mass = peak_storage[p_count].mass;
					intensity_t& intensity = peak_storage[p_count].intensity;
					is >> mass >> intensity;
				
					if (mass <0 || intensity==0)   // the peak probably got rounded down
						continue;

					p_count++;
				}

				ssf->num_peaks = p_count;
				bs.num_peaks =   p_count;
				bs.peaks = peak_storage;

				peak_storage +=  p_count;

				basic_spectra.push_back(bs);
			}
		}
	}

	return true;
}




//H293CoCl2-b1-total-try-500ug-a-2D34-120205-LTQ2-01.mzXML 10078 1 R.PVCDLAAD.A 0.048
void ann_mgf_and_create_mgf(Config *config, char *annotations_file, char *org_mgf_list_file,
							  char *out_dir_name, char *file_prefix, bool output_only_ann_spectra)
{
	FileManager fm;
	FileSet     all_spec_fs;
	vector< vector<int> >    annotation_idxs;
	vector<mzXML_annotation> annotations;  // use this structure should work the same

	cout << "Read anns..." << endl;

	read_mzXML_annotations(org_mgf_list_file,annotations_file, annotation_idxs, 
						   annotations, 50000);

	cout << "Read: " << annotations.size() << endl;

	config->set_need_to_normalize(0);

	// these are the lists and file names that are used
	if (output_only_ann_spectra)
	{
		fm.init_from_list_file(config,org_mgf_list_file,annotation_idxs);
	}
	else
		fm.init_from_list_file(config,org_mgf_list_file);


	char mgf_list_file[256];
	sprintf(mgf_list_file,"%s/%s_ann_list.txt",out_dir_name,file_prefix);

	all_spec_fs.select_all_files(fm,true);
//	all_spec_fs.sort_according_to_m_over_z();
	const int total_spectra = all_spec_fs.get_total_spectra();
	const vector<SingleSpectrumFile *>& all_ssf = all_spec_fs.get_ssf_pointers();

	cout << "MGF ssfs: " << all_ssf.size() << endl;


	// create mgf files

	FILE *list_stream = fopen(mgf_list_file,"w");
	if (! list_stream)
	{
		cout << "Error: couldn't open list stream!: " << mgf_list_file << endl;
		exit(1);
	}



	BasicSpecReader bsr;

	int max_size = 10000;

	
	int file_idx = 0;
	int s_idx;
	for (s_idx=0; s_idx<total_spectra; s_idx += max_size)
	{
		FileSet fs;
		int i;

		int end_idx = s_idx + max_size -1;
		if (end_idx>=total_spectra)
			end_idx = total_spectra-1;

		
		fs.init_from_another_fs(all_spec_fs,s_idx,end_idx);

		char mgf_name[256];
		sprintf(mgf_name,"%s%s_%d.mgf",out_dir_name,file_prefix,++file_idx);

		fstream mgf_stream(mgf_name,ios::out);
		if (! mgf_stream.good())
		{
			cout << "Error: couldn't open mgf for writing: " << mgf_name << endl;
			exit(1);
		}
		
		// read spectra
		static QCPeak peaks[5000];
		const vector<SingleSpectrumFile *>& ssfs = fs.get_ssf_pointers();
		for (i=0; i<ssfs.size(); i++)
		{
			BasicSpectrum bs;
			MGF_single *ssf = (MGF_single *)ssfs[i];

		//	cout << ssf->idx_in_file << " " << ssf->m_over_z << endl;

	
			int mgf_file_idx = ssf->file_idx;
			int ms2_scan_number = ssf->scan_number;
			int ann_idx = annotation_idxs[mgf_file_idx][ms2_scan_number];

			
			if (ann_idx >=0)
			{
				ssf->peptide.parseFromString(config, annotations[ann_idx].pep);
			//	cout << mgf_file_idx << "," << ssf->scan_number << " => " << ann_idx << " " <<
			//		annotations[ann_idx].pep << endl;
			}
			else
			{
				if (output_only_ann_spectra)
					continue;
			} 

			bs.peaks = peaks;
			bs.ssf = ssf;
			bs.ssf->charge = annotations[ann_idx].charge;
		
			ostringstream oss;
			oss << file_prefix << "." << ssf->file_idx << "." << ssf->scan_number;
			ssf->single_name = oss.str();

			

			bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		
			if (ssf->idx_in_file<0)
			{
				cout << "Error: no scan number read from MGF!!!" << endl;
				exit(1);
			}

			bs.output_to_mgf(mgf_stream,config);

		//	cout << bs.ssf->org_pm_with_19 - ssf->peptide.get_mass()- 19.0183 << endl;
			if (fabs(bs.ssf->org_pm_with_19 - ssf->peptide.get_mass()- MASS_OHHH)>5)
			{
				cout << "PM doesn't match peptide, pm_w_19 = " << bs.ssf->org_pm_with_19 << endl;
				cout << "Peptide mass                      = " << ssf->peptide.get_mass()+MASS_OHHH << endl;
				cout << ssf->file_idx << " " << ssf->scan_number << endl;
				exit(0);
			}
		}

		mgf_stream.close();

		fprintf(list_stream,"%s\n",mgf_name);
		fflush(list_stream);
	}
}






// checks that the anntoated spectra have a correct m_over_z
void print_specs(Config *config, char *list_name)
{
	FileManager fm;
	FileSet fs;

	vector<string> file_list;

	read_paths_into_list(list_name, file_list);
	fm.init_from_list(config,file_list);
	fs.select_all_files(fm);


	
	BasicSpecReader bsr;
	int i;
		// read spectra
	const vector<SingleSpectrumFile *>& ssfs = fs.get_ssf_pointers();
	for (i=0; i<ssfs.size(); i++)
	{
		static QCPeak peaks[5000];
		BasicSpectrum bs;
		MZXML_single *ssf = (MZXML_single *)ssfs[i];


		cout << ssf->scan_number << " " << ssf->num_peaks << " " << ssf->m_over_z << " " << ssf->num_peaks << endl;

	
		int mzXML_file_idx = ssfs[i]->file_idx;
		
		bs.peaks = peaks;
		bs.ssf = ssf;
	
		bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);

	//	cout << ssfs[i]->single_name << endl;
	//	cout << "m/z: " << ssfs[i]->m_over_z << endl;
	//	bs.print_peaks();
	}
	cout << endl;
		
}




