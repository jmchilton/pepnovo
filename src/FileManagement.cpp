#include "FileManagement.h"
#include "auxfun.h"

void read_paths_into_list(const char *list_file, vector<string>& list)
{
	ifstream fs(list_file,ios::in);
	if (! fs.good())
	{
		cout << "Error: couldn't open file for reading: " << list_file << endl;
		exit(1);
	}

	char buff[512];
	list.clear();

	while ( fs.getline(buff,512) )
	{
		string path(buff);
		if (path.length()>3)
		{
			list.push_back(path);
		}
	}
}


// parses the type from the file name, -1 if not recognized
int get_file_extension_type(const string& fname)
{
	int length = fname.length();

	int last_pos = length-1;
	while (last_pos>0 && (
		   fname[last_pos] == '\n' || fname[last_pos] == '\r' || 
		   fname[last_pos] == '\t' || fname[last_pos] == '\f') )
		last_pos--;

	if (fname[last_pos-2]=='d' && fname[last_pos-1]=='t' && 
		fname[last_pos  ]=='a')
		return DTA;

	if (fname[last_pos-2]=='m' && fname[last_pos-1]=='g' && 
		fname[last_pos  ]=='f')
		return MGF;

	if (fname[last_pos-4] == 'm' &&
	    fname[last_pos-3] == 'z' &&
		fname[last_pos-2]=='X' && fname[last_pos-1]=='M' && 
		fname[last_pos  ]=='L')
		return MZXML;

	if (fname[last_pos-2]=='d' && fname[last_pos-1]=='a' && 
		fname[last_pos  ]=='t')
		return DAT;

	if (fname[last_pos-2]=='m' && fname[last_pos-1]=='s' && 
		fname[last_pos  ]=='2')
		return MS2;

	if (fname[last_pos-2]=='t' && fname[last_pos-1]=='x' && 
		fname[last_pos  ]=='t')
		return TXT;


	return -1;
}

int SingleSpectrumFile::get_scan() const
{
	if (type == MGF)
	{
		MGF_single *mgf_single = (MGF_single *)this;
		if (mgf_single->scan_number>0)
			return mgf_single->scan_number;
		return mgf_single->idx_in_file;
	}
	if (type == MZXML)
		return ((MZXML_single *)this)->scan_number;
	if (type == DAT)
		return ((DAT_single *)this)->scan_number;

	return 0;
}


/****************************************************************
quickly extracts the charge, sequence, and pm_with_19 from dta
*****************************************************************/
bool DTA_file::scan_dta(const string& fname, const Config *config)
{
	char buff[1024];
	peptide.clear();

	ifstream fs(fname.c_str(),ios::in);
	if (! fs.good())
	{
		cout << "Error: couldn't open "<< fname << endl;
		exit(1);
	}
	
	fs.getline(buff,1024);

	while (buff[0] =='#') 
	{
		if (! strncmp(buff,"#SEQ",4))
			peptide.parseFromString(config,buff+5);


		if (! fs.getline(buff,1024))
		{
			cout << "Error scanning " << fname << endl;
			exit(1);
		}
		continue;
	}
	
	istringstream is(buff);

	org_pm_with_19 = -1;
	charge = -1;
	is >> org_pm_with_19 >> charge;
//	charge = 0; // bad fix

	if (org_pm_with_19<0 || charge <0)
	{
		cout << "Error: couldn't find pm and charge in DTA file, pm: "
			 << org_pm_with_19 << " , charge: " << charge << endl;
		exit(1);
	}

	if (charge>0)
	{
		m_over_z = (org_pm_with_19 + charge -1 ) / (mass_t)charge;
	}

	if (peptide.get_num_aas()>0 && charge>0)
	{
		mass_t diff = org_pm_with_19 - peptide.get_mass() - MASS_OHHH;
		if (fabs(diff)>6.0)
		{
			// try and correct charge!
			for (charge=1; charge<=3; charge++)
			{
				org_pm_with_19 =  m_over_z * charge + MASS_PROTON * (1 - charge);
				diff = org_pm_with_19 - peptide.get_mass() - MASS_OHHH;
				if (diff<5)
					break;
			}

			if (diff>7)
			{
				cout << "Error: sequence mass doesn't add up: " << fname << " (diff: "
					 << diff << ")" << endl;
				cout << "Pepitde: " << peptide.as_string(config) << endl;
				cout << "Mass Cys = " << config->get_session_tables().get_aa2mass(Cys) << endl;
			}
		//	exit(1);
		}
	}

	num_peaks=0;
	while (! fs.eof())
	{
		fs.getline(buff,64);
		float mass,inten;
		if (sscanf(buff,"%f %f",&mass,&inten)==2)
		{
			if (mass<0 || inten<=0)
				continue;
			num_peaks++;
		}
	}

	fs.close();
	return true;
}


/****************************************************************
quickly extracts the charge, sequence, and pm_with_19 from MGF location
*****************************************************************/
bool MGF_single::scan_mgf_single(FILE *stream, const Config *config)
{
	char buff[1024];
	mass_t exp_pm=-1;
	
	peptide.clear();

	
	while (fgets(buff, 256, stream))
	{
		file_pos = ftell(stream);
		if (strncmp(buff,"BEGIN IONS",10) )
		{
			
			continue;
		}
		break;
	}


	charge=0;
	m_over_z=-1;

		// read header info and first peak
	while (1)
	{
	

		if( ! fgets(buff, 1024, stream))
			return false;

		if (! strncmp(buff,"END IONS",7))
			return false;

		if (! strncmp(buff,"TITLE=",6) )
		{
			int len = strlen(buff)-1;
			buff[len]='\0';
			if (buff[len-1]=='\r' || buff[len-1]=='\n' )
				buff[len-1]='\0';
			single_name = buff + 6;

			// see if title includes scan number information.
			// this works only if the title ends with: .xxxx.yyyy.d.dta
			// e.g., MyMSMSData.2000.2000.2.dta
			len = single_name.length();
			if (len>7 &&
				single_name[len-1] == 'a' && single_name[len-2] == 't' && single_name[len-3]=='d' && 
				single_name[len-6] == '.' && single_name[len-4]== '.')
			{
				int pos = len -7;
				int num_dots =0;
				while (pos>0 && num_dots<2)
				{
					pos--;
					if (single_name[pos] == '.')
						num_dots++;
				}

				if (num_dots == 2)
				{
					string scan_str = single_name.substr(pos+1,len-7-pos);
					int i;
					for (i=0; i<scan_str.length(); i++)
						if (scan_str[i] == '.')
						{
							scan_str[i]=' ';
							break;
						}
					
					istringstream iss(scan_str);
					int scan1=-1,scan2=-1;
					iss >> scan1 >> scan2;
					if (scan1<=scan2 && scan1>0)
					{
						scan_number = scan1;
					//	cout << scan1 << "\t" << scan2 << "\t" << single_name << endl;
					}
				}
			}


			continue;
		}
		else
		if (! strncmp(buff,"SEQ=",4) )
		{
			string seq_string(buff+4);
			peptide.parseFromString(config,seq_string);
			peptide.calc_mass(config);
			continue;		
		}
		else
		if (! strncmp(buff,"SCAN=",5) )
		{
			if (sscanf(buff+5,"%d",&scan_number) != 1)
			{
				cout << "Error: couldn't read scan number!" << endl;
				exit(1);
			}
		continue;
		}
				else
		if (! strncmp(buff,"SCANS=",6) )
		{
			if (sscanf(buff+6,"%d",&scan_number) != 1)
			{
				cout << "Error: couldn't read scan number!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buff,"RT=",3) )
		{
			if (sscanf(buff+3,"%f",&retention_time) != 1)
			{
				cout << "Error: couldn't read retention_time!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buff,"RTINSECONDS=",12) )
		{
			if (sscanf(buff+12,"%f",&retention_time) != 1)
			{
				cout << "Error: couldn't read retention_time!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buff,"CLUSTER_SIZE=",13) )
		{
			if (sscanf(buff+13,"%d",&cluster_size) != 1)
			{
				cout << "Error: couldn't read cluster size!" << endl;
				exit(1);
			}
			continue;
		}
		else	
		if ( ! strncmp(buff,"CHARGE=",6))
		{
			if (sscanf(buff,"CHARGE=%d",&charge) != 1)
			{
				cout <<  "Error: couldn't read charge!" << endl;
				return false;
			}
		}
		else
		if (! strncmp(buff,"PEPMASS=",8))
		{
			istringstream is(buff+8);
			is >> m_over_z;
			
			if (m_over_z < 0)
			{
				cout << "Error: reading pepmass:" << m_over_z << endl;
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
				first_peak_mass = p.mass;
				break;
			}
		}
	}

	if (charge<0 || m_over_z<0)
		return false;

	org_pm_with_19 = m_over_z * charge + MASS_PROTON * (1 - charge);
	if (org_pm_with_19<0)
	{
		org_pm_with_19=-1;
	}


	num_peaks=0;
	while ( fgets(buff, 1024, stream) )
	{
		if (! strncmp(buff,"END IONS",8) )
			break;
		num_peaks++;
	}

	if (peptide.get_num_aas()>0 && charge>0)
	{
		mass_t diff = org_pm_with_19 - peptide.get_mass() - MASS_OHHH;
		if (fabs(diff)>7.0)
		{
			// try and correct charge!
			int c;
			for (c=1; c<=4; c++)
			{
				mass_t new_pm_with_19 =  m_over_z * c + MASS_PROTON * (1 - c);
				diff = fabs(new_pm_with_19 - peptide.get_mass() - MASS_OHHH);
				if (diff<10)
				{
					charge = c;
					org_pm_with_19 = new_pm_with_19;
					return true;
				}
			}
		

			cout << "Warning: MGF sequence mass doesn't add up: " << single_name << " (diff: "
				 << diff << ")" <<  endl;
			cout << "m/z: " << m_over_z << " org_pm_with_19: " << org_pm_with_19 << endl;
			cout << "Pepitde: " << peptide.as_string(config) << " (" << peptide.get_mass() << ")" << endl;
			cout << "Mass Cys = " << config->get_session_tables().get_aa2mass(Cys) << endl;
		//	org_pm_with_19 = -1;

			ind_peptide_mass_ok = false;
			return false;
		}
	}

	return true;
}




/****************************************************************
quickly extracts the charge, sequence, and pm_with_19 from MS2 location
*****************************************************************/
bool MS2_single::scan_ms2_single(FILE *stream, const Config *config)
{
	char buff[128];
	mass_t exp_pm=-1;
	
	peptide.clear();

	
	while (fgets(buff, 256, stream))
	{
		if (! strncmp(buff,":",1) )
			break;
	}

	if (feof(stream))
		return true;

	buff[strlen(buff)-1]='\0';
	single_name = buff+1;

	charge=-1;
	org_pm_with_19=-1;

	// read header info and first peak
	if( ! fgets(buff, 1024, stream))
		return false;

	
	if (! strncmp(buff,"SEQ=",4) )
	{
		string seq_string(buff+4);
		peptide.parseFromString(config,seq_string);
		if( ! fgets(buff, 1024, stream))
			return false;		
	}
	
	istringstream is(buff);
	is >> org_pm_with_19 >> charge;
	if (org_pm_with_19<0 || org_pm_with_19>10000)
	{
		cout << "Error: pm_with_19 in MS2 not in range: " << org_pm_with_19 << endl;
		exit(1);
	}

	if (charge<=0 || charge>10)
	{
		cout << "Error: bad charge state in MS2: " << charge << endl;
		exit(1);
	}


	m_over_z = (org_pm_with_19 -1 + charge) / charge;


	num_peaks=0;
	while ( fgets(buff, 128, stream) )
	{
		if (! strncmp(buff,":",1) || strlen(buff)< 3 )
			break;
		num_peaks++;
	}

	if (peptide.get_num_aas()>0 && charge>0)
	{
		mass_t diff = org_pm_with_19 - peptide.get_mass() - MASS_OHHH;
		if (fabs(diff)>7.0)
		{
			cout << "Error: sequence mass doesn't add up: " << single_name << " (diff: "
				 << diff << ")" <<  endl;
			cout << "m/z: " << m_over_z << " org_pm_with_19: " << org_pm_with_19 << endl;
			cout << "Pepitde: " << peptide.as_string(config) << " (" << peptide.get_mass() << ")" << endl;
			cout << "Mass Cys = " << config->get_session_tables().get_aa2mass(Cys) << endl;
			exit(1);
		}
	}

	return true;
}


/*******************************************************************
Special function to overcome mzMXL parsing problems.
Reads the peak lists into a local buff of floats.
********************************************************************/
void FileManager::init_and_read_single_mzXML(
						const Config *config, 
						const char * file_name, 
						int file_idx,
						mass_t min_m_over_z, 
						mass_t max_m_over_z)
{

	min_charge =9999;
	max_charge =0;

	num_spectra.resize(100,0);
    
	mzxml_files.resize(1);

	mzxml_files[0].mzxml_name = file_name;

	int mzxml_file_idx = (file_idx>=0 ? file_idx : 0);

	mzxml_files[0].extract_peak_lists_from_mzXML(config,mzxml_files[0].mzxml_name,mzxml_file_idx,0,10000);

	min_spec_mass = mzxml_files[0].min_spec_mass;
	max_spec_mass = mzxml_files[0].max_spec_mass;
	min_charge = mzxml_files[0].min_charge;
	max_charge = mzxml_files[0].max_charge;

	count_num_spectra();
//	print_summary_stats();
}


void FileManager::init_from_list(const Config *config, const vector<string>& list, bool quick_flag,
								 int file_idx)
{
	int i,num_dta=0,num_mgf=0,num_mzxml=0,num_dat=0,num_ms2=0;

	min_charge =9999;
	max_charge =0;

//	cout << "init from list: " << list.size() << endl;

	num_spectra.resize(100,0);
    
    for (i=0; i<list.size(); i++)
    {
		if (list[i][0] == '#')
			continue;

		int last_pos = list[i].length()-1;
	
		if (list[i][last_pos] == '\n' || list[i][last_pos] == '\r' || 
			list[i][last_pos] == '\t' || list[i][last_pos] == '\f')
			last_pos--;

		if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='t' && 
			list[i][last_pos  ]=='a')
		{
			num_dta++;
		}
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='g' && 
			list[i][last_pos  ]=='f')
		{
			num_mgf++;
		}
		else if (list[i][last_pos-4] == 'm' &&
				 list[i][last_pos-3] == 'z' &&
			list[i][last_pos-2]=='X' && list[i][last_pos-1]=='M' && 
			list[i][last_pos  ]=='L')
		{
			num_mzxml++;
		}
		else if (
			list[i][last_pos-2]=='d' && list[i][last_pos-1]=='a' && 
			list[i][last_pos  ]=='t')
		{
			num_dat++;
		}
		else if (
			list[i][last_pos-2]=='m' && list[i][last_pos-1]=='s' && 
			list[i][last_pos  ]=='2')
		{
			num_ms2++;
		}
		else
		{
			cout << "Error::: couldn't recognize file type for: " << list[i] << endl;
			exit(1);
		}
    }

	dta_files.resize(num_dta);
	mgf_files.resize(num_mgf);
	mzxml_files.resize(num_mzxml);
	dat_files.resize(num_dat);
	ms2_files.resize(num_ms2);

	int mgf_c=0;
	int dta_c=0;
	int mzxml_c=0;
	int dat_c=0;
	int ms2_c=0;

	for (i=0; i<list.size(); i++)
    {
		if (list[i][0] == '#')
			continue;

		int last_pos = list[i].length()-1;
		if (list[i][last_pos] == '\n' || list[i][last_pos] == '\r' || 
			list[i][last_pos] == '\t' || list[i][last_pos] == '\f')
			last_pos--;

		if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='t' && 
			list[i][last_pos  ]=='a')
		{
			dta_files[dta_c].single_name = list[i];

			if (quick_flag)
			{
				dta_files[dta_c].scan_dta(list[i],config);
			}
			else
				dta_files[dta_c].initial_read(config);

			dta_files[dta_c].type = DTA;
			dta_files[dta_c].file_idx =  dta_c;

			if (dta_files[dta_c].org_pm_with_19 < min_spec_mass)
				min_spec_mass = dta_files[dta_c].org_pm_with_19;

			if (dta_files[dta_c].org_pm_with_19 > max_spec_mass)
				max_spec_mass = dta_files[dta_c].org_pm_with_19;

			if (dta_files[dta_c].charge < min_charge)
				min_charge = dta_files[dta_c].charge;

			if (dta_files[dta_c].charge > max_charge)
				max_charge = dta_files[dta_c].charge;

			num_spectra[dta_files[dta_c].charge]++;

			dta_c++;
		}	
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='g' && 
			list[i][last_pos  ]=='f')
		{
			mgf_files[mgf_c].mgf_name =list[i];

//			cout << "DD1 reading idx " << mgf_c << endl;
			int mgf_file_idx = (file_idx>=0 ? mgf_c + file_idx : mgf_c);

			mgf_files[mgf_c].initial_read(config,mgf_file_idx, quick_flag);

			if (mgf_files[mgf_c].min_spec_mass<min_spec_mass)
				min_spec_mass = mgf_files[mgf_c].min_spec_mass;

			if (mgf_files[mgf_c].max_spec_mass>max_spec_mass)
				max_spec_mass = mgf_files[mgf_c].max_spec_mass;

			if (mgf_files[mgf_c].min_charge < min_charge)
				min_charge = mgf_files[mgf_c].min_charge;

			if (mgf_files[mgf_c].max_charge > max_charge)
				max_charge = mgf_files[mgf_c].max_charge;

			int c;
			for (c=0; c<=max_charge; c++)
				num_spectra[c] += mgf_files[mgf_c].num_spectra[c];

			mgf_c++;	
		}
		else if (list[i][last_pos-4] == 'm' &&
				 list[i][last_pos-3] == 'z' &&
				list[i][last_pos-2]=='X' && list[i][last_pos-1]=='M' && 
				list[i][last_pos  ]=='L')
		{
			mzxml_files[mzxml_c].mzxml_name =list[i];

			int mzxml_file_idx = (file_idx>=0 ? mzxml_c + file_idx : mzxml_c);

			mzxml_files[mzxml_c].initial_read(config,mzxml_c);

			if (mzxml_files[mzxml_c].min_spec_mass<min_spec_mass)
				min_spec_mass = mzxml_files[mzxml_c].min_spec_mass;

			if (mzxml_files[mzxml_c].max_spec_mass>max_spec_mass)
				max_spec_mass = mzxml_files[mzxml_c].max_spec_mass;

			if (mzxml_files[mzxml_c].min_charge < min_charge)
				min_charge = mzxml_files[mzxml_c].min_charge;

			if (mzxml_files[mzxml_c].max_charge > max_charge)
				max_charge = mzxml_files[mzxml_c].max_charge;

	//		int c;
	//		for (c=0; c<=max_charge; c++)
	//			num_spectra[c] += mzxml_files[mzxml_c].num_spectra[c];

			mzxml_c++;	
		}
		else if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='a' && 
				 list[i][last_pos  ]=='t')
		{
			dat_files[dat_c].dat_name =list[i];

			int dat_file_idx = (file_idx>=0 ? dat_c + file_idx : dat_c);

			dat_files[dat_c].initial_read(config,dat_file_idx);

			if (dat_files[dat_c].min_spec_mass<min_spec_mass)
				min_spec_mass = dat_files[dat_c].min_spec_mass;

			if (dat_files[dat_c].max_spec_mass>max_spec_mass)
				max_spec_mass = dat_files[dat_c].max_spec_mass;

			if (dat_files[dat_c].min_charge < min_charge)
				min_charge = dat_files[dat_c].min_charge;

			if (dat_files[dat_c].max_charge > max_charge)
				max_charge = dat_files[dat_c].max_charge;

			int c;
			for (c=0; c<=max_charge; c++)
				num_spectra[c] += dat_files[dat_c].num_spectra[c];

			dat_c++;	
		}
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='s' && 
			list[i][last_pos  ]=='2')
		{
			ms2_files[ms2_c].ms2_name =list[i];

			ms2_files[ms2_c].initial_read(config,ms2_c, quick_flag);

			if (ms2_files[ms2_c].min_spec_mass<min_spec_mass)
				min_spec_mass = ms2_files[ms2_c].min_spec_mass;

			if (ms2_files[ms2_c].max_spec_mass>max_spec_mass)
				max_spec_mass = ms2_files[ms2_c].max_spec_mass;

			if (ms2_files[ms2_c].min_charge < min_charge)
				min_charge = ms2_files[ms2_c].min_charge;

			if (mgf_files[ms2_c].max_charge > max_charge)
				max_charge = ms2_files[ms2_c].max_charge;

			int c;
			for (c=0; c<=max_charge; c++)
				num_spectra[c] += ms2_files[ms2_c].num_spectra[c];

			ms2_c++;	
		}
		else
		{
			cout << "Error :: couldn't recognize file type for: " << list[i] << endl;
			exit(1);
		}
    }

	count_num_spectra();
//	print_summary_stats();
}


/*****************************************************************
 returns how many spectra are present in the list file
 also samples m_over_z values to generate an approximate
 histogram in case the set of spectra needs to be spilt.
******************************************************************/
int FileManager::count_num_spectra(const Config *config, const char* list_file,
						  vector<mass_t>& mass_histogram) const
{
	vector<string> list;
	read_paths_into_list(list_file,list);

	int i,num_spectra_read=0;
	vector<mass_t> masses;
	int set_size = 1000;    // number of spectra that are examined for mass collection before fulsing

	masses.clear();
	mass_histogram.clear();

	for (i=0; i<list.size(); i++)
    {
		if (list[i][0] == '#')
			continue;

		int last_pos = list[i].length()-1;
		if (list[i][last_pos] == '\n' || list[i][last_pos] == '\r' || 
			list[i][last_pos] == '\t' || list[i][last_pos] == '\f')
			last_pos--;

		if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='t' && 
			list[i][last_pos  ]=='a')
		{
			DTA_file dta_file;
			masses.push_back(dta_file.m_over_z);
			num_spectra_read++;
		}	
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='g' && 
			list[i][last_pos  ]=='f')
		{
			MGF_file mgf_file;
			mgf_file.mgf_name =list[i];

			mgf_file.initial_read(config,0,true);

			int j;
			for (j=0; j<mgf_file.single_spectra.size(); j++)
				masses.push_back(mgf_file.single_spectra[j].m_over_z);

			num_spectra_read+=mgf_file.single_spectra.size();
		}
		else if (list[i][last_pos-4] == 'm' &&
				 list[i][last_pos-3] == 'z' &&
				list[i][last_pos-2]=='X' && list[i][last_pos-1]=='M' && 
				list[i][last_pos  ]=='L')
		{
			MZXML_file mzxml_file;
			mzxml_file.mzxml_name =list[i];

			mzxml_file.initial_read(config,0);

			int j;
			for (j=0; j<mzxml_file.single_spectra.size(); j++)
				masses.push_back(mzxml_file.single_spectra[j].m_over_z);

			num_spectra_read+=mzxml_file.single_spectra.size();
		}
		else if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='a' && 
				list[i][last_pos  ]=='t')
		{
			DAT_file dat_file;
			dat_file.dat_name =list[i];

			dat_file.initial_read(config,0);

			int j;
			for (j=0; j<dat_file.single_spectra.size(); j++)
				masses.push_back(dat_file.single_spectra[j].m_over_z);

			num_spectra_read+=dat_file.single_spectra.size();
		}
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='s' && 
				list[i][last_pos  ]=='2')
		{
			MS2_file ms2_file;
			ms2_file.ms2_name =list[i];

			ms2_file.initial_read(config,0);

			int j;
			for (j=0; j<ms2_file.single_spectra.size(); j++)
				masses.push_back(ms2_file.single_spectra[j].m_over_z);

			num_spectra_read+=ms2_file.single_spectra.size();
		}
		else
		{
			cout << "Error: couldn't recognize file type for: " << list[i] << endl;
			exit(1);
		}


		// if enough spectra were read, sample their masses and clear buffer
		if (masses.size()>set_size)
		{
			// sample masses per set size
			int num_sample = (int)((masses.size()/set_size)*10 + 0.5);
			int j;

			for (j=0; j<num_sample; j++)
			{
				int idx = (int)(myRandom() * masses.size());
				mass_histogram.push_back(masses[idx]);
			}
		
			masses.clear();
		}

		if (mass_histogram.size() == 100000)
		{
			vector<mass_t> tmp_his;
			int j;

			for (j=0; j<set_size; j++)
				if (myRandom()>=0.5)
					tmp_his.push_back(mass_histogram[j]);

			mass_histogram = tmp_his;
			tmp_his.clear();
			set_size *= 2;
		}
    }


	return num_spectra_read;
}


// counts how many total spectra are available from each charge
void FileManager::count_num_spectra()
{
	num_spectra.clear();
	num_spectra.resize(20,0);
	min_charge =99;
	max_charge =-1;


	int i;
	for (i=0; i<dta_files.size(); i++)
	{
		int charge =dta_files[i].charge;
		num_spectra[charge]++;
		if (charge<min_charge)
			min_charge = charge;
		if (charge>max_charge)
			max_charge = charge;
	}

	for (i=0; i<mgf_files.size(); i++)
	{
		int j;
		for (j=0; j<mgf_files[i].single_spectra.size(); j++)
		{
			int charge =mgf_files[i].single_spectra[j].charge;
			num_spectra[charge]++;
			if (charge<min_charge)
				min_charge = charge;
			if (charge>max_charge)
				max_charge = charge;
		}
	}

	for (i=0; i<mzxml_files.size(); i++)
	{
		int j;
		for (j=0; j<mzxml_files[i].single_spectra.size(); j++)
		{
			int charge =mzxml_files[i].single_spectra[j].charge;
			num_spectra[charge]++;
			if (charge<min_charge)
				min_charge = charge;
			if (charge>max_charge)
				max_charge = charge;
		}
	}

	for (i=0; i<dat_files.size(); i++)
	{
		int j;
		for (j=0; j<dat_files[i].single_spectra.size(); j++)
		{
			int charge =dat_files[i].single_spectra[j].charge;
			num_spectra[charge]++;
			if (charge<min_charge)
				min_charge = charge;
			if (charge>max_charge)
				max_charge = charge;
		}
	}

	total_num_spectra=0;
	for (i=0; i<num_spectra.size(); i++)
		total_num_spectra += num_spectra[i];

}



// Inits the FileManager using mass levels (for very large
// collections of spectra). This initialization uses a quick scan
void FileManager::init_from_list_file(const Config *config, const char* list_file,
		mass_t min_m_over_z, mass_t max_m_over_z)
{
	int i;
	vector<string> list;
	read_paths_into_list(list_file,list);

	num_spectra.resize(20,0);

	for (i=0; i<list.size(); i++)
    {
		if (list[i][0] == '#')
			continue;

		int last_pos = list[i].length()-1;
		if (list[i][last_pos] == '\n' || list[i][last_pos] == '\r' || 
			list[i][last_pos] == '\t' || list[i][last_pos] == '\f')
			last_pos--;

		if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='t' && 
			list[i][last_pos  ]=='a')
		{
			DTA_file dta_file;

			dta_file.scan_dta(list[i],config);

			if (dta_file.m_over_z>=min_m_over_z && dta_file.m_over_z<max_m_over_z)
				dta_files.push_back(dta_file);
			
		}	
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='g' && 
			list[i][last_pos  ]=='f')
		{
			MGF_file mgf_file;
			mgf_file.mgf_name =list[i];

			mgf_file.initial_read(config,mgf_files.size(),true);

			// change the single spectrum pointers in the mgf file record
			// to include only those that have a mass that is in the permitted range

			vector<MGF_single> good_singles;
			int j;
			for (j=0; j<mgf_file.single_spectra.size(); j++)
			{
				if (mgf_file.single_spectra[j].m_over_z >= min_m_over_z && 
					mgf_file.single_spectra[j].m_over_z <  max_m_over_z)
				{
					good_singles.push_back(mgf_file.single_spectra[j]);
				}
			}

			mgf_file.single_spectra = good_singles;
			mgf_files.push_back(mgf_file);
		}
		else if (list[i][last_pos-4] == 'm' &&
				 list[i][last_pos-3] == 'z' &&
				list[i][last_pos-2]=='X' && list[i][last_pos-1]=='M' && 
				list[i][last_pos  ]=='L')
		{
			MZXML_file mzxml;
			mzxml.mzxml_name =list[i];

			mzxml.initial_read(config,mzxml_files.size());

			// change the single spectrum pointers in the mgf file record
			// to include only those that have a mass that is in the permitted range

			vector<MZXML_single> good_singles;
			int j;
			for (j=0; j<mzxml.single_spectra.size(); j++)
			{
				if (mzxml.single_spectra[j].m_over_z >= min_m_over_z && 
					mzxml.single_spectra[j].m_over_z <  max_m_over_z)
				{
					good_singles.push_back(mzxml.single_spectra[j]);
				}
			}

			mzxml.single_spectra = good_singles;
			mzxml_files.push_back(mzxml);
		}
		else if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='a' && 
				list[i][last_pos  ]=='t')
		{
			DAT_file dat;
			dat.dat_name =list[i];

			dat.initial_read(config,dat_files.size());

			// change the single spectrum pointers in the mgf file record
			// to include only those that have a mass that is in the permitted range

			vector<DAT_single> good_singles;
			int j;
			for (j=0; j<dat.single_spectra.size(); j++)
			{
				if (dat.single_spectra[j].m_over_z >= min_m_over_z && 
					dat.single_spectra[j].m_over_z <  max_m_over_z)
				{
					good_singles.push_back(dat.single_spectra[j]);
				}
			}

			dat.single_spectra = good_singles;
			dat_files.push_back(dat);
		}
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='s' && 
			list[i][last_pos  ]=='2')
		{
			MS2_file ms2_file;
			ms2_file.ms2_name =list[i];

			ms2_file.initial_read(config,ms2_files.size(),true);

			// change the single spectrum pointers in the mgf file record
			// to include only those that have a mass that is in the permitted range

			vector<MS2_single> good_singles;
			int j;
			for (j=0; j<ms2_file.single_spectra.size(); j++)
			{
				if (ms2_file.single_spectra[j].m_over_z >= min_m_over_z && 
					ms2_file.single_spectra[j].m_over_z <  max_m_over_z)
				{
					good_singles.push_back(ms2_file.single_spectra[j]);
				}
			}

			ms2_file.single_spectra = good_singles;
			ms2_files.push_back(ms2_file);
		}
		else
		{
			cout << "Error: couldn't recognize file type for:: " << list[i] << endl;
			exit(1);
		}
	}

	count_num_spectra();
}


// Inits the FileManager using mass levels (for very large
// collections of spectra). This initialization uses a quick scan
// Only for files whose idx is true in file indicators
void FileManager::init_from_list_file(const Config *config, const char* list_file,
		const vector<bool>& file_indicators)
{
	int i;
	vector<string> list;
	read_paths_into_list(list_file,list);

	mass_t min_m_over_z=0;
	mass_t max_m_over_z=100000;

	num_spectra.resize(20,0);

	for (i=0; i<list.size(); i++)
    {
		if (list[i][0] == '#')
			continue;

		bool read_file = file_indicators[i];
		
		int last_pos = list[i].length()-1;
		if (list[i][last_pos] == '\n' || list[i][last_pos] == '\r' || 
			list[i][last_pos] == '\t' || list[i][last_pos] == '\f')
			last_pos--;

		if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='t' && 
			list[i][last_pos  ]=='a')
		{
			DTA_file dta_file;

			dta_file.scan_dta(list[i],config);

			if (dta_file.m_over_z>=min_m_over_z && dta_file.m_over_z<max_m_over_z)
				dta_files.push_back(dta_file);
			
		}	
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='g' && 
			list[i][last_pos  ]=='f')
		{
			MGF_file mgf_file;
			mgf_file.mgf_name =list[i];

			if (read_file)
			{
				mgf_file.initial_read(config,mgf_files.size(),true);

				// change the single spectrum pointers in the mgf file record
				// to include only those that have a mass that is in the permitted range

				vector<MGF_single> good_singles;
				int j;
				for (j=0; j<mgf_file.single_spectra.size(); j++)
				{
					if (mgf_file.single_spectra[j].m_over_z >= min_m_over_z && 
						mgf_file.single_spectra[j].m_over_z <  max_m_over_z)
					{
						good_singles.push_back(mgf_file.single_spectra[j]);
					}
				}

				mgf_file.single_spectra = good_singles;
			}
			else
				mgf_file.single_spectra.clear();

			mgf_files.push_back(mgf_file);
		}
		else if (list[i][last_pos-4] == 'm' &&
				 list[i][last_pos-3] == 'z' &&
				list[i][last_pos-2]=='X' && list[i][last_pos-1]=='M' && 
				list[i][last_pos  ]=='L')
		{
			MZXML_file mzxml;
			mzxml.mzxml_name =list[i];

			if (read_file)
			{

				mzxml.initial_read(config,mzxml_files.size());

				// change the single spectrum pointers in the mgf file record
				// to include only those that have a mass that is in the permitted range

				vector<MZXML_single> good_singles;
				int j;
				for (j=0; j<mzxml.single_spectra.size(); j++)
				{
					if (mzxml.single_spectra[j].m_over_z >= min_m_over_z && 
						mzxml.single_spectra[j].m_over_z <  max_m_over_z)
					{
						good_singles.push_back(mzxml.single_spectra[j]);
					}
				}

				mzxml.single_spectra = good_singles;
			}
			else
				mzxml.single_spectra.clear();

			mzxml_files.push_back(mzxml);
		}
		else if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='a' && 
				list[i][last_pos  ]=='t')
		{
			DAT_file dat;
			dat.dat_name =list[i];

			if (read_file)
			{

				dat.initial_read(config,dat_files.size());

				// change the single spectrum pointers in the mgf file record
				// to include only those that have a mass that is in the permitted range

				vector<DAT_single> good_singles;
				int j;
				for (j=0; j<dat.single_spectra.size(); j++)
				{
					if (dat.single_spectra[j].m_over_z >= min_m_over_z && 
						dat.single_spectra[j].m_over_z <  max_m_over_z)
					{
						good_singles.push_back(dat.single_spectra[j]);
					}
				}

				dat.single_spectra = good_singles;
			}
			else
				dat.single_spectra.clear();

			dat_files.push_back(dat);
		}
		else if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='s' && 
			list[i][last_pos  ]=='2')
		{
			MS2_file ms2_file;
			ms2_file.ms2_name =list[i];

			ms2_file.initial_read(config,ms2_files.size(),true);

			// change the single spectrum pointers in the mgf file record
			// to include only those that have a mass that is in the permitted range

			vector<MS2_single> good_singles;
			int j;
			for (j=0; j<ms2_file.single_spectra.size(); j++)
			{
				if (ms2_file.single_spectra[j].m_over_z >= min_m_over_z && 
					ms2_file.single_spectra[j].m_over_z <  max_m_over_z)
				{
					good_singles.push_back(ms2_file.single_spectra[j]);
				}
			}

			ms2_file.single_spectra = good_singles;
			ms2_files.push_back(ms2_file);
		}
		else
		{
			cout << "Error: couldn't recognize file type for:: " << list[i] << endl;
			exit(1);
		}
	}

	count_num_spectra();
}


// only keeps ssfs of mzXML singles that have an annotation
void FileManager::init_from_list_file(const Config *config, const char* list_file,
		const vector< vector<int> >& annotation_idxs)
{
	int i;
	vector<string> list;
	read_paths_into_list(list_file,list);

	mgf_files.clear();
	int total_dat_read=0;

	for (i=0; i<list.size(); i++)
    {
		if (list[i][0] == '#')
			continue;

		int last_pos = list[i].length()-1;

		if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='g' && 
				list[i][last_pos  ]=='f')
		{
			MGF_file mgf;
			mgf.mgf_name =list[i];

			bool has_anns=false;
			const vector<int>& anns = annotation_idxs[i];
			const int ann_size = anns.size();
			int j;
			for (j=0; j<ann_size; j++)
				if (anns[j]>=0)
					break;
			if (j<anns.size())
				has_anns = true;

			if (has_anns)
			{
//				cout << "DD2 reading idx " << mgf_files.size() << endl;

				mgf.initial_read(config,mgf_files.size(),true);

				cout << i <<  " " << mgf.mgf_name << " ";

				// change the single spectrum pointers in the mgf file record
				// to include only those that have a mass that is in the permitted range

				vector<MGF_single> good_singles;
				int j;
				for (j=0; j<mgf.single_spectra.size(); j++)
				{
					int mgf_file_idx = mgf.single_spectra[j].file_idx;
					int scan_number =  mgf.single_spectra[j].scan_number;

				

				//	cout << mgf_file_idx <<" " << scan_number << endl;

					if (mgf_file_idx < annotation_idxs.size() &&
						scan_number < annotation_idxs[mgf_file_idx].size() && 
						annotation_idxs[mgf_file_idx][scan_number]>=0) 
					{
						mgf.single_spectra[j].ann_idx=annotation_idxs[mgf_file_idx][scan_number];
						good_singles.push_back(mgf.single_spectra[j]);
					}
				}

				cout << good_singles.size() << " ..." << endl;
				mgf.single_spectra = good_singles;
			}
			mgf_files.push_back(mgf);
		}
		else if (list[i][last_pos-4] == 'm' &&
				 list[i][last_pos-3] == 'z' &&
				list[i][last_pos-2]=='X' && list[i][last_pos-1]=='M' && 
				list[i][last_pos  ]=='L')
		{
			MZXML_file mzxml;
			mzxml.mzxml_name =list[i];

			bool has_anns=false;
			const vector<int>& anns = annotation_idxs[i];
			const int ann_size = anns.size();
			int j;
			for (j=0; j<ann_size; j++)
				if (anns[j]>=0)
					break;
			if (j<anns.size())
				has_anns = true;

			if (has_anns)
			{
				mzxml.initial_read(config,i);

				vector<MZXML_single> good_singles;
				int j;
				for (j=0; j<mzxml.single_spectra.size(); j++)
				{
					int scan_num  = mzxml.single_spectra[j].scan_number;
					if (scan_num < ann_size &&  annotation_idxs[i][scan_num]>=0)
					{
						int ann_val = annotation_idxs[i][scan_num];
						mzxml.single_spectra[j].ann_idx=annotation_idxs[i][scan_num];
						good_singles.push_back(mzxml.single_spectra[j]);
					}
				}

				mzxml.single_spectra = good_singles;
			}

			mzxml_files.push_back(mzxml);
		}
		else if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='a' && 
				list[i][last_pos  ]=='t')
		{
			DAT_file dat;
			dat.dat_name =list[i];

			dat.initial_read(config,dat_files.size());

			cout << dat.dat_name << " ";

			// change the single spectrum pointers in the mgf file record
			// to include only those that have a mass that is in the permitted range

			vector<DAT_single> good_singles;
			int j;
			for (j=0; j<dat.single_spectra.size(); j++)
			{
				int mzxml_file_idx = dat.single_spectra[j].mzxml_file_idx;
				int scan_number = dat.single_spectra[j].scan_number;
				if (mzxml_file_idx < annotation_idxs.size() &&
					scan_number < annotation_idxs[mzxml_file_idx].size() && 
					annotation_idxs[mzxml_file_idx][scan_number]>=0) 
				{
					dat.single_spectra[j].ann_idx = annotation_idxs[mzxml_file_idx][scan_number];
					good_singles.push_back(dat.single_spectra[j]);
				}
			}

			cout << good_singles.size() << " ..." << endl;
			total_dat_read += good_singles.size();

			dat.single_spectra = good_singles;
			dat_files.push_back(dat);
		}
		else
		{
			cout << "Error: couldn't recognize file type for: " << list[i] << endl;
			exit(1);
		}
	}

	count_num_spectra();

	if (total_dat_read>0)
	{
		cout << "Read " << total_dat_read << " DAT spectra" << endl;
	}
}


// adds the annotations to the ssfs
void FileManager::init_from_list_file_and_add_annotations(const Config *config, const char* list_file,
		const vector< vector<int> >& annotation_idxs, vector<mzXML_annotation>& annotations,
		bool read_only_annotated )
{
	int i;
	vector<string> list;
	read_paths_into_list(list_file,list);

	mgf_files.clear();
	int total_dat_read=0;

	for (i=0; i<list.size(); i++)
    {
		if (list[i][0] == '#')
			continue;

		int last_pos = list[i].length()-1;

		if (list[i][last_pos-2]=='m' && list[i][last_pos-1]=='g' && 
				list[i][last_pos  ]=='f')
		{
			MGF_file mgf;
			mgf.mgf_name =list[i];

//			cout << "DD3 reading mgf idx : " << mgf_files.size() << endl;

			mgf.initial_read(config,mgf_files.size(),true);

			cout << mgf.mgf_name << " ";

			// change the single spectrum pointers in the mgf file record
			// to include only those that have a mass that is in the permitted range

			vector<MGF_single> good_singles;
			int j;
			for (j=0; j<mgf.single_spectra.size(); j++)
			{
				int mgf_file_idx = mgf.single_spectra[j].file_idx;
				int scan_number =  mgf.single_spectra[j].idx_in_file;

			

			//	cout << mgf_file_idx <<" " << scan_number << endl;

				if (mgf_file_idx < annotation_idxs.size() &&
					scan_number < annotation_idxs[mgf_file_idx].size() && 
					annotation_idxs[mgf_file_idx][scan_number]>=0) 
				{
					mgf.single_spectra[j].peptide.parseFromString(config,
						annotations[annotation_idxs[mgf_file_idx][scan_number]].pep);
					mgf.single_spectra[j].charge = annotations[annotation_idxs[mgf_file_idx][scan_number]].charge;
					good_singles.push_back(mgf.single_spectra[j]);
				}
				else
					if (! read_only_annotated)
						good_singles.push_back(mgf.single_spectra[j]);
			}

			cout << good_singles.size() << " ..." << endl;
			mgf.single_spectra = good_singles;
			mgf_files.push_back(mgf);
		}
		else if (list[i][last_pos-4] == 'm' &&
				 list[i][last_pos-3] == 'z' &&
				list[i][last_pos-2]=='X' && list[i][last_pos-1]=='M' && 
				list[i][last_pos  ]=='L')
		{
			MZXML_file mzxml;
			mzxml.mzxml_name =list[i];

			cout << i << " " << list[i] << endl;

			if (i>= annotation_idxs.size() && ! read_only_annotated)
			{
				mzxml_files.push_back(mzxml);
				continue;
			}

			bool has_anns=false;
			const vector<int>& anns = annotation_idxs[i];
			const int ann_size = anns.size();
			int j;
			for (j=0; j<ann_size; j++)
				if (anns[j]>=0)
					break;
			if (j<anns.size())
				has_anns = true;

			if (has_anns)
			{
				mzxml.initial_read(config,i);

				vector<MZXML_single> good_singles;
				int j;
				for (j=0; j<mzxml.single_spectra.size(); j++)
				{
					int scan_num  = mzxml.single_spectra[j].scan_number;
					if (scan_num < ann_size &&  annotation_idxs[i][scan_num]>=0)
					{
						int ann_val = annotation_idxs[i][scan_num];
						mzxml.single_spectra[j].peptide.parseFromString(config,
							annotations[ann_val].pep);
						mzxml.single_spectra[j].charge = annotations[ann_val].charge;
						good_singles.push_back(mzxml.single_spectra[j]);
					}
					else
						if (! read_only_annotated)
							good_singles.push_back(mzxml.single_spectra[j]);
				}

				mzxml.single_spectra = good_singles;
			}

			mzxml_files.push_back(mzxml);
		}
		else if (list[i][last_pos-2]=='d' && list[i][last_pos-1]=='a' && 
				list[i][last_pos  ]=='t')
		{
			DAT_file dat;
			dat.dat_name =list[i];

			dat.initial_read(config,dat_files.size());

			cout << dat.dat_name << " ";

			// change the single spectrum pointers in the mgf file record
			// to include only those that have a mass that is in the permitted range

			vector<DAT_single> good_singles;
			int j;
			for (j=0; j<dat.single_spectra.size(); j++)
			{
				int mzxml_file_idx = dat.single_spectra[j].mzxml_file_idx;
				int scan_number = dat.single_spectra[j].scan_number;
				if (mzxml_file_idx < annotation_idxs.size() &&
					scan_number < annotation_idxs[mzxml_file_idx].size() && 
					annotation_idxs[mzxml_file_idx][scan_number]>=0) 
				{
					int ann_val = annotation_idxs[mzxml_file_idx][scan_number];
					dat.single_spectra[j].peptide.parseFromString(config,
						annotations[ann_val].pep);
					good_singles.push_back(dat.single_spectra[j]);
				}
				else
					if (! read_only_annotated)
						good_singles.push_back(dat.single_spectra[j]);
			}

			cout << good_singles.size() << " ..." << endl;
			total_dat_read += good_singles.size();

			dat.single_spectra = good_singles;
			dat_files.push_back(dat);
		}
		else
		{
			cout << "Error: couldn't recognize file type for: " << list[i] << endl;
			exit(1);
		}
	}

	count_num_spectra();

	if (total_dat_read>0)
	{
		cout << "Read " << total_dat_read << " DAT spectra" << endl;
	}
}







void DTA_file::initial_read(const Config *config)
{

	if (single_name.length() == 0)
	{
		cout << "Error: must first copy name to DTA_file!" << endl;
		exit(1);
	}

	Spectrum s;

	s.read_from_dta(config,single_name.c_str());
	s.init_spectrum();

	this->type = DTA;
	this->single_name = single_name;
	this->peptide = s.getPeptide();
	this->charge = s.getCharge();
	this->org_pm_with_19 = s.get_org_pm_with_19();
	this->pm_with_19 = s.get_corrected_pm_with_19();
	this->num_peaks = s.getNumPeaks();
	this->m_over_z = s.get_m_over_z();	
}



void MGF_file::initial_read(const Config *config, int file_idx, bool quick_flag)
{
	if (mgf_name.length() == 0)
	{
		printf("Error: must first copy name to MGF_file!\n");
		exit(1);
	}

	FILE *stream=fopen(mgf_name.c_str(),"r");
	if (! stream)
	{
		cout << "Error: couldn't open mgf file for reading: " << mgf_name << endl;
		exit(1);
	}

//	cout << "TTT Reading mgf idx: " << file_idx << endl;

	this->single_spectra.clear();
	this->single_spectra.reserve(5000);

	this->num_spectra.resize(100,0);
	
	Spectrum s;
	int counter=0;
	while (1)
	{
	
		MGF_single mg;
		
		long pos = ftell(stream);  // remember stream pointer

		mg.type = MGF;
		mg.file_pos = -1;
		mg.file_idx = file_idx;
		mg.idx_in_file = counter++;

		if (quick_flag)
		{
			
			if (! mg.scan_mgf_single(stream,config) )
			{
				static char tmp_buff[256];
				if (! fgets(tmp_buff,256,stream)) // move forwards, there was a problem with that position
					break;						  // reached eof
				
				continue;
			}

			// there was something wrong with this spectrum
			if (mg.org_pm_with_19<0)
				continue;

			if (mg.org_pm_with_19 < min_spec_mass)
				min_spec_mass = mg.org_pm_with_19;

			if (mg.org_pm_with_19 > max_spec_mass)
				max_spec_mass = mg.org_pm_with_19;

			if (mg.charge <min_charge)
				min_charge = mg.charge;

			if (mg.charge >max_charge)
				max_charge = mg.charge;

			num_spectra[mg.charge]++;
		}
		else
		{
			mg.file_pos = ftell(stream);

			if (! s.read_from_MGF_stream(config,stream) )
				break;

			s.initializeSpectrum(static_cast<void*>(&mg));

			mg.peptide = s.getPeptide();
			mg.single_name = s.getTitle();
			mg.org_pm_with_19 = s.get_org_pm_with_19();
			mg.charge = s.getCharge();
			mg.pm_with_19 = (s.get_corrected_pm_with_19()>0 ?s.get_corrected_pm_with_19() : s.get_org_pm_with_19());
			mg.num_peaks = s.getNumPeaks();
			
			if (mg.org_pm_with_19 < min_spec_mass)
				min_spec_mass = mg.org_pm_with_19;

			if (mg.org_pm_with_19 > max_spec_mass)
				max_spec_mass = mg.org_pm_with_19;

			if (mg.charge <min_charge)
				min_charge = mg.charge;

			if (mg.charge >max_charge)
				max_charge = mg.charge;

			num_spectra[mg.charge]++;
		}
	
		if (mg.file_pos>=0)
		{
			single_spectra.push_back(mg);
		}
		
	}
	fclose(stream);
	
}




void DAT_file::initial_read(const Config *config, int file_idx)
{
	const int bytes_to_read = sizeof(mass_t) + 4 * sizeof(int) + 2 *sizeof(float);

	if (dat_name.length() == 0)
	{
		printf("Error: must first copy name to DAT_file!\n");
		exit(1);
	}

	ifstream dat_stream(dat_name.c_str(),ios::in | ios::binary);

	if (! dat_stream.is_open() || ! dat_stream.good())
	{
		cout << "Error: couldn't open dat file for reading: " << dat_name << endl;
		exit(1);
	}

	this->single_spectra.clear();
	this->single_spectra.reserve(20000);
	this->num_spectra.resize(100,0);
	

	int counter=0;
	while (1)
	{
	
		char header_buff[bytes_to_read];
		long pos = dat_stream.tellg();
		DAT_single dat;
		

		dat.type = DAT;
		dat.file_pos = pos;
		dat.file_idx = file_idx;

		dat_stream.read(header_buff,bytes_to_read);
		if (dat_stream.gcount() != bytes_to_read)
			break;
		
		char *b_pos = header_buff;
		
		dat.m_over_z = *(mass_t *)b_pos;
		b_pos += sizeof(mass_t);
		
		dat.charge = *(int *)b_pos;
		b_pos += sizeof(int);

		dat.mzxml_file_idx = *(int *)b_pos;
		b_pos += sizeof(int);

		dat.scan_number = *(int *)b_pos;
		b_pos += sizeof(int);

		dat.num_peaks = *(int *)b_pos;
		b_pos += sizeof(int);

		dat.retention_time = *(float *)b_pos;
		b_pos += sizeof(float);

		dat.precursor_intensity = *(float *)b_pos;
		b_pos += sizeof(float);

		// sanity checks
		if (dat.m_over_z<0.0 || dat.m_over_z>10000.0 || 
			dat.mzxml_file_idx<0 || dat.mzxml_file_idx>100000 ||
			dat.num_peaks<0  || dat.num_peaks>100000 ||
			dat.scan_number<0 || dat.scan_number>1000000)
		{
			cout << "Error in DAT file " << file_idx << " pos " << pos << "  #" << 
				single_spectra.size() << endl;
			
			cout << "m/z = " << dat.m_over_z << endl;
			cout << "mzxml idx = " << dat.mzxml_file_idx << endl;
			cout << "scan =  " << dat.scan_number << endl;
			cout << "num_peaks = " << dat.num_peaks << endl;
			exit(1);
		}

		dat_stream.seekg(pos  + bytes_to_read + 2*sizeof(float)*dat.num_peaks);

		single_spectra.push_back(dat);
		counter++;
	}

	dat_stream.close();
}



void MS2_file::initial_read(const Config *config, int file_idx, bool quick_flag)
{
	if (ms2_name.length() == 0)
	{
		printf("Error: must first copy name to MS2_file!\n");
		exit(1);
	}

	FILE *stream=fopen(ms2_name.c_str(),"r");
	if (! stream)
	{
		cout << "Error: couldn't open mgf file for reading: " << ms2_name << endl;
		exit(1);
	}

	this->single_spectra.clear();
	this->single_spectra.reserve(5000);

	this->num_spectra.resize(100,0);
	
	Spectrum s;
	int counter=0;
	while (1)
	{
	
		MS2_single ms2;
		
		long pos = ftell(stream);  // remember stream pointer

		ms2.type = MS2;
		ms2.file_pos = pos;
		ms2.file_idx = file_idx;
		ms2.idx_in_file = counter++;

		if (quick_flag)
		{
			
			if (! ms2.scan_ms2_single(stream,config) )
				break;

			if (ms2.org_pm_with_19 < min_spec_mass)
				min_spec_mass = ms2.org_pm_with_19;

			if (ms2.org_pm_with_19 > max_spec_mass)
				max_spec_mass = ms2.org_pm_with_19;

			if (ms2.charge <min_charge)
				min_charge = ms2.charge;

			if (ms2.charge >max_charge)
				max_charge = ms2.charge;

			num_spectra[ms2.charge]++;
		}
		else
		{
			cout << "Reading MS2 not supported in this mode!" << endl;
			exit(0);
		/*	if (! s.read_from_MS2_stream(config,stream) )
				break;

			s.init_spectrum();

			ms2.peptide = s.get_peptide();
			ms2.single_name = s.get_file_name();
			ms2.org_pm_with_19 = s.get_org_pm_with_19();
			ms2.charge = s.get_charge();
			ms2.pm_with_19 = s.get_corrected_pm_with_19();
			ms2.num_peaks = s.get_num_peaks();
			
			if (ms2.org_pm_with_19 < min_spec_mass)
				min_spec_mass = ms2.org_pm_with_19;

			if (ms2.org_pm_with_19 > max_spec_mass)
				max_spec_mass = ms2.org_pm_with_19;

			if (ms2.charge <min_charge)
				min_charge = ms2.charge;

			if (ms2.charge >max_charge)
				max_charge = ms2.charge;

			num_spectra[ms2.charge]++;*/
		}
	
		single_spectra.push_back(ms2);
		
	}
	fclose(stream);
	
}




// reads the summary tsv file and stores stats
// extract all information from the tsv file
void PKL_dir::initial_read(const Config *config, int dir_idx, const string& path, 
						   const string& tsv_file , mass_t min_m_over_z, mass_t max_m_over_z)
{
	this->dir_path = path;
	this->tsv_path = tsv_file;

	FILE *tsv_stream = fopen(tsv_file.c_str(),"r");

	if (! tsv_stream)
	{
		cout << "Error: couldn't open tsv file for reading: " << tsv_file << endl;
		exit(1);
	}

	char line_buff[2048];
	fgets(line_buff,2048,tsv_stream); // first row is headers

	while (fgets(line_buff,2048,tsv_stream))
	{
		int np=-1;
		PKL_single pkl;

		pkl.type = PKL;
		pkl.file_idx = dir_idx;
		pkl.m_over_z = -1;
		pkl.single_name = "";
		pkl.num_peaks=-1;

		istringstream is(line_buff);

		mass_t mz_sel,mz_cen,parent_MH;
		int tmp_charge, ns;
		float max_inten;

		is >> pkl.single_name >> pkl.scan_number >> pkl.retention_time >> np >> 
			pkl.num_peaks >> mz_sel >> pkl.m_over_z >> mz_cen >> parent_MH >> tmp_charge >> ns >>
			max_inten >> pkl.precursor_intensity;

		if (pkl.num_peaks<=5)
			continue;

		if (pkl.m_over_z < min_m_over_z || pkl.m_over_z>max_m_over_z)
			continue;


		// set charge from file name
		pkl.charge = 0;
		if (pkl.single_name.length()>4)
		{
			char charge_sym = pkl.single_name[pkl.single_name.length()-5];
			if (charge_sym<'0' || charge_sym>'9')
			{
				cout << "Error: couldn't extract charge from file name " << pkl.single_name << endl;
				exit(1);
			}
			pkl.charge = int(charge_sym - '0');
		}

		single_spectra.push_back(pkl);
			
	}

	fclose(tsv_stream);
}



// reads a list with dirs and paths to tsv files
void FileManager::init_from_pkl_dir_list(const Config *config, const char *list_file, 
										 mass_t min_m_over_z, mass_t max_m_over_z)
{
	FILE *list_stream=fopen(list_file,"r");

	if (! list_stream)
	{
		cout << "Error: couldn't open pkl dir list for reading: " << list_file << endl;
		exit(1);
	}


	char line_buff[1024];
	int n=0;
	while (fgets(line_buff,1024,list_stream))
		n++;

	fclose(list_stream);

	pkl_dirs.resize(n);

	list_stream=fopen(list_file,"r");


	cout << "Read pkl files (idx path #spectra):" << endl;
	int pkl_dir_idx=0;
	while (fgets(line_buff,1024,list_stream))
	{
		istringstream is(line_buff);

		string pkl_dir = "";
		string tsv_file = "";

		is >> pkl_dir >> tsv_file;

		if (pkl_dir.length()>2 && tsv_file.length() >2)
		{
			pkl_dirs[pkl_dir_idx].initial_read(config,pkl_dir_idx,pkl_dir,tsv_file,
											   min_m_over_z, max_m_over_z);
			
			if (pkl_dirs[pkl_dir_idx].single_spectra.size()>0)
			{
				cout << pkl_dir_idx << "\t" << pkl_dir << "\t" << pkl_dirs[pkl_dir_idx].single_spectra.size() << endl;

				pkl_dir_idx++;
			}
		}
	}	
}


void FileManager::init_from_single_pkl_dir(const Config *config, const string& pkl_dir_path, 
					const string& tsv_file, int pkl_dir_idx, mass_t min_m_over_z, mass_t max_m_over_z)
{

	pkl_dirs.resize(1);
	pkl_dirs[0].initial_read(config, pkl_dir_idx, pkl_dir_path, tsv_file, 
							 min_m_over_z, max_m_over_z);

}



/*************************************************************************
Select all file pointers from the FileManager
**************************************************************************/
void FileSet::select_all_files(const FileManager& fm, bool remove_duplicates)
{
	int i;

	this->ssf_pointers.clear();

	for (i=0; i<fm.dta_files.size(); i++)
	{
		ssf_pointers.push_back((struct DTA_file *)&fm.dta_files[i]);

		// check that this doesn't correpond to the immediate previous spec
		if (remove_duplicates && ssf_pointers.size()>1)
		{
			const SingleSpectrumFile * current_ssf = ssf_pointers[ssf_pointers.size()-1];
			const SingleSpectrumFile * previous_ssf = ssf_pointers[ssf_pointers.size()-2];

			if ( ( previous_ssf->num_peaks == current_ssf->num_peaks ) &&
				 ( fabs(previous_ssf->m_over_z-current_ssf->m_over_z) < 0.005) )
			{
				ssf_pointers.pop_back();
			}
		}
	}

	for (i=0; i<fm.mgf_files.size(); i++)
	{
		int j;
		for (j=0; j<fm.mgf_files[i].single_spectra.size(); j++)
		{
			const SingleSpectrumFile *ssf = &fm.mgf_files[i].single_spectra[j];

			if (ssf->num_peaks<5)
				continue;

			this->ssf_pointers.push_back((struct MGF_single *)ssf);

			// check that this doesn't correpond to the immediate previous spec
			if (remove_duplicates && ssf_pointers.size()>1)
			{
				const SingleSpectrumFile * current_ssf = ssf_pointers[ssf_pointers.size()-1];
				const SingleSpectrumFile * previous_ssf = ssf_pointers[ssf_pointers.size()-2];

				if ( ( previous_ssf->num_peaks == current_ssf->num_peaks ) &&
				 ( fabs(previous_ssf->m_over_z-current_ssf->m_over_z) < 0.005) )
				{
					ssf_pointers.pop_back();
					continue;
				}
			}

//			cout << setprecision(9)<< ssf->m_over_z << " " << ssf->num_peaks << endl;
		}
	}

	for (i=0; i<fm.mzxml_files.size(); i++)
	{
		int j;
		for (j=0; j<fm.mzxml_files[i].single_spectra.size(); j++)
		{
			const SingleSpectrumFile *ssf = &fm.mzxml_files[i].single_spectra[j];

			if (ssf->num_peaks<5)
				continue;

			this->ssf_pointers.push_back((struct MZXML_single *)ssf);

			// check that this doesn't correpond to the immediate previous spec
			if (remove_duplicates && ssf_pointers.size()>1)
			{
				const SingleSpectrumFile * current_ssf = ssf_pointers[ssf_pointers.size()-1];
				const SingleSpectrumFile * previous_ssf = ssf_pointers[ssf_pointers.size()-2];

				if ( ( previous_ssf->num_peaks == current_ssf->num_peaks ) &&
				 ( fabs(previous_ssf->m_over_z-current_ssf->m_over_z) < 0.005) )
				{
					ssf_pointers.pop_back();
					continue;
				}
			}

//			cout << setprecision(9)<< ssf->m_over_z << " " << ssf->num_peaks << endl;
		}
	}

	for (i=0; i<fm.dat_files.size(); i++)
	{
		int j;
		for (j=0; j<fm.dat_files[i].single_spectra.size(); j++)
		{
			const SingleSpectrumFile *ssf = &fm.dat_files[i].single_spectra[j];

		//	if (ssf->num_peaks<5)
		//		continue;

			this->ssf_pointers.push_back((struct DAT_single *)ssf);

			// check that this doesn't correpond to the immediate previous spec
			if (remove_duplicates && ssf_pointers.size()>1)
			{
				const SingleSpectrumFile * current_ssf = ssf_pointers[ssf_pointers.size()-1];
				const SingleSpectrumFile * previous_ssf = ssf_pointers[ssf_pointers.size()-2];

				if ( ( previous_ssf->num_peaks == current_ssf->num_peaks ) &&
				 ( fabs(previous_ssf->m_over_z-current_ssf->m_over_z) < 0.005) )
				{
					ssf_pointers.pop_back();
					continue;
				}
			}

//			cout << setprecision(9)<< ssf->m_over_z << " " << ssf->num_peaks << endl;
		}
	}

	for (i=0; i<fm.pkl_dirs.size(); i++)
	{
		int j;
		for (j=0; j<fm.pkl_dirs[i].single_spectra.size(); j++)
		{
			const SingleSpectrumFile *ssf = &fm.pkl_dirs[i].single_spectra[j];

			if (ssf->num_peaks<5)
				continue;

			this->ssf_pointers.push_back((struct PKL_single *)ssf);

			// check that this doesn't correpond to the immediate previous spec
			if (remove_duplicates && ssf_pointers.size()>1)
			{
				const SingleSpectrumFile * current_ssf = ssf_pointers[ssf_pointers.size()-1];
				const SingleSpectrumFile * previous_ssf = ssf_pointers[ssf_pointers.size()-2];

				if ( ( previous_ssf->num_peaks == current_ssf->num_peaks ) &&
				 ( fabs(previous_ssf->m_over_z-current_ssf->m_over_z) < 0.001) )
				{
					ssf_pointers.pop_back();
					continue;
				}
			}

//			cout << setprecision(9)<< ssf->m_over_z << " " << ssf->num_peaks << endl;
		}
	}


	this->next_ssf_pointer=0;
}


/**********************************************************************
Selects dat ssfs that have an mzxml file idx of at most max_mzxml idx
***********************************************************************/
void FileSet::filter_dat_spectra_by_mzxml_idx(int max_mzxml_idx)
{
	int i;
	vector<SingleSpectrumFile *> new_ssfs;
	
	for (i=0; i<ssf_pointers.size(); i++)
	{
		DAT_single *dat_ssf = (DAT_single *)ssf_pointers[i];
		if (dat_ssf->mzxml_file_idx > max_mzxml_idx)
			continue;
		new_ssfs.push_back(ssf_pointers[i]);
	}
	ssf_pointers = new_ssfs;
}



/*************************************************************************
Select the file pointers from the FileManager
**************************************************************************/
void FileSet::select_files(const FileManager& fm, 
						   mass_t min_pm_with_19, 
						   mass_t max_pm_with_19, 
						   score_t min_sqs,       
						   score_t max_sqs, 
						   int charge,            
						   bool only_unassigned)
{
	int i;

	this->ssf_pointers.clear();

		
	if (charge !=0)
	{
		min_charge = charge;
		max_charge = charge;
	}

	for (i=0; i<fm.dta_files.size(); i++)
	{
		if (fm.dta_files[i].org_pm_with_19 < min_pm_with_19 ||
			fm.dta_files[i].org_pm_with_19 > max_pm_with_19)
			continue;

		if (fm.dta_files[i].sqs>-1 && 
			(fm.dta_files[i].sqs<min_sqs || fm.dta_files[i].sqs>max_sqs) )
			continue;

		if (charge>0 && fm.dta_files[i].charge != charge)
			continue;

		if (only_unassigned && fm.dta_files[i].assigned_cluster>=0)
			continue;

		if (fm.dta_files[i].charge<min_charge)
			min_charge =  fm.dta_files[i].charge;

		if (fm.dta_files[i].charge>max_charge)
			max_charge =  fm.dta_files[i].charge;

		ssf_pointers.push_back((struct DTA_file *)&fm.dta_files[i]);
	}


	for (i=0; i<fm.mgf_files.size(); i++)
	{
		int j;
		for (j=0; j<fm.mgf_files[i].single_spectra.size(); j++)
		{
			if (fm.mgf_files[i].single_spectra[j].org_pm_with_19 < min_pm_with_19 ||
				fm.mgf_files[i].single_spectra[j].org_pm_with_19 > max_pm_with_19)
				continue;

			if ( fm.mgf_files[i].single_spectra[j].sqs > -1 &&
				(fm.mgf_files[i].single_spectra[j].sqs<min_sqs || fm.mgf_files[i].single_spectra[j].sqs>max_sqs))
				continue;

			if (charge>0 && fm.mgf_files[i].single_spectra[j].charge != charge)
				continue;

			if (only_unassigned && fm.mgf_files[i].single_spectra[j].assigned_cluster>=0)
				continue;

			if (fm.mgf_files[i].single_spectra[j].charge>max_charge)
				max_charge = fm.mgf_files[i].single_spectra[j].charge;
			
			if (fm.mgf_files[i].single_spectra[j].charge<min_charge)
				min_charge = fm.mgf_files[i].single_spectra[j].charge;

			this->ssf_pointers.push_back((struct MGF_single *)&fm.mgf_files[i].single_spectra[j]);
		}
	}
	this->next_ssf_pointer=0;
}



/*************************************************************************
Select the file pointers from the FileManager
**************************************************************************/
void FileSet::select_files_in_mz_range(
						   const FileManager& fm, 
						   mass_t min_mz, 
						   mass_t max_mz, 
						   int charge)
{
	int i;

	this->ssf_pointers.clear();

		
	if (charge !=0)
	{
		min_charge = charge;
		max_charge = charge;
	}


	for (i=0; i<fm.mgf_files.size(); i++)
	{
		int j;
		for (j=0; j<fm.mgf_files[i].single_spectra.size(); j++)
		{
			if (fm.mgf_files[i].single_spectra[j].m_over_z < min_mz ||
				fm.mgf_files[i].single_spectra[j].m_over_z > max_mz)
				continue;

			if (charge>0 && fm.mgf_files[i].single_spectra[j].charge != charge)
				continue;

			if (fm.mgf_files[i].single_spectra[j].charge>max_charge)
				max_charge = fm.mgf_files[i].single_spectra[j].charge;
			
			if (fm.mgf_files[i].single_spectra[j].charge<min_charge)
				min_charge = fm.mgf_files[i].single_spectra[j].charge;

			this->ssf_pointers.push_back((struct MGF_single *)&fm.mgf_files[i].single_spectra[j]);
		}
	}
	this->next_ssf_pointer=0;
}


// Randomly select n ssf pointers
void FileSet::randomly_reduce_ssfs(int n)
{
	if (n>= ssf_pointers.size())
		return;

	vector<int> idxs;
	chooseKFromN(n,ssf_pointers.size(),idxs);

	vector<SingleSpectrumFile *> red_ssfs;
	red_ssfs.resize(n,NULL);

	int i;
	for (i=0; i<n; i++)
		red_ssfs[i] = ssf_pointers[idxs[i]];

	ssf_pointers = red_ssfs;
}


void FileSet::remove_spectra_with_PTMs()
{
	int i;
	for (i=0; i<ssf_pointers.size(); i++)
	{
		if (ssf_pointers[i]->peptide.get_num_aas()>0)
		{
			vector<int> aas;
			aas=ssf_pointers[i]->peptide.get_amino_acids();
			int j;
			for (j=0; j<aas.size(); j++)
				if (aas[j]>Val)
					break;
			
			if (j<aas.size())
			{
				ssf_pointers[i]=ssf_pointers[ssf_pointers.size()-1];
				ssf_pointers.pop_back();
				i--;
			}
		}
	}
}

struct mz_pair {
	bool operator< (const mz_pair& other) const
	{
		return (m_over_z < other.m_over_z);
	}

	mass_t m_over_z;
	SingleSpectrumFile *ssf;
};

void FileSet::sort_according_to_m_over_z()
{
	vector<mz_pair> pairs;
		
	pairs.resize(ssf_pointers.size());
	int i;
	for (i=0; i<ssf_pointers.size(); i++)
	{
		pairs[i].m_over_z=ssf_pointers[i]->m_over_z;
		pairs[i].ssf = ssf_pointers[i];
	}

	sort(pairs.begin(),pairs.end());
	
	for (i=0; i<pairs.size(); i++)
		ssf_pointers[i]=pairs[i].ssf;
}


// removes all ssf without a peptides
void FileSet::keep_only_spectra_with_peptides()
{
	int i;
	vector<SingleSpectrumFile *> keep_ssfs;

	keep_ssfs.clear();

	for (i=0; i<ssf_pointers.size(); i++)
		if (ssf_pointers[i]->peptide.get_num_aas()>3)
			keep_ssfs.push_back(ssf_pointers[i]);

	ssf_pointers = keep_ssfs;
}



bool comp_ssf_pointers(const SingleSpectrumFile *ss1,
			           const SingleSpectrumFile *ss2)
{
	return ( (ss1->file_idx < ss2->file_idx) ||
		     (ss1->file_idx == ss2->file_idx && ss1->file_pos < ss2->file_pos));
}

// copies the ssf pointers from another FileSet
// sorts them according to the file order (not m_over_z);
void FileSet::init_from_another_fs(const FileSet& other_fs, 
						int start_ssf_idx, int end_ssf_idx)
{
	int i;

	reset_pointers();

	ssf_pointers.clear();
	ssf_pointers.reserve(end_ssf_idx-start_ssf_idx+1);

	int not_added=0;
	for (i=start_ssf_idx; i<=end_ssf_idx; i++)
		if (other_fs.ssf_pointers[i]->assigned_cluster<0)
		{
			ssf_pointers.push_back(other_fs.ssf_pointers[i]);
		}
	//	else
	//		not_added++;

//	if (not_added>0)
//		cout << "Not added: " << not_added << endl;

	// the comp_ssf_pointers is used so the ssfs are addressed according to 
	// the order in which they appear in the list file
	// if this causes problems can use regular sort according to pointer value
	// this will cause the files to be accessed in an unpredictable order (however
	// the spectra in the files will be accessed in the correct order)
	sort(ssf_pointers.begin(),ssf_pointers.end(),comp_ssf_pointers);

}


/************************************************************
// reads the next spectrum into spec
// returns false if no more spectra are available
// returns the mgf_pointer through mp, if file = -1 this means
// this was a dta
*************************************************************/
bool FileSet::get_next_spectrum(const FileManager& fm, 
								Config *config, 
								Spectrum *spec, 
								SingleSpectrumFile **ssfp,
								bool perform_init_spectrum, 
								bool set_charge_to_zero)
{
	SingleSpectrumFile *ssf;
	if (next_ssf_pointer < ssf_pointers.size())
	{
		ssf = ssf_pointers[next_ssf_pointer++];
		if (ssfp)
			*ssfp=ssf;

		if (ssf->type == DTA)
		{
			if (set_charge_to_zero)
			{
				spec->read_from_dta(config,ssf->single_name.c_str(),0);
			}
			else
				spec->read_from_dta(config,ssf->single_name.c_str());

			if (perform_init_spectrum)
				spec->initializeSpectrum(ssf);

			return true;
		}
		else if (ssf->type == MGF)
		{
			if (this->current_mgf_file_idx != ssf->file_idx)
			{
				if (mgf_stream)
					fclose(mgf_stream);
				
				mgf_stream=fopen(fm.mgf_files[ssf->file_idx].mgf_name.c_str(),"r");
				if (! mgf_stream)
				{
					cout << "Error: Couldn't open file for reading: " << 
						fm.mgf_files[ssf->file_idx].mgf_name << endl;
					exit(1);
				}
			
				this->current_mgf_file_idx = ssf->file_idx;
			}

			if ( fseek(mgf_stream,ssf->file_pos,0) )
			{
				cout << "Error: could not skip in file!" << endl;
				exit(1);
			}

			spec->read_from_MGF_stream(config,mgf_stream);
			if (perform_init_spectrum)
				spec->initializeSpectrum(ssf);

			return true;
		}
		else
		{
			cout << "Error: invalid file type:" << ssf->type <<  endl;
			exit(1);
		}
		
	}
	else
		return false;
}



inline bool comp_ssf(const SingleSpectrumFile* ssf1, const SingleSpectrumFile* ssf2)
{
	return (   ssf1->file_idx < ssf2->file_idx || 
		    ( (ssf1->file_idx == ssf2->file_idx) && (ssf1->file_pos<ssf2->file_pos) ) );
}








// creates an mgf file with the desred number of spectra per charge as designated
// in the spectra_per_charges vector(strating from charge 0). so the maximum
// number of from each charge in the outputted in the mgf file will be at most
// the numbers desginated in the vector
void FileSet::create_mgf_file(const FileManager& fm, Config *config, const char *file_name,
						 vector<int> spectra_per_charges)
{
	int max_charge = spectra_per_charges.size()-1;
	vector<int> spectra_written;

	spectra_written.resize(max_charge+1,0);

	ofstream ofs(file_name,ios::out);
	if (! ofs.good())
	{
		cout << "Error: couldn't open file for writing: " << file_name << endl;
		exit(1);
	}

	this->reset_pointers();

	while (1 )
	{
		Spectrum s;
		if (! this->get_next_spectrum(fm,config,&s))
			break;

	//	s.set_charge(0);

	/*	int charge = s.get_charge();
		cout << "Before: " << s.get_m_over_z() << endl;
		s.set_m_over_z((s.get_true_mass_with_19()+(charge)*MASS_PROTON)/(charge+1));
		s.set_charge(charge+1);
		cout << "After : " << s.get_m_over_z() << endl;	
		int spec_charge = s.get_charge()-1;*/

		int spec_charge = s.getCharge();
		if (spec_charge>max_charge)
			continue;

		if (spectra_written[spec_charge]>=spectra_per_charges[spec_charge])
			continue;

		spectra_written[spec_charge]++;

		s.output_as_MGF(ofs);
	}

	ofs.close();

}





/*************************************************************
// iterates over the files and makes a fasta out of the seq
//  puts 10 in a row, calls it TRUE_X 
  uses the org_aa
**************************************************************/
void FileSet::make_fasta_from_file_seqs(const FileManager& fm, Config *config,
										int inc, ostream& os)
{
	int seq_counter=0;
	int line_counter=0;

	Spectrum s;

	while ( 1 )
	{
		if ( ! this->get_next_spectrum(fm,config,&s))
			break;

		if (s.getPeptide().get_num_aas()>0)
		{
			if (seq_counter == 0)
				os << ">TRUE_" << line_counter*inc << endl;
	
			Peptide p = s.getPeptide();
			p.convert_to_org(config);

			os << p.as_string(config);	
			seq_counter++;
			if (seq_counter == inc)
			{
				os << endl;
				seq_counter=0;
				line_counter++;
			}
		}
	}
	os << endl;
}






// concatonates several mgf files into one large file
// slow and stupid...
void concat_mgf_files(Config *config, const char *mgf_file_list, 
					  const char *big_mgf_file)
{
	int i;
	ofstream ofs(big_mgf_file,ios::out);
	if (! ofs.good() )
	{
		cout << "Error: couldn't open file for writing: " << big_mgf_file << endl;
		exit(1);
	}

	vector<string> list;
	read_paths_into_list(mgf_file_list,list);
	
	for (i=0; i<list.size(); i++)
	{
		Spectrum s;
		
		FILE *mgf_stream = fopen(list[i].c_str(),"r");
		if (mgf_stream)
		{
			cout << "Error: couldn't open file for reading: " << list[i] << endl;
			exit(1);
		}

		while ( s.read_from_MGF_stream(config,mgf_stream) )
		{
		//	s.init_spectrum(); // so we filter peaks
			s.output_as_MGF(ofs);
		}
		fclose(mgf_stream);
	}
	ofs.close();
}


/*********************************************************************
Reads the lines in an mgf file to find the number of BEGIN IONS
and the number of lines.
**********************************************************************/
void examine_mgf_file(char *mgf_file, int *max_num_spectra, int *max_num_lines)
{
	char buff[256];
	int num_lines=0;
	int num_begin_ions=0;
	FILE *mgf_stream = fopen(mgf_file,"r");
	if (! mgf_stream)
	{
		cout << "Error: couldn't open file for reading: " << mgf_file << endl;
		exit(1);
	}

	while (fgets(buff,256,mgf_stream))
	{
		if (strstr(buff,"BEGIN IONS"))
			num_begin_ions++;
		num_lines++;
	}

	*max_num_lines=num_lines;
	*max_num_spectra = num_begin_ions;
}



void FileManager::print_summary_stats() const
{
	int i;

	int total_mgf = 0, total_mzxml=0, total_ms2=0, total_dat =0;
	
	for (i=0; i<mgf_files.size(); i++)
		total_mgf += mgf_files[i].total_num_spectra;

	for (i=0; i<mzxml_files.size(); i++)
		total_mzxml += mzxml_files[i].total_num_spectra;

	for (i=0; i<dat_files.size(); i++)
		total_dat += dat_files[i].total_num_spectra;

	for (i=0; i<ms2_files.size(); i++)
		total_ms2 += ms2_files[i].total_num_spectra;

	cout << "Min mass:    " << setw(8) << fixed << setprecision(2) << min_spec_mass << endl;
	cout << "Max mass:    " << setw(8) << fixed << setprecision(2) << max_spec_mass << endl;
	cout << "DTA files:   " << this->dta_files.size() << " spectra." << endl;
	cout << "MGF_files:   " << this->mgf_files.size() << " files, " << total_mgf
		 << " spectra." << endl;
	cout << "MZXML_files: " << this->mzxml_files.size() << " files, " << total_mzxml
		 << " spectra." << endl;
	cout << "DAT_files:   " << this->dat_files.size() << " files, " << total_dat <<
		    " spectra," << endl;
	cout << "MS2_files:   " << this->ms2_files.size() << " files, " << total_ms2 <<
		    " spectra." << endl;
	cout << "Total: " << total_mgf+dta_files.size()+total_mzxml + total_dat + total_ms2 << " spectra." << endl << endl;

}


void SingleSpectrumFile::print_ssf_stats(const Config *config, ostream& os, bool print_endl) const
{
	if (type == MGF)
	{
		MGF_single *ssf = (MGF_single *)this;
		os << ">> " << ssf->file_idx << " " << ssf->idx_in_file << " " << ssf->single_name;
	}
	else if (type == MZXML)
	{
		MZXML_single *ssf = (MZXML_single *)this;
		os << ">> " << ssf->file_idx << " " << ssf->scan_number;
	}
	else if (type == DAT)
	{
		DAT_single *ssf = (DAT_single *)this;
		os << ">> " << ssf->mzxml_file_idx << " " << ssf->scan_number;
	}	
	else
		os << ">> " << this->file_idx << " " << this->single_name;

	if (this->peptide.get_num_aas()>1)
	{
		os << " " << this->peptide.as_string(config);
		os << " " << peptide.get_mass() + MASS_OHHH;
	}

	if (print_endl)
		os << endl;
	
}





void FileSet::print_file_stats() const
{
	int i;
	cout << this->ssf_pointers.size() << " total spectra." << endl;


	for (i=0; i<ssf_pointers.size() ;i++)
	{

		cout << setw(5) << left << i << setw(8) << setprecision(3) <<
			ssf_pointers[i]->org_pm_with_19 << "  " << ssf_pointers[i]->charge <<
			ssf_pointers[i]->single_name << endl;
	}
}



void FileSet::print_summary() const
{
	int i;
	int max_charge = 20;
	vector<int> charge_counts;

	charge_counts.resize(max_charge+1,0);

	for (i=0; i<ssf_pointers.size(); i++)
	{
		int c = ssf_pointers[i]->charge;
		if (c>max_charge)
			c=max_charge;

		charge_counts[c]++;
	}
	
	cout << "Read total " << ssf_pointers.size() <<  endl;
	int c;
	for (c=0; c<=max_charge; c++)
	{
		if (charge_counts[c]>0)
			cout << "Charge " << c << " " << charge_counts[c] << " (" <<
				(double)charge_counts[c]/ssf_pointers.size() << ")" << endl;
	}
}


/********************************************************************************
	Reads files and outputs a list of the scans m/z
*********************************************************************************/
void extractMZFromFiles(Config *config, char *file_list, char *output_file)
{
	FileManager fm;
	FileSet		fs;
	vector<string> list;

	ofstream out(output_file,ios::out);

	read_paths_into_list(file_list,list);
	fm.init_from_list(config,list);
	fs.select_all_files(fm);


	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	int i;
	for (i=0; i<all_ssf.size(); i++)
	{
		MZXML_single *ssf = (MZXML_single *)all_ssf[i];

		out << ssf->file_idx << "\t" << ssf->scan_number << "\t" << ssf->m_over_z << endl;
	}
}









	
