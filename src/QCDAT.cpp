#include "QuickClustering.h"
#include "base64.h"

DAT_FileBuff::~DAT_FileBuff()
{
	if (buff && pos>buff)
		flush_buff();

	if (buff)
		delete [] buff; 
}

void DAT_FileBuff::init(string& _path,  int dat_buff_size)
{
	if (! ind_was_initialized)
	{
		buff = new char[dat_buff_size];
		if (! buff)
		{
			cout << "Error: couldn't allocate memory for DAT file buff!" << endl;
			exit(1);
		}
	}

	ind_was_initialized = 1;
	path = _path;
	max_pos = buff + dat_buff_size - 100;
	pos = buff;
	ind_first_write = 1; 
}


// copies the files to the DAT file
void DAT_FileBuff::add_spec_to_DAT_file(
							mass_t m_over_z, 
							int charge, 
							int mzxml_file_idx,
							int scan_number, 
							float retention_time, 
							float precursor_intensity, 
							int num_peaks, 
							char *peak_buff)
{
	const int spec_bytes = sizeof(mass_t) + 4 * sizeof(int) + 2 * sizeof(float) * num_peaks;

	if (num_peaks<=2)
	{
		//cout << mzxml_file_idx << " " << scan_number << " " << " p: " << num_peaks << endl;
		return;
	}

	if (mzxml_file_idx<0 || scan_number<0)
	{
		cout <<"Error: bad file idx or scan number! " << mzxml_file_idx << ", " << scan_number << endl;
		exit(1);
	}

	if (pos + spec_bytes >= max_pos)
		flush_buff();

	mass_t *m_ptr = (mass_t *)pos;
	*m_ptr++ = m_over_z;

	int *i_ptr = (int *)m_ptr;

	*i_ptr++ = charge;
	*i_ptr++ = mzxml_file_idx;
	*i_ptr++ = scan_number;
	*i_ptr++ = num_peaks;

	float *f_ptr = (float *)i_ptr;

	*f_ptr++ = retention_time;
	*f_ptr++ = precursor_intensity;
	
	pos = (char *)f_ptr;
	memcpy(pos,peak_buff,2 * sizeof(float) * num_peaks);
	pos += 2 * sizeof(float) * num_peaks;

	this->counter++;

}



void DAT_FileBuff::flush_buff()
{
	FILE* dat_stream;

	if (pos == buff)
		return;

	if (! ind_was_initialized)
	{
		cout << "Error: must first initialize the DAT_FileBuff!" << endl;
		exit(1);
	}

	if (ind_first_write)
	{
		dat_stream=fopen(path.c_str(),"wb");
		ind_first_write = 0;
	}
	else
		dat_stream=fopen(path.c_str(),"ab");

	if (! dat_stream)
	{
		cout << "Error: couldn't open DAT file for writing: " << path.c_str() << endl;
		exit(1);
	}


	fwrite(buff,1,pos-buff,dat_stream);
	fclose(dat_stream);

	pos=buff;
}



void DAT_Converter::init_DAT_Converter(mass_t _max_m_over_z, mass_t _mass_increment, int _dat_buff_size)
{
//	name = _name;
//	out_dir = _out_dir;
	max_m_over_z = _max_m_over_z;
	mass_increment = _mass_increment;
	dat_buff_size = _dat_buff_size;

	max_dat_file_idx = (int)(max_m_over_z / mass_increment) + 1;

	dat_buffs.resize(max_dat_file_idx+1);

	ind_was_initialized =1;
}




// creates a file with the list of DAT files
void DAT_Converter::create_list_file() const
{
	FILE *list_stream;

	ostringstream oss;
	oss << batch;
	string list_path = out_dir + "/" + name + "_" + oss.str() + "_list.txt";

	list_stream = fopen(list_path.c_str(),"w");
	if (! list_stream)
	{
		cout << "Error: couldn't open list file for writing: " << list_path.c_str() << endl;
		exit(1);
	}

	int i;
	for (i=0; i<dat_buffs.size(); i++)
	{
		if (dat_buffs[i].buff)
		{
			fprintf(list_stream,"%s\n",dat_buffs[i].path.c_str());
		}
	}

	fclose(list_stream);
}








int DAT_Converter::convert_single_non_MZXML_file_to_DAT(Config* config, string file, 
							mass_t min_m_over_z, mass_t max_m_over_z, int file_idx)
{
	static QCPeak *peak_buff=NULL;    
	static float  *dat_peak_buff=NULL;

	if (! peak_buff)
	{
		peak_buff     = new QCPeak[20000];
		dat_peak_buff = new  float[40000];
	}

	if (! peak_buff || ! dat_peak_buff)
	{
		cout << "Error: couldn't allocate memory for DAT conversion!" << endl;
		exit(1);
	}

	if (! ind_was_initialized)
	{
		cout << "Error: must initialize DAT_Converter!" << endl;
		exit(0);
	}

	vector<string> list;
	list.push_back(file);

	BasicSpecReader bsr;
	FileManager fm;
	FileSet all_spec_fs;


	
	fm.init_from_list(config,list,true,file_idx);

	all_spec_fs.select_all_files(fm);

	const int total_spectra = all_spec_fs.get_total_spectra();
	const vector<SingleSpectrumFile *>& all_ssf = all_spec_fs.get_ssf_pointers();

	cout << file_idx << " " << file << " ... " << endl;

	int num_spectra_extracted=0;
	int i;
	for (i=0; i<all_ssf.size(); i++)
	{

		SingleSpectrumFile *ssf = all_ssf[i];
		// no small spectra
		if (ssf->num_peaks<5 || ssf->m_over_z<=min_m_over_z || ssf->m_over_z > max_m_over_z)
			continue;

		int j,pos=0;
		int num_spec_peaks = bsr.read_basic_spec(config,fm,ssf,peak_buff,true);

		if (num_spec_peaks<5)
			continue;


		mass_t m_over_z = ssf->m_over_z;
		int scan_number=-1, file_idx=-1;
		float precursor_intensity=0;
		
		if (ssf->type == DTA)
		{
			file_idx = 0;
			scan_number=-1;
		}
		else if (ssf->type == MGF)
		{
			MGF_single *mgf_single = (MGF_single *)ssf;
			file_idx    = mgf_single->file_idx;
			scan_number = (mgf_single->scan_number > 0 ? mgf_single->scan_number : mgf_single->idx_in_file);	
		//	cout << "s " << mgf_single->scan_number << "\tidx " << mgf_single->idx_in_file << endl;  
		}
		else if (ssf->type == MZXML)
		{
			MZXML_single *mzxml_single = (MZXML_single *)ssf;
			file_idx = mzxml_single->file_idx;
			scan_number = mzxml_single->scan_number;
			precursor_intensity = mzxml_single->precursor_intensity;
		}
		else if (ssf->type == DAT)
		{
			DAT_single *dat_single = (DAT_single *)ssf;
			file_idx = dat_single->mzxml_file_idx;
			scan_number = dat_single->scan_number;
		}

		// copy peaks
		
		for (j=0; j<num_spec_peaks; j++)
		{
			dat_peak_buff[pos++]=(float)peak_buff[j].mass;
			dat_peak_buff[pos++]=(float)peak_buff[j].intensity;
		}

		
		// add spectrum
		int DAT_file_idx =  (int)(m_over_z/mass_increment);
		if (DAT_file_idx > max_dat_file_idx)
			DAT_file_idx = max_dat_file_idx;

		if (! dat_buffs[DAT_file_idx].ind_was_initialized)
		{
			ostringstream os,os_batch;
			os << DAT_file_idx;
			os_batch << batch;
			string path = out_dir + "/" + name + "_" + os_batch.str() + "_" + os.str() + ".dat";
			dat_buffs[DAT_file_idx].init(path,dat_buff_size);
		}

		int charge = ssf->charge;
		if (charge<0 || charge>1000)
			charge=0;

	
		dat_buffs[DAT_file_idx].add_spec_to_DAT_file(
			m_over_z,
			charge,
			file_idx,
			scan_number,ssf->retention_time,
			precursor_intensity,
			num_spec_peaks,
			(char *)dat_peak_buff);

		num_spectra_extracted++;
	}
	cout << num_spectra_extracted << " spectra..." << endl;

	return num_spectra_extracted;
}

int DAT_Converter::convert_PKL_dir_to_DAT(Config* config, char *file_list, int file_start_idx,
							mass_t min_m_over_z, mass_t max_m_over_z)
{
	FILE *list_stream=fopen(file_list,"r");
	if (! list_stream)
	{
		cout << "Error: couldn't open pkl dir list for reading: " << file_list << endl;
		exit(1);
	}

	char line_buff[1024];

	cout << "Read pkl files (idx path #spectra):" << endl;

	int total_extracted=0;
	int pkl_dir_idx=0;
	while (fgets(line_buff,1024,list_stream))
	{
		istringstream is(line_buff);

		string pkl_dir_path = "";
		string tsv_file = "";

		is >> pkl_dir_path >> tsv_file;

		if (pkl_dir_path.length()>2 && tsv_file.length() >2)
		{
			FileManager fm;

			fm.init_from_single_pkl_dir(config, pkl_dir_path, tsv_file, pkl_dir_idx,
										min_m_over_z, max_m_over_z);

			FileSet fs;
			fs.select_all_files(fm);
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
			
			cout << pkl_dir_idx << "\t" << pkl_dir_path << "\t" << all_ssf.size() << endl;


			static QCPeak *peak_buff=NULL;    
			static float  *dat_peak_buff=NULL;

			if (! peak_buff)
			{
				peak_buff     = new QCPeak[20000];
				dat_peak_buff = new  float[40000];
			}

			if (! peak_buff || ! dat_peak_buff)
			{
				cout << "Error: couldn't allocate memory for DAT conversion!" << endl;
				exit(1);
			}

			if (! ind_was_initialized)
			{
				cout << "Error: must initialize DAT_Converter!" << endl;
				exit(0);
			}


			BasicSpecReader bsr;

			int num_spectra_extracted=0;
			int i;
			for (i=0; i<all_ssf.size(); i++)
			{

				PKL_single *ssf = (PKL_single *)all_ssf[i];
				// no small spectra
				if (ssf->num_peaks<5 || ssf->m_over_z<=min_m_over_z || ssf->m_over_z > max_m_over_z)
					continue;

				int j,pos=0;
				int num_spec_peaks = bsr.read_basic_spec(config,fm,ssf,peak_buff,true);

				if (num_spec_peaks<5)
					continue;


				mass_t m_over_z = ssf->m_over_z;
				int scan_number=ssf->scan_number;
				float precursor_intensity= (ssf->precursor_intensity >0 ? ssf->precursor_intensity : 0);
		
				// copy peaks
		
				for (j=0; j<num_spec_peaks; j++)
				{
					dat_peak_buff[pos++]=(float)peak_buff[j].mass;
					dat_peak_buff[pos++]=(float)peak_buff[j].intensity;
				}

		
				// add spectrum
				int DAT_file_idx =  (int)(m_over_z/mass_increment);
				if (DAT_file_idx > max_dat_file_idx)
					DAT_file_idx = max_dat_file_idx;

				if (! dat_buffs[DAT_file_idx].ind_was_initialized)
				{
					ostringstream os,os_batch;
					os << DAT_file_idx;
					os_batch << batch;
					string path = out_dir + "/" + name + "_" + os_batch.str() + "_" + os.str() + ".dat";
					dat_buffs[DAT_file_idx].init(path,dat_buff_size);
				}

				int charge = ssf->charge;
				if (charge<0 || charge>1000)
					charge=0;

	
				dat_buffs[DAT_file_idx].add_spec_to_DAT_file(
					m_over_z,
					charge,
					pkl_dir_idx + file_start_idx,
					scan_number,ssf->retention_time,
					precursor_intensity,
					num_spec_peaks,
					(char *)dat_peak_buff);

				num_spectra_extracted++;
			}
			cout << "(" << num_spectra_extracted << ")" << endl;

			total_extracted += num_spectra_extracted;

			pkl_dir_idx++;
		}
	}

	return total_extracted;
}

void DAT_Converter::convert_files_to_DAT_on_the_fly(Config* config, char *file_list, 
							char * _out_dir, char * _name, int _batch, 
							mass_t min_m_over_z, 
							mass_t max_m_over_z, 
							int file_start_idx,
							bool ind_is_pkl_dir)
{
	name = _name;
	out_dir = _out_dir;
	batch = _batch;

	int num_spectra_extracted=0;

	name = _name;
	out_dir = _out_dir;

	if (! ind_was_initialized)
	{
		cout << "Error: must initialize DAT_Converter!" << endl;
		exit(0);
	}


	if (ind_is_pkl_dir)
	{
		num_spectra_extracted = convert_PKL_dir_to_DAT(config,file_list, file_start_idx, min_m_over_z,max_m_over_z);
	}
	else
	{
		vector<string> list;
		read_paths_into_list(file_list,list);

		cout << endl << endl <<"Extracting spectra and writing dat files for " << list.size() << " Files. " << endl << endl;

		int i;
		for (i=0; i<list.size(); i++)
		{
			int file_type = get_file_extension_type(list[i]);

			if (file_type == MZXML)
			{
				num_spectra_extracted += parse_single_MZXML_file(config,list[i],file_start_idx+i,
					min_m_over_z,max_m_over_z);
			}
			else
			{
				num_spectra_extracted += convert_single_non_MZXML_file_to_DAT(config,list[i],
											min_m_over_z,max_m_over_z,file_start_idx+i);
			}
		}
	}
	

	int d;
	for (d=0; d<dat_buffs.size(); d++)
	{
		if (dat_buffs[d].ind_was_initialized && dat_buffs[d].pos > dat_buffs[d].buff)
			dat_buffs[d].flush_buff();
	}
	
	
	cout << "Total spectra extracted and converted to DAT: " << num_spectra_extracted << endl;

	create_list_file();
}





int DAT_Converter::parse_single_MZXML_file(Config *config, 
										   string& mzxml_name, 
										   int file_idx,
										   mass_t min_m_over_z, 
										   mass_t max_m_over_z)
{
	static char* Buffer = NULL;
    static char* DecodedPeakBuffer = NULL;
    static float* Peaks = NULL;
	static float* FilteredPeaks = NULL;
	static int PeakBufferSize = 0;
    int Trail;
    static char* PrecursorStr;
    int FloatIndex;
    char* ByteOrderStr;
    int ByteOrderLittle = 1;

    int BytesToRead;
    int BufferStartPos = 0;
    int BytesRead;
    int BufferEnd = 0;
    FILE* MZXMLFile;
    int ParseState = 0;
    int FilePos = 0;
    

    // allocate
	if (! Buffer)
		Buffer = (char*)calloc(XML_BUFFER_SIZE + 1, sizeof(char));

    MZXMLFile = fopen(mzxml_name.c_str(), "rb");
    if (!MZXMLFile)
    {
        cout << "Error: Can't open MZXML file " <<  mzxml_name << endl;
        exit(1);
    }

    printf("File idx %d , '%s'...\n",file_idx, mzxml_name.c_str());

	int idx_in_file = 0;
	int spec_counter =0;
	char *scan_start_ptr = NULL;
    while (1)
    {
		char* ScanStr;
		char* ScanNumberStr;
		int ScanNumber;
		char* MSLevelStr;
		int MSLevel;
		char *PrecursorStr;
		mass_t PrecursorMZ;
		char *retentionTimeStr;
		float retentionTime;
		char *precursorIntensityStr;
		float precursorIntensity;
		char* PeakStr;
		char* PeakBuffer;
		int  BufferPos;

        // Read more data, to fill up the buffer:
	
		if ( ! scan_start_ptr || 
			( (Buffer + BufferEnd - scan_start_ptr) < XML_BUFFER_HALF_SIZE) )
		{
			// try shunt half of the buffer
			if (scan_start_ptr)
			{
				if (BufferEnd - XML_BUFFER_HALF_SIZE>0)
				{
					memmove(Buffer, Buffer + XML_BUFFER_HALF_SIZE, BufferEnd - XML_BUFFER_HALF_SIZE);
					BufferEnd -= XML_BUFFER_HALF_SIZE;
					scan_start_ptr -= XML_BUFFER_HALF_SIZE;

//					cout << "MOVED!" << endl;
				}
			}
			else
				scan_start_ptr = Buffer;

			BytesToRead = XML_BUFFER_SIZE - BufferEnd;
			BytesRead = fread(Buffer + BufferEnd, sizeof(char), BytesToRead, MZXMLFile);


			if (BytesRead<5)
				break;

			BufferEnd += BytesRead;
			Buffer[BufferEnd] = '\0';


			FilePos += BytesRead;
		}

        // Look for a new <scan tag opening:
		// this scan cannot be done with strstr since there might be NULL termination
		const char *last_pos = Buffer + BufferEnd - 5;
		char *pos = scan_start_ptr;

		while (++pos<last_pos)
		{
			if (*pos != '<')
				continue;

			if (*(pos+1)=='s' && *(pos+2)=='c' && *(pos+3)=='a' && *(pos+4)=='n')
				break;
		}
		ScanStr =  (pos<last_pos) ? pos : NULL;

        if (ScanStr)
        {
			// if this is the case, read over more data to avoid the case where
			// the spectrum's record is not all in the buffer
			if (scan_start_ptr - Buffer > XML_BUFFER_HALF_SIZE)
			{
				scan_start_ptr = ScanStr-2;
				continue;
			}

            BufferPos = ScanStr - Buffer;
        }
        else
        {
            BufferPos = 0;
        }

        if (!ScanStr )
        {
			scan_start_ptr = Buffer + BufferEnd-5;
            continue;
        }

        ScanNumberStr = strstr(ScanStr, "num=");
        if (!ScanNumberStr)
        {
          //  printf("** Warning: mzXML parser encountered a scan with no scan number!  File %s Pos %d\n", 
		  //		   mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

            ScanNumber = -1;
        }
        else
        {
            ScanNumber = ParseIntFromXML(ScanNumberStr);
			scan_start_ptr = ScanNumberStr;
        }

		retentionTimeStr = strstr(ScanStr,"retentionTime=\"PT");
		if (! retentionTimeStr)
		{
			retentionTime = -1;
		}
		else
		{
			retentionTime = ParseMassFromXML(retentionTimeStr);
		}


		char *PeakCountStr = strstr(ScanStr, "peaksCount=\"");
		if (!PeakCountStr)
		{
			cout << "Warning: bad parsing peaks in mzxml. " << endl;
			cout << "Scan: " << ScanNumber << " Pos: " << FilePos << endl;
			scan_start_ptr += 50;
			continue;
		}
		int PeakCount = ParseIntFromXML(PeakCountStr);
		
        MSLevelStr = strstr(ScanStr, "msLevel=");
        if (!MSLevelStr)
        {
            printf("** Warning: mzXML parser encountered a scan with no scan level!  File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

			scan_start_ptr += 50;
			continue; 
        }
        else
        {
            MSLevel = ParseIntFromXML(MSLevelStr);
        }

		if (MSLevel<=1)
		{
			scan_start_ptr += 50;
			continue;
		}


		precursorIntensityStr = strstr(ScanStr,"precursorIntensity=");
		if (! precursorIntensityStr)
		{
			cout << "Warning: no precursor intensity found for scan " << ScanNumber << " pos: " << BufferPos << endl;
			scan_start_ptr += 50;
			continue;
		}
		else
		{
			precursorIntensity = ParseMassFromXML(precursorIntensityStr);
		}


		PrecursorStr = strstr(ScanStr, "<precursorMz");
		if (PrecursorStr)
		{
			PrecursorStr = strstr(PrecursorStr, ">");
			PrecursorMZ = ParseMassFromXML(PrecursorStr);
		}

		if (!PrecursorStr && MSLevel > 1)
		{
			printf("Warning: mzXML parser encountered a scan with no m/z: File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

			scan_start_ptr += 50;
			continue;
		}
		
		//  check that this is a good spectrum to output
		if (MSLevel>1 && PeakCount>2 && PrecursorMZ>= min_m_over_z && PrecursorMZ <= max_m_over_z)
		{

			// read peaks

			PeakStr = strstr(PrecursorStr, "<peaks");
			if (PeakStr)
			{
				// Get byte order:
				ByteOrderStr = strstr(PeakStr, "byteOrder=\"");
				if (ByteOrderStr)
				{
					ByteOrderStr += 11;
					if (!strncmp(ByteOrderStr, "network", 7))
					{
						ByteOrderLittle = 0;
					}
					if (!strncmp(ByteOrderStr, "big", 3))
					{
						ByteOrderLittle = 0;
					}
					if (!strncmp(ByteOrderStr, "little", 6))
					{
						ByteOrderLittle = 1;
					}
				}
				PeakStr = strstr(PeakStr, ">");
			}
			if (!PeakStr)
			{
				cout << "Warning couldn't find peaks tag for scan: " << ScanNumber << " Buffer pos: " << BufferPos << "  skipping..." << endl;
				scan_start_ptr += 50;
				continue;
			}

			PeakStr++;
			PeakBuffer = PeakStr;

			if (PeakBufferSize < PeakCount)
			{
				if (DecodedPeakBuffer)
				{
					char *dbf = DecodedPeakBuffer;
					free(DecodedPeakBuffer);
					DecodedPeakBuffer = NULL;
					free(Peaks);
					Peaks = NULL;
					free(FilteredPeaks);
					FilteredPeaks=NULL;
				}
				PeakBufferSize = (int)(PeakCount*1.5);
				DecodedPeakBuffer = (char*)calloc(PeakBufferSize * 8 + 8, 1);
				Peaks = (float*)calloc(PeakBufferSize * 2, sizeof(float));
				FilteredPeaks = (float*)calloc(PeakBufferSize * 2, sizeof(float));
			}
			
			Trail = (PeakCount % 3);

			int pbf = PeakBufferSize;

			if (!(PeakCount % 3))
			{
				PeakBuffer[PeakCount * 32/3] = '\0';
			}
			else
			{
				PeakBuffer[(PeakCount * 32/3) + Trail + 1] = '\0';
			}

			b64_decode_mio( DecodedPeakBuffer, PeakBuffer);
			for (FloatIndex = 0; FloatIndex < (2 * PeakCount); FloatIndex++)
			{
		#ifdef BYTEORDER_LITTLE_ENDIAN
				if (!ByteOrderLittle)
				{
					char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
					DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
					DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
					ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
					DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
					DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
				}
				memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
		#else
				if (ByteOrderLittle)
				{
					char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
					DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
					DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
					ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
					DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
					DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
				}
				memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
		#endif
			}


			// add spectrum
			int DAT_file_idx =  (int)(PrecursorMZ/mass_increment);
			if (DAT_file_idx > max_dat_file_idx)
				DAT_file_idx = max_dat_file_idx;

			if (! dat_buffs[DAT_file_idx].ind_was_initialized)
			{
				ostringstream os, os_batch;
				os << DAT_file_idx;
				os_batch << batch;
				string path = out_dir + "/" + name + "_" + os_batch.str() + "_" + os.str() + ".dat";
				dat_buffs[DAT_file_idx].init(path,dat_buff_size);
			}

			// join and filter peaks
			int new_num_peaks = join_and_filter_peak_list(config,PrecursorMZ,Peaks,
														  PeakCount,FilteredPeaks);
			if (new_num_peaks>3)
			{
				dat_buffs[DAT_file_idx].add_spec_to_DAT_file(PrecursorMZ, 0, file_idx, ScanNumber, retentionTime,  
									 precursorIntensity, new_num_peaks, (char *)FilteredPeaks);
				spec_counter++;
			}
			scan_start_ptr = PeakStr + 8 * PeakCount;
		}
		else
			scan_start_ptr = ScanStr +50;
	}

	fclose(MZXMLFile);
	cout << spec_counter << " spectra..." << endl << endl;

	return spec_counter;
}


void create_annotated_mgf_from_dat(Config *config, 
                                   char *dat_list,
                                   char *mzxml_list,
                                   char *anns_file,
                                   char *output_file)
{
        vector<string> list, mzxml_names;

        read_paths_into_list(dat_list,list);

        vector< vector<int> >    annotation_idxs;
        vector<mzXML_annotation> annotations;

        read_mzXML_annotations(mzxml_list,anns_file, annotation_idxs, annotations);

        FileManager fm;
        fm.init_from_list_file_and_add_annotations(config,dat_list,annotation_idxs,annotations,true);

        FileSet fs;
        fs.select_all_files(fm);
        const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();


        BasicSpecReader bsr;
        int i;
        ofstream mgf_stream(output_file);

        QCPeak peaks[5000];
        for (i=0; i<all_ssf.size(); i++)
        {
                int num_spec_peaks = bsr.read_basic_spec(config,fm,all_ssf[i],peaks);
                BasicSpectrum bs;
                bs.num_peaks = num_spec_peaks;
                bs.peaks = peaks;
                bs.ssf = all_ssf[i];
                bs.ssf->charge = 2;


                DAT_single *dat_single = (DAT_single *)all_ssf[i];
                char single_name[256];
                sprintf(single_name,"spec_%d_%d",dat_single->mzxml_file_idx,dat_single->scan_number);

                bs.ssf->single_name = single_name;

                bs.output_to_mgf(mgf_stream,config);
        }

        mgf_stream.close();
}







int parse_single_MZXML_file_print_peaks(Config *config, 
										   string& mzxml_name, 
										   int the_scan)
										   
{
	static char* Buffer = NULL;
    static char* DecodedPeakBuffer = NULL;
    static float* Peaks = NULL;
	static float* FilteredPeaks = NULL;
	static int PeakBufferSize = 0;
    int Trail;
    static char* PrecursorStr;
    int FloatIndex;
    char* ByteOrderStr;
    int ByteOrderLittle = 1;

    int BytesToRead;
    int BufferStartPos = 0;
    int BytesRead;
    int BufferEnd = 0;
    FILE* MZXMLFile;
    int ParseState = 0;
    int FilePos = 0;
    

    // allocate
	if (! Buffer)
		Buffer = (char*)calloc(XML_BUFFER_SIZE + 1, sizeof(char));

    MZXMLFile = fopen(mzxml_name.c_str(), "rb");
    if (!MZXMLFile)
    {
        cout << "Error: Can't open MZXML file " <<  mzxml_name << endl;
        exit(1);
    }

  
	int idx_in_file = 0;
	int spec_counter =0;
	char *scan_start_ptr = NULL;
    while (1)
    {
		char* ScanStr;
		char* ScanNumberStr;
		int ScanNumber;
		char* MSLevelStr;
		int MSLevel;
		char *PrecursorStr;
		mass_t PrecursorMZ;
		char *retentionTimeStr;
		float retentionTime;
		char *precursorIntensityStr;
		float precursorIntensity;
		char* PeakStr;
		char* PeakBuffer;
		int  BufferPos;

        // Read more data, to fill up the buffer:
	
		if ( ! scan_start_ptr || 
			( (Buffer + BufferEnd - scan_start_ptr) < XML_BUFFER_HALF_SIZE) )
		{
			// try shunt half of the buffer
			if (scan_start_ptr)
			{
				if (BufferEnd - XML_BUFFER_HALF_SIZE>0)
				{
					memmove(Buffer, Buffer + XML_BUFFER_HALF_SIZE, BufferEnd - XML_BUFFER_HALF_SIZE);
					BufferEnd -= XML_BUFFER_HALF_SIZE;
					scan_start_ptr -= XML_BUFFER_HALF_SIZE;

//					cout << "MOVED!" << endl;
				}
			}
			else
				scan_start_ptr = Buffer;

			BytesToRead = XML_BUFFER_SIZE - BufferEnd;
			BytesRead = fread(Buffer + BufferEnd, sizeof(char), BytesToRead, MZXMLFile);


			if (BytesRead<5)
				break;

			BufferEnd += BytesRead;
			Buffer[BufferEnd] = '\0';

			
			// remove all '\0' from buffer (these cause parsing errors.
			// replace them with ' '
		//	const char *end_ptr = Buffer + BufferEnd - 2;
		//	for (char *ptr=Buffer; ptr<end_ptr; ptr++)
		//		if (*ptr == '\0')
		//			*ptr = ' ';
			
			FilePos += BytesRead;
		}

        // Look for a new <scan tag opening:
		// this scan cannot be done with strstr since there might be NULL termination
		const char *last_pos = Buffer + BufferEnd - 5;
		char *pos = scan_start_ptr;

		while (++pos<last_pos)
		{
			if (*pos != '<')
				continue;

			if (*(pos+1)=='s' && *(pos+2)=='c' && *(pos+3)=='a' && *(pos+4)=='n')
				break;
		}
		ScanStr =  (pos<last_pos) ? pos : NULL;

        if (ScanStr)
        {
			// if this is the case, read over more data to avoid the case where
			// the spectrum's record is not all in the buffer
			if (scan_start_ptr - Buffer > XML_BUFFER_HALF_SIZE)
			{
				scan_start_ptr = ScanStr-2;
				continue;
			}

            BufferPos = ScanStr - Buffer;
        }
        else
        {
            BufferPos = 0;
        }

        if (!ScanStr )
        {
			scan_start_ptr = Buffer + BufferEnd-5;
            continue;
        }

        ScanNumberStr = strstr(ScanStr, "num=");
        if (!ScanNumberStr)
        {
          //  printf("** Warning: mzXML parser encountered a scan with no scan number!  File %s Pos %d\n", 
		  //		   mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

            ScanNumber = -1;
        }
        else
        {
            ScanNumber = ParseIntFromXML(ScanNumberStr);
			scan_start_ptr = ScanNumberStr;
        }

		retentionTimeStr = strstr(ScanStr,"retentionTime=\"PT");
		if (! retentionTimeStr)
		{
		//	printf("Error: mzXML parser encountered a scan with no retnetion time: File %s Pos %d\n", 
		//		mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);
		//	cout <<"SCANSTR:"<<endl << ScanStr << endl;
         //   exit(1);
			retentionTime = -1;
		}
		else
		{
			retentionTime = ParseMassFromXML(retentionTimeStr);
		}


		char *PeakCountStr = strstr(ScanStr, "peaksCount=\"");
		if (!PeakCountStr)
		{
			cout << "Warning: bad parsing peaks in mzxml. " << endl;
			cout << "Scan: " << ScanNumber << " Pos: " << FilePos << endl;
			scan_start_ptr += 50;
			continue;
		//	exit(1);
		}
		int PeakCount = ParseIntFromXML(PeakCountStr);
		
        MSLevelStr = strstr(ScanStr, "msLevel=");
        if (!MSLevelStr)
        {
            printf("** Warning: mzXML parser encountered a scan with no scan level!  File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

			scan_start_ptr += 50;
			continue;
            
        }
        else
        {
            MSLevel = ParseIntFromXML(MSLevelStr);
        }

		if (MSLevel<=1)
		{
			scan_start_ptr += 50;
			continue;
		}

		if (ScanNumber != the_scan)
		{
			scan_start_ptr += 150;
			continue;
		}



		precursorIntensityStr = strstr(ScanStr,"precursorIntensity=");
		if (! precursorIntensityStr)
		{
			cout << "Warning: no precursor intensity found for scan " << ScanNumber << " pos: " << BufferPos << endl;
			scan_start_ptr += 50;
			continue;
		}
		else
		{
			precursorIntensity = ParseMassFromXML(precursorIntensityStr);
		}


		PrecursorStr = strstr(ScanStr, "<precursorMz");
		if (PrecursorStr)
		{
			PrecursorStr = strstr(PrecursorStr, ">");
			PrecursorMZ = ParseMassFromXML(PrecursorStr);
		}

		if (!PrecursorStr && MSLevel > 1)
		{
			printf("Warning: mzXML parser encountered a scan with no m/z: File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

			scan_start_ptr += 50;
			continue;
		}
		
		//  check that this is a good spectrum to output
		if (MSLevel>1 && PeakCount>2)
		{

			// read peaks

			PeakStr = strstr(PrecursorStr, "<peaks");
			if (PeakStr)
			{
				// Get byte order:
				ByteOrderStr = strstr(PeakStr, "byteOrder=\"");
				if (ByteOrderStr)
				{
					ByteOrderStr += 11;
					if (!strncmp(ByteOrderStr, "network", 7))
					{
						ByteOrderLittle = 0;
					}
					if (!strncmp(ByteOrderStr, "big", 3))
					{
						ByteOrderLittle = 0;
					}
					if (!strncmp(ByteOrderStr, "little", 6))
					{
						ByteOrderLittle = 1;
					}
				}
				PeakStr = strstr(PeakStr, ">");
			}
			if (!PeakStr)
			{
				cout << "Warning couldn't find peaks tag for scan: " << ScanNumber << " Buffer pos: " << BufferPos << "  skipping..." << endl;
				scan_start_ptr += 50;
				continue;
				
			}

			PeakStr++;
			PeakBuffer = PeakStr;

			if (PeakBufferSize < PeakCount)
			{
				if (DecodedPeakBuffer)
				{
					char *dbf = DecodedPeakBuffer;
					free(DecodedPeakBuffer);
					DecodedPeakBuffer = NULL;
					free(Peaks);
					Peaks = NULL;
					free(FilteredPeaks);
					FilteredPeaks=NULL;
				}
				PeakBufferSize = (int)(PeakCount*1.5);
				DecodedPeakBuffer = (char*)calloc(PeakBufferSize * 8 + 8, 1);
				Peaks = (float*)calloc(PeakBufferSize * 2, sizeof(float));
				FilteredPeaks = (float*)calloc(PeakBufferSize * 2, sizeof(float));
			}
			
			Trail = (PeakCount % 3);
			if (!(PeakCount % 3))
			{
				PeakBuffer[PeakCount * 32/3] = '\0';
			}
			else
			{
				PeakBuffer[(PeakCount * 32/3) + Trail + 1] = '\0';
			}
	
		//	cout << "dd " << spec_counter << " " << ScanNumber << endl;

			b64_decode_mio( DecodedPeakBuffer, PeakBuffer);
			for (FloatIndex = 0; FloatIndex < (2 * PeakCount); FloatIndex++)
			{
		#ifdef BYTEORDER_LITTLE_ENDIAN
				if (!ByteOrderLittle)
				{
					char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
					DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
					DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
					ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
					DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
					DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
				}
				memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
		#else
				if (ByteOrderLittle)
				{
					char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
					DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
					DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
					ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
					DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
					DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
				}
				memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
		#endif
			}


			// add spectrum
			int i;
			for (i=0; i<FloatIndex; i+=2)
			{
				cout << fixed << setprecision(3) << Peaks[i] << " " << Peaks[i+1] << endl;
			}

			exit(0);
			// join and filter peaks
		}
		else
			scan_start_ptr = ScanStr +50;

	}

	fclose(MZXMLFile);
	cout << spec_counter << " spectra..." << endl << endl;

	return spec_counter;
}




int DAT_Converter::parse_annotated_spectra_from_single_MZXML(
								Config *config, 
								string& mzxml_name, 
								int file_idx,
								map<mzXML_annotation,int>& ann_map)
{
	static char* Buffer = NULL;
    static char* DecodedPeakBuffer = NULL;
    static float* Peaks = NULL;
	static float* FilteredPeaks = NULL;
	static int PeakBufferSize = 0;
    int Trail;
    static char* PrecursorStr;
    int FloatIndex;
    char* ByteOrderStr;
    int ByteOrderLittle = 1;

    int BytesToRead;
    int BufferStartPos = 0;
    int BytesRead;
    int BufferEnd = 0;
    FILE* MZXMLFile;
    int ParseState = 0;
    int FilePos = 0;
    

    // allocate
	if (! Buffer)
		Buffer = (char*)calloc(XML_BUFFER_SIZE + 1, sizeof(char));

    MZXMLFile = fopen(mzxml_name.c_str(), "rb");
    if (!MZXMLFile)
    {
        cout << "Error: Can't open MZXML file " <<  mzxml_name << endl;
        exit(1);
    }

    printf("File idx %d , '%s'...\n",file_idx, mzxml_name.c_str());

	int idx_in_file = 0;
	int spec_counter =0;
	char *scan_start_ptr = NULL;
    while (1)
    {
		char* ScanStr;
		char* ScanNumberStr;
		int ScanNumber;
		char* MSLevelStr;
		int MSLevel;
		char *PrecursorStr;
		mass_t PrecursorMZ;
		char *retentionTimeStr;
		float retentionTime;
		char *precursorIntensityStr;
		float precursorIntensity;
		char* PeakStr;
		char* PeakBuffer;
		int  BufferPos;

        // Read more data, to fill up the buffer:
	
		if ( ! scan_start_ptr || 
			( (Buffer + BufferEnd - scan_start_ptr) < XML_BUFFER_HALF_SIZE) )
		{
			// try shunt half of the buffer
			if (scan_start_ptr)
			{
				if (BufferEnd - XML_BUFFER_HALF_SIZE>0)
				{
					memmove(Buffer, Buffer + XML_BUFFER_HALF_SIZE, BufferEnd - XML_BUFFER_HALF_SIZE);
					BufferEnd -= XML_BUFFER_HALF_SIZE;
					scan_start_ptr -= XML_BUFFER_HALF_SIZE;

//					cout << "MOVED!" << endl;
				}
			}
			else
				scan_start_ptr = Buffer;

			BytesToRead = XML_BUFFER_SIZE - BufferEnd;
			BytesRead = fread(Buffer + BufferEnd, sizeof(char), BytesToRead, MZXMLFile);


			if (BytesRead<5)
				break;

			BufferEnd += BytesRead;
			Buffer[BufferEnd] = '\0';


			FilePos += BytesRead;
		}

        // Look for a new <scan tag opening:
		// this scan cannot be done with strstr since there might be NULL termination
		const char *last_pos = Buffer + BufferEnd - 5;
		char *pos = scan_start_ptr;

		while (++pos<last_pos)
		{
			if (*pos != '<')
				continue;

			if (*(pos+1)=='s' && *(pos+2)=='c' && *(pos+3)=='a' && *(pos+4)=='n')
				break;
		}
		ScanStr =  (pos<last_pos) ? pos : NULL;

        if (ScanStr)
        {
			// if this is the case, read over more data to avoid the case where
			// the spectrum's record is not all in the buffer
			if (scan_start_ptr - Buffer > XML_BUFFER_HALF_SIZE)
			{
				scan_start_ptr = ScanStr-2;
				continue;
			}

            BufferPos = ScanStr - Buffer;
        }
        else
        {
            BufferPos = 0;
        }

        if (!ScanStr )
        {
			scan_start_ptr = Buffer + BufferEnd-5;
            continue;
        }

        ScanNumberStr = strstr(ScanStr, "num=");
        if (!ScanNumberStr)
        {
          //  printf("** Warning: mzXML parser encountered a scan with no scan number!  File %s Pos %d\n", 
		  //		   mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

            ScanNumber = -1;
        }
        else
        {
            ScanNumber = ParseIntFromXML(ScanNumberStr);
			scan_start_ptr = ScanNumberStr;
        }

		retentionTimeStr = strstr(ScanStr,"retentionTime=\"PT");
		if (! retentionTimeStr)
		{
			retentionTime = -1;
		}
		else
		{
			retentionTime = ParseMassFromXML(retentionTimeStr);
		}


		char *PeakCountStr = strstr(ScanStr, "peaksCount=\"");
		if (!PeakCountStr)
		{
			cout << "Warning: bad parsing peaks in mzxml. " << endl;
			cout << "Scan: " << ScanNumber << " Pos: " << FilePos << endl;
			scan_start_ptr += 50;
			continue;
		}
		int PeakCount = ParseIntFromXML(PeakCountStr);
		
        MSLevelStr = strstr(ScanStr, "msLevel=");
        if (!MSLevelStr)
        {
            printf("** Warning: mzXML parser encountered a scan with no scan level!  File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

			scan_start_ptr += 50;
			continue; 
        }
        else
        {
            MSLevel = ParseIntFromXML(MSLevelStr);
        }

		if (MSLevel<=1)
		{
			scan_start_ptr += 50;
			continue;
		}


		precursorIntensityStr = strstr(ScanStr,"precursorIntensity=");
		if (! precursorIntensityStr)
		{
			cout << "Warning: no precursor intensity found for scan " << ScanNumber << " pos: " << BufferPos << endl;
			scan_start_ptr += 50;
			continue;
		}
		else
		{
			precursorIntensity = ParseMassFromXML(precursorIntensityStr);
		}


		PrecursorStr = strstr(ScanStr, "<precursorMz");
		if (PrecursorStr)
		{
			PrecursorStr = strstr(PrecursorStr, ">");
			PrecursorMZ = ParseMassFromXML(PrecursorStr);
		}

		if (!PrecursorStr && MSLevel > 1)
		{
			printf("Warning: mzXML parser encountered a scan with no m/z: File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + BufferPos - BufferEnd);

			scan_start_ptr += 50;
			continue;
		}
		

		// check if this spectrum is in the list of annotated spectra
		// if not, don't output it.
		bool is_annotated_spectrum = false;
		map<mzXML_annotation,int>::const_iterator it;
		mzXML_annotation ann_pos;
		ann_pos.mzXML_file_idx = file_idx;
		ann_pos.scan = ScanNumber;

		it = ann_map.find(ann_pos);
		if (it != ann_map.end())
			is_annotated_spectrum = true;

		//  check that this is a good spectrum to output
		if (MSLevel>1 && PeakCount>2 && is_annotated_spectrum)
		{

			// read peaks

			PeakStr = strstr(PrecursorStr, "<peaks");
			if (PeakStr)
			{
				// Get byte order:
				ByteOrderStr = strstr(PeakStr, "byteOrder=\"");
				if (ByteOrderStr)
				{
					ByteOrderStr += 11;
					if (!strncmp(ByteOrderStr, "network", 7))
					{
						ByteOrderLittle = 0;
					}
					if (!strncmp(ByteOrderStr, "big", 3))
					{
						ByteOrderLittle = 0;
					}
					if (!strncmp(ByteOrderStr, "little", 6))
					{
						ByteOrderLittle = 1;
					}
				}
				PeakStr = strstr(PeakStr, ">");
			}
			if (!PeakStr)
			{
				cout << "Warning couldn't find peaks tag for scan: " << ScanNumber << " Buffer pos: " << BufferPos << "  skipping..." << endl;
				scan_start_ptr += 50;
				continue;
			}

			PeakStr++;
			PeakBuffer = PeakStr;

			if (PeakBufferSize < PeakCount)
			{
				if (DecodedPeakBuffer)
				{
					char *dbf = DecodedPeakBuffer;
					free(DecodedPeakBuffer);
					DecodedPeakBuffer = NULL;
					free(Peaks);
					Peaks = NULL;
					free(FilteredPeaks);
					FilteredPeaks=NULL;
				}
				PeakBufferSize = (int)(PeakCount*1.5);
				DecodedPeakBuffer = (char*)calloc(PeakBufferSize * 8 + 8, 1);
				Peaks = (float*)calloc(PeakBufferSize * 2, sizeof(float));
				FilteredPeaks = (float*)calloc(PeakBufferSize * 2, sizeof(float));
			}
			
			Trail = (PeakCount % 3);
			if (!(PeakCount % 3))
			{
				PeakBuffer[PeakCount * 32/3] = '\0';
			}
			else
			{
				PeakBuffer[(PeakCount * 32/3) + Trail + 1] = '\0';
			}

			b64_decode_mio( DecodedPeakBuffer, PeakBuffer);
			for (FloatIndex = 0; FloatIndex < (2 * PeakCount); FloatIndex++)
			{
		#ifdef BYTEORDER_LITTLE_ENDIAN
				if (!ByteOrderLittle)
				{
					char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
					DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
					DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
					ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
					DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
					DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
				}
				memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
		#else
				if (ByteOrderLittle)
				{
					char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
					DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
					DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
					ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
					DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
					DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
				}
				memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
		#endif
			}


			// add spectrum
			int DAT_file_idx =  (int)(PrecursorMZ/mass_increment);
			if (DAT_file_idx > max_dat_file_idx)
				DAT_file_idx = max_dat_file_idx;

			if (! dat_buffs[DAT_file_idx].ind_was_initialized)
			{
				ostringstream os, os_batch;
				os << DAT_file_idx; 
				os_batch << batch;
				string path = out_dir + "/" + name + "_" + os_batch.str() + "_" + os.str() + ".dat";
				dat_buffs[DAT_file_idx].init(path,dat_buff_size);
			}

			// join and filter peaks
			int new_num_peaks = join_and_filter_peak_list(config,PrecursorMZ,Peaks,
														  PeakCount,FilteredPeaks);
			if (new_num_peaks>3)
			{
				dat_buffs[DAT_file_idx].add_spec_to_DAT_file(PrecursorMZ, 0, file_idx, ScanNumber, retentionTime,  
									 precursorIntensity, new_num_peaks, (char *)FilteredPeaks);
				spec_counter++;
			}
			scan_start_ptr = PeakStr + 8 * PeakCount;
		}
		else
			scan_start_ptr = ScanStr +50;
	}

	fclose(MZXMLFile);
	cout << spec_counter << " spectra..." << endl << endl;

	return spec_counter;
}





/*******************************************************************************
Create DAT files only for annoated spectra
********************************************************************************/
void DAT_Converter::create_dat_files_for_anns(Config *config, char *mzXML_list, char *anns_file,
								char *_out_dir, char *_name)
{
	map<mzXML_annotation,int> ann_map;
//	read_mzXML_annotations_to_map(anns_file,ann_map);

	vector<string> list;
	read_paths_into_list(mzXML_list,list);

	init_DAT_Converter((float)2000.0,(float)20.0,524288);

	name	= _name;
	out_dir = _out_dir;
	batch	= 0;

	int num_spectra_parsed=0;
	int f;
	for (f=0; f<list.size(); f++)
	{

		num_spectra_parsed += parse_annotated_spectra_from_single_MZXML(
								config, 
								list[f],
								f,
								ann_map);
	}

	int d;
	for (d=0; d<dat_buffs.size(); d++)
	{
		if (dat_buffs[d].ind_was_initialized && dat_buffs[d].pos > dat_buffs[d].buff)
			dat_buffs[d].flush_buff();
	}
	
	
	cout << "Total spectra extracted and converted to DAT: " << num_spectra_parsed << endl;

	create_list_file();
}


