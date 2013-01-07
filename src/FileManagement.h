#ifndef __FILEMANAGEMENT_H__
#define __FILEMANAGEMENT_H__

#include "Config.h"
#include "Spectrum.h"
#include "BasicDataStructs.h"

#include "includes.h"

#define DTA   1
#define MGF   2
#define MZXML 3
#define DAT   4
#define PKL   5
#define MS2   6
#define TXT   7

// parses the type from the file name, -1 if not recognized
int get_file_extension_type(const string& fname);

struct mzXML_annotation {
	mzXML_annotation() : charge(0), mzXML_file_idx(-1), scan(-1), retention_time(-1.0) {};

	bool operator< (const mzXML_annotation& other) const
	{
		return (mzXML_file_idx<other.mzXML_file_idx ||
			    (mzXML_file_idx == other.mzXML_file_idx &&
				 scan < other.scan));
	}

	int charge;
	int mzXML_file_idx;
	int scan;
	float retention_time;
	string pep;
};






void read_paths_into_list(const char *list_file, vector<string>& list);

///////////////////////////////////////////////////////////////////////////
//
//     INPUT FILE MANAGEMENT
//

//  class for all types of files
struct SingleSpectrumFile {
	SingleSpectrumFile() : 	 org_pm_with_19(-1), pm_with_19(-1), m_over_z(-1), 
		charge(0), type(MGF), file_idx(-1), file_pos(0),retention_time(-1),
		precursor_intensity(0),
		num_peaks(0), assigned_cluster(-1), ann_idx(-1), sqs(-1), ind_peptide_mass_ok(true) { };

	bool operator< (const SingleSpectrumFile& other) const
	{
		return (file_idx<other.file_idx || (file_idx == other.file_idx &&
				file_pos < other.file_pos) ) ;
	}

	int get_scan() const;
	

	void print_ssf_stats(const Config *config, ostream& os = cout, bool print_endl = true) const;
	
	mass_t org_pm_with_19;
	mass_t pm_with_19;
	mass_t m_over_z;

	int charge;

	int type;  // DTA =1 , MGF =2, MZXML = 3, DAT=4

	int file_idx; // the number of the MGF file (if known) / -1 default for DTA

	long file_pos; // position pointer in stream

	float retention_time;

	float precursor_intensity;

	int num_peaks;
 
	int assigned_cluster; // if >-1 means that this spectrum was already assigned
						  // to a cluster.

	int ann_idx;

	score_t sqs;

	string  single_name;

	Peptide peptide;

	bool ind_peptide_mass_ok;
};




struct DTA_file : public SingleSpectrumFile {

	
	// makes initial read of dta and stores the stats
	void initial_read(const Config *config);

	/****************************************************************
	quickly extracts the charge, sequence, and pm_with_19 from dta, num peaks.
	*****************************************************************/
	bool scan_dta(const string& file_name, const Config *config);
};


struct  MGF_single : public SingleSpectrumFile {
	MGF_single() : scan_number(-1), cluster_size(1), idx_in_file(-1), first_peak_mass(-1) {};

	int scan_number;        // if created from mzXML
	int cluster_size;       // if represents a cluster
	int idx_in_file;
	mass_t first_peak_mass;
	

	/****************************************************************
	quickly extracts the charge, sequence, and pm_with_19 from MGF location
	*****************************************************************/
	bool scan_mgf_single(FILE *stream, const Config *config);
};


struct  MZXML_single : public SingleSpectrumFile {

	MZXML_single() : scan_number(-1), MS_level(-1),
					 retention_time(-1), peak_buff_start_idx(-1) {};

	int scan_number;
	int MS_level;
	float retention_time;   //retentionTime="PT575.021S"
	int	  peak_buff_start_idx;  // in case peak lists were stored in the mzXML file


	/****************************************************************
	quickly extracts the charge...
	*****************************************************************/
	bool scan_mzxml_single(FILE *stream, const Config *config)
	{
		cout << "Error: scan_mzxml_single not implemented!" << endl;
		exit(1);
	}
};

struct  DAT_single : public SingleSpectrumFile {
	DAT_single() : scan_number(-1), mzxml_file_idx(-1)  {};

	int scan_number;
	int mzxml_file_idx;
	


	/****************************************************************
	quickly extracts the charge, sequence, and pm_with_19 from DAT location
	*****************************************************************/
	bool scan_dat_single(FILE *stream, const Config *config)
	{
		cout << "Error: scan_DAT_single not implemented!" << endl;
		exit(1);
	}
};


struct PKL_single : public SingleSpectrumFile {
	PKL_single() : scan_number(-1),  retention_time(-1) {};

	int   scan_number;
	float retention_time;

	bool scan_pkl_single(const string& file, const Config& config)
	{
		cout << "Error: scan_pkl_single not implemented!" << endl;
		exit(1);
	}
};


struct  MS2_single : public SingleSpectrumFile {
	MS2_single() :  idx_in_file(-1) {};

	int idx_in_file;		
	

	/****************************************************************
	quickly extracts the charge, sequence, and pm_with_19 from MGF location
	*****************************************************************/
	bool scan_ms2_single(FILE *stream, const Config *config);
};






struct MGF_file {
	MGF_file() :  total_num_spectra(0), min_spec_mass(1E7), max_spec_mass(0),
					min_charge(9999), max_charge(0) { }

	// makes initial read of mgf and stores the stats

	void initial_read(const Config *config, int file_idx , bool quick_flag = false);

	int total_num_spectra;
	mass_t min_spec_mass, max_spec_mass;
	int min_charge, max_charge;
	string mgf_name;

	vector<MGF_single> single_spectra;

	vector<int> num_spectra; // number of spectra per charge
};

struct MZXML_file {
	MZXML_file() :  total_num_spectra(0), min_spec_mass(1E7), max_spec_mass(0),
					min_charge(9999), max_charge(0), file_peak_buff_pos(0) { }

	// makes initial read of mgf and stores the stats

	void initial_read(const Config *config, int file_idx );

	int extract_peak_lists_from_mzXML(const Config *config, 
								  string& mzxml_name, 
								  int file_idx,
								  mass_t min_m_over_z, 
								  mass_t max_m_over_z);

	int    total_num_spectra;
	mass_t min_spec_mass, max_spec_mass;
	int    min_charge, max_charge;

	int file_peak_buff_pos;

	string mzxml_name;

	vector<MZXML_single> single_spectra;

	vector<int> num_spectra; // number of spectra per charge

	vector<float> file_peak_buff; // for storing all peak info of a single file
};


struct DAT_file {
	DAT_file() :  total_num_spectra(0), min_spec_mass(1E7), max_spec_mass(0),
					min_charge(9999), max_charge(0) { }

	// makes initial read of DAT and stores the stats
	void initial_read(const Config *config, int file_idx );

	int    total_num_spectra;
	mass_t min_spec_mass, max_spec_mass;
	int    min_charge, max_charge;
	string dat_name;

	vector<DAT_single> single_spectra;

	vector<int> num_spectra; // number of spectra per charge
};



struct PKL_dir {
	PKL_dir() :  total_num_spectra(0), min_spec_mass(1E7), max_spec_mass(0),
				 min_charge(9999), max_charge(0) { }

	// reads the summary tsv file and stores stats
	void initial_read(const Config *config, int dir_idx, const string& path, const string& tsv_file,
		mass_t min_m_over_z = 0, mass_t max_m_over_z = 9999999);

	int    total_num_spectra;
	mass_t min_spec_mass, max_spec_mass;
	int    min_charge, max_charge;

	string dir_path;
	string tsv_path;

	vector<PKL_single> single_spectra;

	vector<int> num_spectra; // number of spectra per charge
};



struct MS2_file {
	MS2_file() :  total_num_spectra(0), min_spec_mass(1E7), max_spec_mass(0),
					min_charge(9999), max_charge(0) { }

	// makes initial read of mgf and stores the stats

	void initial_read(const Config *config, int file_idx , bool quick_flag = false);

	int total_num_spectra;
	mass_t min_spec_mass, max_spec_mass;
	int min_charge, max_charge;
	string ms2_name;

	vector<MS2_single> single_spectra;

	vector<int> num_spectra; // number of spectra per charge
};





struct FileSet;

// contains the data for the entire file

struct FileManager {
friend struct FileSet;
public:
	FileManager() : config(NULL), min_spec_mass(9999999), max_spec_mass(0),
						min_charge(9999), max_charge(0), total_num_spectra(0) {};


	// returns how many spectra are present in the list file
	// also samples m_over_z values to generate an approximate
	// histogram in case the set of spectra needs to be spilt
	int count_num_spectra(const Config *config, const char* list_file,
						  vector<mass_t>& mass_histogram) const;



	// Inits the FileManager using mass levels (for very large
	// collections of spectra). This initialization uses a quick scan
	void init_from_list_file(const Config *config, const char* list_file,
		mass_t min_m_over_z, mass_t max_m_over_z);

	void init_from_list_file(const Config *config, const char* list_file,
		const vector<bool>& file_indicators);

	// only keeps ssfs of mzXML singles that have an annotation
	void init_from_list_file(const Config *config, const char* list_file,
		const vector< vector<int> >& annotation_idxs);

	// only keeps ssfs of mzXML singles that have an annotation
	void init_from_list_file_and_add_annotations(const Config *config, const char* list_file,
		const vector< vector<int> >& annotation_idxs, vector<mzXML_annotation>& annotations,
		bool read_only_annotated = false);


	void init_from_dat_list_extract_only_annotated(const Config *config, char* dat_list_file,
		char *ann_file);




	// if quick init is used, files are not completely scaneed, the correct_pm is not
	// calculated, and neither is the SQS
	
	void init_from_list_file(const Config *config, const char* list_file, 
		bool quick_flag = true)
	{
		vector<string> list;
		read_paths_into_list(list_file,list);
		init_from_list(config,list,quick_flag);
	}

	void init_from_list(const Config *config, const vector<string>& list,  
		bool quick_flag = true, int file_idx = -1);

	
	void init_from_mgf(const Config *config, const char * mgf_name, 
		bool quick_flag = true)
	{
		vector<string> list;
		list.push_back(string(mgf_name));
		init_from_list(config,list,quick_flag);
	}

	void init_from_file(const Config *config, const char * file_name, 
		bool quick_flag = true)
	{
		if (get_file_extension_type(file_name) != MZXML)
		{
			vector<string> list;
			list.push_back(string(file_name));
			init_from_list(config,list,quick_flag);
		}
		else
			init_and_read_single_mzXML(config,file_name);
	}

	// reads a list with dirs and paths to tsv files
	void init_from_pkl_dir_list(const Config *config, const char *list,
			mass_t min_m_over_z = 0, mass_t max_m_over_z = 99999999);

	void init_from_single_pkl_dir(const Config *config, const string& pkl_dir_path, 
		const string& tsv_file, int pkl_dir_idx, mass_t min_m_over_z, mass_t max_m_over_z);


	void init_and_read_single_mzXML(const Config *config, 
									const char * file_name,
									int file_idx=0,
									mass_t min_m_over_z=0,
									mass_t max_m_over_z=100000);

	void copy_mzxml_peak_buff_ptr(float **ptr) const  { *ptr=NULL;
										if (mzxml_files.size()>0)
											*ptr=(float *)&mzxml_files[0].file_peak_buff[0]; }



	const SingleSpectrumFile * get_dta_ssf(int idx) const { return &dta_files[idx]; }
	const SingleSpectrumFile * get_mgf_ssf(int mgf_file, int idx) const
	{ return &mgf_files[mgf_file].single_spectra[idx]; }


	const MGF_file&   get_mgf_file(int file_idx) const   { return mgf_files[file_idx]; }
	const MZXML_file& get_mzxml_file(int file_idx) const { return mzxml_files[file_idx]; }
	const DAT_file&   get_dat_file(int file_idx) const { return dat_files[file_idx]; }
	const PKL_dir&    get_pkl_dir(int dir_idx) const   { return pkl_dirs[dir_idx]; }
	const MS2_file&   get_ms2_file(int file_idx) const { return ms2_files[file_idx]; }


	mass_t get_min_spec_mass() const { return min_spec_mass; }
	mass_t get_max_spec_mass() const { return max_spec_mass; }

	int get_min_charge() const { return min_charge; }
	int get_max_charge() const { return max_charge; }

	const vector<int>& get_spectra_counts() const { return num_spectra; }
	int get_num_spectra(int charge) const { return num_spectra[charge]; }
	void count_num_spectra();
	
	const string& get_list_name() const { return list_name; }
	void set_list_name(const string& name) { list_name=name; }


	void print_summary_stats() const;

	
private:

	Config *config;

	string list_name;

	vector<DTA_file>   dta_files;
	vector<MGF_file>   mgf_files;
	vector<MZXML_file> mzxml_files;
	vector<DAT_file>   dat_files;
	vector<MS2_file>   ms2_files;
	vector<PKL_dir >   pkl_dirs;

	vector<int> num_spectra; // number of spectra for each charge
	int total_num_spectra;

	mass_t min_spec_mass, max_spec_mass;

	int  min_charge, max_charge;
};


// This data structure is used to access the files in the FileManager
struct FileSet {
	FileSet() : mgf_stream(NULL), current_mgf_file_idx(-1), 
				next_ssf_pointer(0), min_mass(-1), max_mass(-1),
				min_charge(9999), max_charge(0) {};

	~FileSet() { if (mgf_stream ) fclose(mgf_stream); }

	// reads the next spectrum into spec
	// returns false if no more spectra are available
	bool get_next_spectrum(const FileManager& fm, Config *config, Spectrum *spec, 
						   SingleSpectrumFile **ssf = NULL, bool perform_init_spectrum=true,
						   bool set_charge_to_zero = false); 


	// selects all SSFs from the fm
	// if the remove duplicates is true, checks that previos SSF soesn't have
	// same m_over_z and number of peaks, if so, ignores it
	void select_all_files(const FileManager& fm, bool remove_duplicates=false);



	// select only spectra with these conditions
	// charge = 0 is all charges
	// reads each spec and calculates corrected pm and SQS
	void select_files(const FileManager& fm, mass_t min_pm_with_19, mass_t max_pm_with_19,
		              score_t min_sqs, score_t max_sqs, 
					  int charge = 0,
					  bool only_unassigned = true);

	void select_files_in_mz_range( const FileManager& fm, 
						   mass_t min_mz, 
						   mass_t max_mz, 
						   int charge=0);

	void randomly_reduce_ssfs(int n);

	void remove_spectra_with_PTMs();

	void filter_dat_spectra_by_mzxml_idx(int max_mzxml_idx);



	void sort_according_to_m_over_z();


	// removes all ssf without a peptides
	void keep_only_spectra_with_peptides();



	// copies the ssf pointers from another FileSet
	void init_from_another_fs(const FileSet& other_fs, int start_ssf_idx, int end_ssf_idx);



	void reset_pointers() { current_mgf_file_idx = -1;
							next_ssf_pointer = 0;  
							if (mgf_stream) fclose(mgf_stream); 
							}



	

	// creates an mgf file with the desred number of spectra per charge as designated
	// in the spectra_per_charges vector(strating from charge 0). so the maximum
	// number of from each charge in the outputted in the mgf file will be at most
	// the numbers desginated in the vector
	void create_mgf_file(const FileManager& fm, Config *config, const char *file_name,
						 vector<int> spectra_per_charges);




	
	// iterates over the files and makes a fasta out of the seq
	//  puts 10 in a row, calls it TRUE_X 
	void make_fasta_from_file_seqs(const FileManager& fm, Config *config, 
		int inc = 10, ostream& os=cout);

	int get_total_spectra() const { return ssf_pointers.size(); }

	void print_file_stats() const;


	void print_summary() const;

	const vector<SingleSpectrumFile *>& get_ssf_pointers() const { return ssf_pointers; }
	vector<SingleSpectrumFile *>& get_non_const_ssf_pointers() { return ssf_pointers; }

	int get_min_charge() const { return min_charge; }
	int get_max_charge() const { return max_charge; }

private:

	vector<SingleSpectrumFile *>  ssf_pointers;

	FILE *mgf_stream;          // the current MGF file being scanned (its open stream)
	
	int current_mgf_file_idx;  // the file index of the current mgf that is open

	int next_ssf_pointer;      // the idx of the next single spectrum file pointer that will be used
	
	mass_t min_mass, max_mass; // for all the selected spectra

	int min_charge,max_charge;

        
};





int ParseIntFromXML(char* AttributeString);

mass_t ParseMassFromXML(char* AttributeString);

void examine_mgf_file(char *mgf_file, int *max_num_spectra, int *max_num_lines);


void make_training_mgf(Config *config, char *list_file, 
					   int target_num_spectra, char *out_mgf_name);

// concatonates several mgf files into one large file
void concat_mgf_files(const Config *config, const string& mgf_file_list, 
					  const string& big_mgf_file);



void extractMZFromFiles(Config *config, char *file_list, char *output_file);






#endif


