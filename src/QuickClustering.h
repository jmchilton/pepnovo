#ifndef __QUICKCLUSTERING_H__
#define __QUICKCLUSTERING_H__

#include "FileManagement.h"
#include "includes.h"


#define XML_BUFFER_SIZE 2097152
#define XML_BUFFER_HALF_SIZE 1048576

#define LARGISH_CLUSTER_SIZE 6
#define NUM_CLUSTERS_PER_FILE 10000
#define DAT_BUFF_SIZE 10000000
#define DAT_FILE_INCREMENT 25.0


#define NUM_TOP_CLUSTER_PEAKS 5

#define MAX_SIZE_FOR_SQS_REP 4


// values at which the cluster should be remade
const int cluster_reset_values[]={2,3,4,6,9,15,30,50,100,200};
const int num_cluster_reset_values = sizeof(cluster_reset_values)/sizeof(int);


struct MassRankPair {
	
	bool operator< (const MassRankPair& other) const
	{
		return (rank<other.rank);
	}

	mass_t mass;
	int    rank;
};


struct QCPeak {
	QCPeak() : mass(-1), intensity(0), scaled_intensity(0), adjusted_inten(0),
			   num_occurences(1), max_num_occurences(0), source_spec_idx(-1) {};

	mass_t mass;
	intensity_t intensity;
	intensity_t scaled_intensity;
	float   adjusted_inten;  // holds the adjusted value of the peak intensity for
	int	    num_occurences;
	int		max_num_occurences;
	int	    source_spec_idx;
};


struct BasicSpectrum {
public:
	BasicSpectrum() : peaks(NULL), prm_peaks(NULL), num_peaks(-1), num_prm_peaks(-1),
		ssf(NULL), squared_adjusted_inten(-1),	signal_level(0) {};

	void output_to_mgf(ostream& mgf, const Config *config, const char *seq=NULL) const;

	void output_to_mgf(FILE* mgf_stream, const Config *config, const char *seq=NULL) const;

	// returns number of peaks that match the given masses
	int get_number_of_matching_peaks(mass_t tolerance, const vector<mass_t>& masses) const;

	float calc_signal_level();

	void  calc_peak_isotope_levels(mass_t tolerance, vector<float>& iso_levels) const;
	void  mark_all_possible_isotope_peaks(mass_t tolerance, vector<bool>& iso_inds) const;

	bool  select_strong_peak_idxs(const vector<float>& iso_levels, vector<bool>& indicators) const;
	
	void print_peaks() const;

	QCPeak *peaks;
	QCPeak *prm_peaks;

	int     num_peaks;
	int		num_prm_peaks;
	SingleSpectrumFile *ssf;
	float   squared_adjusted_inten;
	float	signal_level;
};





// Similar to FileSet but has minimum overhead and works with BasicSpectra
class BasicSpecReader {
public:

	BasicSpecReader() : max_peak_list_size(1000), 
						mgf_stream(NULL), current_mgf_file_idx(-1),
						mzxml_stream(NULL), current_mzxml_file_idx(-1),
						current_dat_file_idx(-1), current_ms2_file_idx(-1) {};

	~BasicSpecReader() { if (mgf_stream) fclose(mgf_stream);
						 if (mzxml_stream) fclose(mzxml_stream);}
						 

	// Reads the important info from the single spectrum
	// pretty much does what the get_next_spectrum() does with the FileSet
	// but without much of the overhead. Returns number of peaks read (after
	// joining close adjacent peaks)
	int read_basic_spec(const Config *config, 
						const FileManager& fm, 
						SingleSpectrumFile *ssf, 
						QCPeak* peaks,
						bool override_file_idx = false,
						bool no_processing = false);

private:

	int max_peak_list_size;

	FILE *mgf_stream;          // the current MGF file being scanned (its open stream)
	int current_mgf_file_idx;  // the file index of the current mgf that is open

	FILE *mzxml_stream;          // the current MZXML file being scanned (its open stream)
	int current_mzxml_file_idx;  // the file index of the current MZXML that is open

	ifstream dat_file;
	int current_dat_file_idx;

	FILE *ms2_stream;
	int current_ms2_file_idx;

	vector<QCPeak> peak_list;      // used for temporary storage of a spectrum's peak list

	// these functions just extract the peak list from the spectrum file, return the actual
	// number of peaks (after joining)
	int get_peak_list_from_DTA(const char* dta_name);
	int get_peak_list_from_MGF(FILE *mgf_stream);
	int get_peak_list_from_MZXML(FILE *mzxml_stream);
	int get_peak_list_from_DAT(ifstream& dat_file, QCPeak *peaks);
	int get_peak_list_from_MS2(FILE *ms2_stream);
	int get_peak_list_from_PKL(const string& pkl_path);
};



struct PeakListPointer {
	QCPeak *peaks;
	int num_peaks;
};

struct CutProb {
	bool operator< (const CutProb& other) const
	{
		return (mass < other.mass);
	}

	float mass;
	float prob;
};

class ClusterSpectrum {
	friend class QCOutputter;
public:
	ClusterSpectrum() : tmp_cluster_idx(-1), sim_matrix_row_start(NULL), charge(0), config(NULL), m_over_z(-1), retention_time(-1),
		                tolerance(-1), 
						maximum_peaks_vector_size(0), maximum_good_peaks_to_output(0),
						num_spectra_in_cluster(0), best_sqs_spec_idx(-1) {};


//	void create_new_cluster(Config *config, BasicSpectrum& bs,  int cluster_idx);
	
	// adds the spectrum to this cluster. If the total number of spectra 
	// is one of the "reset" values (2,3,4,6,9,15,30,50,100), the cluster is
	// recreated.
//	void add_spectrum_to_cluster(BasicSpectrum& bs, 
//								 const vector<int>& spec_top_idxs,
//								 float top_x_masses[NUM_TOP_CLUSTER_PEAKS]);

	// tries to add the cluster
	// succeeds only if the similarity of the two originals to the new consensus
	// is above the sim_tresh (returns true if it made the addition, false otherwise)
//	bool add_cluster(ClusterSpectrum& cs, float sim_thresh);

	void filter_peaks_with_slidinig_window();
	
	const vector<BasicSpectrum>& get_basic_spectra() const { return basic_spectra; }

	// sets the m_over_z as the average of the m_over_z of the basic_spectra
//	void set_cluster_m_over_z();

	// checks that all basic spectra have the same charge, otherwise cahrge =0
//	void set_charge();
//	void set_charge(int c) { charge=c; }

	void set_title( const string& new_title) { title = new_title; }
	int  get_num_basic_spectra() const { return basic_spectra.size(); }
	int  get_tmp_cluster_idx() const { return tmp_cluster_idx; }
	void set_tmp_cluster_idx(int idx) { tmp_cluster_idx=-1; }

	unsigned char *get_sim_matrix_row_start() const { return sim_matrix_row_start; }
	void set_sim_matrix_row_start(unsigned char *pos) { sim_matrix_row_start = pos; }


	float get_retention_time() const { return retention_time; }
	void  set_retention_time(float r) { retention_time = r; }

	mass_t get_m_over_z() const { return m_over_z; }

	const string& get_title() const { return title; }
	BasicSpectrum& get_basic_spectrum(int idx) { return basic_spectra[idx]; }
	void set_basic_spectra(vector<BasicSpectrum>& new_spectra) { basic_spectra = new_spectra; }

	const vector<int>& get_top_ranked_idxs() const { return top_ranked_peak_idxs; }

	void set_top_ranked_idxs(const vector<int>& top_idxs) { top_ranked_peak_idxs = top_idxs; }

	void set_top_masses(float top_x_masses[NUM_TOP_CLUSTER_PEAKS]) 
	{
		int i;
		for (i=0; i<NUM_TOP_CLUSTER_PEAKS; i++)
			top_peak_masses[i]=top_x_masses[i];
	}


	float* get_top_peak_masses() const { return (float *)top_peak_masses; }

	bool find_match_in_top_masses(float top_x_masses[NUM_TOP_CLUSTER_PEAKS]) const
	{
		int a_idx=0;
		int b_idx=0;
		while (1)
		{
			float diff = top_x_masses[a_idx] - top_peak_masses[b_idx];
			if (diff<0.65 && diff>-0.65)
				return true;

			if (diff<0)
			{
				if (++a_idx==NUM_TOP_CLUSTER_PEAKS)
					return false;
				continue;
			}

			if (++b_idx==NUM_TOP_CLUSTER_PEAKS)
				return false;
		}
	}

	const QCPeak * get_peaks_pointer() const { return &peaks[0]; }
	int   get_num_peaks() const { return peaks.size(); }

	const string& get_peptide_str() const { return peptide_str; }
	void set_peptide_str(string& str) { peptide_str = str; }

	Config *get_config() { return config; }
	void    set_config(Config *c) { config = c; }

	int     get_charge() const { return charge; }

	// joins the spectra to form a new cluster spectrum
	void create_cluster_by_binning_basic_spectra();

	


	static void init_statics(Config *config)
	{
//		init_min_num_occurences(config->get_tolerance());
	}

//	static void set_num_top_peaks_per_1000_da(int val) { num_top_peaks_per_1000_da = val; }
//	static int get_num_top_peaks_per_1000_da()  { return num_top_peaks_per_1000_da; }


	// sets the consensus spectrum to be the basic spectrum with the maximal similarity
	// to other spectra. Returns the idx of the most similar spectrum
	int select_max_similarity_spectrum_as_consensus();

	int select_max_sqs_spectrum_as_consensus();

	// creates the title (file name) for the cluster
	void make_title(string& name, int batch_idx, int cluster_idx);

	// writes the spectrum to the output file in the mgf format
	void write_spectrum_to_mgf(ostream& mgf, bool write_peak_count=false, bool ind_write_charge=false) const;

	void write_spectrum_to_pkl_single(string& file_name) const;


	// finds how many spectra have a peptide sequence that doesn't match the majority assignment
//	int  get_num_misassigned_spectra(mass_t pm_tolerance= 0.08,
//		 int *mismatched_with_pep = NULL) const;

	// returns true if there is a mjority annotation (with 75% of the annotated spectra)
	// if so, returns its string and mass in the variables
	bool has_majority_annotation(string& pep_str, mass_t& pep_mass) const;

	int  get_num_basic_spectra_with_peptide() const;

//	void print_cluster_peptides() const;
	
//	void print_cluster_similarities();

	void print_explained_intensity_stats(Peptide& pep) const;


private:
	int	   tmp_cluster_idx; // the idx given when the clustering is done, -1 means invalid cluster
	unsigned char  *sim_matrix_row_start;  // the address of this clusters entry in the sim_matrix
	int    charge;

	Config *config;  

	mass_t      m_over_z;
	float       retention_time; // the average time for the spectra in the cluster
	mass_t      tolerance;

	int		    maximum_peaks_vector_size;
	int		    maximum_good_peaks_to_output;
	int			num_spectra_in_cluster;

	int			best_sqs_spec_idx; // holds the basic spectrum idx for the spectrum with the
								   // highest sqs score, otherwise holds -1 (if the cluster is
								   // bigger than MAX_SIZE_FOR_SQS_REP

	string title;
	string peptide_str; // if the cluster as an assigned peptide

	float top_peak_masses[NUM_TOP_CLUSTER_PEAKS]; // the added spectrum needs to match one of these 
								             

	vector<QCPeak>        peaks;
	vector<BasicSpectrum> basic_spectra;
	vector<int>           top_ranked_peak_idxs;  // idxs of top ranked peaks (sorted list by idx)


	void create_consensus_by_binning_basic_spectra();
	
	// adds the given list to the clusters current list.
	// if the new cluster is larger than X, peak scores are weighted according
	// to the percentage in which the peaks appear
	// list of peaks is then filtered to remove excess peaks
	bool add_peak_list(const QCPeak *second_peaks, 
					   int num_second_peaks, 
					   mass_t tolerance,
					   int num_basic_spectra_added, 
					   bool need_to_scale = true);

	void create_consensus_sepctrum_from_peak_list_pointers(
				vector<PeakListPointer>& plp, 
				int total_num_peaks);

	// recursively merges the peak lists from the various spectra
	// performs a merge until the pointers list has only one entry
	void merge_peak_lists(vector<QCPeak>& org_peaks,
						  vector<QCPeak>& new_peaks,
						  vector<PeakListPointer>& pointers);


	// joins adjacent peaks, markes invaldiated peaks by assigning their mass to -1
	// and weeds them out
	void join_merged_peak_lists(PeakListPointer& plp,
								PeakListPointer& alt_plp,
								int num_merged_spectra,
								mass_t tolerance);


	// selects the consensus peaks - those that appear more than the expected cutoff
	// also takes some of the stronger peaks that didn't make the cutoff
	// writes the selected peaks into the peaks of the cluster spectrum
	void select_consensus_peaks(PeakListPointer& plp, 
								PeakListPointer& alt_plp,
								int num_org_spectra);


	// static
	static void increase_tmp_storage_size(int num_peaks);
//	static void init_min_num_occurences(mass_t tolerance);
	

	static vector<QCPeak> tmp_peak_area1;  // used in peak merging
	static vector<QCPeak> tmp_peak_area2;
    
	static vector<int>   min_num_occurences;
//	static int   num_top_peaks_per_1000_da;
	static float large_num_occurence_ratios;
	static int   large_cluster_size;

};




bool mark_top_peaks_with_sliding_window(const QCPeak *peaks, int num_peaks, mass_t window_size, 
					    int num_peaks_per_window, vector<bool>& indicators);


/************************************************************
// Functions for the dot product similarity distance
*************************************************************/
void select_top_peak_idxs(QCPeak *peaks, int num_peaks, 
						  mass_t m_over_z, mass_t tolerance, 
						  vector<int>& top_ranked_peak_idxs,  
						  float top_x_masses[NUM_TOP_CLUSTER_PEAKS]=NULL,
						  int top_peaks_per_100da = 20,
						  Config *config = NULL);


// sets the adjusted intensity of the peaks
void set_adjusted_inten(QCPeak *peaks, int num_peaks);

// Sets the adusted intensity of the peaks that are in top_ranked_idxs
// the spectrum is reduced only to these peaks, so the intensity given to each peak
// is the ration I/I_total where I_total is for all the peaks in the top_ranked_idxs
void set_top_ranked_peak_inten(QCPeak *peaks, int num_peaks, const vector<int>& top_ranked_idxs, 
							   vector<float>& top_ranked_peak_inten);


float calc_sum_adjusted_inten_squared(const QCPeak *peaks, int num_peaks);

float calc_sum_adjusted_inten_squared(const QCPeak *peaks, int num_peaks, 
									  const vector<int>& top_ranked_peak_idxs);




float calc_selected_dot_prod(mass_t tolerance, 
							 const QCPeak *pa, int na, const vector<int>& peak_idxs_a,
	  					     const QCPeak* pb, int nb, const vector<int>& peak_idxs_b,
							 bool verbose = false);


void collect_dot_product_stats(char *list_file);


void dot_prod_exp();



// Outputs clusters 
class QCOutputter {
public:
	QCOutputter() : init_flag(false), name("clusters"), 
		dir("."), batch_idx(0), spectra_counter(0), total_spectra_counter(0), file_counter(0) {};

	~QCOutputter();
	
	void init(string _name , string _dir, int _batch_idx=-1,
			  mass_t min_m_over_z=POS_INF, mass_t max_m_over_z=NEG_INF,
			  float min_similarity=-1, int min_cluster_size=0);

	void output_basic_spectrum_to_mgf(BasicSpectrum &bs, Config *config);

//	void output_cluster_spectrum(ClusterSpectrum& cs, bool ind_write_charge=false);

//	void output_cluster_spectrum_as_single_pkl(ClusterSpectrum& cs);

	void output_cluster_anns(ClusterSpectrum& cs);


private:
	bool   init_flag;
	
	string name;
	string dir;
	string mgf_name;
	string cluster_file_name;
	string batch_str;

	int    batch_idx;
	int    spectra_counter;
	int    total_spectra_counter;
	int	   file_counter;
	
	fstream mgf_stream;           // mgf containing cluster consensus spectra
	fstream cluster_file_stream;  // file listing foreach cluster what single spectra belong to it 
	fstream summary_stream;               // for each consensus spectrum lists how many single spectra belong to it
	fstream file_list_stream;             // a list of all the mgf files created
	fstream anns_stream;                  // stream for annotations 
};





class DAT_FileBuff {
	friend class DAT_Converter;
public:

	DAT_FileBuff()  : buff(NULL), max_pos(NULL), pos(NULL), counter(0),
					  ind_first_write(1), ind_was_initialized(0) {};
	~DAT_FileBuff();

	void init(string& _path, int buff_size); // initializes ans allocates buffer memory

	// copies the files to the DAT file
	void add_spec_to_DAT_file(mass_t m_over_z, int charge, int mzxml_file_idx,
							  int scan_number, float retention_time, 
							  float precursor_intensity, int num_peaks, char *peak_buff);

	void flush_buff(); // writes buff to file
private:

	string path;
	char *buff;
	char *max_pos;
	char *pos;

	int  counter;
	int  ind_first_write;      // is this the first write (if not, append)
	int  ind_was_initialized;  
};


class DAT_Converter {
public:

	DAT_Converter() : max_m_over_z((mass_t)2000.0), mass_increment((mass_t)DAT_FILE_INCREMENT),
		dat_buff_size(DAT_BUFF_SIZE), max_dat_file_idx(-1), ind_was_initialized(0),
		batch(0) {};

	~DAT_Converter()
	{
		int d;
		for (d=0; d<dat_buffs.size(); d++)
			if (dat_buffs[d].ind_was_initialized && dat_buffs[d].pos > dat_buffs[d].buff)
				dat_buffs[d].flush_buff();
	}

	void init_DAT_Converter(mass_t _max_m_over_z, mass_t _mass_increment, 
			  int dat_buff_size);

	void convert_files_to_DAT_on_the_fly(Config* config, char *file_list, 
							char * _out_dir, char * _name, int _batch, 
							mass_t min_m_over_z, mass_t max_m_over_z, int file_start_idx,
							bool ind_is_pkl_dir);



	int convert_PKL_dir_to_DAT(Config* config, char *file_list, int file_start_idx,
							mass_t min_m_over_z, mass_t max_m_over_z);


	int convert_single_non_MZXML_file_to_DAT(Config* config, string file, 
							mass_t min_m_over_z, mass_t max_m_over_z, int file_idx);

	int parse_annotated_spectra_from_single_MZXML(
								Config *config, 
								string& mzxml_name, 
								int file_idx,
								map<mzXML_annotation,int>& ann_map);

	void create_dat_files_for_anns(Config *config, char *mzXML_list, char *anns_file,
								char *_out_dir, char *_name);



	// creates a file with the list of DAT files
	void create_list_file() const;

private:
	string out_dir;
	string name;

	mass_t max_m_over_z;
	mass_t mass_increment;
	int    dat_buff_size;
	int    max_dat_file_idx;
	int    batch;

	int    ind_was_initialized;

	vector<DAT_FileBuff> dat_buffs;

	int parse_single_MZXML_file(Config *config, string& mzxml_name, int file_idx,
								mass_t min_m_over_z, mass_t max_m_over_z);
};

/////////////////////////////////////////////////////////////////////////////
// For mzXML parsing
struct MassInten {
	float mass, intensity;
};


int join_and_filter_peak_list(const Config *config, 
							  mass_t m_over_z, 
							  float *org_peaks, 
							  int num_org_peaks, 
							  float *new_peaks);
/////////////////////////////////////////////////////////////////////////////




/***************************************************************************
	This function creates clusters from a list of files containing spectra
	(possibly different file types).
	The cluster spectra are outputted as mgf files in the output dir (x spectra
	per file). In addition, for each cluster file there is a map file that holds
	the indices (position in list, and idx in file) of the original spectra
	that are part of the cluster.
****************************************************************************/
void cluster_full_dataset(Config *config,
							  char *list_file,
							  const string& out_dir,
							  const string& clust_name,
							  int batch_idx,
							  int max_spec_per_slice = 20000,
							  mass_t min_m_over_z = 0,
							  mass_t max_m_over_z = 10000,
							  float  min_similarity = 0.8, 
							  int min_cluster_size = 2,
							  int max_cluster_size = 1000,
							  bool verbose = false,
							  int  max_small_cluster_size = -1,
							  int  k_value = 20,
							  void *pmcsqs = NULL,
							  bool ind_pkl_mode = false,
							  float filter_prob = 0.075,
							  bool assign_charges = false,
							  int max_mzxml_idx = -1,
							  char *good_anns_file = NULL);


/**************************************************************************
Given a query spectrum finds the x number of spectra that have the
highest similarity to the query
***************************************************************************/
void find_spectra_with_max_similarity_to_query(
							  char *model_name,
							  char *dat_list,
							  char *query_spectrum,
							  int num_top_spectra);


int cluster_spec_in_file_set(Config *config, const FileManager& fm, FileSet& cluster_fs,
							  mass_t tolerance, 
							  QCPeak *basic_peaks,
							  vector<ClusterSpectrum>& clusters, 
							  float min_similarity,
							  int   max_small_cluster_size,
							  int	max_cluster_size,
							  int   num_top_peaks_per_1000_da,
							  bool verbose,
							  void *pmcqsqs = NULL,
							  float filter_prob = 0.075,
							  bool	assign_charges = false,
							  map<mzXML_annotation,int> *ann_map_ptr = NULL);



int add_additional_spectra_to_existing_clusters(Config *config, const FileManager& fm, 
							  FileSet& additional_fs, mass_t tolerance, QCPeak *basic_peaks, 
							  vector<ClusterSpectrum>& clusters, 
							  float min_similarity,
							  int max_cluster_size = 100000,
							  void *pmcqsqs = NULL,
							  float filter_prob = 0.075,
							  bool  ind_assign_charges = false,
							  map<mzXML_annotation,int> *ann_map_ptr = NULL,
							  bool verbose = false);



// reading mzXML, mgf files,  annotations




bool read_mgf_file_into_basic_spectra(Config *config, char *mgf_file, QCPeak *basic_peaks,
									  vector<BasicSpectrum>& basic_spectra);

void ann_mgf_and_create_mgf(Config *config, char *annotations_file, char *mgf_list_file,
							  char *out_dir_name, char *file_prefix, bool output_only_ann_spectra = false);


void read_mzXML_annotations(char *mzXML_list,char *ann_file, vector< vector<int> >& annotation_idxs, 
							vector<mzXML_annotation>& annotations, int max_ann_size = 50000);

//void read_mzXML_annotations_to_map(char *ann_file, map<mzXML_annotation,int>& ann_map);

void read_mzXML_annotations_limited(char *mzXML_list, 
							char *ann_file, 
							vector< vector<int> >& annotation_idxs, 
							vector<mzXML_annotation>& annotations);

void extract_spectra_stats_from_mzXML(Config *config, char *annotations_file, 
									  char *mzXML_list_file, char *stat_file);

void read_annotated_dataset_into_clusters(Config *config, FileManager& fm,
				const vector<SingleSpectrumFile *>& all_ssf, char *ann_mgf, 
				int num_specs_per_cluster, int max_num_clusters, vector<ClusterSpectrum>& clusters);

void benchmark_similarity_measures(Config *config, char *ann_mgf, 
								   int num_specs_per_cluster);

void benchmark_inter_similarity_vs_outer_similarity(Config *config, char *ann_mgf, 
								   int num_specs_per_cluster);

void benchmark_similarity_to_consensus(Config *config, char *ann_mgf, 
								   int num_specs_per_cluster);

void benchmark_signal_to_noise(void *model, char *ann_mgf, 
								   int num_specs_per_cluster);

void benchmark_signal(void *model, char *ann_mgf, 
								   int num_specs_per_cluster);
  


void print_specs_from_mgf(Config *config,char *mgf_name);



void qc_exp();

void qc_ann_exp(char *name, bool verbose = false);

void benchmark_large_clusters(void *model_ptr, char *ann_mgf, 
								   int num_specs_per_cluster);

void create_spectra_mgf_file(void *model_ptr, char *ann_mgf, 
								   int num_specs_per_cluster);

void create_file_with_rt_scores(void *model_ptr, char *ann_mgf, 
								   int num_specs_per_cluster);

// checks that the anntoated spectra have a correct m_over_z
void print_specs(Config *config, char *list_name);

void read_ms2(Config *config, char *ms2_list);

void ann_mgf_and_create_mgf_with_sim_masses(Config *config, 
											char *annotations_file, 
											char *org_mgf_list_file,
											char *good_peptide_list,
											char *out_dir_name, 
											char *file_prefix);

void make_specified_benchmark_clustering_dataset(
				Config *config, 
				char *ann_mgf_list,
				mass_t min_m_over_z, 
				mass_t max_m_over_z,
				char *out_dir_name, 
				char *file_prefix,
				int num_clusters,
				int cluster_size,
				int number_non_assigned_per_assigned);

void make_benchmark_clustering_dataset(Config *config, char *ann_mgf_list,
				mass_t min_m_over_z, mass_t max_m_over_z,bool only_confident_anns,
				char *out_dir_name, char *file_prefix);

void benchmark_clustering_performance(Config *config,
							  char *list_file,
							  int  k_value);

void benchmark_k_value(Config *config, char *list_file);

void benchmark_heuristic_filtering(Config *config, char *dat_list);

void benchmark_top7_and_sim_thresh(Config *config, char *dat_list, char *ann_file);

void benchmark_retention_thresh(Config *config, char *dat_list, char *ann_file);

void find_pair_similarities(Config *config, char *mgf_file, char *pair_file);

void extractAnnoatedScansFromFiles(Config *config, char *file_list, char *anns, 
								   char *output_file, bool extract_no_process);

void extract_annotate_scans_from_dat(Config *config, char *dat_list, char *anns_file, 
									 char *out_mgf_file);

void create_annotated_mgf_from_dat(Config *config, 
								   char *dat_list,
								   char *mzxml_list,
								   char *anns_file,
								   char *output_file);

void convert_dat_to_mgf(Config *config, 
						char *dat_list, 
						char *out_name, 
						char *out_dir,
						char *anns_file = NULL);

void find_best_similar_pairs(char *model_name, char *mgf1, char *mgf2, int num_pairs);
void find_self_similarity(char *model_name, char *mgf1, char *mgf2);
void find_similar_pairs_ditrib(char *model_name, char *mgf1, char *mgf2);
void find_homeometric_similarity_distrib(char *model_name, char *mgf1, char *mgf2);
void find_self_similarity_ranges(char *model_name, char *mgf1);
void peptide_distances();
void find_matches_similarity_distrib(char *model_name, char *mgf1, char *mgf2);

/*void create_annotated_mgf_from_mgf(Config *config, 
								   char *mgf_list,
								   char *anns_file,
								   char *output_file);*/


void test_sims();



void check_fs_for_missing_anns_and_test_sqs(Config *config,
											const vector<SingleSpectrumFile *>& ssfs,
											const FileManager& fm,
											void *pmcqsqs,
											char *anns_file);


#endif

