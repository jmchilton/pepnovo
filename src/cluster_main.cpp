#include "FileManagement.h"
#include "QuickClustering.h"
#include "AdvancedScoreModel.h"
#include "includes.h"

void print_help()
{

printf("\n\nMass spectra clustering tool. Created by Ari Frank (arf@cs.ucsd.edu).\n");
printf("Copyright 2007, The Regents of the University of California\nAll Rights Reserved.\n\n");
 
printf("program prameters:\n");
printf("------------------\n");
printf("-list <list_file>  - list of full path to input files\nor\n");
printf("-dat_list <dat_list_file> - use existing dat files\nor\n"); 
printf("-pkl_dir_list <pkl_list_file> - use a list of pkl cpick_in directories.\n\n");
 
printf("optional parameters:\n");
printf("--------------------\n");
printf("-no_normalize        - do not normalize the peaks intensities.\n");
printf("-make_dat_only       - only makes dat files (from the input files)\n");
printf("-use_spectrum_charge - only cluster spectra with the same charge\n");
printf("-tolerance    <xx>  - peak tolerance in Da. (default xx=0.5 Da.)\n");
printf("-slice_width  <xx>  - the width of the m/z clustering window (default xx=2.5 Da.)\n");
printf("-similarity   <xx>  - the minimal similarity for clustering (default xx=0.5, values that can be used 0.1-0.999)\n");
printf("-sim_peaks    <xx>  - the number of peaks to use for similarity calculations per 1000 Da. of mass (default xx=15)\n");
printf("-min_m_over_z <xx>  - minimal m/z to cluster (default xx=0)\n");
printf("-max_m_over_z <xx>  - maximial m/z to cluster (default xx=10000)\n");
printf("-min_size     <xx>  - the minimal cluster size (default xx=1)\n");
printf("-max_size     <xx>  - the maximal cluster size (default xx=10000)\n");
printf("-tmp_dir     <xx>  - directory where the temporary files are written (default xx=\".\\tmp\")\n");
printf("-out_dir      <xx>  - directory where the cluster files are written (default xx=\".\\out\")\n");
printf("-name         <xx>  - name to be given cluster files (default xx=\"Cluster\")\n");
printf("-batch        <xx>  - the index given to this batch of clusters (default xx=0)\n");
printf("-file_idx_start <xx>    - the number to add to the file indices (used if dat files are created from split lists).\n");
printf("-specs_per_slice <xx>   - number of expected spectra per slice (default xx = 25000 is suitable for up to 5 million spectra)\n");
printf("-filter_model_name <xx> - name of model PMC and SQS model files (should reside in \".\\Models\", currently only CID_IT_TRYP has these models)\n");
printf("-model_dir       <path> - directory where model files are kept (default ./Models)\n\n");
printf("-min_filter_prob   <xx> - the minimal probability score to keep a spectrum ([0-1], default xx=0.1)\n");
printf("-assign_charges         - if filtering is performed, then a charge is typically asigned to each cluster, using this flag adds the CHARGE field to the outputted MGF\n");
printf("-output_peak_density <xx> - the number of peaks per 100 Da. that will be in the outputted clusters' spectra, default xx=8.\n");
printf("-mass_exclude_range <xx yy> - mass ranges to exclude for similarity computation (e.g.,for iTraq) use \"-mass_exclude_range 113.5 118.5\"\n");
printf("-use_spectrum_header	- sets the TITLE= field in the consensus spectrum to the TITLE of the first spectrum in the cluster.\n");

exit(1); 
}


int main(int argc, char **argv)
{
	Config default_config;
	
	int i;
	
	char list_file[256];
	char tmp_dir[256];
	char out_dir[256]; 
	char name[256];
	char dat_list[256];
	char anns_file[256];
	char filter_model_name[256];
	char model_dir[256];
	char sim_query_file[256];
	char ptm_line[256];
	char out_mgf_file[256];
	char good_anns_file[256];

	
	mass_t tolerance = 0.5;
	mass_t slice_width = 2.5;
	float  similarity = 0.5;
	mass_t min_m_over_z = 0;
	mass_t max_m_over_z = 10000;
	int    min_size = 1;
	int	   max_size = 10000;
	int	   batch=0;
	int	   max_small_cluster_size = -1;
	int    specs_per_slice = 25000;
	int    file_idx_start = 0;
	int	   output_peak_density=8;

	float filter_prob = 0.1;
	int max_mzxml_idx = -1;
	int k_value = 15; 

	bool need_to_normalize=true;
	bool got_model_name = false;
	bool got_model_dir  = false;
	bool got_list_file  = false;
	bool has_only_mzxml_files = true;
	bool skip_dat_creation  = false;
	bool make_dat_only      = false;
	bool perform_filtering  = false;
	bool perform_sim_query  = false;
	bool extract            = false;
	bool extract_no_process = false;
	bool got_ptms           = false;
	bool ind_is_pkl_dir     = false;
	bool got_good_anns_file = false;
	bool dat_to_mgf         = false;
	bool assign_charges     = false;
	bool use_spectrum_charge = false;
	bool use_spectrum_header = false;

	AdvancedScoreModel model;
	void *pmcsqs_ptr = NULL;

	vector<mass_t> min_exclude_masses, max_exclude_masses;

	min_exclude_masses.clear();
	max_exclude_masses.clear();
	
	sprintf(tmp_dir,"./tmp");
	sprintf(out_dir,"./out");
	sprintf(name,"Cluster"); 

	
	// read command line arguments
	i=1;
	while (i<argc)
	{
		if (!strcmp(argv[i],"-list"))
		{
			if (++i == argc)
				print_help();

			strcpy(list_file,argv[i]);
			got_list_file=true;
		}
		else if (! strcmp(argv[i],"-pkl_dir_list"))
		{
			if (++i == argc)
				print_help();

			strcpy(list_file,argv[i]);
			got_list_file=true;
			ind_is_pkl_dir=true;
		}
		else if (!strcmp(argv[i],"-tolerance"))
		{
			if (++i == argc)
				print_help();

			tolerance = atof(argv[i]);
			if (tolerance<0 || tolerance>1.0)
			{
				printf("Error: illegal value for tolerance (should be between 0-1.0 Da)\n");
				exit(1);
			}
		}
		else if (!strcmp(argv[i],"-slice_width"))
		{
			if (++i == argc)
				print_help();

			slice_width = atof(argv[i]); 
		}
		else if (!strcmp(argv[i],"-similarity"))
		{
			if (++i == argc)
				print_help();

			similarity = atof(argv[i]);
			if (similarity<0.1 || similarity>1.0)
			{
				printf("Error: illegal value for similarity : %f (should be between 0.1-0.999 Da)\n",similarity);
				exit(1);
			}
		}
		else if (!strcmp(argv[i],"-min_m_over_z"))
		{
			if (++i== argc)
				print_help();

			min_m_over_z = atof(argv[i]); 
		}
		else if (!strcmp(argv[i],"-max_m_over_z"))
		{
			if (++i== argc)
				print_help(); 

			max_m_over_z = atof(argv[i]);  
		}
		else if (!strcmp(argv[i],"-min_size"))
		{
			if (++i== argc)
				print_help();

			min_size = atoi(argv[i]);
			if (min_size<1)
			{
				printf("Error: illegal value for min_size (should be  >=1)\n");
				exit(1);
			}
		}
		else if (!strcmp(argv[i],"-max_size"))
		{
			if (++i== argc)
				print_help();

			max_size = atoi(argv[i]);
			if (max_size<1)
			{
				printf("Error: illegal value for max_size (should be  >=1)\n");
				exit(1);
			}
		}
		else if (!strcmp(argv[i],"-max_mzxml_idx"))
		{
			if (++i== argc)
				print_help();

			max_mzxml_idx = atoi(argv[i]);
		}
		else if (!strcmp(argv[i],"-sim_peaks"))
		{
			if (++i== argc)
				print_help();

			k_value = atoi(argv[i]);
			if (k_value<1)
			{
				printf("Error: illegal value for sim_peaks (should be  >=1)\n");
				exit(1);
			}
		}
		else if (!strcmp(argv[i],"-batch"))
		{
			if (++i== argc)
				print_help();

			batch = atoi(argv[i]);
		 	if (batch<0)
			{
				printf("Error: illegal value for batch (should be  >=0)\n");
				exit(1);
			}
		}
		else if (!strcmp(argv[i],"-extract"))
		{
			if (++i== argc)
			{
				printf("\nERROR: use:  -extract <ann_file> <out_file> -list <mzXML list>\n");
				print_help();
			}

			strcpy(anns_file,argv[i]);

			if (++i== argc)
			{
				printf("\nERROR: use:  -extract <ann_file> <out_file> -list <mzXML list>\n");
				print_help();
			}

			strcpy(out_mgf_file,argv[i]); 

	
			extract=true;
		}
		else if (!strcmp(argv[i],"-dat_to_mgf"))
		{
			if (++i== argc)
			{
				printf("\nERROR: use:  -dat_to_mgf <ann_file> <out_dir> <name> -dat_list <DAT list>\n");
				print_help();
			}

			strcpy(anns_file,argv[i]);

			if (++i== argc)
			{
				printf("\nERROR: use:  -dat_to_mgf <ann_file> <out_dir> <name> -dat_list <DAT list>\n");
				print_help();
			}

			strcpy(out_dir,argv[i]); 

			if (++i== argc)
			{
				printf("\nERROR: use:  -dat_to_mgf <ann_file> <out_dir> <name> -dat_list <DAT list>\n");
				print_help();
			}

			strcpy(name,argv[i]); 

			dat_to_mgf =true;
		}
		else if (! strcmp(argv[i],"-extract_no_process"))
		{
			extract_no_process = true;
		}
		else if (!strcmp(argv[i],"-good_anns_file"))
		{
			if (++i== argc)
			{
				printf("\nERROR: must supply annotations file\n");
				print_help();
			}

			strcpy(good_anns_file,argv[i]);
			got_good_anns_file=true;
		}
		else if (!strcmp(argv[i],"-ptms"))
		{
			if (++i== argc)
				print_help();

			strcpy(ptm_line,argv[i]);
			got_ptms=true;
		}
		else if (!strcmp(argv[i],"-tmp_dir"))
		{
			if (++i== argc)
				print_help();

			strcpy(tmp_dir,argv[i]);
		}
		else if (!strcmp(argv[i],"-out_dir"))
		{
			if (++i== argc)
				print_help();

			strcpy(out_dir,argv[i]);
		}
		else if (!strcmp(argv[i],"-name"))
		{
			if (++i== argc)
				print_help();

			strcpy(name,argv[i]);
		}
		else if (! strcmp(argv[i],"-make_dat_only"))
		{
			make_dat_only = true;
		}
		else if (!strcmp(argv[i],"-no_normalize"))
		{
			need_to_normalize=false;
		}
		else if (!strcmp(argv[i],"-dat_list"))
		{
			if (++i== argc)
				print_help();

			strcpy(dat_list,argv[i]);
			skip_dat_creation= true; 
		}
		else if (!strcmp(argv[i],"-max_small_cluster_size"))
		{
			if (++i== argc)
				print_help();

			max_small_cluster_size = atoi(argv[i]);
			if (max_small_cluster_size<2)
			{
				printf("Error: max_small_cluster_size must be at least 2!\n");
				exit(1);
			}
		} 
		else if (!strcmp(argv[i],"-specs_per_slice"))
		{
			if (++i== argc)
				print_help();

			specs_per_slice = atoi(argv[i]);
		}
		else if (!strcmp(argv[i],"-file_idx_start"))
		{
			if (++i== argc)
				print_help();

			file_idx_start = atoi(argv[i]);
		}
		else if (! strcmp(argv[i],"-filter_model_name")) 
		{
			if (++i == argc)
				print_help();
			strcpy(filter_model_name,argv[i]);

			got_model_name = true;
		}
		else if (! strcmp(argv[i],"-model_dir"))
		{
			if (++i == argc)
				print_help();

			strcpy(model_dir,argv[i]);

			got_model_dir = true;
		}
		else if (! strcmp(argv[i],"-min_filter_prob"))
		{
			if (++i == argc)
				print_help();

			filter_prob = atof(argv[i]);
			if (filter_prob>1.0)
			{
				printf("Error: filter probability should be in range [0-1]\n");
				exit(1);
			}
		}
		else if (! strcmp(argv[i],"-assign_charges"))
		{
			assign_charges = true;
		}
		else if (! strcmp(argv[i],"-use_spectrum_charge"))
		{
			use_spectrum_charge=true;
		}
		else if (! strcmp(argv[i],"-use_spectrum_header"))
		{
			use_spectrum_header=true;
		}
		else if (! strcmp(argv[i],"-mass_exclude_range"))
		{
			mass_t min_exclude_mass=NEG_INF;
			mass_t max_exclude_mass=NEG_INF;

			if (++i == argc)
				print_help();

			min_exclude_mass = atof(argv[i]);

			if (++i == argc)
				print_help();

			max_exclude_mass = atof(argv[i]);

			if (min_exclude_mass == NEG_INF || max_exclude_mass == NEG_INF ||
				min_exclude_mass > max_exclude_mass)
			{
				printf("Error: the flag \"-mass_exclude_range\" should be followd by two numbers (e.g., \"-exclude_mas_range 113.5 118.5\")\n");
				exit(1);
			}

			min_exclude_masses.push_back(min_exclude_mass);
			max_exclude_masses.push_back(max_exclude_mass);
		}
		else if (! strcmp(argv[i],"-output_peak_density"))
		{
			if (++i == argc)
				print_help();

			output_peak_density = atoi(argv[i]);

			if (output_peak_density<1)
			{
				printf("Error: output peak density should be at least 1!\n");
				exit(1);
			}
		}
		else if (! strcmp(argv[i],"-sim_query")) 
		{
			if (++i == argc)
				print_help();

			strcpy(sim_query_file,argv[i]);

			perform_sim_query = true;
		}
		else
		{
			printf("Unkown command line option: %s\n\n",argv[i]);
			print_help();
			exit(0);
		}
		i++;
	}

	if (! got_list_file && ! skip_dat_creation)
	{
		print_help();
		exit(0);
	}

	// make sure output directories exist
	string tmp_test = string(tmp_dir) + "/test.txt";
	FILE *tmp_stream = fopen(tmp_test.c_str(),"w");
	if (! tmp_stream)
	{
		printf("Error: could not write to %s, please make sure directory exists!\n", tmp_dir);
		exit(1);

	}
	unlink(tmp_test.c_str());


	string out_test = string(out_dir) + "/test.txt";
	FILE *out_stream = fopen(out_test.c_str(),"w");
	if (! out_stream)
	{
		printf("Error: could not write to %s, please make sure directory exists!\n", out_dir);
		exit(1);
	}
	unlink(out_test.c_str());

	time_t start_time = time(NULL);

	// init model if needed
	if (got_model_dir)
	{
		model.get_config()->set_resource_dir(model_dir);
	}

	if (got_model_name)
	{
		model.read_model(filter_model_name,true);

		if (! model.get_ind_pmcsqs_was_intialized())
		{
			printf("Error: PMC and SQS model files were not found for this model!\n");
			exit(1);
		}
		pmcsqs_ptr=(void *)model.get_pmcsqs_ptr();
	}
	else
	{
		if (got_model_dir)
			default_config.set_resource_dir(model_dir);
		default_config.init_with_defaults(); 
	}


	Config *config = (pmcsqs_ptr ? model.get_config() : &default_config);
	
	config->set_need_to_normalize(false);
	config->setTolerances(tolerance);
	config->setPrecursorMassTolerance(slice_width*0.5);
	config->set_max_number_peaks_per_local_window(output_peak_density*2);

	for (i=0; i<min_exclude_masses.size(); i++)
		config->add_exclude_range(min_exclude_masses[i],max_exclude_masses[i]);
	
	if (! got_ptms)
		config->apply_selected_PTMs("M+16:C+57:Q-17");

	config->set_use_spectrum_charge(0);	
	if (use_spectrum_charge)
		config->set_use_spectrum_charge(1);
	
	
	ClusterSpectrum::init_statics(config);

	if (extract)
	{
		if (got_ptms)
			config->apply_selected_PTMs(ptm_line);

		if (skip_dat_creation)
		{	// we got a dat list
			extract_annotate_scans_from_dat(config,dat_list, anns_file, out_mgf_file);
		}
		else
			extractAnnoatedScansFromFiles(config, list_file, anns_file,
									  out_mgf_file, extract_no_process);
		exit(0);
	}

	if (dat_to_mgf)
	{
		if (got_ptms)
			config->apply_selected_PTMs(ptm_line);
		convert_dat_to_mgf(config, dat_list, name, out_dir, anns_file);
		exit(0);
	}
	
	
	if (! skip_dat_creation)
	{
		DAT_Converter dat;

		dat.init_DAT_Converter((float)2000.0,(float)20.0,524288);

		dat.convert_files_to_DAT_on_the_fly(config, list_file, tmp_dir, name, batch,
				min_m_over_z, max_m_over_z, file_idx_start, ind_is_pkl_dir);
		

		ostringstream oss;
		oss << batch;
		string batch_str = oss.str(); 

		sprintf(dat_list,tmp_dir); 
		strcat(dat_list,"/");
		strcat(dat_list,name);
		strcat(dat_list,"_");
		strcat(dat_list,batch_str.c_str());
		strcat(dat_list,"_list.txt");

		if (make_dat_only)
		{
			cout << "Made dat files, list is at:" << endl;
			cout << dat_list << endl;
			exit(0);
		}
	}


	// regular clustering
	if (1)
	{
		
		char *good_anns_ptr = NULL;
		if (got_good_anns_file)
			good_anns_ptr = (char *)good_anns_file;


		cluster_full_dataset(config, 
							 dat_list , 
							 out_dir , 
							 name ,
							 batch, 
							 specs_per_slice, 
							 min_m_over_z, 
							 max_m_over_z,  
							 similarity, 
							 min_size, 
							 max_size, 
							 false, 
							 max_small_cluster_size, 
							 k_value, 
							 pmcsqs_ptr,
							 ind_is_pkl_dir, 
							 filter_prob,
							 assign_charges,
							 max_mzxml_idx,
							 good_anns_ptr);
	}

	time_t current_time = time(NULL);
	double total_time = current_time - start_time;
	cout << endl << "Total running time: " << setprecision(2) << fixed << setw(8) << total_time << endl;

	return 0;
}
