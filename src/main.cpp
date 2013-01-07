#include "AdvancedScoreModel.h"
#include "FileManagement.h"
#include "Fragmentation.h"
#include "DeNovoDp.h"
#include "QuickClustering.h"
#include "EdgeModel.h"
#include "FragmentSelection.h"
#include "MSBlast.h"
#include "PMCSQS.h"
#include "PeakRankStats.h"
#include "PeptideRankScorer.h"
#include "DeNovoSolutions.h"
#include "auxfun.h"
#include "includes.h"



int main()
{
	Config *config=NULL;  
	FileManager fm;
	FileSet fs;
	AdvancedScoreModel model; 
	EdgeModel edge_model;
	PeptideRankScorer drs;
	 
	

	rand_seed(112233);


//	train_all(); 

	model.read_model("CID_IT_TRYP");
	config = model.get_config(); 
	config->apply_selected_PTMs("C+57:M+16:Q-17");

//	fm.init_from_file(config,"C:\\Work\\msms5\\DnvScore\\all_ds\\HEK_98_3_unique_30.mgf");
	fs.select_files_in_mz_range(fm,300,2000,3);
	fs.randomly_reduce_ssfs(100);
	
	vector<int> cc(4,0);
	cc[3]=100;
	fs.create_mgf_file(fm,config,"HEK_4_30e.mgf",cc);
	exit(0);

//	model.get_config()->apply_selected_PTMs("C+57:M+16:Q-17"); 


	PMCSQS_Scorer *pmcsqs = (PMCSQS_Scorer *)model.get_pmcsqs_ptr();

	fm.init_from_list_file(config,"C:\\Work\\msms5\\PepNovoHQ\\train3.txt");

	pmcsqs->benchmark_pm_selection(config,fm,0.3);

	exit(0);

/*	benchmark_shew(model,"C:\\Work\\msms5\\PepNovoHQ\\Shew_test_10.mgf");
	exit(0);

		drs.set_model_type(0);
	drs.read_denovo_rank_scorer_model("C:\\Work\\msms5\\PepNovoHQ\\Models\\DBSCORE\\DBSCORE_rank_model.txt");
	drs.give_de_novo_and_peak_match_examples("C:\\Work\\msms5\\DnvScore\\all_db_hits",
									 "C:\\Work\\msms5\\DnvScore\\seq_freqs\\sequences",
									 "C:\\Work\\msms5\\DnvScore\\dnv_full_parts",
									 "C:\\Work\\msms5\\DnvScore\\dicty2_all.txt",
									 2,2);
	exit(0);
	
	drs.read_denovo_rank_scorer_model("C:\\Work\\msms5\\PepNovoHQ\\Models\\DBSCORE\\DBSCORE_rank_model.txt");
	drs.rescore_inspect_results("C:\\Work\\msms5\\DnvScore\\inspect_res\\H293-40ul-08.mzXML",
								"C:\\Work\\msms5\\DnvScore\\inspect_res\\H293-40ul-08.txt",
								"C:\\Work\\msms5\\DnvScore\\inspect_res\\H293-40ul-08_new.txt");

	exit(0);

//	test_denovo_integrity(model,"C:\\Work\\msms5\\DnvScore\\all_ds\\Dicty_98_2_unique_8.mgf", 20000, 8); 

//	benchmark_ranking_on_denovo("C:\\Work\\msms5\\PepNovoHQ\\Models\\DNVSCORE5\\DNVSCORE5_rank_model.txt",
//		"C:\\Work\\msms5\\DnvScore\\test\\DNVSCORE4_test_10.mgf",400,10); // 

//	benchmark_ranking_on_full_denovo("C:\\Work\\msms5\\PepNovoHQ\\Models\\DNVFULL\\DNVFULL_rank_model.txt",
//		"C:\\Work\\msms5\\DnvScore\\test\\FULL_test_10.mgf",1000);
//	exit(0); 



	fm.init_from_list_file(config,//"C:\\Work\\msms5\\DnvScore\\short2_train_mgf_list.txt");
	 "C:\\Work\\msms5\\DnvScore\\comp2_train_mgf_list.txt");
	// "C:\\Work\\msms5\\NewScore\\lists\\Shew_98_3_unique_mgf_list.txt");
	fs.select_files(fm,0,2500,-1,-1,2);
	
	find_special_PTM_frags_using_offset_counts("S",fm,fs.get_ssf_pointers(),&model,2);

	exit(0); 
	
	drs.read_denovo_rank_scorer_model("C:\\Work\\msms5\\PepNovoHQ\\Models\\DNVSC_RANK\\LTQ_DNVRANK_model.txt");
	drs.test_model("C:\\Work\\msms5\\DnvScore\\test_sets\\LTQ_DNVRANK_test_10.mgf",2000);

	drs.train_denovo_partition_model("C:\\Work\\msms5\\DnvScore\\all_db_hits",
									 "C:\\Work\\msms5\\DnvScore\\seq_freqs\\sequences",
									 "C:\\Work\\msms5\\DnvScore\\comp_all_parts",
								//	 "C:\\Work\\msms5\\DnvScore\\short2_train_mgf_list.txt",
									 "C:\\Work\\msms5\\DnvScore\\comp2_train_mgf_list.txt",
									 2,
									 1,
									 30000,
									 5); 




//	model.read_model("ETD");
//	config = model.get_config();
//	config->apply_selected_PTMs("M+16:Q-17:N+1:C+57");

//	fm.init_from_list_file(config,"C:\\Work\\msms5\\PepNovoHQ\\ETD2\\ETD_unique_train.txt");
//	model.full_train_model("ETD",fm,0.5);
	
//	model.train_pmc_rank_models("C:\\Work\\msms5\\PepNovoHQ\\ETD2\\ETD_all_train.txt");
//	model.write_model();
	exit(0); 

//	train_all();

//	create_training_sets();
//	exit(0);
//	generate_size_reports();
//	test_sims();
//	data_set_stats();

//	convert_list_to_trianing_peptide_file(
//		"C:\\Work\\msms5\\NewScore\\lists\\Dicty_98_3_unique_mgf_list.txt",
//		"C:\\Work\\msms5\\NewScore\\tps\\Dicty_98_3_unique_tps.txt");

//	proline_cleavage_reports("b",2); 
//	exit(0);
//	center_cleavage_reports("y",3);
//	n_terminal_cleavage_reports("y",2);
//	c_terminal_cleavage_reports("y",-2);


//	find_best_similar_pairs("LTQ_LOW_TRYP",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\HEK_pos2.mgf",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\Shew_pos2.mgf",8);

//	find_self_similarity("LTQ_LOW_TRYP",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\HEK_pos2.mgf",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\Shew_pos2.mgf");

//	find_similar_pairs_ditrib("LTQ_LOW_TRYP",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\HEK_pos2.mgf",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\Shew_pos2.mgf");

//	find_homeometric_similarity_distrib("LTQ_LOW_TRYP",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\HEK_pos2.mgf",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\Shew_pos2.mgf");

//	find_self_similarity_ranges("LTQ_LOW_TRYP",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\SHEW18.mgf");

//	peptide_distances();

//	find_matches_similarity_distrib("LTQ_LOW_TRYP",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\HEK_pos2.mgf",
//		"C:\\Work\\msms5\\PepNovoHQ\\pairs\\Shew_pos2.mgf");

//	match_sim_exp();  

//	exit(0);

//	edge_model.train_all_edge_models("C:\\Work\\clust_exp\\LTQ_train2_ann_list.txt","LTQ",2); 
//	saa.train_saa_models("C:\\Work\\clust_exp\\LTQ_train2_ann_list.txt","LTQ",2); 
//	saa.train_saancd_models("C:\\Work\\clust_exp\\LTQ_train2_ann_list.txt","LTQ",2); 
//	daa.train_daa_models("C:\\Work\\clust_exp\\LTQ_train2_ann_list.txt","LTQ",2,0.25); 
//	daa.train_daancd_model("C:\\Work\\clust_exp\\LTQ_train2_ann_list.txt","LTQ",2); 

	
//	dot_prod_exp();   
//	qc_exp();  
//	qc_ann_exp("64068",true);
//	exit(0); 

	config = model.get_config();
	config->init_with_defaults();
	config->apply_selected_PTMs("M+16:Q-17:N+1:C+57");

	fm.init_from_file(config,"C:\\Work\\msms5\\PepNovoHQ\\ETD\\train_etd.mgf");
	model.full_train_model("ETD",fm,0.5);
	exit(0); 

	if (1)
	{
		model.read_model("LTQ_LOW_TRYP");

	//	model.get_config()->apply_selected_PTMs("C+57:M+16");
	//	model.get_config()->init_with_defaults();
		model.get_config()->apply_selected_PTMs("M+16:C+57:Q-17");

	//	model.test_pmc("C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\sqs_train_1.mgf",1);
	//	model.compute_sqs_cum_stats_for_ided("C:\\Work\\msms5\\NewScore\\lists\\all_HEK_mgf_list.txt");
	//	model.compute_sqs_cum_stats_for_crap("C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\crap_list.txt");
	//	model.write_model();

		model.compute_sqs_cum_stats_for_ided("C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40good\\H40good_mgf_list.txt");

	///	model.benchmark_sqs("C:\\Work\\msms5\\PepNovoHQ\\small_list.txt",
	//						"C:\\Work\\msms5\\PepNovoHQ\\small_anns.txt");

	//	model.benchmark_sqs("C:\\Work\\msms5\\PepNovoHQ\\tmp\\H40ul_0_list.txt",
	//						"C:\\Work\\msms5\\PepNovoHQ\\H40ul55_missed.txt");
	//						"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul98_anns.txt");
	exit(0);


	DAT_Converter dat;
	dat.create_dat_files_for_anns(model.get_config(),
								  "C:\\Work\\Data\\Briggs\\HEK293\\40ul_list.txt",
								//  "C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul98_anns.txt",
								  "C:\\Work\\msms5\\PepNovoHQ\\H40ul55_missed.txt",
								  "C:\\Work\\msms5\\PepNovoHQ\\tmp\\",
								  "H4055");
//	
		
	//	model.train_pmc_rank_models( 
		//	"C:\\Work\\msms5\\NewScore\\lists\\HEK_98_1_unique_mgf_list.txt",0);
	//		"C:\\Work\\msms5\\NewScore\\lists\\all_unique_mgf_list.txt",0);
	//		"C:\\Work\\msms5\\NewScore\\lists\\all_HEK_mgf_list.txt",0);

	//	model.write_model();
	//	make_before_and_after_matrices(model.get_config(),"C:\\Work\\msms5\\lists\\mgf10.txt",3,"y");
		exit(0);

		FileManager fm;
		fm.init_from_list_file(model.get_config(),"C:\\Work\\msms5\\lists\\LTQ_train_list.txt");

		model.full_train_model("LTQ_IT_TRYP",fm,0.45);

		model.write_model(); 

	//	model.train_pmc("C:\\Work\\msms5\\lists\\pos_sqs_list.txt");

		vector< vector<float> > weights;
		weights.resize(4);
		weights[1].resize(3,0);
		weights[2].resize(3,0);
		weights[3].resize(3,0);
		weights[1][0] = 0.1; weights[1][1] = 0.1;  weights[1][2] = 0.4;
		weights[2][0] = 0.6; weights[2][1] = 0.75; weights[2][2] = 0.5;
		weights[3][0] = 0.3; weights[3][1] = 0.15; weights[3][2] = 0.1;
	
	//	model.train_sqs("C:\\Work\\msms5\\lists\\pos_sqs_list.txt",
	//					"C:\\Work\\msms5\\lists\\neg_sqs_list.txt",&weights);

	//	model.train_sqs("C:\\Work\\msms5\\NewScore\\lists\\all_unique_mgf_list.txt",
	//					"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\crap_list.txt",&weights);
	//	config = model.get_config();

	//	config->set_tolerance(0.5);
//
	//	find_pair_similarities(config,"C:/Work/clust_exp/Results/Shew_bm/ShewBM40_0_1.mgf",
	//		"C:/Work/clust_exp/Results/Shew_bm/ShewBM40_pairs.txt");

		exit(0);
	}

	//	make_y_vectors("C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\sqs10_train_2.mgf",&model);
	//	create_training_files(config);
	//	exit(0);

	if (0)
	{
		PMCSQS_Scorer sqs;


	//	exit(0); 



	//	create_training_files(config);
		exit(0);
	}

	if (1)
	{
	//	fm.init_from_file(model.get_config(),"C:\\Work\\msms5\\PepNovoHQ\\orbi_ann.mgf");
	//	create_MSB_query_for_file_list(fm,&model);

		vector< vector<int> >    annotation_idxs;
		vector<mzXML_annotation> annotations;
		read_mzXML_annotations("C:/Work/Data/Briggs/HEK293_mzxml_list.txt", 
					"C:/Work/ClusterAnn/mzxml_anns3.txt", annotation_idxs, annotations, 35000);

	//	read_mzXML_annotations("C:/Work/ClusterAnn/H40ul_mgf_list.txt", 
	//				"C:/Work/ClusterAnn/mgf_anns.txt", annotation_idxs, annotations, 35000);


		cout << "Read annotations: " << annotations.size() << endl;

		fm.init_from_list_file_and_add_annotations(config,"C:/Work/Data/Briggs/HEK293_mzxml_list.txt",
			annotation_idxs, annotations,true);

	//	fm.init_from_list_file_and_add_annotations(config,"C:/Work/ClusterAnn/H40ul_mgf_list.txt",annotation_idxs,
	//		annotations,true);


		FileSet all_spec_fs;
		all_spec_fs.select_all_files(fm,true);

	//	config->set_need_to_normalize(0);
	//	all_spec_fs.create_MGF_file(fm,config,"C:/Work/ClusterAnn/mgf_spectra.mgf");
	//	exit(0);

		ofstream mgf_stream("C:/Work/ClusterAnn/mzxml_spectra3.mgf",ios::out);
		BasicSpecReader bsr;
		const vector<SingleSpectrumFile *>& ssfs = all_spec_fs.get_ssf_pointers();
		int i;
		for (i=0; i<ssfs.size(); i++)
		{
			static QCPeak peaks[5000];
			BasicSpectrum bs;
			MZXML_single *ssf = (MZXML_single *)ssfs[i];

			bs.peaks = peaks;
			bs.ssf = ssf;
						
			ostringstream oss;
			oss << ssf->file_idx << " " << ssf->scan_number;
			ssf->single_name = oss.str();

			bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
			
			if (ssf->scan_number<0)
			{
				cout << "Error: no scan number read from mzXML!!!" << endl;
				exit(1);
			}

			cout << "scan: " << ssf->scan_number << " " << bs.num_peaks << endl;

			bs.output_to_mgf(mgf_stream,config);
		//	bs.output_to_mgf(cout,&config);
		}
		//all_spec_fs.create_MGF_file(fm,config,"C:/Work/ClusterAnn/mzxml_spectra.mgf");
	//	extractMZFromFiles(model.get_config(),,"C:/Work/Data/Briggs/HEK293/H29340ul_mz.txt");
	

	//
		exit(0);
	}


	if (1) 
	{ 
		model.read_model("LTQ_LOW_TRYP");  
		config = model.get_config();
		config->apply_selected_PTMs("C+57 M+16");
		config->set_tolerances(0.5);
		config->set_pm_tolerance(2.5);
		config->set_digest_type(TRYPSIN_DIGEST);
	
		config->set_max_number_peaks_per_local_window(15);

	//	fm.init_from_list_file(config,"C:\\Work\\clust_exp\\LTQ_train2_ann_list.txt");
	//	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\LTQ_train_list.txt");
	//	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\orbi_train.txt");
	//	edge_model.train_all_edge_models(fm,&model);
	//	tm.train_models(fm,&model); 

	//	fm.init_from_list_file(config,"C:\\Work\\clust_exp\\ShewMGF\\BM2000_ann_list.txt");
	//	make_frag_rank_histogram(fm,config);
	//	exit(0);

	//	benchmark_k_value(config,"C:\\Work\\clust_exp\\ShewMGF\\BM2000_ann_list.txt");

	//	make_benchmark_clustering_dataset(config, "C:\\Work\\clust_exp\\ShewMGF\\AnnsPlus_ann_list.txt",
	//			600, 750, false, "C:\\Work\\clust_exp\\ShewMGF\\", "BMNEW"); 

	//	exit(0);
		benchmark_clustering_performance(config,
			"C:\\Work\\clust_exp\\ShewMGF\\BM2000_ann_list.txt",15);

 
	//	print_dataset_spectra_by_stats(config,"C:\\Work\\clust_exp\\ann_mgf\\Sings_1.mgf");

	//	benchmark_top7_and_sim_thresh(config,"C:\\Work\\clust_exp\\tmp\\H293_40ul_list.txt",
	//		"C:\\Work\\clust_exp\\Results\\BM40ul\\BM40ul_anns.txt");

	//	benchmark_heuristic_filtering(config,"C:\\Work\\clust_exp\\tmp\\H293_40ul_list.txt");
	//	benchmark_retention_thresh(config,"C:\\Work\\clust_exp\\tmp\\H293_40ul_list.txt",
	//		"C:\\Work\\clust_exp\\Results\\BM40ul\\BM40ul_anns.txt");

		exit(0);

		FileManager fm;

	//	fm.init_from_mgf(config,"C:\\Work\\clust_exp\\ShewMGF\\OnlyAnn_1.mgf");


	//	make_specified_benchmark_clustering_dataset(config,"C:\\Work\\clust_exp\\ShewMGF\\both_list.txt",400,1000,
	//		"C:\\Work\\clust_exp\\ShewMGF\\","BM3000",3000,10,0);

		make_benchmark_clustering_dataset(config, "C:\\Work\\clust_exp\\ShewMGF\\AnnsPlus_ann_list.txt",
				800, 1200, true, "C:\\Work\\clust_exp\\ShewMGF\\", "AnnOnly"); 
		exit(0);

		ann_mgf_and_create_mgf_with_sim_masses(config,"K:\\Work\\Data\\Shewenella\\FT_anns.txt",
			"K:\\Work\\Data\\Shewenella\\FT_mgf_list.txt",
			"K:\\Work\\Data\\Shewenella\\FT_peptides.txt",
			"C:\\Work\\clust_exp\\ShewMGF\\",
			"AnnsPlus");

		exit(0);

		ann_mgf_and_create_mgf(config,"C:\\Work\\Data\\FT_mgf\\FT_single_anns.txt", 
			"C:\\Work\\Data\\FT_mgf\\FT_single_mgf.txt",
			"C:\\Work\\clust_exp\\ShewMGF\\",
			"Single",true);

		exit(0);

		ann_mgf_and_create_mgf(config,"C:\\Work\\Data\\FT_mgf\\FT_anns.txt", 
			"C:\\Work\\Data\\FT_mgf\\FT_mgf_list.txt",
			"C:\\Work\\clust_exp\\ShewMGF\\",
			"OnlyAnn",true);

		exit(0);

	//	create_16O_18O_dataset("C:\\Work\\msms5\\lists\\p19_list.txt",config);
	//	exit(0);

	//	config->set_need_to_estimate_pm(0);

		model.clone_charge_model(2,1);
		model.clone_charge_model(2,3);
		model.clone_charge_model(2,4);
		model.clone_charge_model(2,5);

	//	dataset_eval(&model,"C:\\Work\\msms5\\lists\\CAD_376.txt",0.05);
	//	dataset_eval(&model,"C:\\Work\\msms5\\lists\\ann_qtof_list.txt",0.1);
	//	dataset_eval(&model,"C:\\Work\\msms5\\lists\\list280_mgf.txt",0.6);

		

		vector<int>   set_sizes;  
		vector<float> probs;
		denovo_sequencing_and_aa_probs(&model,"C:\\Work\\msms5\\lists\\m280_list.txt",
			set_sizes,probs,2); 

//		denovo_sequencing_and_aa_probs(&model,"C:\\Work\\clust_exp\\LTQ_train2_ann_list.txt",
//			set_sizes,probs,2); 
//
//		output_denovo_results(&model,"C:\\Work\\msms5\\lists\\LTQ-FT_mgf_list.txt");

	//	denovo_sequencing_and_aa_probs(&model,"C:\\Work\\msms5\\lists\\LTQ-FT_mgf_list.txt",
	//		set_sizes,probs, 2);
		exit(0);

	//	print_specs(model.get_config(), "C:\\Work\\msms5\\lists\\one_mzxml.txt");
	//	check_m_over_z(&model,"C:\\Work\\msms5\\lists\\CoCl345sann_ann_list.txt");
	//	calc_parent_mass_tolerance_distribution(&model, "C:\\Work\\msms5\\lists\\ann_mgf_list.txt" , 0.6, 0.98);
	//	calc_tolerance_distribution(&model,"C:\\Work\\msms5\\lists\\ann_mgf_list.txt", 0.6, 0.95);

	//	perfrom_inital_evalutation("C:\\Work\\msms5\\lists\\ann_mgf_list.txt",0.5,2,0.05);

	//	denovo_sequencing_results(&model,"C:\\Work\\msms5\\lists\\CAD_376.txt" ,0.0075);
	//	denovo_sequencing_results(&model,"C:\\Work\\msms5\\lists\\ann_qtof_list.txt",0.1);
	//	denovo_sequencing_results(&model,"C:\\Work\\msms5\\lists\\ann_orbi_list.txt",0.008);

	//	perfrom_inital_evalutation("C:\\Work\\msms5\\lists\\ann_qtof_list.txt",0.2,2,0.05);
	//	perfrom_inital_evalutation("C:\\Work\\msms5\\lists\\ann_orbi_list.txt",0.025,2,0.05);
	//	calc_tolerance_distribution(&model,"C:\\Work\\msms5\\lists\\ann_qtof_list.txt",0.1,0.95);
	//	calc_tolerance_distribution(&model,"C:\\Work\\msms5\\lists\\ann_orbi_list.txt",0.1,0.95);
	//	calc_tolerance_distribution(&model,"C:\\Work\\msms5\\lists\\CAD_376.txt",0.1,0.95);
	//	calc_tolerance_distribution(&model,"C:\\Work\\Data\\Omics04\\omics_ann_list.txt",0.75,0.95);

	//	calc_parent_mass_tolerance_distribution(&model,"C:\\Work\\msms5\\lists\\ann_qtof_list.txt",0.1,0.975);
	//	calc_parent_mass_tolerance_distribution(&model,"C:\\Work\\msms5\\lists\\ann_orbi_list.txt",0.1,0.975);
	//	calc_parent_mass_tolerance_distribution(&model,"C:\\Work\\msms5\\lists\\CAD_376.txt",0.1,0.975);
		exit(0);

	//	fm.init_from_mgf(config,"C:\\Work\\msms4\\PepNovo\\test\\m280.mgf");
	//	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\good_list2.txt");
	//	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\p215.txt");
	//	fm.init_from_list_file(config,"C:\\Work\\Data\\Omics04\\omics_ann_list.txt");
	//	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\omics_mgf_list.txt");
	//	fm.init_from_mgf(config,"C:\\Work\\msms5\\PepNovoHQ\\Omics04Spectra.mgf");

	//	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\omics02_dta.txt");
	//	FileSet fs;
	//	fs.select_all_files(fm);
	//	fs.sort_according_to_m_over_z();
	//	fs.create_MGF_file(fm,config,"Omics02Spectra.mgf");


	//	exit(0);

	//	collect_denovo_statistics(fm,&model);
	//	denovo_histograms(fm,&model);
	//	config->set_tolerance(0.0075);
	//	random_check_homemorphic(config,50000,25);

	//	create_spectrum_clusters(config,"C:\\Work\\msms5\\lists\\Drosophila_list.txt",".","cc",0,5E6,0,1E6);
	//	create_spectrum_clusters(config,"C:\\Work\\msms5\\lists\\Dros_short.txt",".","cc",0,1E6,0,1E6);
//		create_spectrum_clusters(config,"C:\\Work\\msms5\\lists\\all_clust.txt","clust_out","ikkb",0,5E6,0,1E6);


	//	create_spectrum_clusters(config,"C:\\Work\\msms5\\lists\\omics_mgf_list.txt","clust_out",
	//		"Omics04b",0,5E6,0,1E6);

	//	create_spectrum_clusters(config,"C:\\Work\\msms5\\lists\\omics02_mgf_list.txt","clust_out",
	//		"Omics02",0,5E6,0,1E6);
		

	//	create_spectrum_clusters(config,"C:\\Work\\msms5\\PepNovoHQ\\clust_out2\\h293_dat_list.txt",
	//		"clust_out2","H293b_2nd_digest_abd3",0,5E6,1938.76,1E6,2);

	//	create_spectrum_clusters(config,"C:\\Work\\msms5\\PepNovoHQ\\clust_out\\h29s_list.txt","clust_out",
	//		"xxxx",0,5E6,835.397,1E6,2);

		exit(0);



		DAT_Converter dat;
		dat.init_DAT_Converter(2000,25,1048576);

		exit(0);
	}


	
//	config->add_selected_PTMs("C+57 M+16 S+80 T+80 Y+80 N+1 Q+1 K+42 D+16 K+16 P+16 N+16");
//	config->set_tolerances(0.0075);
//	config->set_pm_tolerance(0.011);




	
//	fdb.create_db_from_fasta("C:\\Work\\msms5\\PepNovoHQ\\DB\\contaminants.fasta",config);
//  fdb.create_db_from_fasta("C:\\Work\\msms5\\PepNovoHQ\\DB\\Homo_sapiens.NCBI35.dec.pep.fa",config);
//	fdb.create_db_from_fasta("C:\\Work\\msms5\\PepNovoHQ\\DB\\fa50mb.fa",config,true,5,6);
//	fdb.create_db_from_fasta("C:\\Work\\msms5\\PepNovoHQ\\DB\\homo_pos.fa",config);
//	fdb.print_protein_names();
//	fdb.write_FastaDB("C:\\Work\\msms5\\PepNovoHQ\\DB\\homo.dat");
//	fdb.read_FastaDB("C:\\Work\\msms5\\PepNovoHQ\\DB\\homo.dat",config);
//	fdb.read_FastaDB("C:\\Work\\msms5\\PepNovoHQ\\DB\\qqq.dat",config);
//	fdb.read_FastaDB("C:\\Work\\msms5\\PepNovoHQ\\DB\\fa500k.dat",config);
//	fdb.print_protein_names();
//	fdb.write_FastaDB("C:\\Work\\msms5\\PepNovoHQ\\DB\\fa5--0mb.dat");

//  exit(0);

//	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\CAD_seq_list.txt");
//	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\CAD_u_list.txt");
	fm.init_from_list_file(config,"C:\\Work\\msms5\\lists\\CAD_376.txt");

	

//	make_bin_histograms(fm,config);
//	calc_avg_rand(fm,config);
//	explore_fragment_set(fm,config);
//	show_occurences(fm,config,"p-25.0");
//	find_internal_offset_counts(fm,config);
//	make_frag_rank_histogram(fm,config);
//	exit(0);


//	internal_fragment_annotation(fm,model.get_config());
//	find_internal_offset_counts(fm,model.get_config());
//	exit(0);

	FileSet fs;
	fs.select_all_files(fm);
	fstream ofs("true_seq.txt",ios::out);
	fs.make_fasta_from_file_seqs(fm,config,45,ofs);

//	dbsm.train_regression_models(fm,fdb,12,&model);

//	analyze_cad_spec(fm,config);
//	make_rank_histograms(fm,&model);
//	dbsm.read_model("DBS.txt");
//	CAD_histograms(fm,&model,fdb,&dbsm);
//	rand_db_stats(fm,&model,&dbsm);
//	db_search_stats(fm,&model,&dbsm);
//	neg_db_search_stats(fm,&model,&dbsm);
//	CAD_edge_stats(fm,&model);
//	CAD_denovo_histograms(fm,&model,fdb);
//	CAD_edge_stats(fm,&model);
//	collect_denovo_statistics(fm,&model);
	exit(0);

//	int *arr;
//	int t_size;
//	read_fasta("cshort.fasta",&arr,&t_size,&config);
//	read_fasta("C:\\Work\\msms4\\PepNovo\\Homo_sapiens.NCBI35.dec.pep.fa",&arr,&t_size,config);
	read_fasta("Homo_sapiens.NCBI35.dec.pep.fa",&arr,&t_size,&config);
	read_fasta("Homo_sapiens.NCBI35.dec.pep.fa",&arr,&t_size,&config);
	read_fasta("Homo_sapiens.NCBI35.dec.pep.fa",&arr,&t_size,&config);
	read_fasta("Homo_sapiens.NCBI35.dec.pep.fa",&arr,&t_size,&config);
	read_fasta("Homo_sapiens.NCBI35.dec.pep.fa",&arr,&t_size,&config);
	read_fasta("Homo_sapiens.NCBI35.dec.pep.fa",&arr,&t_size,&config); 
//	homeomorphic_levels(&config,arr,t_size,1000.0,1001,"hp_res.txt");

	config->set_tolerance(0.0075);
//	full_exp(config,0,arr,t_size,"hp_res_dis_a0.txt");
//	full_exp(config,1,arr,t_size,"hp_res_dis_a1.txt");
//	full_exp(config,2,arr,t_size,"hp_res_dis_a2.txt");
//	full_exp(config,3,arr,t_size,"hp_res_dis_a3.txt");
//	full_exp(config,4,arr,t_size,"hp_res_dis_a4.txt");
	

//	homeomorphic_exp3(&config,100);
//	exit(0);
//	config.print_supported_PTMs();

//	config.print_session_aas();
	string p;
//	ifstream fs("D:\\msms4\\PepNovo\\test\\C25A19_IP_01.mgf",ios::in);
//	ifstream fs("D:\\msms4\\PepNovo\\test\\m280.mgf",ios::in);
//	fm.init_from_mgf(&config,"D:\\msms4\\PepNovo\\test\\m280.mgf");
//	fm.init_from_mgf(&config,"D:\\msms4\\PepNovo\\mgf_2600.2.mgf");
//	fm.init_from_list_file(&config,"D:\\msms4\\lists\\unique_good2.txt");
	
//	fm.init_from_list_file(&config,"D:\\msms4\\lists\\short2.txt");
//	fm.init_from_list_file(&config,"D:\\Data2\\ikkb_unmod_list.txt");

//	rand_seed(1111);
	rand_seed(1426347);


//	model.read_model("ESI_RANKS");
//	model.get_config()->add_selected_PTMs("C+57");

//	SpectrumScore sqs;
//	sqs.learn_score_params(&model,"D:\\msms4\\lists\\sqs2_short.txt",
//		"D:\\msms4\\lists\\sqs_neg1.txt");	
//	exit(0);


//	model.set_model_name(string("ESI2"));


//	random_peak_match_exp(model.get_config(),fm,800,1200,10000000);
//	model.print_joint_scores();

//	me_reg_exp();
//	me_exp();
//	exit(0);





//	fm.init_from_mgf(model.get_config(),"c:\\work\\msms4\\PepNovo\\test\\m280.mgf");
//	fm.init_from_mgf(model.get_config(),"c:\\work\\msms4\\PepNovo\\hpep.mgf");
	fm.init_from_list_file(config,"c:\\Work\\msms4\\lists\\efast2.txt");
//	fm.init_from_list_file(&config,"D:\\msms4\\lists\\charge2.txt");
//m.init_from_list_file(model.get_config(),"C:\\Work\\msms4\\PepNovo\\lll.txt");
//m.init_from_list_file(model.get_config(),"D:\\msms4\\lists\\l1.txt");
	


*/


	return 0;
}