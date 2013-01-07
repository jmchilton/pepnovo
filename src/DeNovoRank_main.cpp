#include "AdvancedScoreModel.h"
#include "PeakRankModel.h"
#include "PeptideRankScorer.h"
#include "CumulativeSeqProb.h"
#include "PMCSQS.h"
#include "auxfun.h"
#include "includes.h"


/*
 
  -model_type 1 -charge 2 -size_idx 1 -mgf_list C:\Work\msms5\PepNovoHQ\hek2.txt -report_dir C:\Work\msms5\DnvScore\reports -max_num_training_samples 100000 -max_weight_ratio 5 -num_rounds 4000 -name DNV 
  -name CID_IT_TRYP -PTMs c+57 -test C:\Work\msms5\PepNovoHQ\m280.mgf 2000 -rerank_depth 1000
  
*/
  
int main(int argc, char **argv) 
{
	AdvancedScoreModel model;
	PeptideRankScorer drs;
 
	int charge=-1, size_idx=-1, num_rounds = 10000, max_num_training_samples = 200000;
	float max_weight_ratio = 1000.0;
	int model_type=0;
	int length_limit = 30;
	int min_length = 6;
	int max_length = 30;
	int num_test_cases=0;
	int rerank_depth=0;
	int new_mgf_size=0;
	int obs_rank =0;
	char name[64];
	char mgf_list[512],report_dir[512], test_file_path[512],new_mgf_name[512],peak_benchmark_file[256];

	char *test_path   = NULL;
	char *report_path = NULL;

	bool got_name=false;
	bool got_mgf_list=false;
	bool got_report_dir=false;
	bool got_stop_signal_file=false;
	bool got_test = false;
	bool got_make_mgf = false;
	bool got_peak_benchmark = false;
	bool got_peak_hist = false;
	bool got_examples = false;
	
	rand_seed(112233);   

	int i;
	i=1;

	while (i<argc)
	{
		if (! strcmp(argv[i],"-charge"))
		{
			if (sscanf(argv[++i],"%d",&charge) != 1)
				error("bad charge");
		}
		else if (!strcmp(argv[i],"-size_idx"))
		{
			if (sscanf(argv[++i],"%d",&size_idx) != 1)
				error("size idx");
		}
		else if (!strcmp(argv[i],"-max_num_training_samples"))
		{
			if (sscanf(argv[++i],"%d",&max_num_training_samples) != 1)
				error("max_num_training_samples");
		}
		else if (!strcmp(argv[i],"-num_rounds"))
		{
			if (sscanf(argv[++i],"%d",&num_rounds) != 1)
				error("num_rounds");
		}
		else if (!strcmp(argv[i],"-name")) 
		{
			if (sscanf(argv[++i],"%s",name) != 1)
				error("name");

			got_name = true;
		}
		else if (!strcmp(argv[i],"-mgf_list"))
		{
			if (sscanf(argv[++i],"%s",mgf_list) != 1)
				error("mgf_list");

			got_mgf_list = true; 
		}
		else if (! strcmp(argv[i],"-make_mgf"))
		{
			if (sscanf(argv[++i],"%s",new_mgf_name) != 1)
				error("make_mgf1");

			if (sscanf(argv[++i],"%d",&new_mgf_size) != 1)
				error("make_mgf2");

			got_make_mgf = true;
		}
		else if (!strcmp(argv[i],"-report_dir"))
		{
			if (sscanf(argv[++i],"%s",report_dir) != 1)
				error("report_dir");

			report_path = report_dir;
			got_report_dir = true;
		}
		else if (!strcmp(argv[i],"-max_weight_ratio"))
		{
			if (sscanf(argv[++i],"%f",&max_weight_ratio) != 1)
				error("max weight ratio");

		}
		else if (!strcmp(argv[i],"-model_type"))
		{
			if (sscanf(argv[++i],"%d",&model_type) != 1)
				error("model type");
		}
		else if (! strcmp(argv[i],"-length_limit"))
		{
			if (sscanf(argv[++i],"%d",&length_limit) != 1)
				error("length limit");
		}
		else if (! strcmp(argv[i],"-min_length"))
		{
			if (sscanf(argv[++i],"%d",&min_length) != 1)
				error("min_length");
		}
		else if (! strcmp(argv[i],"-max_length"))
		{
			if (sscanf(argv[++i],"%d",&max_length) != 1)
				error("max_length");
		}
		else if (!strcmp(argv[i],"-rerank_depth"))
		{
			if (sscanf(argv[++i],"%d",&rerank_depth) != 1)
				error("rerank_depth");
		}
		else if (! strcmp(argv[i],"-test"))
		{
		
			if (sscanf(argv[++i],"%s",test_file_path) != 1)
				error("no test file");

			test_path = test_file_path;

			if (sscanf(argv[++i],"%d",&num_test_cases) != 1)
				error("no num cases");


			got_test=true;;
		}
		else if (! strcmp(argv[i],"-examples"))
		{
			if (sscanf(argv[++i],"%s",test_file_path) != 1)
				error("no test file");

			test_path = test_file_path;
			got_examples=true;
		}
		else if (! strcmp(argv[i],"-peak_benchmark_file"))
		{
			if (sscanf(argv[++i],"%s",peak_benchmark_file) != 1)
				error("benchmark_file");

			got_peak_benchmark = true;
		}

		else if (! strcmp(argv[i],"-peak_hist"))
		{
			if (sscanf(argv[++i],"%s",peak_benchmark_file) != 1)
				error("peak_hist");

			if (sscanf(argv[++i],"%d",&obs_rank) != 1)
				error("peak_hist obs rank");

			got_peak_hist = true;
		}
		else
		{
			cout << "Error: don't recognize option " << argv[i] << endl;
			exit(1);
		}
		i++;
	}


	if (got_make_mgf)
	{
		model.read_model(name);
		model.get_config()->apply_selected_PTMs("C+57:M+16:Q-17");
 
		create_bench_mgf(model.get_config(),mgf_list,new_mgf_name,new_mgf_size,false);
		exit(0);
	}

	// 
	if (got_test)
	{
		model.read_model(name);
		model.get_config()->apply_selected_PTMs("C+57");
		model.read_rank_models(name);
		model.get_config()->set_use_spectrum_charge(1);

		if (! model.get_rank_model_ptr(1))
		{
			cout << "Didn't get model pointer!" << endl;
			exit(1);
		}

		if (! test_path)
			error("must supply path to test file!");
		if (num_test_cases<1)
			error("must supply num cases!");
		if (rerank_depth<1 || rerank_depth>100000)
			error("must supply rerank depth!");

	//	model.test_pmc(test_path,2);
	//	exit(0); 
	
		benchmark_ranking_on_full_denovo(&model,
									  test_path, 
									  num_test_cases,
									  rerank_depth,
									  report_path,
									  min_length,
									  max_length);



		exit(0);
	}

	if (got_examples)
	{
		model.read_model(name);
		model.get_config()->apply_selected_PTMs("C+57");
		model.read_rank_models(name);
		model.get_config()->set_use_spectrum_charge(1);

		if (! model.get_rank_model_ptr(1))
		{
			cout << "Didn't get model pointer!" << endl;
			exit(1);
		}

		if (! test_path)
			error("must supply path to test file!");

		make_ranking_examples(&model, test_path);

		exit(0);
	}

	if (got_peak_benchmark)
	{
		model.read_model(name);
		model.get_config()->apply_selected_PTMs("C+57:M+16:Q-17");
		model.read_rank_models(name);
		run_peak_benchmark(&model, peak_benchmark_file);
		exit(0);
	}


	if (got_peak_hist)
	{
		model.read_model(name);
		model.get_config()->apply_selected_PTMs("C+57:M+16:Q-17");
		model.read_rank_models(name);
		make_peak_hist_for_obs_rank(&model, peak_benchmark_file, obs_rank);
		exit(0);
	}

	
	if (! got_report_dir)
		error("must supply report_dir!");

	if (charge<0 || size_idx<0)
		error("must suupply charge and size_idx!");

	if (model_type<0 || model_type>3)
		error("model type must be 0 (db), 1 (de novo), 2 (de novo for complete sequences), or 3 (tags) only!");


	cout << "STARTING TRAINING..." << endl;
	cout << "Model type: " << model_type << endl; 

// -model_type 1 -charge 2 -size_idx 1 -mgf_list C:\Work\msms5\DnvScore\dicty2_all.txt -report_dir C:\Work\msms5\DnvScore\reports -max_num_training_samples 100000 -max_weight_ratio 5 -num_rounds 2000 -name DNV

/*	const string db_hits = "C:\\Work\\msms5\\DnvScore\\new_db_hits";
	const string sequences = "C:\\Work\\msms5\\DnvScore\\seq_freqs\\sequences_mqscore";
	const string full_dnv = "C:\\Work\\msms5\\DnvScore\\new_dnv_hits"; */
	

	const string db_hits =   "/scratch/amfrank/DnvRank/db_small_dbh";
//	const string db_hits = 	"/scratch/amfrank/DnvRank/all_db_hits";
	const string sequences = "/scratch/amfrank/DnvRank/seq_freqs/sequences_mqscore";
	const string full_dnv =  "/scratch/amfrank/DnvRank/dnv_dbh"; 

	if (model_type == 0)
	{
		drs.set_type(0); 
		drs.train_partition_model_for_complete_sequences(
								db_hits,
								sequences, 
								full_dnv,
							//	"C:\\Work\\msms5\\DnvScore\\short2_train_mgf_list.txt",
							//	"C:\\Work\\msms5\\DnvScore\\comp2_train_mgf_list.txt",
								 mgf_list,
								 report_dir,
								 name,
								 charge,
								 size_idx,
								 num_rounds,
								 max_weight_ratio,
								 max_num_training_samples); 
	}
	else if (model_type == 1)
	{
		drs.set_type(1); 
		drs.train_partial_denovo_partition_model(
								 mgf_list,
								 report_dir,
								 name,
								 charge,
								 size_idx,
								 num_rounds,
								 max_weight_ratio,
								 max_num_training_samples,
								 length_limit,
								 NULL); 
	}
	else if (model_type == 2)
	{
		drs.set_type(2); 
		drs.train_partition_model_for_complete_sequences(
								db_hits,
								sequences, 
								full_dnv,
								 mgf_list,
								 report_dir,
								 name,
								 charge,
								 size_idx,
								 num_rounds,
								 max_weight_ratio,
								 max_num_training_samples,
								 0.2,
								 0.8,
								 0,
								 NULL,
								 rerank_depth); 
	}
	if (model_type == 3)
	{
		drs.set_model_length(min_length);
		drs.set_type(3); 
		drs.train_partial_denovo_partition_model(
								 mgf_list,
								 report_dir,
								 name,
								 charge,
								 size_idx,
								 num_rounds,
								 max_weight_ratio,
								 max_num_training_samples,
								 min_length,
								 NULL);
	}


	drs.write_denovo_rank_scorer_model(name);   


//	rank_model.train_all_partition_models(sample_prefix_path, charge,size_idx, mobility, frag,
//		report_dir, num_rounds, test_set_ptr, test_length, stop_signal_file_ptr, max_weight_ratio);

//	rank_model.write_peak_rank_model(name,model_out_dir);
	

	return 0;
}


