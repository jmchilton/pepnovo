#include "AdvancedScoreModel.h"
#include "PeakRankModel.h"
#include "PMCSQS.h"
#include "auxfun.h"
#include "includes.h"




/*
 -charge 3 -size_idx 2 -mobility 3 -frag 8 -num_rounds 100000 -name F5 -sample_prefix_path C:\\Work\\msms5\\NewScore\\sams/tr -report_dir C:\\Work\\msms5\\NewScore\\report -model_out_dir C:\\Work\\msms5\\NewScore\\model_out -stop_signal_file C:\\Work\\msms5\\NewScore\\signal\\STOP_F5 -max_weight_ratio 5 -test_set C:\\Work\\msms5\\NewScore\\sams\\ts_3_2_3.txt -test_lengh -1

-feature_type 5 -charge 2 -size_idx 1 -mobility 1 -num_frags 4 -num_rounds 1000 -name F5 -sample_prefix_path C:\\Work\\msms5\\NewScore\\sams/tr -report_dir C:\\Work\\msms5\\NewScore\\report -model_out_dir C:\\Work\\msms5\\NewScore\\model_out -stop_signal_file C:\\Work\\msms5\\NewScore\\signal\\STOP_F5 -max_weight_ratio 5 -test_set C:\\Work\\msms5\\NewScore\\sams\\ts_2_1_1.txt -test_lengh -1
*/


int main(int argc, char **argv) 
{
	int charge=-1,size_idx=-1,mobility=-1,frag=-1,test_length=-1, num_rounds = 10000;
	int num_frags=-1;
	float max_weight_ratio = 1000.0;
	char name[64],score_model[64];
	char sample_prefix_path[256],report_dir[256],stop_signal_file[256],
		 model_out_dir[256],test_set[256];

	bool got_name=false;
	bool got_print_idxs=false;
	bool got_sample_prefix_path=false;
	bool got_report_dir=false;
	bool got_stop_signal_file=false;
	bool got_model_out_dir=false;
	bool got_test_set=false;
	bool got_score_model=false;
	bool got_benchmark=false;

	int feature_type = 5;

	rand_seed(112233); 

	int i;
	i=1;

	while (i<argc)
	{
		if (!strcmp(argv[i],"-print_idxs"))
		{
			got_print_idxs=true;
			break;
		}
		else if (! strcmp(argv[i],"-charge"))
		{
			if (sscanf(argv[++i],"%d",&charge) != 1)
				error("bad charge");
		}
		else if (!strcmp(argv[i],"-size_idx"))
		{
			if (sscanf(argv[++i],"%d",&size_idx) != 1)
				error("size idx");
		}
		else if (!strcmp(argv[i],"-mobility"))
		{
			if (sscanf(argv[++i],"%d",&mobility) != 1)
				error("mobility");
		}
		else if (!strcmp(argv[i],"-frag"))
		{
			if (sscanf(argv[++i],"%d",&frag) != 1)
				error("frag type idx");
		}
		else if (!strcmp(argv[i],"-num_frags"))
		{
			if (sscanf(argv[++i],"%d",&num_frags) != 1)
				error("num frags");
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
		else if (!strcmp(argv[i],"-sample_prefix_path"))
		{
			if (sscanf(argv[++i],"%s",sample_prefix_path) != 1)
				error("sample_prefix_path");

			got_sample_prefix_path = true;
		}
		else if (!strcmp(argv[i],"-report_dir"))
		{
			if (sscanf(argv[++i],"%s",report_dir) != 1)
				error("report_dir");

			got_report_dir = true;
		}
		else if (!strcmp(argv[i],"-stop_signal_file"))
		{
			if (sscanf(argv[++i],"%s",stop_signal_file) != 1)
				error("stop_signal_file");

			got_stop_signal_file = true;
		}
		else if (!strcmp(argv[i],"-model_out_dir"))
		{
			if (sscanf(argv[++i],"%s",model_out_dir) != 1)
				error("model_out_dir");

			got_model_out_dir = true;
		}
		else if (!strcmp(argv[i],"-test_set"))
		{
			if (sscanf(argv[++i],"%s",test_set) != 1)
				error("test_set");

			got_test_set = true;
		}
		else if (!strcmp(argv[i],"-test_length"))
		{
			if (sscanf(argv[++i],"%d",&test_length) != 1)
				error("test length");
		}
		else if (!strcmp(argv[i],"-score_model"))
		{
			if (sscanf(argv[++i],"%s",score_model) != 1)
				error("score model");

			got_score_model = true;
		}
		else if (!strcmp(argv[i],"-max_weight_ratio"))
		{
			if (sscanf(argv[++i],"%f",&max_weight_ratio) != 1)
				error("max weight ratio");

		}
		else if (!strcmp(argv[i],"-feature_type"))
		{
			if (sscanf(argv[++i],"%d",&feature_type) != 1)
				error("feature_type");

		}

		i++;
	}

	AdvancedScoreModel model;
	PeakRankModel    rank_model;

	if (! got_score_model)
		strcpy(score_model,"CID_IT_TRYP");


	model.read_model(score_model);
	model.get_config()->apply_selected_PTMs("M+16:C+57:Q-17");
	Config *config = model.get_config();




	if (! got_name && ! got_print_idxs)
	{
		cout << "Error: must supply rank model name! (-name XXX)" << endl;
		exit(1);
	}

	rank_model.init_peak_rank_model_with_defaults(config,name,feature_type);

	if (got_print_idxs)
	{
		rank_model.list_all_model_idxs();
		return 0;
	}
	
	if (! got_report_dir)
		error("must supply report_dir!");

	if (! got_model_out_dir)
		error("must supply model out dir!");

	char *test_set_ptr = NULL;
	if (got_test_set)
		test_set_ptr = test_set;

	char *stop_signal_file_ptr=NULL;
	if (got_stop_signal_file)
		stop_signal_file_ptr = stop_signal_file;

	cout << "STARTING TRAINING..." << endl;

	if (num_frags>0)
	{
		rank_model.train_all_combined_partition_models(feature_type, sample_prefix_path, charge, size_idx, mobility, num_frags,
			report_dir, num_rounds, test_set_ptr, test_length, stop_signal_file_ptr, max_weight_ratio);
	}
	else
	{
		rank_model.train_all_partition_models(feature_type, sample_prefix_path, charge,size_idx, mobility, frag,
			report_dir, num_rounds, test_set_ptr, test_length, stop_signal_file_ptr, max_weight_ratio);
	}

	rank_model.write_peak_rank_model(name,model_out_dir);
	

	return 0;
}


