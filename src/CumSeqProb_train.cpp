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
	
	int charge=-1, size_idx=-1, num_rounds = 10000, max_num_training_samples = 200000;
	float max_weight_ratio = 1000.0;
	int tag_length=0;
	int num_test_cases=0;
	char name[64];
	char mgf_list[512],report_dir[512], test_file_path[512];

	char *test_path   = NULL;
	char *report_path = NULL;

	bool got_name=false;
	bool got_mgf_list=false;
	bool got_report_dir=false;
	bool got_stop_signal_file=false;
	bool got_test = false;

	
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
		else if (!strcmp(argv[i],"-model")) 
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
		else if (! strcmp(argv[i],"-tag_length"))
		{
			if (sscanf(argv[++i],"%d",&tag_length) != 1)
				error("tag length");
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
		else
		{
			cout << "Error: don't recognize option " << argv[i] << endl;
			exit(1);
		}
		i++;
	}

	if (tag_length<0 || tag_length>9)
	{
		cout << "Error: bad tag length: " << tag_length << endl;
		exit(1);
	}



	// 
	if (got_test)
	{
		model.read_model(name);
		model.get_config()->apply_selected_PTMs("C+57");
		model.read_rank_models(name);
		exit(0);
	}


	if (! got_report_dir)
		error("must supply report_dir!");

	if (charge<0 || size_idx<0)
		error("must suupply charge and size_idx!");


	
	model.read_model(name);
	model.get_config()->apply_selected_PTMs("C+57:M+16:Q-17");
	model.read_rank_models(name);
	Config *config= model.get_config();

	cout << "STARTING TRAINING..." << endl;
	cout << "tag length: " << tag_length << " charge " << charge << " size_idx " << size_idx << endl;



	CumulativeSeqProbModel csp;

	FileManager fm;

	fm.init_from_list_file(config,mgf_list);
	csp.train_seq_prob_models(fm,&model,tag_length,charge,size_idx);

	

	return 0;
}


