#include "PeakRankModel.h"
#include "AllScoreModels.h"

void train_all()
{
	PeakRankModel rm,rm2;
	AllScoreModels model;
	vector<FileManager> fms;

	model.read_model("LTQ_LOW_TRYP"); 
	Config *config= model.get_config();

	char report_dir[]={"C:\\Work\\msms5\\NewScore\\report\\"};
	char singal_dir[]={"C:\\Work\\msms5\\NewScore\\signal"};

	config->apply_selected_PTMs("C+57:M+16:Q-17");

	int frag_set_type = 4;

	rm.init_peak_rank_model_with_defaults(config,"COMB3dnv",frag_set_type);

//	rm.read_peak_rank_model(config,"RM/RM");
//	rm.write_peak_rank_model("RM/RM2");

//	exit(0); 

	rm.train_all_combined_partition_models(frag_set_type,"C:\\Work\\msms5\\NewScore\\sams\\tr",2,2,2,4,
			report_dir, 200, "C:\\Work\\msms5\\NewScore\\sams\\ts_2_2_2.txt",-1,NULL,5.0);


//	rm.train_all_combined_partition_models(frag_set_type,"C:\\Work\\msms5\\NewScore\\sams\\tr",2,1,1,4,
//			report_dir, 2000, "C:\\Work\\msms5\\NewScore\\sams\\ts_2_1_1.txt",-1,NULL,5.0);

//	rm.train_all_combined_partition_models(frag_set_type,"C:\\Work\\msms5\\NewScore\\sams\\trz",3,2,3,3,
//			report_dir, 200, "C:\\Work\\msms5\\NewScore\\sams\\ts_2_1_1.txt",-1,NULL,5.0);

//	rm.train_all_partition_models(frag_set_type,"C:\\Work\\msms5\\NewScore\\sams\\tr",2,1,1,0,
//			report_dir, 200, "C:\\Work\\msms5\\NewScore\\sams\\ts_2_1_1.txt",-1,NULL,5.0);

//	rm.train_all_combined_partition_models(frag_set_type,"C:\\Work\\msms5\\NewScore\\sams\\trz",3,2,3,4,
//					report_dir, 50, "C:\\Work\\msms5\\NewScore\\sams\\tsz_3_2_3.txt",-1,NULL);


//	rm.write_peak_rank_model("RM/RM6");

//	rm2.read_peak_rank_model(config,"RM/RM");
	
}