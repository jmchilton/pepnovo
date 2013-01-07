#include "AllScoreModels.h"
#include "PeptideRankScorer.h"
#include "CumulativeSeqProb.h"
#include "FragmentSelection.h"






// reads a model and all relevant files
// the model files are assumed to be in the resource_dir
// all this model's files are assumed to have a name <model_name>_XXXXX.txt
// the main model file is <model_name>.txt
void AllScoreModels::read_model(const char* name, bool silent_ind)
{
	char file[256];

	model_name = name;

	if (config_.get_resource_dir().length()<2)
	{
		config_.set_resource_dir("Models");
	}

	config_.set_model_name(string(name));

	strcpy(file,config_.get_resource_dir().c_str());
	strcat(file,"/");
	strcat(file,name); 
	strcat(file,".txt");   

	fstream fs(file,ios::in);
	if (! fs.good() )  
	{
		cout << "Error: couldn't open model file: " << file << endl;
		exit(1);
	}

	while (! fs.eof())
	{
		char buff[1024];
		fs.getline(buff,1024);
		if (fs.gcount()<4)
			continue;

		char arg[128];
		if (sscanf(buff,"#CONFIG_FILE %s",arg) == 1)
		{
			config_.read_config(arg);
			config_.set_model_name(string(model_name));
			continue;
		}

		if (! strncmp("#CONF",buff,5))
		{
			string path = config_.get_resource_dir() + "/" + string(buff);
			config_.parse_config_parameter((char *)path.c_str());
			continue;
		}

		if (sscanf(buff,"#BREAK_SCORE_MODEL %s",arg) ==1)
		{
			prmNodeScoreModel_.read_score_model(&config_, arg,silent_ind);
			continue;
		}

		if (sscanf(buff,"#EDGE_MODEL %s",arg) ==1)
		{
			edgeModel_.read_edge_models(&config_,arg,silent_ind);
			continue;
		}

		if (sscanf(buff,"#SQS_MODEL %s",arg) == 1)
		{
			pmcsqs.read_sqs_models(&config_,arg);
			continue;
		}
		
		if (sscanf(buff,"#PMCR_MODEL %s",arg) == 1)
		{
			pmcsqs.read_pmc_rank_models(&config_,arg);
			continue;
		}

		if (sscanf(buff,"#COMP_ASSIGNER %s",arg) == 1)
		{
			compAssigner_.read_and_init_from_tables(&config_,arg);
			continue;
		}

		if (sscanf(buff,"#AAP_MODEL %s",arg) == 1)
		{
			amino_acid_probs.read_amino_acid_prob_models(&config_,arg);
			continue;
		}

	}

	// try default comp assigned
	if (! compAssigner_.get_ind_was_initialized() &&
		! compAssigner_.read_and_init_from_tables(&config_,"LTQ_COMP/IT_TRYP"))
	{
		cout << endl << "Warning: no peptide composition assigner was found (e.g., LTQ_COMP/IT_TRYP)!" << endl;
		cout <<         "This may cause problems when trying to run PepNovo!" << endl << endl;
	}

	// check if some of the defaults need to be changed
	if (config_.get_max_edge_length() != 2)
		config_.calc_aa_combo_masses();

}



// writes a model and all relevant files
// the model files are assumed to be in the resource_dir
// all this model's files are assumed to have a name <model_name>_XXXXX.txt
// the main model file is <model_name>.txt
void AllScoreModels::write_model()
{
	string model_file;

	model_file = config_.get_resource_dir() + "/" + model_name + ".txt";

	fstream os(model_file.c_str(),ios::out);
	if ( ! os.good())
	{
		cout << "Error writing model to " << model_file << endl;
		exit(1);
	}


	string config_file = config_.get_resource_dir() + "/" + model_name + "_config.txt";
	config_.set_config_file(config_file);
	config_.set_model_name(model_name);
	os << "#CONFIG_FILE " << model_name + "_config.txt" << endl;
	config_.write_config();


	if (pmcsqs.getIndInitializedPmcr())
	{
		os << "#PMCR_MODEL " << model_name + "_PMCR.txt" << endl;
		string path = config_.get_resource_dir() + "/" + model_name + "_PMCR.txt";
		pmcsqs.write_pmc_rank_models(path.c_str());
	}

	if (pmcsqs.getIndInitializedSqs())
	{
		os << "#SQS_MODEL " << model_name + "_SQS.txt" << endl;
		string path = config_.get_resource_dir() + "/" + model_name + "_SQS.txt";
		pmcsqs.write_sqs_models(path.c_str());
	}

	if (compAssigner_.get_ind_was_initialized())
	{
		os << "#COMP_ASSIGNER " << compAssigner_.get_model_name() << endl;
	}
	else
	{
		cout << "Warning: no peptide composition assigner was written" << endl;
		cout << "You might need to add the line: \"#COMP_ASSIGNER LTQ_COMP/IT_TRYP\" to the main model file " <<
			model_name << ".txt to get everything to work!" << endl;
	}

	if (prmNodeScoreModel_.get_ind_was_initialized())
	{
		os << "#BREAK_SCORE_MODEL " << model_name << endl;
		prmNodeScoreModel_.write_score_model(model_name.c_str());
	}

	
	if (edgeModel_.get_ind_was_initialized())
	{
		os << "#EDGE_MODEL " << model_name << endl;
		edgeModel_.write_edge_models(model_name.c_str());
	}

	if (amino_acid_probs.get_ind_initialized())
	{
		os << "#AAP_MODEL " << model_name << endl;
		string path = config_.get_resource_dir() + "/" + model_name + "_AAP.txt";
		amino_acid_probs.write_amino_acid_prob_models(path.c_str());
	}
}








void AllScoreModels::read_rank_models(const char *name, bool silent_ind, bool use_large_db_model)
{
	const string suffixes[]={"DB","DNVPART","DNVCOMP"};
	
	peptideRankModels_[0]=NULL;
	peptideRankModels_[1]=NULL;
	peptideRankModels_[2]=NULL;

	int i;
	for (i=0; i<10; i++)
		tagRankModels_[i]=NULL;

	const int num_models = sizeof(peptideRankModels_)/sizeof(PeptideRankScorer*);

	for (i=0; i<num_models; i++)
	{
		const string suffix = suffixes[i];
		PeptideRankScorer *rank_model = peptideRankModels_[i];
	
		string rank_name = string(name) + "_" + suffix;
		string path = config_.get_resource_dir() + "/" + rank_name + "/" + suffix + "_rank_model.txt";

		ifstream fs(path.c_str());
		if (! fs.is_open() || ! fs.good())
		{
			if (! silent_ind)
				cout << "No " << path << endl;
			continue;
		}

		fs.close();
		peptideRankModels_[i] = new PeptideRankScorer;
		peptideRankModels_[i]->set_type(i);
		peptideRankModels_[i]->set_model(this);

		// add hack for DB models that will tell it to change name for large db
		peptideRankModels_[i]->read_denovo_rank_scorer_model(path.c_str(),suffix,silent_ind, (i==0 && use_large_db_model));
	
		if (! silent_ind)
			cout << "Read " << path << endl;

	}

	// read tag models
	for (i=3; i<10; i++)
	{
		PeptideRankScorer *rank_model = tagRankModels_[i];
	
		
		ostringstream oss;
		oss << "TAG" << i;
		string suffix = oss.str();
		string rank_name = string(name) + "_" + suffix;
		string path = config_.get_resource_dir() + "/" + rank_name + "/" + suffix + "_rank_model.txt";

		ifstream fs(path.c_str());
		if (! fs.is_open() || ! fs.good())
		{
			if (! silent_ind)
				cout << "No " << path << endl;
			continue;
		}

		fs.close();
		tagRankModels_[i] = new PeptideRankScorer;
		tagRankModels_[i]->set_type(3);
		tagRankModels_[i]->set_model(this);
		tagRankModels_[i]->read_denovo_rank_scorer_model(path.c_str(),suffix,silent_ind);
		 
		if (! silent_ind)
			cout << "Read " << path << endl;

	}

//	cout << endl;
}


void AllScoreModels::read_cum_seq_prob_models(const char *name, bool silent_ind)
{
	const int max_tag_model = 9;
	string dir_name = string(name)+"_CSP";
	int i;
	for (i=0; i<=max_tag_model; i++)
	{
		CumulativeSeqProbModel *csp_model = (CumulativeSeqProbModel *)cumulativeSequenceModels_[i];
		char name_buff[64];

		sprintf(name_buff,"%s_CSP_%d.txt",name,i);
		
		string path = config_.get_resource_dir() + "/" +dir_name + "/" + string(name_buff);

		ifstream ifs(path.c_str());
		if (! ifs.is_open() || ! ifs.good())
		{
			if (! silent_ind)
				cout << "No " << path << endl;

			if (ifs.is_open())
				ifs.close();
			continue;
		}

		
		csp_model = new CumulativeSeqProbModel;
		csp_model->read_model(&config_,ifs);
		ifs.close();

		cumulativeSequenceModels_[i] = csp_model;
	
		if (! silent_ind)
			cout << "Read " << path << endl;

	}
}






	


void AllScoreModels::score_graph_edges(PrmGraph& prm) const
{
	edgeModel_.score_graph_edges(prm);
}



int AllScoreModels::get_max_score_model_charge() const
{
	return prmNodeScoreModel_.get_max_score_model_charge();
}

























