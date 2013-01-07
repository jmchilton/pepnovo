/******************************************************************
Parser for the website's interface.
Controls everything that involves running de novo from the website
*******************************************************************/


#include "AllScoreModels.h"
#include "PeptideRankScorer.h" 
#include "FileManagement.h"
#include "DeNovoDp.h"
#include "DeNovoSolutions.h"
#include "auxfun.h"
#include "includes.h"
#include "MSBlast.h"


void print_error(char *line) 
{
	cout << "Error: coudln't parse input file line:\n" << line << endl;
	exit(1);
}



int main(int argc, char **argv) 
{ 

	char *input_file = argv[1];

	FILE *inp_file = fopen(input_file,"r");
	if (! inp_file)
	{
		cout << "Error: couldn't open input file: " << input_file << endl;
		exit(1);
	}

	vector<string> files;
	vector<string> ptm_lines;
	
	char model_name[64]; 
	char model_dir[128];
	char output_file[256];
	char out_dir[256];
	int  protease_digest=TRYPSIN_DIGEST;
	bool got_out_dir = false;

	mass_t tolerance=0.5;
	mass_t pm_tolerance=2.5;

	int num_solutions=10;
	int tag_length = 0;
	int file_start_idx = 0;

	float filter_prob = 0.0;
	bool got_filter_spectra = false;
	bool ind_pmcsqs = false;
	bool got_output_file = false;
	bool use_spectrum_charge = false;
	bool use_spectrum_mz = false;


	strcpy(model_dir,"Models");
	strcpy(model_name,"CID_IT_TRYP");

	char buff[256];
	while (fgets(buff,256,inp_file))
	{
		if (buff[0]=='#')
			continue;

		if (! strncmp(buff,"spectra,",8))
		{
			char file_path[256];
			if (sscanf(buff,"spectra,%s",file_path) != 1)
				print_error(buff);

			files.push_back(file_path);
			continue;
		}

		if (! strncmp(buff,"output,",7))
		{
			if (sscanf(buff,"output,%s",output_file) != 1)
				print_error(buff);

			got_output_file = true;

			continue;
		}

		if (! strncmp(buff,"protease,",9))
		{
			char prot_name[64];
			if (sscanf(buff,"protease,%s",prot_name) != 1)
				print_error(buff);

			if (! strcmp(prot_name,"Trypsin"))
			{
				protease_digest=TRYPSIN_DIGEST;
			}
			else
				protease_digest = NON_SPECIFIC_DIGEST;

			continue;
		}

		if (! strncmp(buff,"modelDir,",7))
		{
			if (sscanf(buff,"modelDir,%s",model_dir) != 1)
				print_error(buff);

			continue;
		}

		if (! strncmp(buff,"instrument,",11))
		{
			char instrument_name[64];
			if ( sscanf(buff,"instrument,%s",instrument_name) != 1)
				print_error(buff);

			if (! strcmp(instrument_name,"ESI-ION-TRAP"))
			{
				strcpy(model_name,"CID_IT_TRYP");
			}
			else if (! strcmp(instrument_name,"FT-HYBRID"))
			{
				strcpy(model_name,"CID_IT_TRYP");
			}
			else
			{
				cout << "Warning: model for " << instrument_name << " not supported. Using " << model_name << endl;
			}

			continue;
		}

		if (! strncmp(buff,"mod,",4))
		{
			ptm_lines.push_back(buff);
			continue;
		}

		if (! strncmp(buff,"useSpectrumCharge",16))
		{
			use_spectrum_charge = true;
		}

		if (! strncmp(buff,"useSpectrumMZ",13))
		{
			use_spectrum_mz = true;
		}

		if (! strncmp(buff,"PMTolerance,",10))
		{
			float val;
			if (sscanf(buff,"PMTolerance,%f",&val) != 1)
				print_error(buff);
			pm_tolerance = val;
			continue;
		}

		if (! strncmp(buff,"IonTolerance,",10))
		{
			float val;
			if (sscanf(buff,"IonTolerance,%f",&val) != 1)
				print_error(buff);
			tolerance = val;

			continue;
		}

		if (! strncmp(buff,"numSolutions,",11))
		{
			int val;
			if (sscanf(buff,"numSolutions,%d",&val) != 1)
				print_error(buff);
			num_solutions = val;
			continue;
		}

		if (! strncmp(buff,"tagLength,",8))
		{
			int val;
			if (sscanf(buff,"tagLength,%d",&val) != 1)
				print_error(buff);
			tag_length = val;
			continue;
		}


		if (! strncmp(buff,"filterSpectra,",12))
		{
			float prob;
			if (sscanf(buff,"filterSpectra,%f",&prob) != 1)
				print_error(buff);
			filter_prob = prob;
			got_filter_spectra = true;
		}

		if (! strncmp(buff,"mgfOutDir,",9))
		{

			if (sscanf(buff,"mgfOutDir,%s",out_dir) != 1)
				print_error(buff);
			got_out_dir=true;
			continue;
		}


		if (! strncmp(buff,"pmcsqs",6))
		{
			ind_pmcsqs = true;
			continue;
		}


		if (! strncmp(buff,"fileStartIdx,",12) )
		{
			if (sscanf(buff,"fileStartIdx,%d",&file_start_idx) != 1)
				print_error(buff);
			continue;
		}
	}

	time_t start_time,last_time;
	start_time = time(NULL);
	last_time = start_time;


	// Read model and set model parameters
	AllScoreModels model;
	PMCSQS_Scorer *pmcsqs = (PMCSQS_Scorer *)model.get_pmcsqs_ptr();
	Config *config = model.get_config();
	
	config->set_resource_dir(string(model_dir));

	cout << "Resource dir: " << config->get_resource_dir() << endl;

	model.read_model(model_name,true);
	config->apply_site_input_PTMs(ptm_lines);
	model.read_rank_models(model_name,true);
	
	config->set_tolerance(tolerance);
	config->setPrecursorMassTolerance(pm_tolerance);
	config->set_digest_type(protease_digest);

	if (use_spectrum_charge)
		config->set_use_spectrum_charge(1);

	if (use_spectrum_mz)
		config->set_use_spectrum_mz(1);
	

	ofstream  out_file_stream;
	if (got_output_file)
	{
		out_file_stream.open(output_file);
		if (! out_file_stream.is_open() || ! out_file_stream.good())
		{
			cout << "Error: couldn't open output file for writing: " << output_file << endl;
			exit(1);
		}
	}
	
	
	//////////////////////////////////////////////////////////////////
	// read pmc sqs models
	if (config->get_need_to_estimate_pm() || got_filter_spectra)
	{
		if (! model.get_ind_pmcsqs_was_intialized() )
		{
			cout << "Error: could not find PMC and SQS models for " << config->get_model_name() << endl;
			cout << "Cannot perform precursor mass correction and charge determiniation!" << endl;
			exit(1);
		}
	}



	///////////////////////////////////////////////////////////////////
	// FILTER SPECTRA
	if (got_filter_spectra)
	{
		int num_written =0;
		int num_read = 0;
		PMCSQS_Scorer *pmcsqs = (PMCSQS_Scorer *)model.get_pmcsqs_ptr();

		pmcsqs->output_filtered_spectra_to_mgfs(config, files, out_dir, filter_prob, num_written, num_read);
		
		time_t curr_time = time(NULL);
		double elapsed_time = (curr_time - start_time);
		cout << "Processed " << files.size() << " (" << num_read << " spectra)." << endl;
		cout << "Wrote " << num_written << " spectra to mgfs in " << out_dir << endl;
		cout << "Elapsed time " << fixed << elapsed_time << " seconds." << endl;
		return 0;
	}


	//////////////////////////////////////////////////////////////////
	// PMCSQS
	if (ind_pmcsqs)
	{
		perform_pmcsqs_on_list_of_files(model, files, file_start_idx, ( got_output_file ? out_file_stream : cout));
		return 0;
	}
	
	//////////////////////////////////////////////////////////////////
	// DENOVO AND TAGS

	config->set_filter_flag(0);
	if (tag_length<=0)
	{
		new_perform_denovo_on_list_of_files(model, files, file_start_idx, num_solutions, 7, 16, 99999.0, true, filter_prob, false, false,
			( got_output_file ? out_file_stream : cout));
	}
	else
	{
		perform_tags_on_list_of_files(model, files, file_start_idx, num_solutions, tag_length, true, 
			filter_prob, false, false, ( got_output_file ? out_file_stream : cout));	
		
	}

	return 0;
}





