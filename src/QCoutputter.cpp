#include "QuickClustering.h"

void QCOutputter::init(string _name , string _dir, int _batch_idx,
			  mass_t min_m_over_z, mass_t max_m_over_z,
			  float min_similarity, int min_cluster_size)
{
	batch_idx = _batch_idx;
	

	if (batch_idx>=0)
	{

		ostringstream oss;
		oss << batch_idx;
		batch_str = oss.str();

		dir = _dir;
		name = _name;
		string sum_name = dir + "/" + name + "_" + batch_str + "_sum.txt";
		string list_name = dir + "/" + name + "_" + batch_str + "_list.txt";
		string param_name = dir + "/" + name + "_" + batch_str + "_params.txt";

		summary_stream.open(sum_name.c_str(),ios::out);               
		file_list_stream.open(list_name.c_str(),ios::out);
		fstream param_stream(param_name.c_str(),ios::out);

		
		if (! summary_stream.is_open() || 
			! file_list_stream.is_open() ||
			! param_stream.is_open() )
		{
			cout << "Error: couldn't open outputter file streams!" << endl;
			exit(1);
		}

		param_stream << "batch:     " << batch_idx << endl;
		param_stream << "out dir:   " << dir << endl;
		param_stream << "min m/z:   " << min_m_over_z << endl;
		param_stream << "max m/z:   " << max_m_over_z << endl;
		param_stream << "min similarity: " << min_similarity << endl;
		param_stream << "min cluster size: " << min_cluster_size << endl;

		param_stream.close();
	}
	else
	{
		dir = _dir;
		name = _name;
		cout << "Sending output MGF files called " << name << " to dir: " << dir << endl;
	}
}

QCOutputter::~QCOutputter()
{
	if (mgf_stream.is_open())
		mgf_stream.close();

	if (cluster_file_stream.is_open())
		cluster_file_stream.close();

	if (summary_stream.is_open())
		summary_stream.close();

	if (file_list_stream.is_open())
		file_list_stream.close();
}


void QCOutputter::output_basic_spectrum_to_mgf(BasicSpectrum &bs, Config *config)
{
	// check if a new file should be opened
	if (total_spectra_counter == 0 ||
		spectra_counter == NUM_CLUSTERS_PER_FILE)
	{
		if (mgf_stream.is_open())
		{
			mgf_stream.close();
			cluster_file_stream.close();
		}

		file_counter++;
		ostringstream oss;
		oss << file_counter;
		mgf_name = name +  "_" + oss.str() + ".mgf";
		string mgf_path = dir + "/" + mgf_name;
		mgf_stream.open(mgf_path.c_str(),ios::out);

		string list_name = dir + "/" + name +  "_" + oss.str() + "_list.txt";
		cluster_file_stream.open(list_name.c_str(),ios::out);

		spectra_counter=0;
	}

	const SingleSpectrumFile * ssf = bs.ssf;
	if (ssf->type == MZXML)
	{
		MZXML_single* mzxml_ssf = (MZXML_single *)ssf;

		cluster_file_stream << mzxml_ssf->file_idx << "\t" << mzxml_ssf->scan_number << "\t" <<
			mzxml_ssf->m_over_z << "\t" << mzxml_ssf->charge ;
	}
	else if (ssf->type == DAT)
	{
		DAT_single* dat_ssf = (DAT_single *)ssf;
		cluster_file_stream << dat_ssf->mzxml_file_idx << "\t" << dat_ssf->scan_number << "\t" <<
				dat_ssf->m_over_z << "\t" << dat_ssf->charge;
	}
	else if (ssf->type == MGF)
	{
		MGF_single* mgf_ssf = (MGF_single *)ssf;
		cluster_file_stream << mgf_ssf->file_idx << "\t" << mgf_ssf->idx_in_file << "\t" <<
			mgf_ssf->m_over_z << "\t" << mgf_ssf->charge;
	}
	else if (ssf->type == DTA)
	{
		cluster_file_stream << ssf->single_name << "\t" << ssf->m_over_z << "\t" << ssf->charge;
	}

	cluster_file_stream << "\t" << ssf->sqs << endl;
	bs.output_to_mgf(mgf_stream,config);

	spectra_counter++;
	total_spectra_counter++;

}



