#include "FileManagement.h"
#include "QuickClustering.h"
#include "includes.h"

void print_help()
{
	printf("To run:\n");
	printf("mzXML2MGF.exe <mzXML>\n");
	printf("\nThe program writes an mgf file to the same location.\n");
}

int main(int argc, char **argv)
{
	Config config;
	FileManager fm;
	FileSet     fs;
 
	char mzxml_name[256];
	char mgf_name[256];
	
	config.init_with_defaults();
	config.set_need_to_normalize(0);


	// read command line arguments
	if (argc<2)
	{
		print_help();
		exit(0);
	}

	strcpy(mzxml_name,argv[1]);
	vector<string> list;
	list.push_back(mzxml_name);

	fm.init_from_list(&config,list);
	fs.select_all_files(fm);

	
	strncpy(mgf_name,mzxml_name,strlen(mzxml_name)-5);
	strcat(mgf_name,"mgf");

	printf("Converting to:\n%s\n",mgf_name);
	
	fstream mgf_stream;
	mgf_stream.open(mgf_name,ios::out);

	// read spectra
	int i;
	const vector<SingleSpectrumFile *>& ssfs = fs.get_ssf_pointers();
	BasicSpecReader bsr;

	for (i=0; i<ssfs.size(); i++)
	{
		static QCPeak peaks[5000];
		BasicSpectrum bs;
		MZXML_single *ssf = (MZXML_single *)ssfs[i];

		bs.peaks = peaks;
		bs.ssf = ssf;
					
		ostringstream oss;
		oss << "scan_" << ssf->scan_number;
		ssf->single_name = oss.str();

		bs.num_peaks = bsr.read_basic_spec(&config,fm,ssf,peaks);
		
		if (ssf->scan_number<0)
		{
			cout << "Error: no scan number read from mzXML!!!" << endl;
			exit(1);
		}

		cout << "scan: " << ssf->scan_number << " " << bs.num_peaks << endl;

		bs.output_to_mgf(mgf_stream,&config);
	//	bs.output_to_mgf(cout,&config);
	}


	return 0;
}
