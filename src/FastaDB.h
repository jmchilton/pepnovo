#ifndef __FASTA_DB_H__
#define __FASTA_DB_H__

#include "Config.h"
#include "includes.h"


struct PeptideLocation {
	int aa_loc; // in the whole sequence database
	int protein_idx;
	int loc_in_protein;
};



struct TagListPointer {
	TagListPointer() : num_locations(0), list_start_idx(-1) {};

	int num_locations;
	int list_start_idx;
};

typedef map< int , TagListPointer, less<int> > INT2TLP_MAP;

class FastaDB {
public:

	FastaDB();
	~FastaDB();

	// creates all relevant data structures from the fasta file
	// includes the sequences (stored as aa - ints), protein names,
	// and tag hashes
	void create_db_from_fasta(char *file, Config *con,
		bool create_tags = true, int min_length=3, int max_length=6);

	int get_num_cands_with_mass(float mass, float tolerance) const;

	int get_total_seq_length() const { return all_aa_seqs.size(); }

	int get_max_tag_length() const { return max_tag_length; }
	int get_min_tag_length() const { return min_tag_length; }

	void print_stats() const;
	void print_protein_names() const;

	int calc_tag_index(const vector<int>& tag_aas) const 
	{	int i,idx=0;
		const vector<int>& org_aa = config->get_org_aa();
		for (i=0; i<tag_aas.size()-1; i++)
		{
			idx += aa_codes[org_aa[tag_aas[i]]];
			idx *= mult_val;
		}
		idx+=aa_codes[org_aa[tag_aas[i]]];
		return idx;
	}
	

	// gets number and list of tag locations
	int get_tag_locations(const vector<int>& tag, int **list) const
	{
		const int tag_length = tag.size();
		if (tag_length<min_tag_length || tag_length>max_tag_length)
		{
			cout << "Tag length " << tag_length << " not supported!" << endl;
			exit(0);
		}
		const int idx = calc_tag_index(tag);

		INT2TLP_MAP::const_iterator iter = tag_maps[tag_length].find(idx);

		if (iter == tag_maps[tag_length].end())
			return 0;

		*list = (int *)(&tag_locations[tag_length][(*iter).second.list_start_idx]);
		return (*iter).second.num_locations;
	}

	// returns the protein idx for an amino acid location, -1 if invalid locaiton is given
	int get_protein_number(int loc, int *idx_in_protein = NULL) const
	{
		if (loc<0 || loc>all_aa_seqs.size())
			return -1;

		INT_MAP::const_iterator it;
		it = aa_seq_starts.lower_bound(loc);
		if (it != aa_seq_starts.end())
		{
			if (idx_in_protein)
			{
				if ((*it).first == loc)
				{
					*idx_in_protein = 0;
					return (*it).second;
				}
				else
				{
					it--;
					*idx_in_protein = loc - (*it).first;
					return (*it).second;
				}
			}
			else
			{
				return  (*it).first == loc ? (*it).second : (*it).second -1;
			}
			
		}

		return -1;
	}

	// gets pointer to sequnce location
	const int * get_aa_seq_pointer(int loc=0) const
	{
		return &all_aa_seqs[loc];
	}

	void print_aas_at_loc(int loc_idx, int num_aas = 10) const;


	Config *get_config() { return config; }

	void read_FastaDB(const char *file_name, Config *con);
	void write_FastaDB(const char *file_name) const;

private:
	Config *config;
	char *fasta_file;
	int min_tag_length, max_tag_length;

	int mult_val;  // the number in which we multiply the tag indices

	vector<int> aa_codes;        // for each amino acid holds an integer code
	INT_MAP aa_seq_starts;      // the starting idx of each sequence
	vector< int> protein_name_starts; // the start position of each protein name

	vector<int> all_aa_seqs; // holds the concatenated aa sequences for all proteins (seq_starts points to
						 // positions in this array)

	vector<char> all_protein_names; // holds concatenated char sequences of all protein names
	                                // (protein_name_starts points to positions in this array).

	vector<INT2TLP_MAP> tag_maps; // holds a map to each tag locations

	vector<vector<int> > tag_locations; // holds for each tag length a list of locations
										// which are organized according to the order of tags
										// (the entries of tag_hashes point to places in this list).
};



#endif


