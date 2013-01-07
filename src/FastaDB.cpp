#include "FastaDB.h"



FastaDB::FastaDB()
{
	int a,i;

	config = NULL;
	fasta_file=NULL;
	a=0;

	aa_codes.resize(Val+1);
	for (i=0; i<=Val; i++)
		aa_codes[i]=0;

	for (i=Xle; i<=Val; i++)
		aa_codes[i]=a++;

	aa_codes[Ile]=Xle; // 
	aa_codes[Leu]=Xle;

	mult_val=a;
}

FastaDB::~FastaDB()
{
	if (fasta_file)
		delete [] fasta_file;
}


struct TagLoc {
	bool operator< (const TagLoc& other) const
	{
		if (idx<other.idx)
			return true;
		if (idx>other.idx)
			return false;
		if (loc<other.loc)
			return true;
		return false;
	}
	int idx;
	int loc;
};

/*******************************************************************
  creates all relevant data structures from the fasta file
  includes the sequences (stored as aa - ints), protein names,
  and tag hashes (direct hash table)
********************************************************************/
void FastaDB::create_db_from_fasta(char *file_name, Config *con,
		bool create_tags, int min_length, int max_length)
{
	char buff[1024];
	int file_size;
	int seq_p=0;

	this->config = con;
	const vector<int>& char2aa = config->get_char2aa();

	ifstream file (file_name, ios::in|ios::ate);
	if (file.is_open())
	{
		file_size = file.tellg();
		file.seekg (0, ios::beg);
	}
	else
	{
		cout << "Error: reading!"<< file_name << endl;
		exit(1);
	}
	
	min_tag_length = min_length;
	max_tag_length = max_length;

	if (create_tags && (min_length <3 || max_length>6 || min_tag_length>max_tag_length))
	{
		printf("Tag length must be 3-6 !\n");
		exit(1);
	}

	fasta_file = new char[strlen(file_name)+1];
	strcpy(fasta_file,file_name);

	aa_seq_starts.clear();     
	protein_name_starts.clear(); 
	all_aa_seqs.clear();
	all_protein_names.clear();


	// add sequence terminatng symbol -1
	// before first sequence
	all_aa_seqs.push_back(-1);

	file.getline(buff,1024);
	while(1)
	{
		if (file.eof())
			break;

		if (file.gcount()>0 && buff[0] != '>')
		{
			file.getline(buff,1024);
			if (file.gcount() <= 0)
				break;

			continue;
		}

		// push the protein name
		int len=strlen(buff);
		int name_start_idx=all_protein_names.size();
		int i;
		for (i=1; i<len; i++)
		{
			if (buff[i] != '\n')
				all_protein_names.push_back(buff[i]);
		}
		all_protein_names.push_back('\0');

		aa_seq_starts.insert(INT_MAP::value_type(all_aa_seqs.size(),protein_name_starts.size()));
		protein_name_starts.push_back(name_start_idx);
		

		// read protein sequence
		while ( 1)
		{
			file.getline(buff,1024);
			if (file.gcount()<=0 || buff[0] == '>')
				break;

			const int len=strlen(buff);
			for (i=0; i<len; i++)
				if (buff[i]>= 'A' && buff[i]<'Z')
					all_aa_seqs.push_back(char2aa[buff[i]]);
		}
		// add sequence terminatng symbol -1
		all_aa_seqs.push_back(-1);
	}

	// creates tags using a vector of vector, then transforms it into a tag_hash
	// and a sequence of tag locations
	if (create_tags)
	{
		tag_maps.resize(max_tag_length+1);
		tag_locations.resize(max_tag_length+1);

		int tag_length;

		for (tag_length=min_tag_length; tag_length<=max_tag_length; tag_length++)
		{
			vector<TagLoc> tag_locs;
			vector<int> tag;
			const int max_loc = all_aa_seqs.size()-tag_length;
			const int max_tag_idx = static_cast<int>(pow(static_cast<double>(mult_val),
														 static_cast<double>(tag_length)));
			int i;

			tag_locs.reserve(all_aa_seqs.size());
			tag.resize(tag_length);

			for (i=0; i<max_loc; i++)
			{
				int j;
				for (j=0; j<tag_length; j++)
				{
					if (all_aa_seqs[i+j]<Ala)
						break;

					tag[j]=all_aa_seqs[i+j];
				}
				if (j<tag_length)
					continue;

				int idx=calc_tag_index(tag);
				
				TagLoc tl;
				tl.idx = idx;
				tl.loc = i;
				tag_locs.push_back(tl);	
			}
			sort(tag_locs.begin(),tag_locs.end());
		
			// transfer tags to map
			int total_locs=tag_locs.size();

			tag_locations[tag_length].clear();
			tag_locations[tag_length].reserve(total_locs);
			

			tag_maps[tag_length].clear();

			int idx=-1;
			int start=-1;
			int num_locs=0;
			for (i=0; i<tag_locs.size(); i++)
			{
				if (tag_locs[i].idx != idx)
				{
					if (num_locs>0)
					{
						TagListPointer tlp;
						tlp.list_start_idx = start;
						tlp.num_locations = num_locs;

						tag_maps[tag_length].insert(INT2TLP_MAP::value_type(idx,tlp));
					}

					idx=tag_locs[i].idx;
					start=i;
					num_locs=1;
					
					tag_locations[tag_length].push_back(tag_locs[i].loc);
				}
				else
				{
					tag_locations[tag_length].push_back(tag_locs[i].loc);
					num_locs++;
				}
			}

		
		}
	}
}


/*************************************************************
	reads all info from dat file.
**************************************************************/
void FastaDB::read_FastaDB(const char *file_name, Config *con)
{
	fstream ifs(file_name, ios::in|ios::binary);
	if (! ifs.good())
	{
		cout << "Error: couldn't open for writing: "<< file_name << endl;
		exit(1);
	}

	config = con;

	// read file name
	int name_len;
	ifs.read(reinterpret_cast<char *>(&name_len),sizeof(int));
	fasta_file = new char[name_len+1];
	ifs.read(reinterpret_cast<char *>(fasta_file),name_len*sizeof(char));
	fasta_file[name_len]='\0';

	// read protein names
	int num_proteins;
	int all_protein_names_length;
	int all_seqs_length;

	ifs.read(reinterpret_cast<char *>(&num_proteins),sizeof(int));
	protein_name_starts.resize(num_proteins);
	ifs.read(reinterpret_cast<char *>(&protein_name_starts[0]),sizeof(int) * num_proteins);

	ifs.read(reinterpret_cast<char *>(&all_protein_names_length), sizeof(int));
	all_protein_names.resize(all_protein_names_length);
	ifs.read(reinterpret_cast<char *>(&all_protein_names[0]),sizeof(char) * all_protein_names_length);

	// read aa_ses_starts map 
	aa_seq_starts.clear();
	int i;
	for (i=0; i<num_proteins; i++)
	{
		int aa_pos, prot_num;
		ifs.read(reinterpret_cast<char *>(&aa_pos),sizeof(int));
		ifs.read(reinterpret_cast<char *>(&prot_num),sizeof(int));
		aa_seq_starts.insert(INT_MAP::value_type(aa_pos,prot_num));
	}

	// write all seqs
	ifs.read(reinterpret_cast<char *>(&all_seqs_length), sizeof(int));
	all_aa_seqs.resize(all_seqs_length);
	ifs.read(reinterpret_cast<char *>(&all_aa_seqs[0]),sizeof(int)*all_seqs_length);

	// read tags
	
	ifs.read(reinterpret_cast<char *>(&min_tag_length),sizeof(int));
	ifs.read(reinterpret_cast<char *>(&max_tag_length),sizeof(int));
	tag_locations.resize(max_tag_length+1);
	tag_maps.resize(max_tag_length+1);

	int t;
	for (t=min_tag_length; t<=max_tag_length; t++)
	{
		int i;
		int num_tag_locations;

		ifs.read(reinterpret_cast<char *>(&num_tag_locations),sizeof(int));
		tag_locations[t].resize(num_tag_locations);
		ifs.read(reinterpret_cast<char *>(&tag_locations[t][0]),sizeof(int)*num_tag_locations);

		int map_size;

		ifs.read(reinterpret_cast<char *>(&map_size),sizeof(int));
		tag_maps[t].clear();

		for (i=0; i<map_size; i++)
		{
			int idx;
			TagListPointer tlp;

			ifs.read(reinterpret_cast<char *>(&idx),sizeof(int));
			ifs.read(reinterpret_cast<char *>(&tlp),sizeof(TagListPointer));
			tag_maps[t].insert(INT2TLP_MAP::value_type(idx,tlp));
		}
	}


	ifs.close();

}

/*************************************************************
	writes all info to dat file.
**************************************************************/
void FastaDB::write_FastaDB(const char *file_name) const
{
	fstream ofs(file_name, ios::out|ios::binary);
	if (! ofs.good() || ! ofs.is_open())
	{
		cout << "Error: couldn't open for writing: "<< file_name << endl;
		exit(1);
	}

	// write file name
	if (! fasta_file)
	{
		cout << "Error: must initialize from a fasta file!" << endl;
		exit(0);
	}
	int name_len = strlen(this->fasta_file);
	ofs.write(reinterpret_cast<const char *>(&name_len),sizeof(int));
	ofs.write(reinterpret_cast<const char *>(fasta_file),name_len*sizeof(char));
	
	// write protein names
	int num_proteins = protein_name_starts.size();
	int all_protein_names_length = all_protein_names.size();
	int all_seqs_length = all_aa_seqs.size();

	ofs.write(reinterpret_cast<char *>(&num_proteins),sizeof(int));
	ofs.write(reinterpret_cast<const char *>(&protein_name_starts[0]),sizeof(int) * num_proteins);
	ofs.write(reinterpret_cast<char *>(&all_protein_names_length), sizeof(int));
	ofs.write(reinterpret_cast<const char *>(&all_protein_names[0]),sizeof(char) * all_protein_names_length);

	// write aa_seq_starts map
	INT_MAP::const_iterator it;
	for (it=aa_seq_starts.begin(); it!= aa_seq_starts.end(); it++)
	{
		ofs.write(reinterpret_cast<const char *>(&it->first),sizeof(int));
		ofs.write(reinterpret_cast<const char *>(&it->second),sizeof(int));
	}

	ofs.write(reinterpret_cast<char *>(&all_seqs_length), sizeof(int));
	const char *seqs=(char *)(&all_aa_seqs[0]);
	ofs.write(seqs,sizeof(int)*all_seqs_length);

	// write tags
	
	ofs.write(reinterpret_cast<const char *>(&min_tag_length),sizeof(int));
	ofs.write(reinterpret_cast<const char *>(&max_tag_length),sizeof(int));
	
	int t;
	for (t=min_tag_length; t<=max_tag_length; t++)
	{
		int num_tag_locations = tag_locations[t].size();

		ofs.write(reinterpret_cast<const char *>(&num_tag_locations),sizeof(int));
		ofs.write(reinterpret_cast<const char *>(&tag_locations[t][0]),sizeof(int)*num_tag_locations);
		
		int map_size = tag_maps[t].size();
		ofs.write(reinterpret_cast<const char *>(&map_size),sizeof(int));

		INT2TLP_MAP::const_iterator it;
		for (it=tag_maps[t].begin(); it!= tag_maps[t].end(); it++)
		{
			ofs.write(reinterpret_cast<const char *>(&it->first),sizeof(int));
			ofs.write(reinterpret_cast<const char *>(&it->second),sizeof(TagListPointer));
		}
	}


	ofs.close();
}


/*
// returns a merged list of the locations of all the tags
// the indices in the tag_loc records are according to the position
// in the tag_idxs array
list_record FastaDB::get_merged_tag_locs(const vector<int *>& tag_seqs)
{
	int i;
	list_record merged;
	vector<list_record> lists;
	lists.resize(tag_seqs.size());
	
	for (i=0; i<tag_seqs.size(); i++)
	{
		int j;

		const vector<int>& locations =get_tag_locations(tag_seqs[i]);
		lists[i].list.resize(locations.size());
		lists[i].size = locations.size();
		lists[i].free_memory = true;

		for (j=0; j<locations.size(); j++)
		{
			lists[i].list[j].loc= locations[j];
			lists[i].list[j].tag_idx = i;
		}
	}
	merge_lists(lists,merged);
	return merged;
}

*/

void FastaDB::print_stats() const
{
	printf("FILE          : %s\n",fasta_file);
	printf("SEQUENCES     : %d\n",protein_name_starts.size());
	printf("AMINO ACIDS   : %d\n",all_aa_seqs.size());
//	printf("TAG LOCATIONS : %d\n",num_tag_locations);
}


void FastaDB::print_protein_names() const
{
	int i;

	for (i=0; i<protein_name_starts.size(); i++)
	{
		cout << i << " " << 
			&all_protein_names[protein_name_starts[i]]<< endl;
	}
}



int FastaDB::get_num_cands_with_mass(float mass, float tolerance) const
{
	const float min_mass = mass - tolerance;
	const float max_mass = mass + tolerance;
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	int i;
	int n=0;

	for (i=0; i<all_aa_seqs.size(); i++)
	{
		int j=0;
		float m=0;
		while (all_aa_seqs[i+j]>0 && m<min_mass)
		{
			m+=aa2mass[all_aa_seqs[i+j]];
			j++;
		}

		if (m>=min_mass && m<=max_mass)
			n++;
	}

	return n;
}



void FastaDB::print_aas_at_loc(int loc_idx, int num_aas) const
{
	int i;

	for (i=0; i<num_aas; i++)
	{
		const int sym = all_aa_seqs[i+loc_idx];
		if (sym<0)
		{
			cout << "$$$";
			break;
		}

		cout<<config->get_aa2label()[sym];
	}
	cout << endl;
}
