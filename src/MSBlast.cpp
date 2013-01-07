#include "MSBlast.h"
#include "DeNovo.h"


struct MSB_line_ptr {
	MSB_line_ptr() : seq(NULL), max_score(0), length(0) {};

	bool operator < (const MSB_line_ptr& other) const
	{
		return max_score>other.max_score;
	}

	char * seq;
	float max_score;
	int length;
};



void extract_dta_from_mgf(char *mgf_file, ScoreModel *model,  vector<int>& file_idxs)
{
	Config *config = model->get_config();
	FILE *mgf_stream = fopen(mgf_file,"r");

	if (! mgf_file)
	{
		printf("Error: couldn't open mgf %s file for reading!\n",mgf_file);
		exit(1);
	}

	sort(file_idxs.begin(),file_idxs.end());
	if (file_idxs.size() == 0)
		return;

	int pos = 0;
	int next_file_idx = file_idxs[0];
	int count=0;
	while(1)
	{
		BasicSpectrum s;
	

		if (! s.read_and_init_from_MGF_stream(mgf_stream,model->get_config()))
		{
			printf("Done reading, read a total of %d dta files.\n",count);
			break;
		}
		count++;

		if (count == next_file_idx)
		{
			s.write_dta(s.get_file_name());
			pos++;
			if (pos == file_idxs.size()) 
				break;

			next_file_idx = file_idxs[pos];
		}

	}

	fclose(mgf_stream);
}



void generate_MSB_candidates_for_dta(char *file_name, ScoreModel *model)
{
	BasicSpectrum s;
	PrmGraph prm;
	Config *config = model->get_config();

	s.read_and_init(file_name,config,ISB);
	model->iso_model.calc_all_peaks_iso_scores(&s);
	const float opt_pm = calc_corrected_pm_with_19(&s,model,1.5,0.6,2);
		
	prm.create_graph(&s,model,opt_pm,0.5,2,true);
		
	if (prm.get_num_nodes()<4)
	{
		printf("only %d nodes in %s\n",prm.get_num_nodes(),s.get_file_name());
		exit(1);
	}
	
	DeNovoPath	path=prm.find_best_denovo_path();

	if (path.seq_length<6)
	{
		printf("Only %d aa in seq for %s\n",path.seq_length,s.get_file_name());
		exit(1);
	}

	path.calc_aa_probs(&prm,&s,model,opt_pm,1);

	vector<MSB_seq> final_seqs;
	generate_msblast_candidates(&prm,model,&s,opt_pm,path,final_seqs,10,
		8, 2);

}


/**********************************************************************
// generates a file with all the sequences from the mgf dta's
// also generates a file with all the denovo results
***********************************************************************/
void generate_mgf_sequences(char *mgf_file, ScoreModel *model, 
			float min_cand_score, int num_candidates, int max_size_bytes, char * name,
			int charge)
{
	char msb_out[256],denovo_out[256];
	Config *config = model->get_config();
	FILE *mgf_stream = fopen(mgf_file,"r");

	vector<MSB_line_ptr> all_res;
	all_res.clear();


	
	if (! mgf_stream)
	{
		printf("Error: couldn't open mgf %s file for reading!\n",mgf_file);
		exit(1);
	}


	if (name)
	{
		sprintf(msb_out,"%s_msb.txt",name);
		sprintf(denovo_out,"%s_denovo.txt",name);
	}
	else
	{
		sprintf(msb_out,"results_msb.txt");
		sprintf(denovo_out,"results_denovo.txt");
	}

	if (! fopen(msb_out,"w") )
	{
		printf("Error: couldn't open %s for writing!\n",msb_out);
		exit(1);
	}

	FILE *den_stream=fopen(denovo_out,"w");
	if (! den_stream )
	{
		printf("Error: couldn't open %s for writing!\n",denovo_out);
		exit(1);
	}


	fprintf(den_stream,"PepNovo2MSB v1.02 - de Novo peptide sequencing for MS-Blast.\n");
	fprintf(den_stream,"Copyright (c) 2005 Ari Frank, Pavel Pevzner and The Regents of the University of California, All Rights Reserved.");
	fprintf(den_stream," This software is not to be used commercially without the explicit permission of the authors.\n");
	fprintf(den_stream,"\nFor more information on the algorithms see: Frank, A. and Pevzner, P. \"PepNovo: De Novo Peptide Sequencing via Probabilistic Network Modeling.\" Analytical Chemistry, 77:964-973, 2005.\n");
	fprintf(den_stream,"\nPlease send comments and bug reports to Ari Frank (arf@cs.ucsd.edu).\n\n");

	fclose(den_stream);

	int count=0;
	while(1)
	{
		int i;
		BasicSpectrum s;
		PrmGraph prm;
        DeNovoPath path;

		if (! s.read_and_init_from_MGF_stream(mgf_stream,model->get_config()))
		{
			printf("Done reading, read a total of %d dta files.\n",count);
			break;
		}

		if (s.get_num_peaks()<5)
			continue;

		count++;

		if (charge >0 && s.get_charge() != charge)
			continue;

		if (s.get_charge()>3)
		{
			printf("%d\t%s\t0\n",count,s.get_file_name());
			continue;
		}


		model->iso_model.calc_all_peaks_iso_scores(&s);
		const float opt_pm = calc_corrected_pm_with_19(&s,model,1.5,0.6,2);
		
		prm.create_graph(&s,model,opt_pm,0.6,2,true);
		
		if (prm.get_num_nodes()<4)
		{
			printf("%d\t%s\t0\n",count,s.get_file_name());
			continue;
		}
	
//		path=prm.find_best_path_with_forbidden_vertices(opt_pm,0);
		path=prm.find_best_denovo_path();

		if (path.seq_length<6)
		{
			printf("%d\t%s\t0\n",count,s.get_file_name());
			continue;
		}
		if (path.seq_length>40)
		{
			printf("%d\t%s\t0\n",count,s.get_file_name());
			continue;
		}

		path.calc_aa_probs(&prm,&s,model,opt_pm,1);

		vector<MSB_seq> final_seqs;
		generate_msblast_candidates(&prm,model,&s,opt_pm,path,final_seqs,num_candidates,
			8, min_cand_score);
		printf("%d\t%s\t%d\n",count,s.get_file_name(),final_seqs.size());

		if (path.seq_length)
			path.append_denovo_results(denovo_out,config,(char *)s.get_file_name());

		// results to screen
		for (i=0; i<final_seqs.size(); i++)
			final_seqs[i].print_short(config);

		float percent = prm.get_percent_explained_intensity(path);

		// results to files
		if (final_seqs.size()>0)
		{
			// extract scan number
			const char *dta_name = s.get_file_name();

			int p=strlen(dta_name)-4;
			int last_dot_loc=-1;

			while (p>0)
			{
				if (dta_name[p]=='.' && dta_name[p+1]=='d' && 
					dta_name[p+2]=='t' && dta_name[p+3]=='a')
				{
					last_dot_loc=p;
					break;
				}
				p--;
			}

			p = last_dot_loc-1;
			int dot_count=0;
			char scan_name[20];

			while (p>0 && dot_count<3)
			{
				if (dta_name[--p] == '.')
					dot_count++;
			}
			
			if (dot_count == 3)
			{
				int q;
				for (q=0; q<last_dot_loc-p-1; q++)
					scan_name[q]=dta_name[1+p+q];
				scan_name[q]='\0';
			}
			else
			{
				char buff[128];

				if (strlen(dta_name)<3)
				{
					sprintf(scan_name,"...");
				}
				else
				{
					int a;
					int start=0;
					strcpy(buff,dta_name);
					for (a=0; a<strlen(buff); a++)
					{
						if (buff[a] != '0' && ( buff[a]<'1' || buff[a]>'9'))
						{
							buff[a]='.';
						}
						else if (start == 0)
							start=a;
					}
						
					if (start == strlen(buff) -1)
					{
						sprintf(scan_name,"...");
					}
					else
						sprintf(scan_name,"%s",buff+start);					
				}
				
			}

		
			char msb_line[1024];
			sprintf(msb_line,"%.2f ",s.get_original_pm_with_19()-1);
			strcat(msb_line,scan_name);

			int k;
			float max_score=0;
			for (k=0; k<final_seqs.size(); k++)
				if (final_seqs[k].expected_score>max_score)
					max_score = final_seqs[k].expected_score;

			if (max_score == 0)
				continue;

			float total_inten = s.get_total_intensity();
			char intensity_line[32];
			sprintf(intensity_line," %f %f ",total_inten,percent);
			strcat(msb_line,intensity_line);

			char seq_score[10];
			sprintf(seq_score," %.1f ",max_score);
			strcat(msb_line,seq_score);

			int total_cand_length=0;
			for (k=0; k<final_seqs.size(); k++)
			{
				char cand_seq[64];
				
				final_seqs[k].write_seq_to_string(config,cand_seq);
				strcat(msb_line,cand_seq);
				total_cand_length+=strlen(cand_seq);
			}

			if (! max_size_bytes)
			{
				write2file(msb_out,msb_line);
			}
			else
			{
				// add to list 
				strcat(msb_line,"\n");
				char *seq_mem=new char[strlen(msb_line)+1];
				MSB_line_ptr mlp;
				mlp.max_score=max_score;
				mlp.seq = seq_mem;
				mlp.length = total_cand_length;

				strcpy(seq_mem,msb_line);

				all_res.push_back(mlp);
			}
		}
	}
	
	if (max_size_bytes>0)
	{
		FILE *stream;
		int out_size=0;
		int i;

		sort(all_res.begin(),all_res.end());

		stream=fopen(msb_out,"w");

		if (! stream)
		{
			printf("Error: couldn't open %s for writing!\n",msb_out);
			exit(1);
		}

		for (i=0; i<all_res.size(); i++)
		{
			if (out_size + all_res[i].length> max_size_bytes-1)
				break;

			fprintf(stream,"%s",all_res[i].seq);
			out_size += all_res[i].length;
		}
		fclose(stream);
	}
}





/**********************************************************************
// generates a file with all the sequences from the mgf dta's
// also generates a file with all the denovo results
***********************************************************************/
void generate_mgf_sequences(FileManager& fm, vector<cluster_info>& clusters, 
							ScoreModel *model, float min_cand_score, int num_candidates, 
							int max_size_bytes, char * name, int charge)
{
	char msb_out[256],denovo_out[256];

	Config *config = model->get_config();

	vector<MSB_line_ptr> all_res;
	all_res.clear();

	if (name)
	{
		sprintf(msb_out,"%s_msb.txt",name);
		sprintf(denovo_out,"%s_denovo.txt",name);
	}
	else
	{
		sprintf(msb_out,"results_msb.txt");
		sprintf(denovo_out,"results_denovo.txt");
	}

	if (! fopen(msb_out,"w") )
	{
		printf("Error: couldn't open %s for writing!\n",msb_out);
		exit(1);
	}

	FILE *den_stream=fopen(denovo_out,"w");
	if (! den_stream )
	{
		printf("Error: couldn't open %s for writing!\n",denovo_out);
		exit(1);
	}


	fprintf(den_stream,"PepNovo2MSB v1.03 - de Novo peptide sequencing for MS-Blast.\n");
	fprintf(den_stream,"Copyright (c) The Regents of the University of California, All Rights Reserved.");
	fprintf(den_stream," This software is not to be used commercially without the explicit permission of the authors.\n");
	fprintf(den_stream,"\nCitation: Frank, A. and Pevzner, P. \"PepNovo: De Novo Peptide Sequencing via Probabilistic Network Modeling.\" Analytical Chemistry, 77:964-973, 2005.\n");
	fprintf(den_stream,"\nPlease send comments and bug reports to Ari Frank (arf@cs.ucsd.edu).\n\n");

	fclose(den_stream);

	int c;
	for (c=0; c<clusters.size(); c++)
	{
		ClusterSpectrum cs;
		PrmGraph prm;
        DeNovoPath path;

		const int cluster_charge = clusters[c].org_spectra_SSF_pointers[0]->charge;
		const int cluster_size = clusters[c].org_spectra_SSF_pointers.size();

		if (clusters[c].total_num_spectra == 1 && 
			clusters[c].org_spectra_SSF_pointers[0]->sqs >=0 &&
			clusters[c].org_spectra_SSF_pointers[0]->sqs < 0.2)
		{
			printf("Cluster %d\t (%d spectra, charge %d) - low spectral quality\n",c,
				clusters[c].total_num_spectra,cluster_charge);
			continue;
		}
		
		if (charge >0 && cluster_charge != charge)
			continue;

		if (cluster_charge>3)
		{
			printf("Cluster %d\t (%d spectra, charge %d)\n",c,
				clusters[c].total_num_spectra,cluster_charge);
			continue;
		}


		cs.full_initialization(clusters[c],fm,model);
	
		
		prm.create_graph_from_cluster(&cs,model,cs.get_corrected_pm_with_19(),0.6,true);
		
		if (prm.get_num_nodes()<4)
		{
			printf("Cluster %d\t (%d spectra, charge %d)\n",c,
				clusters[c].total_num_spectra,cluster_charge);
			continue;
		}

	
		path=prm.find_best_denovo_path();

		if (path.seq_length<6)
		{
			printf("Cluster %d\t (%d spectra, charge %d) - prediction too short\n",c,
				clusters[c].total_num_spectra,cluster_charge);
			continue;
		}
		if (path.seq_length>40)
		{
			printf("Cluster %d\t (%d spectra, charge %d)\n",c,
				clusters[c].total_num_spectra,cluster_charge);
			continue;
		}

	

		path.calc_aa_probs(&prm,&cs,model,cs.get_corrected_pm_with_19(),1);

		vector<MSB_seq> final_seqs;
		generate_msblast_candidates(&prm,model,&cs,cs.get_corrected_pm_with_19(),
			path,final_seqs,num_candidates, 8, min_cand_score);
		printf("Cluster %d\t (%d spectra, charge %d)\t%d\n",c,
			clusters[c].total_num_spectra,cluster_charge,final_seqs.size());

		float percent = prm.get_percent_explained_intensity(path);

		if (path.seq_length)
		{
			char c_name[128];
			sprintf(c_name,"Cluster %d (%d spectra, charge %d)",c,
				clusters[c].total_num_spectra,cluster_charge);
			path.append_denovo_results(denovo_out,config,c_name);
		}

			
		// results to screen
		int q;
		for (q=0; q<final_seqs.size(); q++)
			final_seqs[q].print_short(config);

		// results to files
		if (final_seqs.size()>0)
		{
			// extract scan number
		
			char msb_line[1024];
			sprintf(msb_line,"%.2f %d ",cs.get_original_pm_with_19()-1,c);

			int k;
			float max_score=0;
			for (k=0; k<final_seqs.size(); k++)
				if (final_seqs[k].expected_score>max_score)
					max_score = final_seqs[k].expected_score;

			if (max_score == 0)
				continue;

			float total_inten = cs.get_total_intensity();
			char intensity_line[32];
			sprintf(intensity_line," %f %f ",total_inten,percent);
			strcat(msb_line,intensity_line);

			char seq_score[10];
			sprintf(seq_score," %.1f ",max_score);
			strcat(msb_line,seq_score);

			int total_cand_length=0;
			for (k=0; k<final_seqs.size(); k++)
			{
				char cand_seq[64];
				
				final_seqs[k].write_seq_to_string(config,cand_seq);
				strcat(msb_line,cand_seq);
				total_cand_length+=strlen(cand_seq);
			}

			if (! max_size_bytes)
			{
				write2file(msb_out,msb_line);
			}
			else
			{
				// add to list 
				strcat(msb_line,"\n");
				char *seq_mem=new char[strlen(msb_line)+1];
				MSB_line_ptr mlp;
				mlp.max_score=max_score;
				mlp.seq = seq_mem;
				mlp.length = total_cand_length;

				strcpy(seq_mem,msb_line);

				all_res.push_back(mlp);
			}
		}
	}
	
	if (max_size_bytes>0)
	{
		FILE *stream;
		int out_size=0;
		int i;

		sort(all_res.begin(),all_res.end());

		stream=fopen(msb_out,"w");

		if (! stream)
		{
			printf("Error: couldn't open %s for writing!\n",msb_out);
			exit(1);
		}

		for (i=0; i<all_res.size(); i++)
		{
			if (out_size + all_res[i].length> max_size_bytes-1)
				break;

			fprintf(stream,"%s",all_res[i].seq);
			out_size += all_res[i].length;
		}
		fclose(stream);
	} 
}




/********************************************************************
// returns a vector of alternate sequences along with their probabilities
// uses same mechanism as tag generation.
// Limits program to the first 50 alternate paths..
*********************************************************************/
void PrmGraph::get_alternate_sub_paths(BasicSpectrum *spec, ScoreModel *model, float pm_with_19,
									   int n_idx, int n_aa, int c_idx, int c_aa,
									   vector<DeNovoPath>& all_subs) const
{
	vector<DeNovoPath> tags;
	all_subs.clear();
	tags.clear();
	const int max_num_subs = 20;
	const int max_sub_length = 5;
	const int charge = spec->get_charge();

	DeNovoPath tag_seed;

	tag_seed.score=0;
	tag_seed.seq_length=0;
	tag_seed.breakage_node_idxs[0]=n_idx;
	tag_seed.n_node_idx=n_idx;
	tag_seed.c_node_idx=n_idx;
	tag_seed.n_mass = nodes[n_idx].mass;
	tag_seed.c_mass = nodes[n_idx].mass;

	tags.clear();
	tags.push_back(tag_seed);

	bool stop_search = false;
	while (tags.size())
	{
		int j;
		
		DeNovoPath t;

		if (stop_search)
			break;

		t=tags[tags.size()-1];
		tags.pop_back();

		// try to expand the path
		for (j=0; j<nodes[t.c_node_idx].out_edges.size(); j++)
		{
			const int edge_idx=nodes[t.c_node_idx].out_edges[j];
			const Edge &e = edges[edge_idx];
			
		
			if (e.c_node_idx <= c_idx)
			{
				const float mass_diff = nodes[c_idx].mass - e.c_mass;
				if (mass_diff<90 && is_allowed_prefix_mass(mass_diff))
					continue;

				DeNovoPath new_tag;
				int prev_aa= (t.seq_length>0) ? t.seq[t.seq_length-1] : Gap;
				int next_aa=e.aa[0];

				new_tag=t;
				new_tag.edges.push_back(e);
				new_tag.seq[new_tag.seq_length++]=e.aa[0];
				new_tag.breakage_node_idxs[new_tag.seq_length]=e.c_node_idx;

				if (e.n_aa==2)
				{
					new_tag.edges.push_back(e);
					new_tag.seq[new_tag.seq_length++]=e.aa[1];
					new_tag.breakage_node_idxs[new_tag.seq_length]=e.c_node_idx;

					// make the middle breakage idx the same as the first to signal
					// that this was a double edge, with no node for it in the graph
					new_tag.breakage_node_idxs[new_tag.seq_length-1]=
						new_tag.breakage_node_idxs[new_tag.seq_length-2];
				}

				// add score for previous vertex 
				new_tag.c_node_idx = e.c_node_idx;
				new_tag.score += nodes[t.c_node_idx].get_combo_score(prev_aa,next_aa);
				new_tag.score += e.add_score;
				new_tag.breakage_node_idxs[new_tag.seq_length]=e.c_node_idx;
				new_tag.c_mass = e.c_mass;
		
				// check if tag reaches the right node
				if (new_tag.c_node_idx == c_idx)
				{
					int last_aa=new_tag.seq[new_tag.seq_length-1];
					new_tag.score += nodes[c_idx].get_combo_score(last_aa,Gap);
					
					all_subs.push_back(new_tag);
					if (all_subs.size() >= max_num_subs)
					{
						stop_search = true;
						break;
					}
					
					continue;
						
				}

				// tag is not long enough.. push it back into stack

				if (new_tag.seq_length< max_sub_length)
				{
				//	new_tag.print(config);
					tags.push_back(new_tag);
				}
			}
		}
	}
	sort(all_subs.begin(),all_subs.end());

	// add additional info to final tags
	int i;
	for (i=0; i<all_subs.size(); i++)
	{
		DeNovoPath & tag = all_subs[i];
		int j;

		for (j=0; j<=tag.seq_length; j++)
			tag.breakage_masses[j] = nodes[tag.breakage_node_idxs[j]].breakage.breakage_mass;

		tag.n_mass=tag.breakage_masses[0];
		tag.c_mass=tag.breakage_masses[tag.seq_length];
	}

	// calc average prob of amino acids in each sub path
	for (i=0; i<all_subs.size(); i++)
	{
		int j;
		vector<ME_Regression_Sample> aa_sams;

		if (n_aa> 100 || nodes[n_idx].in_aa_idxs[n_aa]<0)
			n_aa = Gap;

		if (c_aa> 100 || nodes[c_idx].out_aa_idxs[c_aa]<0)
			c_aa = Gap;

		create_aa_samples_from_path(this,spec,model,pm_with_19,all_subs[i],aa_sams,1,
			n_aa,c_aa);

		// assume aa's belong to tag in the second order of ranking
		// and belong to a tag of length 3..
		float prod = 1.0;
		for (j=0; j<all_subs[i].seq_length; j++)
		{
		//	all_subs[i].aa_probs[j]= model->local_aa_models[3]->models[1].p_y_given_x(0,aa_sams[j]);
			all_subs[i].aa_probs[j]= model->denovo_aa[charge]->me.p_y_given_x(0,aa_sams[j]);
			prod *= all_subs[i].aa_probs[j];
		}

		all_subs[i].seq_prob = pow(prod,1.0/all_subs[i].seq_length);

	//	model->get_config()->print_ints(all_subs[i].seq,all_subs[i].seq_length);
	//	printf(" %.3f\n");
	}
	sort(all_subs.begin(), all_subs.end(), comp_DeNovoPathProb);
}


/********************************************************************
// returns a single tag corresponding to the desires sequence
// if the sequence was not found in the spectrum graph, returns
// prob -1 and seq_length 0
********************************************************************/
DeNovoPath PrmGraph::create_tag_from_seq(BasicSpectrum *spec, ScoreModel *model, float pm_with_19,
						int *seq, int seq_length, int n_idx, int n_aa, int c_aa) const
{
	DeNovoPath ret;
	const int charge = spec->get_charge();

	ret.seq_prob=-1;
	ret.seq_length=0;
	ret.n_node_idx = n_idx;
	ret.c_node_idx = n_idx;
	ret.breakage_node_idxs[0] = n_idx;
	ret.n_mass = nodes[n_idx].mass;

	int i;
	for (i=0; i<seq_length; i++)
	{
		int j;
		bool added_edge=false;
		
		// try to expand the path
		for (j=0; j<nodes[ret.c_node_idx].out_edges.size(); j++)
		{
			const int edge_idx=nodes[ret.c_node_idx].out_edges[j];
			const Edge &e = edges[edge_idx];
			

			if (e.n_aa == 1)
			{
				if (e.aa[0] == seq[i] || (e.aa[0] == Leu && seq[i] == Ile) )
				{
					ret.edges.push_back(e);
					ret.seq[ret.seq_length++]=e.aa[0];
					ret.breakage_node_idxs[ret.seq_length]=e.c_node_idx;
					ret.c_node_idx=e.c_node_idx;
					added_edge=true;
					break;
				}
				else
					continue;
			}
			else if (i<seq_length -1)
			{
				if ( (e.aa[0] == seq[i] || (e.aa[0] == Leu && seq[i] == Ile)) &&
					 (e.aa[1] == seq[i+1] || (e.aa[1] == Leu && seq[i+1] == Ile) ) )
				{
					ret.edges.push_back(e);
					ret.edges.push_back(e);
					ret.seq[ret.seq_length++]=e.aa[0];
					ret.breakage_node_idxs[ret.seq_length]=e.n_node_idx;
					ret.seq[ret.seq_length++]=e.aa[1];
					ret.breakage_node_idxs[ret.seq_length]=e.c_node_idx;
					ret.c_node_idx=e.c_node_idx;
					i++;
					added_edge=true;
					break;
				}
			}		
		}
		if (! added_edge)
			break;
	}

	if (ret.seq_length<seq_length)
	{
		ret.seq_length=0;
		ret.seq_prob=-1;
		return ret;
	}

	for (i=0; i<=ret.seq_length; i++)
		ret.breakage_masses[i] = nodes[ret.breakage_node_idxs[i]].breakage.breakage_mass;

	ret.n_mass=ret.breakage_masses[0];
	ret.c_mass=ret.breakage_masses[ret.seq_length];
	ret.n_node_idx=ret.breakage_node_idxs[0];
	ret.c_node_idx=ret.breakage_node_idxs[ret.seq_length];
	

	// this means we got connected to a bad node (probably misaligned)
	if (c_aa != Gap && nodes[ret.c_node_idx].out_aa_idxs[c_aa]<0)
	{
		ret.seq_length=0;
		ret.seq_prob=-1;
		return ret;
	}

	vector<ME_Regression_Sample> aa_sams;

	create_aa_samples_from_path(this,spec,model,pm_with_19,ret,aa_sams,1,
			n_aa,c_aa);

	float prod = 1.0;
	for (i=0; i<ret.seq_length; i++)
	{
	//	all_subs[i].aa_probs[j]= model->local_aa_models[3]->models[1].p_y_given_x(0,aa_sams[j]);
		ret.aa_probs[i]= model->denovo_aa[charge]->me.p_y_given_x(0,aa_sams[i]);
		prod *= ret.aa_probs[i];
	}

	ret.seq_prob = pow(prod,1.0/ret.seq_length);

	return ret;
}


/****************************************************************************
// swaps in a tag instead of a portion of the path,
// updates all relevant info (breakages, masses, edges, probs, etc.)
*****************************************************************************/
bool DeNovoPath::swap_in_tag(int path_pos, int swap_length, 
							 const DeNovoPath& path, const DeNovoPath& tag)
{
	int i;

	seq_length = path.seq_length + tag.seq_length - swap_length;
	edges.resize(seq_length);

	n_mass = path.n_mass;
	c_mass = path.c_mass; 
	n_node_idx = path.n_node_idx;
	c_node_idx = path.c_node_idx;
	has_aa_probs = path.has_aa_probs;
	score = path.score;
	seq_prob= path.seq_prob;
	mirror_ratio = path.mirror_ratio;

	// copy every thing upto the swap position
	for (i=0; i<path_pos; i++)
	{
		seq[i]=path.seq[i];
		breakage_node_idxs[i]=path.breakage_node_idxs[i];
		breakage_masses[i]=path.breakage_masses[i];
		aa_probs[i]=path.aa_probs[i];
		edges[i]=path.edges[i];
	}

	// swap in tag stuff
	for (i=0; i<tag.seq_length; i++)
	{
		seq[path_pos+i]= tag.seq[i];
		breakage_node_idxs[path_pos+i]=tag.breakage_node_idxs[i];
		breakage_masses[path_pos+i]=tag.breakage_masses[i];
		aa_probs[path_pos+i]=tag.aa_probs[i];
		edges[path_pos+i]=tag.edges[i];
	}

	int delta_pos = tag.seq_length - swap_length;

	for (i=path_pos + swap_length; i<path.seq_length; i++)
	{
		seq[delta_pos+i]= path.seq[i];
		breakage_node_idxs[delta_pos+i]=path.breakage_node_idxs[i];
		breakage_masses[delta_pos+i]=path.breakage_masses[i];
		aa_probs[delta_pos+i]=path.aa_probs[i];
		edges[delta_pos+i]=path.edges[i];
	}

	// correct prob (not messing with score)
	float prod = 1.0;
	for (i=0; i<seq_length; i++)
		prod *= aa_probs[i];

	seq_prob = pow(prod,1.0/seq_length);

	return true;
}



/*********************************************************************
Corrects a path according to a defined set of errors (N <=> GG)
**********************************************************************/
void PrmGraph::correct_path(DeNovoPath &path, BasicSpectrum *spec, 
							ScoreModel *model, float pm_with_19) const
{
	// correction combos
	const int num_combos = 11;
	const int s_aa[num_combos]  = {Asn,Gln,Gln,Arg,Arg,Trp,Trp,Trp,Trp,Trp,Trp};
	const int d_aa[num_combos][2] = { 
		{Gly,Gly},{Ala,Gly},{Gly,Ala},{Gly,Val},{Val,Gly},
		{Asp,Ala},{Ala,Asp},{Gly,Glu},{Glu,Gly},{Val,Ser},{Ser,Val} };

	int i;
	

	// make corrections of the type N  => GG
	for (i=0; i<path.seq_length; i++)
	{
		int t;
		DeNovoPath best_swap;
		float best_prob=0;

		for (t=0; t<num_combos; t++)
			if (path.seq[i] == s_aa[t])
			{
				// don't mess around with double edges
				if (path.breakage_node_idxs[i] == path.breakage_node_idxs[i+1])
					continue;
				if (i>0 && path.breakage_node_idxs[i] == path.breakage_node_idxs[i-1])
					continue;

				int n_aa = (i>0) ? path.seq[i-1] : Gap;
				int c_aa = (i<path.seq_length-1) ? path.seq[i+1] : Gap;
				DeNovoPath alt_tag=create_tag_from_seq(spec,model,pm_with_19,(int *)d_aa[t],2, 
					path.breakage_node_idxs[i],n_aa, c_aa);

				if (alt_tag.seq_prob - path.aa_probs[i] > 0.2 && alt_tag.seq_prob > best_prob)
				{
					best_swap = alt_tag;
					best_prob = alt_tag.seq_prob;
				}
			}

		if (best_prob>0)
		{
			printf("\nBest Swap:\n");
			best_swap.print_all_info(config);

			printf("\nBEFORE: ");
			path.print_all_info(config);
			
			DeNovoPath new_path;
			new_path.swap_in_tag(i,1,path,best_swap);
			path= new_path;

			printf("\nAFTER:  ");
			path.print_all_info(config);

		}

	}

	// make corrections of the type GG => N
	for (i=0; i<path.seq_length-1; i++)
	{
		int t;
		DeNovoPath best_swap;
		float best_prob=0;
		float path_d_aa_prob = sqrt(path.aa_probs[i]*path.aa_probs[i+1]);

		for (t=0; t<num_combos; t++)
			if (path.seq[i] == d_aa[t][0] && path.seq[i+1] == d_aa[t][1])
			{

				// don't mess with the middle of double edges
				if (i>0 && path.breakage_node_idxs[i] == path.breakage_node_idxs[i-1])
					continue;
				if (path.breakage_node_idxs[i+1] == path.breakage_node_idxs[i+2])
					continue;

				int n_aa = (i>0) ? path.seq[i-1] : Gap;
				int c_aa = (i<path.seq_length-2) ? path.seq[i+2] : Gap;

				Edge e0=path.edges[0] , e1=path.edges[1], e2=path.edges[2], 
					e3 = path.edges[3] ,e4 = path.edges[4];

				DeNovoPath alt_tag=create_tag_from_seq(spec,model,pm_with_19,(int *)(&s_aa[t]),1,
					path.breakage_node_idxs[i],n_aa, c_aa);
				
				if (alt_tag.seq_prob - path_d_aa_prob > 0.1 &&
					alt_tag.seq_prob > best_prob)
				{
					best_swap = alt_tag;
					best_prob = alt_tag.seq_prob;
				}
			}

		if (best_prob>0)
		{
			printf("\nBEFORE: ");
			path.print_all_info(config);

			DeNovoPath new_path;
			new_path.swap_in_tag(i,2,path,best_swap);
			path=new_path;

			printf("\nAFTER:  ");
			path.print_all_info(config);
		}
	}
}

/************************************************************************
// returns the best value of a match between the two sequences. Assumes
// that prior to the true seq there is a cleavge aa,
// all other aa's before or after will cause a mismatch
*************************************************************************/
float get_msblast_score(int *true_seq, int true_seq_length,
						int *query_seq, int query_seq_length,
						bool print, Config *config)
{
	int i;

	int seq[64];
	int start=10;

	float best_match_score = NEG_INF;
	int best_pos = -1;

	const int miss_aa = -1;

	for (i=0; i<64; i++)
		seq[i]= miss_aa;

	seq[start-1]=B_SYM;

	for (i=0; i<true_seq_length; i++)
		seq[i+start]=true_seq[i];


	const int first_pos = start - 5;
	const int last_pos  = start + true_seq_length - 3;

	// find best match
	for (i=first_pos; i<= last_pos; i++)
	{
		float match_score = 0;
		int j;

		for (j=0; j<query_seq_length; j++)
		{
			const int query_aa = query_seq[j];
			const int true_aa = seq[i+j];

			if (query_aa == true_aa)
			{
				match_score++;
				continue;
			}
			
			// check I/L
			if ((query_aa == Ile || query_aa == Leu) &&
					(true_aa == Ile || true_aa == Leu) )
			{
				match_score++;
				continue;
			}
				
			// check Z/ QK
			if (query_aa == Z_SYM && (true_aa == Lys || true_aa == Gln))
			{
				match_score += 0.5;
				continue;
			}

			// check B / KR
			if (query_aa == B_SYM && (true_aa == Lys || true_aa == Arg))
			{
				match_score += 0.5;
				continue;
			}

			// check X
			if (query_aa == X_SYM)
			{
				continue;
			}

			// else there was a mismatch
			match_score--;
		}

		if (match_score>best_match_score)
		{
			best_match_score = match_score;
			best_pos = i;
		}
	}
	
	if (print && config)
	{
		printf("Match = %.1f\n",best_match_score);
		printf("TRUE SEQ: ");
		for (i=0; i<query_seq_length; i++)
			print_msb_sym(seq[best_pos+i],config);
		printf("\n");
		printf("MSB  SEQ: ");
		for (i=0; i<query_seq_length; i++)
			print_msb_sym(query_seq[i],config);
		printf("\n\n");

	}

	return best_match_score;
}

float get_max_msblast_score(const vector<MSB_seq>& candidates, 
							int *true_seq, int true_seq_length)
{
	float best_score = NEG_INF;
	int i;

	for (i=0; i<candidates.size(); i++)
	{
		float cand_score = get_msblast_score(true_seq,true_seq_length,
			(int *)candidates[i].seq,candidates[i].seq_length);

		if (cand_score>best_score)
			best_score = cand_score;
	}
	return best_score;
}


/*************************************************************
// generates all possible modifications to the current seq
**************************************************************/
void MSB_candidate::generate_all_mods(PrmGraph *prm, BasicSpectrum *spec, 
		ScoreModel *model,float pm_with_19, int core_size)
{
	Config *config = model->get_config();

//	printf(">>> ");

	// check if we reach terminal end
	bool reach_c_term = (current_seq.breakage_idxs[current_seq.seq_length] == prm->get_num_nodes() -1);

	generate_double_single(config);

//	printf(" A");

	generate_alternate_single_cands(config,reach_c_term);

//	printf(" B");

	generate_alt_sub_path(prm,spec,model,pm_with_19);

//	printf(" C");

	treat_low_prob_holes(config);

//	printf(" D");

	replace_with_gap(config);

//	printf(" E");
	
	alter_prefix(config);

//	printf(" F");

	trim_ends(config,core_size);

//	printf(" G\n");
}


/************************************************************************
// Generates the redundant candidates from the de novo sequence
// Each possible redundancy is given an "importance score"
*************************************************************************/
void generate_msblast_candidates(PrmGraph *prm, ScoreModel *model,
			BasicSpectrum *spec, float pm_with_19, const DeNovoPath& path,
			vector<MSB_seq>&final_seqs, int num_candidates, int core_size,
			float min_cand_score)
{
	int i;
	Config *config = prm->get_config();
	vector<MSB_candidate> candidates;
	candidates.clear();

	// create initial candidate

	MSB_seq denovo_seq;

	// add B symbol
	if (path.n_mass<1.0)
	{
		denovo_seq.seq[0]=B_SYM;
		denovo_seq.breakage_idxs[0]=-1;
		denovo_seq.aa_probs[0] = 0.5; // equal prob for K and R
		denovo_seq.marked_aa[0] = 1;

		for (i=0; i<=path.seq_length; i++)
		{
			denovo_seq.seq[i+1]=path.seq[i];
			denovo_seq.breakage_idxs[i+1]=path.breakage_node_idxs[i];
			denovo_seq.aa_probs[i+1]=path.aa_probs[i];
		}
		denovo_seq.seq_length = path.seq_length+1;
	}
	else
	{
		for (i=0; i<=path.seq_length; i++)
		{
			denovo_seq.seq[i]=path.seq[i];
			denovo_seq.breakage_idxs[i]=path.breakage_node_idxs[i];
			denovo_seq.aa_probs[i]=path.aa_probs[i];
		}
		denovo_seq.seq_length = path.seq_length;
	}

	for (i=0; i<denovo_seq.seq_length; i++)
		denovo_seq.marked_aa[i]=0;

	denovo_seq.calc_expected_score(config);
//	printf("Exp score = %.3f\n",denovo_seq.expected_score);

	MSB_candidate denovo_cand(denovo_seq);

	denovo_cand.generate_all_mods(prm,spec,model,pm_with_19,core_size);

	

	candidates.push_back(denovo_cand);

	
	
	while (candidates.size() < (int)(num_candidates * 1.25 + 5))
	{
		int i;

		// find best mod
		float best_mod_delta = NEG_INF;
		int best_mod_candidate = -1;
		for (i=0; i<candidates.size(); i++)
		{
			if (candidates[i].candidate_mods.size()>0 &&
				candidates[i].candidate_mods[0].avg_delta>best_mod_delta)
			{
				best_mod_delta = candidates[i].candidate_mods[0].avg_delta;
				best_mod_candidate = i;
			}
		}

		if (best_mod_delta == NEG_INF) 
			break;

		MSB_alternates alts = candidates[best_mod_candidate].candidate_mods[0];

		pop_heap(candidates[best_mod_candidate].candidate_mods.begin(),
				 candidates[best_mod_candidate].candidate_mods.end());
		candidates[best_mod_candidate].candidate_mods.pop_back();
		
		// add the candidates from the mods
		for (i=0; i<alts.alt_seqs.size(); i++)
		{
			
			if (alts.alt_seqs[i].seq_length>MAX_PEPTIDE_SIZE)
				continue;

			// check that the sequence is not already present
			int j;
			bool present=false;
			
			for (j=0; j<candidates.size(); j++)
			{
				
				if (candidates[j].current_seq.seq_length ==
					alts.alt_seqs[i].seq_length)
				{
					int k;
					bool same=true;

					for (k=0; k<candidates[j].current_seq.seq_length; k++)
					{
						if (candidates[j].current_seq.seq[k] != alts.alt_seqs[i].seq[k])
						{
							same=false;
							break;
						}
					}

					// check to see which has a higher expected score
					if (same)
					{
						present = true;
						if (alts.alt_seqs[i].expected_score > candidates[j].current_seq.expected_score)
						{
							// replace previous with new one
							MSB_candidate new_cand(alts.alt_seqs[i]);
							new_cand.generate_all_mods(prm,spec,model,pm_with_19,core_size);
							candidates[j]=new_cand;
						}
						break;
					}
				}
			}

			if (present)   // the mod seq is already here as one of the previous candidates
				continue; 

			// make new candidate and push it back
			MSB_candidate new_cand(alts.alt_seqs[i]);


			new_cand.generate_all_mods(prm,spec,model,pm_with_19,core_size);
		//	printf("mods:\n");
		//	int k;
		//	for (k=0; k<new_cand.candidate_mods.size(); k++)
		//		new_cand.candidate_mods[k].print(config);

			candidates.push_back(new_cand);
			
		}

	//	break;
	}

//	printf("Done pushing...\n");
	// extract candidates according to order of expected score
	// later on might also look at similarity to other candidates

	int num_5=0; 
	int num_6=0;
	while (1)
	{
		float max_exp_score = -1000;
		int max_cand_idx=-1;
		for (i=0; i<candidates.size(); i++)
		{
			if (candidates[i].current_seq.expected_score > max_exp_score)
			{
				max_exp_score = candidates[i].current_seq.expected_score;
				max_cand_idx = i;
			}
		}


		if (max_cand_idx<0)
			break;

		// check that the seq is good
		if (candidates[max_cand_idx].current_seq.expected_score< min_cand_score || 
			(candidates[max_cand_idx].current_seq.seq_length==5 && num_5++ == 1 ) ||
			(candidates[max_cand_idx].current_seq.seq_length==6 && num_6++ == 3 ) )
		{
			candidates[max_cand_idx].current_seq.expected_score=NEG_INF;
			continue;
		}

		// check that we do not have an equivalent sequence in there already
		bool equiv=false;
		for (i=0; i<final_seqs.size(); i++)
		{
			if (final_seqs[i].equivalent(candidates[max_cand_idx].current_seq))
			{
				candidates[max_cand_idx].current_seq.expected_score=NEG_INF;
				equiv=true;
				break;
			}
		}
		if (equiv)
			continue;

//		printf("CAND %d: ",final_seqs.size()+1);
//		candidates[max_cand_idx].current_seq.print_short(config);
		final_seqs.push_back(candidates[max_cand_idx].current_seq);
		candidates[max_cand_idx].current_seq.expected_score=NEG_INF;

		if (final_seqs.size() == num_candidates)
			break;
	}
}


/**********************************************************
Calcs the expected score for an MSB sequence. Bases calculation
on the probabilities of the predicted amino acids, the 
presence of marked aa's and the presense of problematic
combinations of amino acids. Uses lots of empirically 
derived thresholds.
***********************************************************/
void MSB_seq::calc_expected_score(Config *config)
{
	// Emprical probabilities 
	const float prob_R = 0.75;
	const float prob_W = 0.65;
	const float prob_N = 0.96;
	const float prob_Q = 0.96;
	const float prob_GV = 0.8;
	const float prob_VS = 0.95;
	const float prob_SV = 0.95;
	const float prob_DA = 0.9;
	const float prob_AD = 0.95;
	const float prob_AG = 0.96;
                                        
	const float prob_F_vs_Met_ox = 0.9; // probs for single aa mixup
	const float prob_K_first = 0.45;
	const float prob_K_mid   = 0.2;
	const float prob_K_last  = 1.0;  

	const int ox_met_aa = config->get_oxidized_met_aa();

	int i;

	float prefix_scores[MAX_PEPTIDE_SIZE];
	float suffix_scores[MAX_PEPTIDE_SIZE];
	
	if (seq[0] == X_SYM)
	{
		prefix_scores[0]=0;
	}
	else if (seq[0] == Z_SYM)
	{
		prefix_scores[0] = aa_probs[0]*0.5;
	}
	else
		prefix_scores[0] =( (marked_aa[0]) ? aa_probs[0] : (aa_probs[0] - (1-aa_probs[0])) );

	for (i=1; i<seq_length; i++)
	{
		prefix_scores[i] = prefix_scores[i-1];
		if (seq[i] == X_SYM)
			continue;

		if (seq[i]==Z_SYM)
		{
			prefix_scores[i]+=aa_probs[i]*0.5;
		}
		else
			prefix_scores[i] += ((marked_aa[i]) ? aa_probs[i] : (aa_probs[i] - (1-aa_probs[i])));
	}

	suffix_scores[seq_length]=0;
	for (i=seq_length-1; i>=0; i--)
	{
		suffix_scores[i] = suffix_scores[i+1];
		if (seq[i]==X_SYM)
			continue;

		if (seq[i]==Z_SYM)
		{
			suffix_scores[i] += aa_probs[i] * 0.5;
		}
		else
			suffix_scores[i] += ((marked_aa[i]) ? aa_probs[i] : (aa_probs[i] - (1-aa_probs[i])));
	}

	

	float simple_match_score = prefix_scores[seq_length-1];

	// discount occurences of substituteable single amino acids
	// F <=> M* ,  K <=>Q
	// the score adjustment works as follows:
	// with prob 1-k we have the wrong amino acid, we first need to substract
	// the positive score given (aa_probs[i]), and then add the mismatch
	// penalty (+1).
	if (seq[0] == B_SYM && seq[1] ==Lys)
	{
		simple_match_score -= ((1-prob_K_first)*(aa_probs[1]+1));
	}
	else if (seq[0] == B_SYM && seq[1] ==Gln)
	{
		simple_match_score -= (prob_K_first*(aa_probs[1]+1));
	}

	int start = (seq[0] == B_SYM) ? 2 : 0;
	for (i=start; i<seq_length-1; i++)
	{
		if (marked_aa[i])
				continue;

		if (seq[i] == Lys)
		{
			simple_match_score -= ((1-prob_K_mid)*(aa_probs[i]+1));
		}
		else if (seq[i] == Gln)
		{
			simple_match_score -= (prob_K_mid*(aa_probs[i]+1));
		}
	}

	if (seq[seq_length-1] == Lys)
	{
		simple_match_score -= ((1-prob_K_last)*(aa_probs[seq_length-1]+1));
	}
	else if (seq[seq_length-1] == Gln)
	{
		simple_match_score -= (prob_K_last*(aa_probs[seq_length-1]+1));
	}


	if (ox_met_aa >0)
	{
		for (i=0; i<seq_length; i++)
		{
			if (marked_aa[i])
				continue;

			if (seq[i] == ox_met_aa)
			{
				simple_match_score -= ((prob_F_vs_Met_ox)*(aa_probs[i]+1));
			}
			else if (seq[i] == Phe )
			{
				simple_match_score -= ((1-prob_F_vs_Met_ox)*(aa_probs[i]+1));
			}
		}
	}
	
	float min_score = simple_match_score;

	// look for the worse error that can occur where a double <=> single aa
	// error occurs. Only consider cases where there are unmarked amino acids
	// involved.

	

	// first consider single amino acids that could be replaced by doubles
	// look only in the center, because if it appears on the edge, it can do
	// much damage.
	for (i=1; i< seq_length-1; i++)
	{
		if (marked_aa[i] || seq[i] == X_SYM)
			continue;
		
		float prob = -1;

		switch (seq[i])
		{
			case Arg : prob = prob_R; break;
			case Asn : prob = prob_N; break;
			case Gln : prob = prob_Q; break;
			case Trp : prob = prob_W; break;
		}
		
		if (prob <0)
			continue;

		float pre_match = prefix_scores[i-1] - seq_length+ i;
		float suf_match = suffix_scores[i+1] - i - 1;
		float mismatch_score = (pre_match>suf_match ) ? pre_match : suf_match;
		float exp_score = prob * simple_match_score + (1-prob)*mismatch_score;
		if (exp_score<min_score)
			min_score=exp_score;


	//	printf("%d : pre score %.2f  suf_score: %.2f\n",i,pre_match,suf_match);
	//	printf("mismatach score: %.2f  prob: %.2f   exp_score %.2f\n",
	//		mismatch_score,prob,exp_score);
	}

	// look for double combos that appear within the peptide
	for (i=1; i<seq_length-2; i++)
	{
		if (marked_aa[i] || seq[i]== X_SYM)
			continue;

		float prob=-1;

		if (seq[i]==Gly && seq[i+1] == Val)
		{
			prob = prob_GV;
		}
		else if (seq[i]==Val && seq[i+1]==Ser)
		{
			prob = prob_VS;
		}
		else if (seq[i]==Ser && seq[i+1]==Val)
		{
			prob = prob_SV;
		}
		else if (seq[i]==Asp && seq[i+1]==Ala)
		{
			prob = prob_DA;
		}
		else if (seq[i]==Ala && seq[i+1]==Asp)
		{
			prob = prob_AD;
		}
		else if (seq[i]==Ala && seq[i+1]==Gly)
		{
			prob = prob_AG;
		}
		else 
			continue;

		float pre_match = prefix_scores[i-1] - seq_length+ i;
		float suf_match = suffix_scores[i+1] - i - 1;
		float mismatch_score = (pre_match>suf_match ) ? pre_match : suf_match;
		float exp_score = prob * simple_match_score + (1-prob)*mismatch_score;
		if (exp_score<min_score)
			min_score=exp_score;


	//	printf("%d : pre score %.2f  suf_score: %.2f\n",i,pre_match,suf_match);
	//	printf("mismatach score: %.2f  prob: %.2f   exp_score %.2f\n",
	//		mismatch_score,prob,exp_score);
	}
	expected_score = min_score;
}


/***********************************************************************
// calulates the scores of all the alternates, computes the delta scores
// relative to the original score
************************************************************************/
void MSB_alternates::calaculate_scores(const MSB_seq& org_seq, Config *config)
{
	int i;

	total_delta_score= 0;

	for (i=0; i<alt_seqs.size(); i++)
	{
		alt_seqs[i].calc_expected_score(config);
		total_delta_score += alt_seqs[i].expected_score - org_seq.expected_score;
	}

	avg_delta = total_delta_score / alt_seqs.size();

	num_alternates= alt_seqs.size();
}




/***************************************************************************
// copies the original MSB_seq, and replaces the portion from the original
// with the new amino acids
// replaced new portions get the average prob of the old portion
****************************************************************************/
void MSB_seq::clone_and_replace(const MSB_seq& org, int org_pos, int org_segment_length,
						   int *new_aas, int new_segment_length)
{
	int i;
	int delta_length = new_segment_length - org_segment_length;
	float prob = 1;

	for (i=0; i<org_segment_length; i++)
		prob *= org.aa_probs[org_pos+i];

	prob=pow(prob,1.0/org_segment_length); // average prob of swapped out aa's
	// if swapped in portion is longer than swapped out, decrease probs accordingly
	if (new_segment_length > org_segment_length)
	{
		prob = (prob * org_segment_length) / new_segment_length;
	}
	

	int double_edge_ind[MAX_PEPTIDE_SIZE];

	for (i=0; i<org.seq_length; i++)
		double_edge_ind[i]=0;

	for (i=1; i<org.seq_length; i++)
		if (org.breakage_idxs[i]==org.breakage_idxs[i-1])
			double_edge_ind[i]=1;

	*this=org;

	// update
	for (i=0; i<new_segment_length; i++)
	{
		seq[i+org_pos]=new_aas[i];
		aa_probs[i+org_pos]=prob;
		
		if (new_aas[i] == X_SYM)        // give fixed low prob for X
			aa_probs[i+org_pos] = 0.05;

		marked_aa[i+org_pos]=0;
		breakage_idxs[i+org_pos]=-1; // if breakage idx info is needed, must use different
									 // function that swaps from a path or tag
	}
	breakage_idxs[org_pos]=org.breakage_idxs[org_pos]; // the first brekage is still correct

	// add right portion

	for (i=org_pos+org_segment_length; i<=org.seq_length; i++)
	{
		seq[i+delta_length]=org.seq[i];
		aa_probs[i+delta_length]=org.aa_probs[i];
		marked_aa[i+delta_length]=org.marked_aa[i];
		breakage_idxs[i+delta_length]=org.breakage_idxs[i];
	}

	// invalidate the breakage idx if the swap was performed in the
	// middle of a double edge - this will probably cause problems later
	if (double_edge_ind[org_pos+ org_segment_length])
		breakage_idxs[org_pos + new_segment_length]=-1;

	seq_length = org.seq_length + delta_length;
}

/*************************************************************************
// same as clone and replace, only also copies in breakage_idxs and probs
**************************************************************************/
void MSB_seq::clone_and_swap_in(const MSB_seq& org, int org_pos, int org_segment_length,
						   const DeNovoPath& tag)
{
	int i;
	const int delta_length = tag.seq_length - org_segment_length;
	*this=org;

	// update
	for (i=0; i<tag.seq_length; i++)
	{
		seq[i+org_pos]=tag.seq[i];
		aa_probs[i+org_pos]=tag.aa_probs[i];
		
		marked_aa[i+org_pos]=0;
		breakage_idxs[i+org_pos]= tag.breakage_node_idxs[i];
	}


	// add right portion

	for (i=org_pos+org_segment_length; i<=org.seq_length; i++)
	{
		seq[i+delta_length]=org.seq[i];
		aa_probs[i+delta_length]=org.aa_probs[i];
		marked_aa[i+delta_length]=org.marked_aa[i];
		breakage_idxs[i+delta_length]=org.breakage_idxs[i];
	}
	seq_length = org.seq_length + delta_length;
}



// writes the sequence in MSB format
void MSB_seq::write_seq_to_string(Config *config, char *buff) const
{
	int i;
	strcpy(buff,"-");

	for (i=0; i<seq_length; i++)
	{
		if (seq[i]<100)
		{
			strcat(buff,config->get_aa_label(seq[i]));
		}
		else
		{
			if (seq[i] == X_SYM)
			{
				strcat(buff,"X");
			}
			else if (seq[i] == Z_SYM)
			{
				strcat(buff,"Z");
			}
			else if (seq[i] == B_SYM)
			{
				strcat(buff,"B");
			}
			else
			{
				printf("Error: unknown AA code: %d !\n",seq[i]);
				exit(1);
			}
		}
	}

	// replace all '*' with '.'
	for (i=0; i<strlen(buff); i++)
		if (buff[i]=='*')
			buff[i]='.';

}


void MSB_seq::print(Config *config) const
{
	int i;

	for (i=0; i<seq_length; i++)
		print_msb_sym(seq[i],config);

	printf("  %.2f\n",expected_score);
	
	for (i=0; i<seq_length; i++)
	{
		if (marked_aa[i])
		{
			printf("*");
		}
		else
			printf(".");

		if (seq[i]<100 && strlen(config->get_aa_label(seq[i]))>1)
			printf(" ");
	}
	printf("\n");

	printf("AA_PROBS:");
	for (i=0; i<seq_length; i++)
		printf(" %.2f",aa_probs[i]);
	printf("\n");

	printf("AA_BREAK_IDXS:");
	for (i=0; i<seq_length; i++)
		printf(" %d",breakage_idxs[i]);
	printf("\n");
}


void MSB_seq::print_short(Config *config) const
{
	int i;

	for (i=0; i<seq_length; i++)
		print_msb_sym(seq[i],config);
	printf("  %.2f\n",expected_score);
}


void MSB_alternates::print(Config *config) const
{
	int i;

	printf("AVG_DELTA: %.2f  NUM_CANDS: %d\n",avg_delta,num_alternates);
	for (i=0; i<alt_seqs.size(); i++)
	{
		printf("%d.\n",i+1);
		alt_seqs[i].print(config);
		printf("\n");
	}
}

void MSB_candidate::pop_and_print_all_mods(Config *config)
{

	printf("Total %d candidate mods.\n",candidate_mods.size());
	int count=0;
	while (candidate_mods.size()>0)
	{
		pop_heap(candidate_mods.begin(),candidate_mods.end());
		printf("\nMOD CAND %d:\n",++count);
		printf(  "-----------\n");
		candidate_mods[candidate_mods.size()-1].print(config);
		candidate_mods.pop_back();
	}

}



