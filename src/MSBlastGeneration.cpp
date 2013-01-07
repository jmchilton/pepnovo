#include "MSBlast.h"
#include "DeNovoDp.h"
#include "PrmGraph.h"
#include "AllScoreModels.h"


/*! \fn convert_SeqPath_to_MSB_sequence_set
	\brief converts the de novo generated SeqPath into a set of MSB_query sequences

	The conversion includes creating gaps (XXX), adding prefix and suffix amino acids (B/Z),
	performing common AA replacements (GG->N W->.., Q->DA, ...).
*/
void convert_SeqPath_to_MSB_sequence_set(const PrmGraph& prm, 
										 const SeqPath& path, 
										 vector<MSBSequence>& msb_seqs)
{

}









void create_MSB_query_for_file_list(const FileManager& fm, 
									AllScoreModels *model, 
									int max_query_size)
{
	Config *config = model->get_config();
	FileSet		fs;

	fs.select_all_files(fm);

	while (1)
	{
		Spectrum spec;
		PrmGraph prm;
		SingleSpectrumFile *ssf;
		vector<MSBSequence> msb_sequences;

		if (! fs.get_next_spectrum(fm,config,&spec,&ssf)) 
			break;

		
		model->init_model_for_scoring_spectrum(&spec);
		/*
		AllScoreModels *_model, 
										  Spectrum *spectrum,
										  mass_t _pm_with_19, 
										  int spec_charge, 
										  bool add_all_pepitde_nodes, 
										  bool only_basic_score*/
		prm.create_graph_from_spectrum(model, &spec, spec.get_org_pm_with_19());
		model->score_graph_edges(prm);

		spec.print_expected_by();
		prm.print_with_multi_edges();

//\\	generate_MSB_sequences_from_PrmGraph(prm,msb_sequences,10);
		
		exit(0);
		
	}

}


