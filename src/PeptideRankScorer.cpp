#include "PeptideRankScorer.h"
#include "BasicDataStructs.h"
#include "FragmentSelection.h"
#include "DeNovoSolutions.h"
#include "auxfun.h"


void	*PeptideRankScorer::allScoreModelsPtr_=NULL;


void PeptideRankScorer::fillCompletePeptideRankboostSample(const PeptideSolution& sol,
							   AnnotatedSpectrum& as,
							   const vector<PmcSqsChargeRes>& pmcSqsRes,
							   RankBoostSample& rbs,
							   int size_idx) const
{

	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();
	vector< vector<intensity_t> > ann_intens;
	vector< vector<mass_t> >	  ann_masses;

	Peptide org_peptide = as.getPeptide();
	as.set_peptide(sol.pep);
	as.annotate_spectrum(sol.pm_with_19, sol.charge, true);
	as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

	rbs.clear();
	const int charge=sol.charge;
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	if (size_idx<0)
		size_idx=peak_model->get_size_group(charge,sol.pm_with_19);

	DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];

	if (! part_model || ! part_model->ind_was_initialized)
	{
		cout << "Error: de novo partition model was not initialized, charge " <<
			charge << " size " << size_idx << endl;
		exit(1);
	}

//	cout << "PEP: " << sol.pep.as_string(model->get_config()) << endl;

	PrmGraph prm;
	SeqPath  sol_seq_path;
	if (part_model->use_prm_features)
	{
		prm.create_graph_for_peptide_and_spectrum(allScoreModels,&as,sol.pep.get_mass_with_19(),sol.charge,sol.pep);
		allScoreModels->score_graph_edges(prm);
		sol_seq_path = prm.get_path_from_peptide_prm_graph(sol.pep);
	}

	if (model_type == 0 || model_type == 2)
	{
		if (part_model->use_PTM_peak_features)
			part_model->fill_PTM_peak_features(allScoreModels->get_config(),sol,ann_masses,ann_intens,as,rbs);

		if (part_model->use_tryp_terminal_features)
			part_model->fill_tryp_terminal_features(sol,sol_seq_path,rbs);

		if (part_model->use_ann_peak_features)
			part_model->fill_ann_peak_features(sol,ann_masses,ann_intens,as,rbs);

		if (part_model->use_inten_balance_features)
			part_model->fill_inten_balance_features(prm.get_config(),sol,sol_seq_path,rbs);

		if (part_model->use_peak_offset_features)
			part_model->fill_peak_offset_features(as.getConfig(),sol,ann_masses,ann_intens,rbs);

		if (part_model->use_comp_features)
			part_model->fill_composition_features(sol,as.getConfig(),&compAssigner, sol_seq_path, rbs);

		if (part_model->use_pmc_features)
			part_model->fill_pmcsqs_features(sol,pmcSqsRes,allScoreModels->get_pmcsqs_ptr(),rbs);

		if (part_model->use_ppp_features && peak_model->get_feature_set_type()<=2)
			part_model->fill_peak_prediction_features(sol, ann_intens,peak_model,rbs,size_idx);

		if (part_model->use_prm_features)
			part_model->fill_prm_features(sol, sol_seq_path, model_type, rbs);
			
		if (part_model->use_combined_ppp_features && peak_model->get_feature_set_type()>2)
			part_model->fill_combined_peak_prediction_features(sol, ann_intens,peak_model,rbs,size_idx);
		
		as.set_peptide(org_peptide);

		return;
	}
	else
	{
		cout << "Error: fill_complete_peptide_rbs should only be used with model types 0 or 2, not " <<
			model_type << endl;
		exit(1);
	}
}


/*****************************************************************************************
The main function for rank scoring a complete peptide (co).
******************************************************************************************/
void PeptideRankScorer::fill_complete_peptide_rbs(
							   const PeptideSolution& sol,
							   QCPeak* peaks, 
							   int num_peaks, 
							   AnnotatedSpectrum& as,
							   const vector<PmcSqsChargeRes>& pmcsqs_res,
							   RankBoostSample& rbs,
							   int size_idx) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();
	vector< vector<intensity_t> > ann_intens;
	vector< vector<mass_t> >	  ann_masses;

	Peptide org_peptide = as.getPeptide();
	as.set_peptide(sol.pep);
	as.annotate_spectrum(sol.pm_with_19, sol.charge, true);
	as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

	rbs.clear();
	const int charge=sol.charge;
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	if (size_idx<0)
		size_idx=peak_model->get_size_group(charge,sol.pm_with_19);

	DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];

	if (! part_model || ! part_model->ind_was_initialized)
	{
		cout << "Error: de novo partition model was not initialized, charge " <<
			charge << " size " << size_idx << endl;
		exit(1);
	}

//	cout << "PEP: " << sol.pep.as_string(model->get_config()) << endl;

	PrmGraph prm;
	SeqPath  sol_seq_path;
	if (part_model->use_prm_features)
	{
		prm.create_graph_for_peptide_and_spectrum(allScoreModels,&as,sol.pep.get_mass_with_19(),sol.charge,sol.pep);
		allScoreModels->score_graph_edges(prm);
		sol_seq_path = prm.get_path_from_peptide_prm_graph(sol.pep);
	}

	if (model_type == 0 || model_type == 2)
	{
		if (part_model->use_PTM_peak_features)
			part_model->fill_PTM_peak_features(allScoreModels->get_config(),sol,ann_masses,ann_intens,as,rbs);

		if (part_model->use_tryp_terminal_features)
			part_model->fill_tryp_terminal_features(sol,sol_seq_path,rbs);

		if (part_model->use_ann_peak_features)
			part_model->fill_ann_peak_features(sol,ann_masses,ann_intens,as,rbs);

		if (part_model->use_inten_balance_features)
			part_model->fill_inten_balance_features(prm.get_config(),sol,sol_seq_path,rbs);

		if (part_model->use_peak_offset_features)
			part_model->fill_peak_offset_features(as.getConfig(),sol,ann_masses,ann_intens,rbs);

		if (part_model->use_comp_features)
			part_model->fill_composition_features(sol,as.getConfig(),&compAssigner, sol_seq_path, rbs);

		if (part_model->use_pmc_features)
			part_model->fill_pmcsqs_features(sol,pmcsqs_res,allScoreModels->get_pmcsqs_ptr(),rbs);

		if (part_model->use_ppp_features && peak_model->get_feature_set_type()<=2)
			part_model->fill_peak_prediction_features(sol, ann_intens,peak_model,rbs,size_idx);

		if (part_model->use_prm_features)
			part_model->fill_prm_features(sol, sol_seq_path, model_type, rbs);
			
		if (part_model->use_combined_ppp_features && peak_model->get_feature_set_type()>2)
			part_model->fill_combined_peak_prediction_features(sol, ann_intens,peak_model,rbs,size_idx);
		
		as.set_peptide(org_peptide);

		return;
	}
	else
	{
		cout << "Error: fill_complete_peptide_rbs should only be used with model types 0 or 2, not " <<
			model_type << endl;
		exit(1);
	}
}



/*****************************************************************************************
The main function for rank scoring a partial peptide.
******************************************************************************************/
void PeptideRankScorer::fill_denovo_peptide_rbs(
							   PeptideSolution& sol,
							   const SeqPath& path,
							   QCPeak* peaks, 
							   int num_peaks, 
							   AnnotatedSpectrum& as,
							   const vector<PmcSqsChargeRes>& pmcsqs_res,
							   RankBoostSample& main_rbs,
							   int size_idx) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();
	const PrmGraph *prm = path.prm_ptr;

	vector< vector<intensity_t> > ann_intens;
	vector< vector<mass_t> >	  ann_masses;

	as.set_peptide(sol.pep);
	as.annotate_spectrum(sol.pm_with_19, sol.charge, true);
	as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

	main_rbs.clear();
	const int charge=sol.charge;

	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	if (size_idx<0)
		size_idx=peak_model->get_size_group(charge,sol.pm_with_19);

	DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];

	if (! part_model || ! part_model->ind_was_initialized)
	{
		cout << "Error: de novo partition model was not initialized, charge " <<
			charge << " size " << size_idx << endl;
		exit(1);
	}

	

	if (model_type == 1)
	{
		if (part_model->use_PTM_peak_features)
			part_model->fill_PTM_peak_features(allScoreModels->get_config(),sol,ann_masses,ann_intens,as,main_rbs);

		if (part_model->use_tryp_terminal_features)
			part_model->fill_tryp_terminal_features(sol,path,main_rbs);

		if (part_model->use_ann_peak_features)
			part_model->fill_ann_peak_features(sol,ann_masses,ann_intens,as,main_rbs);

		if (part_model->use_inten_balance_features)
			part_model->fill_inten_balance_features(allScoreModels->get_config(),sol,path,main_rbs);

		if (part_model->use_peak_offset_features)
			part_model->fill_peak_offset_features(as.getConfig(),sol,ann_masses,ann_intens,main_rbs);

		if (part_model->use_comp_features)
			part_model->fill_composition_features(sol,as.getConfig(),&compAssigner, path, main_rbs);

		if (part_model->use_pmc_features && sol.reaches_n_terminal && sol.reaches_c_terminal)
			part_model->fill_pmcsqs_features(sol, pmcsqs_res, allScoreModels->get_pmcsqs_ptr(),main_rbs);

		if (part_model->use_prm_features)
			part_model->fill_prm_features(sol, path, model_type, main_rbs);
			
		if (part_model->use_combined_ppp_features) 
		{
			if (peak_model->get_feature_set_type() != 4)
			{
				cout << "Error: feature type for filling peak predictions should be 4 with this function,";
				cout << " not " << peak_model->get_feature_set_type() << endl;
				exit(1);
			}

			// depending on wether the peptide reaches all the way to the N- and C-terminals
			// we will examine different combos of the missing n and c basic amino acids

		
			part_model->fill_combined_peak_prediction_features(sol, ann_intens,peak_model,main_rbs,size_idx);
		}
		return;

	}
	else
	{
		cout << "Error: fill_denovo_peptide_rbs should only be used with model type 1, not " <<
			model_type << endl;
		exit(1);
	}
}


/*****************************************************************************************
The main function for rank scoring a complete peptide (co).
This function assumes that in case of partial de novo sequences, we don't know the most 
basic amino acids on the N- and C-terminal sides of the predicted (partial) peptide so
we store features for the different posssibilities. The scoring function will choose
the combo that had the highes score. The main set of features which is good for all
combos (these features don't depend on the aa combo), are stored separately in main_rbs.
******************************************************************************************/
void PeptideRankScorer::fillDenovoPeptideRankboostSampleWithCombos(PeptideSolution& sol,
								 const SeqPath& path,
								 AnnotatedSpectrum& as,
								 const vector<PmcSqsChargeRes>& pmcsqs_res,
								 RankBoostSample& main_rbs,
								 vector<RankBoostSample>& peak_prediction_rbs,
							     int size_idx) const
{
	const int missing_aas[]={0,Arg,Lys};
	const int num_missing_aas = sizeof(missing_aas)/sizeof(int);
	const PrmGraph *prm = path.prm_ptr;
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();

	vector< vector<intensity_t> > ann_intens;
	vector< vector<mass_t> >	  ann_masses;

	as.set_peptide(sol.pep);
	as.annotate_spectrum(sol.pm_with_19, sol.charge, true);
	as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

	main_rbs.clear();
	const int charge=sol.charge;

	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	if (size_idx<0)
		size_idx=peak_model->get_size_group(charge,sol.pm_with_19);

	DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];

	if (! part_model || ! part_model->ind_was_initialized)
	{
		cout << "Error: de novo partition model was not initialized, charge " <<
			charge << " size " << size_idx << endl;
		exit(1);
	}

	if (model_type == 1)
	{
		if (part_model->use_PTM_peak_features)
			part_model->fill_PTM_peak_features(allScoreModels->get_config(),sol,ann_masses,ann_intens,as,main_rbs);

		if (part_model->use_tryp_terminal_features)
			part_model->fill_tryp_terminal_features(sol,path,main_rbs);

		if (part_model->use_ann_peak_features)
			part_model->fill_ann_peak_features(sol,ann_masses,ann_intens,as,main_rbs);

		if (part_model->use_inten_balance_features)
			part_model->fill_inten_balance_features(allScoreModels->get_config(),sol,path,main_rbs);

		if (part_model->use_peak_offset_features)
			part_model->fill_peak_offset_features(as.getConfig(),sol,ann_masses,ann_intens,main_rbs);

		if (part_model->use_comp_features)
			part_model->fill_composition_features(sol,as.getConfig(),&compAssigner, path, main_rbs);

		if (part_model->use_pmc_features && sol.reaches_n_terminal && sol.reaches_c_terminal)
			part_model->fill_pmcsqs_features(sol, pmcsqs_res, allScoreModels->get_pmcsqs_ptr(),main_rbs);

		if (part_model->use_prm_features)
			part_model->fill_prm_features(sol, path, model_type, main_rbs);
			
		if (part_model->use_combined_ppp_features) 
		{
			if (peak_model->get_feature_set_type() != 4)
			{
				cout << "Error: feature type for filling peak predictions should be 4 with this function,";
				cout << " not " << peak_model->get_feature_set_type() << endl;
				exit(1);
			}

			// depending on wether the peptide reaches all the way to the N- and C-terminals
			// we will examine different combos of the missing n and c basic amino acids

			vector<int> n_aas,c_aas;
			n_aas.clear();
			c_aas.clear();

			if (sol.reaches_n_terminal && sol.reaches_c_terminal)
			{
				peak_prediction_rbs.resize(1);
				peak_prediction_rbs[0].clear();
				n_aas.push_back(0);
				c_aas.push_back(0);
			}
			else if (sol.reaches_n_terminal && ! sol.reaches_c_terminal)
			{
				peak_prediction_rbs.resize(num_missing_aas);
				int i;
				for (i=0; i<num_missing_aas; i++)
				{
					peak_prediction_rbs[i].clear();
					n_aas.push_back(0);
					c_aas.push_back(missing_aas[i]);
				}
			}
			else if (! sol.reaches_n_terminal && sol.reaches_c_terminal)
			{
				peak_prediction_rbs.resize(num_missing_aas);
				int i;
				for (i=0; i<num_missing_aas; i++)
				{
					peak_prediction_rbs[i].clear();
					n_aas.push_back(missing_aas[i]);
					c_aas.push_back(0);
				}
			}
			else
			{	
				peak_prediction_rbs.resize(num_missing_aas*num_missing_aas);
				int i;
				for (i=0; i<num_missing_aas; i++)
				{
					int j;
					for (j=0; j<num_missing_aas; j++)
					{
						peak_prediction_rbs[i*num_missing_aas+j].clear();
						n_aas.push_back(missing_aas[i]);
						c_aas.push_back(missing_aas[j]);
					}
				}
			}

			int i;
			for (i=0; i<peak_prediction_rbs.size(); i++)
			{
				sol.most_basic_aa_removed_from_n=n_aas[i];
				sol.most_basic_aa_removed_from_c=c_aas[i];
				part_model->fill_combined_peak_prediction_features(sol, ann_intens,peak_model,peak_prediction_rbs[i],size_idx);
			}
		}
		return;

	}
	else
	{
		cout << "Error: fill_denovo_peptide_rbs with combos should only be used with model type 1, not " <<
			model_type << endl;
		exit(1);
	}
}


/*****************************************************************************************
The main function for rank scoring a complete peptide (co).
This function assumes that in case of partial de novo sequences, we don't know the most 
basic amino acids on the N- and C-terminal sides of the predicted (partial) peptide so
we store features for the different posssibilities. The scoring function will choose
the combo that had the highes score. The main set of features which is good for all
combos (these features don't depend on the aa combo), are stored separately in main_rbs.
******************************************************************************************/
void PeptideRankScorer::fill_denovo_peptide_rbs_with_combos(
							   PeptideSolution& sol,
							   const SeqPath& path,
							   QCPeak* peaks, 
							   int num_peaks, 
							   AnnotatedSpectrum& as,
							   const vector<PmcSqsChargeRes>& pmcsqs_res,
							   RankBoostSample& main_rbs,
							   vector<RankBoostSample>& peak_prediction_rbs,
							   int size_idx) const
{
	const int missing_aas[]={0,Arg,Lys};
	const int num_missing_aas = sizeof(missing_aas)/sizeof(int);
	const PrmGraph *prm = path.prm_ptr;

	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();

	vector< vector<intensity_t> > ann_intens;
	vector< vector<mass_t> >	  ann_masses;

	as.set_peptide(sol.pep);
	as.annotate_spectrum(sol.pm_with_19, sol.charge, true);
	as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

	main_rbs.clear();
	const int charge=sol.charge;

	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	if (size_idx<0)
		size_idx=peak_model->get_size_group(charge,sol.pm_with_19);

	DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];

	if (! part_model || ! part_model->ind_was_initialized)
	{
		cout << "Error: de novo partition model was not initialized, charge " <<
			charge << " size " << size_idx << endl;
		exit(1);
	}

	if (model_type == 1)
	{
		if (part_model->use_PTM_peak_features)
			part_model->fill_PTM_peak_features(allScoreModels->get_config(),sol,ann_masses,ann_intens,as,main_rbs);

		if (part_model->use_tryp_terminal_features)
			part_model->fill_tryp_terminal_features(sol,path,main_rbs);

		if (part_model->use_ann_peak_features)
			part_model->fill_ann_peak_features(sol,ann_masses,ann_intens,as,main_rbs);

		if (part_model->use_inten_balance_features)
			part_model->fill_inten_balance_features(allScoreModels->get_config(),sol,path,main_rbs);

		if (part_model->use_peak_offset_features)
			part_model->fill_peak_offset_features(as.getConfig(),sol,ann_masses,ann_intens,main_rbs);

		if (part_model->use_comp_features)
			part_model->fill_composition_features(sol,as.getConfig(),&compAssigner, path, main_rbs);

		if (part_model->use_pmc_features && sol.reaches_n_terminal && sol.reaches_c_terminal)
			part_model->fill_pmcsqs_features(sol, pmcsqs_res, allScoreModels->get_pmcsqs_ptr(),main_rbs);

		if (part_model->use_prm_features)
			part_model->fill_prm_features(sol, path, model_type, main_rbs);
			
		if (part_model->use_combined_ppp_features) 
		{
			if (peak_model->get_feature_set_type() != 4)
			{
				cout << "Error: feature type for filling peak predictions should be 4 with this function,";
				cout << " not " << peak_model->get_feature_set_type() << endl;
				exit(1);
			}

			// depending on wether the peptide reaches all the way to the N- and C-terminals
			// we will examine different combos of the missing n and c basic amino acids

			vector<int> n_aas,c_aas;
			n_aas.clear();
			c_aas.clear();

			if (sol.reaches_n_terminal && sol.reaches_c_terminal)
			{
				peak_prediction_rbs.resize(1);
				peak_prediction_rbs[0].clear();
				n_aas.push_back(0);
				c_aas.push_back(0);
			}
			else if (sol.reaches_n_terminal && ! sol.reaches_c_terminal)
			{
				peak_prediction_rbs.resize(num_missing_aas);
				int i;
				for (i=0; i<num_missing_aas; i++)
				{
					peak_prediction_rbs[i].clear();
					n_aas.push_back(0);
					c_aas.push_back(missing_aas[i]);
				}
			}
			else if (! sol.reaches_n_terminal && sol.reaches_c_terminal)
			{
				peak_prediction_rbs.resize(num_missing_aas);
				int i;
				for (i=0; i<num_missing_aas; i++)
				{
					peak_prediction_rbs[i].clear();
					n_aas.push_back(missing_aas[i]);
					c_aas.push_back(0);
				}
			}
			else
			{	
				peak_prediction_rbs.resize(num_missing_aas*num_missing_aas);
				int i;
				for (i=0; i<num_missing_aas; i++)
				{
					int j;
					for (j=0; j<num_missing_aas; j++)
					{
						peak_prediction_rbs[i*num_missing_aas+j].clear();
						n_aas.push_back(missing_aas[i]);
						c_aas.push_back(missing_aas[j]);
					}
				}
			}

			int i;
			for (i=0; i<peak_prediction_rbs.size(); i++)
			{
				sol.most_basic_aa_removed_from_n=n_aas[i];
				sol.most_basic_aa_removed_from_c=c_aas[i];
				part_model->fill_combined_peak_prediction_features(sol, ann_intens,peak_model,peak_prediction_rbs[i],size_idx);
			}
		}
		return;

	}
	else
	{
		cout << "Error: fill_denovo_peptide_rbs with combos should only be used with model type 1, not " <<
			model_type << endl;
		exit(1);
	}
}



void PeptideRankScorer::fill_tag_rbs(PeptideSolution& sol,
					  const SeqPath& path,
					  QCPeak* peaks, 
					  int num_peaks, 
					  AnnotatedSpectrum& as,
					  RankBoostSample&	rbs,
					  int size_idx) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();

	const PrmGraph *prm = path.prm_ptr;
	vector< vector<intensity_t> > ann_intens;
	vector< vector<mass_t> >	  ann_masses;

	as.set_peptide(sol.pep);
	as.annotate_spectrum(sol.pm_with_19, sol.charge, true);
	as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

	rbs.clear();
	const int charge=sol.charge;

	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	if (size_idx<0)
		size_idx=peak_model->get_size_group(charge,sol.pm_with_19);

	DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];

	if (! part_model || ! part_model->ind_was_initialized)
	{
		cout << "Error: de novo partition model was not initialized, charge " <<
			charge << " size " << size_idx << endl;
		exit(1);
	}

	if (model_type == 3)
	{
		if (part_model->use_PTM_peak_features)
			part_model->fill_PTM_peak_features(allScoreModels->get_config(),sol,ann_masses,ann_intens,as,rbs);

		if (part_model->use_tryp_terminal_features)
			part_model->fill_tryp_terminal_features(sol,path,rbs);

		if (part_model->use_ann_peak_features)
			part_model->fill_ann_peak_features(sol,ann_masses,ann_intens,as,rbs);

		if (part_model->use_peak_offset_features)
			part_model->fill_peak_offset_features(as.getConfig(),sol,ann_masses,ann_intens,rbs);

		if (part_model->use_comp_features)
			part_model->fill_composition_features(sol,as.getConfig(),&compAssigner, path, rbs);

		if (part_model->use_prm_features)
			part_model->fill_prm_features(sol, path, model_type, rbs);
			
		if (part_model->use_combined_ppp_features) 
		{
			if (peak_model->get_feature_set_type() != 4)
			{
				cout << "Error: feature type for filling peak predictions should be 4 with this function,";
				cout << " not " << peak_model->get_feature_set_type() << endl;
				exit(1);
			}

			if (! sol.reaches_c_terminal)
				sol.most_basic_aa_removed_from_c=Lys;

			part_model->fill_combined_peak_prediction_features(sol, ann_intens,peak_model,rbs,size_idx);
		}
		return;
	}
	else
	{
		cout << "Error: fill_tag_rbs should only be used with model type 3, not " <<
			model_type << endl;
		exit(1);
	}
}


void PeptideRankScorer::score_complete_sequences(
								const vector<PeptideSolution>& peptide_sols,
								SingleSpectrumFile *ssf,
								QCPeak* peaks, 
								int num_peaks,
								vector<score_pair>& score_pairs,
								int forced_size_idx) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	AnnotatedSpectrum as;
	vector<PmcSqsChargeRes> pmc_sqs_res;

	score_pairs.clear();
	if (peptide_sols.size()==0)
		return;

	int charge1=0,charge2=0;
	mass_t mz1=0,mz2=0;
	float prob1=0,prob2=0;

	BasicSpectrum bs;
	bs.peaks = peaks;
	bs.num_peaks = num_peaks;
	bs.ssf = ssf;

	Config *config = allScoreModels->get_config();
	score_pairs.clear();

	allScoreModels->get_best_mz_charge(config,bs, &mz1, &charge1, &prob1, 
							  &mz2, &charge2, &prob2, &pmc_sqs_res);

	as.init_from_QCPeaks(config,peaks,num_peaks,ssf);

	
	vector<ScalingFactor> same_scaling_factors;
	same_scaling_factors.resize(10);

	int i;
	for (i=0; i<peptide_sols.size(); i++)
	{
		const PeptideSolution& sol = peptide_sols[i];
		const Peptide& pep = sol.pep;
		const mass_t mass_with_19 = pep.get_mass_with_19();
		const int charge = sol.charge;
		const int size_idx = (forced_size_idx>=0 ? forced_size_idx : 
							  peak_model->get_size_group(charge,mass_with_19) );
		
		
		DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];
		if (! part_model|| ! part_model->ind_was_initialized)
		{
			cout << "Error: no de novo rank model for charge " << charge << " size " << size_idx << endl;
			exit(1);
		}
		if (same_scaling_factors[charge].score_scale == 1.0)  // makes sure the same scale is used for all sequences with this charge
			same_scaling_factors[charge] = part_model->get_scaling_factor(mass_with_19);	
	
		RankBoostSample rbs;
		as.set_peptide(pep);
		as.annotate_spectrum(mass_with_19,charge);
		fill_complete_peptide_rbs(sol, peaks, num_peaks, as, pmc_sqs_res, rbs, size_idx);

		float score=part_model->boost_model.calc_rank_score(rbs);
		score += same_scaling_factors[charge].score_shift;
		score *= same_scaling_factors[charge].score_scale;


		if (pep.get_num_aas()<8)
		{
			score -=1;
			if (pep.get_num_aas()<7)
				score -=1;
			if (pep.get_num_aas()<6)
				score -=5;
		}
		score_pairs.push_back(score_pair(i,score));
	}
}

void PeptideRankScorer::scoreCompleteSequences(const vector<PeptideSolution>& peptideSols,
								AnnotatedSpectrum& as,
								vector<score_pair>& scorePairs,
								int forcedSizeIndex) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	vector<PmcSqsChargeRes> pmc_sqs_res;

	scorePairs.clear();
	if (peptideSols.size()==0)
		return;

	int charge1=0,charge2=0;
	mass_t mz1=0,mz2=0;
	float prob1=0,prob2=0;

	Config *config = allScoreModels->get_config();
	scorePairs.clear();

	allScoreModels->findBestMzAndCharge(config, as, &mz1, &charge1, &prob1, 
							  &mz2, &charge2, &prob2, &pmc_sqs_res);

	vector<ScalingFactor> same_scaling_factors;
	same_scaling_factors.resize(10);

	int i;
	for (i=0; i<peptideSols.size(); i++)
	{
		const PeptideSolution& sol = peptideSols[i];
		const Peptide& pep = sol.pep;
		const mass_t mass_with_19 = pep.get_mass_with_19();
		const int charge = sol.charge;
		const int size_idx = (forcedSizeIndex>=0 ? forcedSizeIndex : 
							  peak_model->get_size_group(charge,mass_with_19) );
		
		
		DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];
		if (! part_model|| ! part_model->ind_was_initialized)
		{
			cout << "Error: no de novo rank model for charge " << charge << " size " << size_idx << endl;
			exit(1);
		}
		if (same_scaling_factors[charge].score_scale == 1.0)  // makes sure the same scale is used for all sequences with this charge
			same_scaling_factors[charge] = part_model->get_scaling_factor(mass_with_19);	
	
		RankBoostSample rbs;
		as.set_peptide(pep);
		as.annotate_spectrum(mass_with_19,charge);
	//	fill_complete_peptide_rbs(sol, as.getPeaks(), as.getNumPeaks(), num_peaks, as, pmc_sqs_res, rbs, size_idx);
		fillCompletePeptideRankboostSample(sol, as, pmc_sqs_res, rbs, size_idx);

	//	rbs.print();
	//	exit(0);

		float score=part_model->boost_model.calc_rank_score(rbs);
		score += same_scaling_factors[charge].score_shift;
		score *= same_scaling_factors[charge].score_scale;


		if (pep.get_num_aas()<8)
		{
			score -=1.0;
			if (pep.get_num_aas()<7)
				score -=1.0;
			if (pep.get_num_aas()<6)
				score -=5.0;
		}
		scorePairs.push_back(score_pair(i,score));
	}

}


void PeptideRankScorer::scoreDenovoSequences(
							  const vector<SeqPath>& seqPaths,
							  AnnotatedSpectrum& as,
							  vector<score_pair>& scorePairs,
							  int forcedSizeIndex,
							  int maximalIndexForRanking) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	vector<PmcSqsChargeRes> pmc_sqs_res;

	int charge1=0,charge2=0;
	mass_t mz1=0,mz2=0;
	float prob1=0,prob2=0;


	Config *config = allScoreModels->get_config();
	scorePairs.clear();

//	model->get_best_mz_charge(config,bs, &mz1, &charge1, &prob1, 
//							  &mz2, &charge2, &prob2, &pmc_sqs_res);
	allScoreModels->findBestMzAndCharge(config, as, &mz1, &charge1, &prob1, 
							  &mz2, &charge2, &prob2, &pmc_sqs_res);

	if (mz1<0.0)
		return;

	const mass_t corr_mass_with_19 = mz1*charge1 - (charge1-1.0)*MASS_PROTON;
	as.set_corrected_pm_with_19(corr_mass_with_19);

	vector<ScalingFactor> scaling_factors;
	scaling_factors.resize(10);

	int i;
	for (i=0; i<seqPaths.size(); i++)
	{
		if (maximalIndexForRanking>0 && i>= maximalIndexForRanking)
		{
			scorePairs.push_back(score_pair(i,-(float)i));
			continue;
		}

		const SeqPath& path = seqPaths[i];
		// create a peptide solution to represent the path
		PeptideSolution sol;

		vector<int> amino_acids;
		path.get_amino_acids(amino_acids);

		sol.charge = path.charge;
		sol.pm_with_19 = path.pm_with_19;
		sol.pep.set_peptide_aas(amino_acids);
		sol.pep.set_n_gap(path.n_term_mass);
		sol.reaches_n_terminal = (sol.pep.get_n_gap()<1.0);
		sol.pep.calc_mass(config);
		sol.reaches_c_terminal = (sol.pep.get_mass_with_19()>sol.pm_with_19 - 7);
		sol.type = -1;

		if (amino_acids.size()<3)
		{
			scorePairs.push_back(score_pair(i,NEG_INF));
			continue;
		}
		
		// choose pm according to if peptide reaches both ends

		const mass_t mass_with_19 =( (sol.reaches_n_terminal && sol.reaches_c_terminal) ?
			sol.pep.get_mass_with_19() : corr_mass_with_19);
	
		const int charge   = path.charge;
		const int size_idx = (forcedSizeIndex>=0 ? forcedSizeIndex : 
							  peak_model->get_size_group(charge,mass_with_19) );

		RankBoostSample		    main_rbs;			// holds the features that are common to all variants of the SeqPath
		vector<RankBoostSample> combo_variant_rbs; // holds features that might depend on the identity of the
												  // missing amino acids on the n and c-terminal sides

		as.set_peptide(sol.pep);
		as.annotate_spectrum(mass_with_19,charge);

//		fill_denovo_peptide_rbs_with_combos(sol, path, peaks, num_peaks, as, pmc_sqs_res, main_rbs, combo_variant_rbs, size_idx);
		fillDenovoPeptideRankboostSampleWithCombos(sol, path, as, pmc_sqs_res, main_rbs, combo_variant_rbs, size_idx);

		// check if we have models in this range at all
		if (charge >= dnv_part_models.size() ||
			size_idx >= dnv_part_models[charge].size())
		{
			scorePairs.push_back(score_pair(i,NEG_INF));
			continue;
		}

		DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];
		if (! part_model || ! part_model->ind_was_initialized)
		{
			cout << "Error: no de novo rank model for charge " << charge << " size " << size_idx << endl;
			exit(1);
		}
		if (scaling_factors[charge].score_scale == 1.0)  // makes sure the same scale is used for all sequences with this charge
			scaling_factors[charge] = part_model->get_scaling_factor(mass_with_19);	

		float max_score = NEG_INF;
		int v_idx;
		for (v_idx = 0; v_idx<combo_variant_rbs.size(); v_idx++)
		{
			RankBoostSample  rbs = main_rbs;
			RankBoostSample& var_rbs = combo_variant_rbs[v_idx];
			int j;
			for (j=0; j<var_rbs.real_active_idxs.size(); j++)
				rbs.add_real_feature(var_rbs.real_active_idxs[j],var_rbs.real_active_values[j]);

			const float score=part_model->boost_model.calc_rank_score(rbs);
			if (score>max_score)
				max_score=score;
		}

		max_score += scaling_factors[charge].score_shift;
		max_score *= scaling_factors[charge].score_scale;
		scorePairs.push_back(score_pair(i,max_score));
	}
}




void PeptideRankScorer::score_denovo_sequences(
								const vector<SeqPath>& seq_paths,
								SingleSpectrumFile *ssf,
								QCPeak* peaks, 
								int num_peaks,
								vector<score_pair>& score_pairs,
								int forced_size_idx,
								int max_idx_for_ranking) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	AnnotatedSpectrum as;
	vector<PmcSqsChargeRes> pmc_sqs_res;

	int charge1=0,charge2=0;
	mass_t mz1=0,mz2=0;
	float prob1=0,prob2=0;

	BasicSpectrum bs;
	bs.peaks = peaks;
	bs.num_peaks = num_peaks;
	bs.ssf = ssf;

	Config *config = allScoreModels->get_config();
	score_pairs.clear();

	allScoreModels->get_best_mz_charge(config,bs, &mz1, &charge1, &prob1, 
							  &mz2, &charge2, &prob2, &pmc_sqs_res);

	const mass_t corr_mass_with_19 = mz1*charge1 - (charge1-1.0);
	as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
	as.set_corrected_pm_with_19(corr_mass_with_19);

	vector<ScalingFactor> scaling_factors;
	scaling_factors.resize(10);

	int i;
	for (i=0; i<seq_paths.size(); i++)
	{
		if (max_idx_for_ranking>0 && i>=max_idx_for_ranking)
		{
			score_pairs.push_back(score_pair(i,-(float)i));
			continue;
		}

		const SeqPath& path = seq_paths[i];
		// create a peptide solution to represent the path
		PeptideSolution sol;

		vector<int> amino_acids;
		path.get_amino_acids(amino_acids);

		sol.charge = path.charge;
		sol.pm_with_19 = path.pm_with_19;
		sol.pep.set_peptide_aas(amino_acids);
		sol.pep.set_n_gap(path.n_term_mass);
		sol.reaches_n_terminal = (sol.pep.get_n_gap()<1.0);
		sol.pep.calc_mass(config);
		sol.reaches_c_terminal = (sol.pep.get_mass_with_19()>sol.pm_with_19 - 7);
		sol.type = -1;

		if (amino_acids.size()<3)
		{
			score_pairs.push_back(score_pair(i,NEG_INF));
			continue;
		}
		
		// choose pm according to if peptide reaches both ends

		const mass_t mass_with_19 =( (sol.reaches_n_terminal && sol.reaches_c_terminal) ?
			sol.pep.get_mass_with_19() : corr_mass_with_19);
	
		const int charge   = path.charge;
		const int size_idx = (forced_size_idx>=0 ? forced_size_idx : 
							  peak_model->get_size_group(charge,mass_with_19) );

		RankBoostSample		    main_rbs;			// holds the features that are common to all variants of the SeqPath
		vector<RankBoostSample> combo_variant_rbs; // holds features that might depend on the identity of the
												  // missing amino acids on the n and c-terminal sides

		as.set_peptide(sol.pep);
		as.annotate_spectrum(mass_with_19,charge);

		fill_denovo_peptide_rbs_with_combos(sol, path, peaks, num_peaks, as, pmc_sqs_res, main_rbs, combo_variant_rbs, size_idx);

		DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];
		if (! part_model || ! part_model->ind_was_initialized)
		{
			cout << "Error: no de novo rank model for charge " << charge << " size " << size_idx << endl;
			exit(1);
		}
		if (scaling_factors[charge].score_scale == 1.0)  // makes sure the same scale is used for all sequences with this charge
			scaling_factors[charge] = part_model->get_scaling_factor(mass_with_19);	

		float max_score = NEG_INF;
		int v_idx;
		for (v_idx = 0; v_idx<combo_variant_rbs.size(); v_idx++)
		{
			RankBoostSample  rbs = main_rbs;
			RankBoostSample& var_rbs = combo_variant_rbs[v_idx];
			int j;
			for (j=0; j<var_rbs.real_active_idxs.size(); j++)
				rbs.add_real_feature(var_rbs.real_active_idxs[j],var_rbs.real_active_values[j]);

			const float score=part_model->boost_model.calc_rank_score(rbs);
			if (score>max_score)
				max_score=score;
		}

		max_score += scaling_factors[charge].score_shift;
		max_score *= scaling_factors[charge].score_scale;
		score_pairs.push_back(score_pair(i,max_score));
	}
}



void PeptideRankScorer::list_feature_differences(const vector<SeqPath>& seq_paths,
								  SingleSpectrumFile *ssf,
								  QCPeak* peaks, 
								  int num_peaks) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	AnnotatedSpectrum as;
	vector<PmcSqsChargeRes> pmc_sqs_res;

	int charge1=0,charge2=0;
	mass_t mz1=0,mz2=0;
	float prob1=0,prob2=0;

	BasicSpectrum bs;
	bs.peaks = peaks;
	bs.num_peaks = num_peaks;
	bs.ssf = ssf;

	Config *config = allScoreModels->get_config();

	allScoreModels->get_best_mz_charge(config,bs, &mz1, &charge1, &prob1, 
							  &mz2, &charge2, &prob2, &pmc_sqs_res);

	const mass_t corr_mass_with_19 = mz1*charge1 - (charge1-1.0);
	as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
	as.set_corrected_pm_with_19(corr_mass_with_19);

	vector<ScalingFactor> scaling_factors;
	scaling_factors.resize(10);

	const mass_t mass_with_19 = seq_paths[0].pm_with_19;
	const int charge   = seq_paths[0].charge;
	const int size_idx = peak_model->get_size_group(charge,mass_with_19);
	

	DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];
	if (! part_model || ! part_model->ind_was_initialized)
	{
		cout << "Error: no de novo rank model for charge " << charge << " size " << size_idx << endl;
		exit(1);
	}
	if (scaling_factors[charge].score_scale == 1.0)  // makes sure the same scale is used for all sequences with this charge
		scaling_factors[charge] = part_model->get_scaling_factor(mass_with_19);	

	vector<RankBoostSample> boost_samples;
	int i;
	for (i=0; i<seq_paths.size(); i++)
	{
	

		const SeqPath& path = seq_paths[i];
		// create a peptide solution to represent the path
		PeptideSolution sol;

		vector<int> amino_acids;
		path.get_amino_acids(amino_acids);

		sol.charge = path.charge;
		sol.pm_with_19 = path.pm_with_19;
		sol.pep.set_peptide_aas(amino_acids);
		sol.pep.set_n_gap(path.n_term_mass);
		sol.reaches_n_terminal = (sol.pep.get_n_gap()<1.0);
		sol.pep.calc_mass(config);
		sol.reaches_c_terminal = (sol.pep.get_mass_with_19()>sol.pm_with_19 - 7);
		sol.type = -1;
		
		// choose pm according to if peptide reaches both ends

		RankBoostSample		    main_rbs;			// holds the features that are common to all variants of the SeqPath
		vector<RankBoostSample> combo_variant_rbs; // holds features that might depend on the identity of the
												  // missing amino acids on the n and c-terminal sides
		as.set_peptide(sol.pep);
		as.annotate_spectrum(mass_with_19,sol.charge);
		fill_denovo_peptide_rbs_with_combos(sol, path, peaks, num_peaks, as, pmc_sqs_res, main_rbs, combo_variant_rbs, size_idx);

		float max_score = NEG_INF;
		int v_idx;
		for (v_idx = 0; v_idx<combo_variant_rbs.size(); v_idx++)
		{
			RankBoostSample  rbs = main_rbs;
			RankBoostSample& var_rbs = combo_variant_rbs[v_idx];
			int j;
			for (j=0; j<var_rbs.real_active_idxs.size(); j++)
				rbs.add_real_feature(var_rbs.real_active_idxs[j],var_rbs.real_active_values[j]);

			boost_samples.push_back(rbs);
		}
	}
	part_model->boost_model.list_feature_vals_according_to_score(boost_samples);
}


void PeptideRankScorer::score_tag_sequences(
								const vector<SeqPath>& seq_paths,
								SingleSpectrumFile *ssf,
								QCPeak* peaks, 
								int num_peaks,
								vector<score_pair>& score_pairs,
								int forced_size_idx) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	AnnotatedSpectrum as;
	vector<PmcSqsChargeRes> pmc_sqs_res;

	int charge1=0,charge2=0;
	mass_t mz1=0,mz2=0;
	float prob1=0,prob2=0;

	BasicSpectrum bs;
	bs.peaks = peaks;
	bs.num_peaks = num_peaks;
	bs.ssf = ssf;

	Config *config = allScoreModels->get_config();
	

	as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
	
	vector<ScalingFactor> scaling_factors;
	scaling_factors.resize(10);
	score_pairs.clear();

	int i;
	for (i=0; i<seq_paths.size(); i++)
	{
		const SeqPath& path = seq_paths[i];
		// create a peptide solution to represent the path
		PeptideSolution sol;

		vector<int> amino_acids;
		path.get_amino_acids(amino_acids);

		sol.charge = path.charge;
		sol.pm_with_19 = path.pm_with_19;
		sol.pep.set_peptide_aas(amino_acids);
		sol.pep.set_n_gap(path.n_term_mass);
		sol.reaches_n_terminal = (sol.pep.get_n_gap()<1.0);
		sol.pep.calc_mass(config);
		sol.reaches_c_terminal = (sol.pep.get_mass_with_19()>sol.pm_with_19 - 7);
		sol.type = -1;

		if (amino_acids.size()<3)
		{
			score_pairs.push_back(score_pair(i,NEG_INF));
			continue;
		}
		
		// choose pm according to if peptide reaches both ends

		const mass_t mass_with_19 =( (sol.reaches_n_terminal && sol.reaches_c_terminal) ?
			sol.pep.get_mass_with_19() : path.pm_with_19);
	
		const int charge   = path.charge;
		const int size_idx = (forced_size_idx>=0 ? forced_size_idx : 
							  peak_model->get_size_group(charge,mass_with_19) );

		RankBoostSample		    rbs;			// holds the features that are common to all variants of the SeqPath
	
		as.set_peptide(sol.pep);
		as.annotate_spectrum(mass_with_19,sol.charge);
		fill_tag_rbs(sol,path,peaks,num_peaks,as,rbs,size_idx);

		if (charge>=dnv_part_models.size() || size_idx>=dnv_part_models[charge].size())
		{
			score_pairs.push_back(score_pair(i,NEG_INF));
			continue;
		}

		DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];
		if (! part_model || ! part_model->ind_was_initialized)
		{
			cout << "Error: no de novo rank model for charge " << charge << " size " << size_idx << endl;
			exit(1);
		}
		if (scaling_factors[charge].score_scale == 1.0)  // makes sure the same scale is used for all sequences with this charge
			scaling_factors[charge] = part_model->get_scaling_factor(mass_with_19);	

		float score = part_model->boost_model.calc_rank_score(rbs);

		score += scaling_factors[charge].score_shift;
		score *= scaling_factors[charge].score_scale;
		score_pairs.push_back(score_pair(i,score));
	}
}



void PeptideRankScorer::read_denovo_rank_scorer_model(const char *path, string type_string, bool silent_ind, bool use_large_db_model)
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();

	ifstream ifs(path);
	if (! ifs.good() || ! ifs.is_open())
	{
		cout << "Error: couldn't open file for reading: " << path << endl;
		exit(1);
	}

	ifs >> dnv_model_name >> model_type;
	if (model_type<0 || model_type>10)
	{
		cout << "Error: bad model type specified!" << endl;
		exit(1);
	}

	// HACK
	if (use_large_db_model)
		dnv_model_name = "DBsn";


	string model_str, peak_model_str, comp_assign_str;
	int max_charge;
	ifs >> max_charge >> model_str >> peak_model_str >> comp_assign_str;

	if (! allScoreModels)
	{
		cout << "Error: must set score model first!" << endl;
		exit(1);
	}

	PeakRankModel*&  peakRankModels = allScoreModels->get_peak_prediction_model_ptr(model_type);
	if (! peakRankModels)
	{
		peakRankModels = new PeakRankModel;
		if (! peakRankModels->read_peak_rank_model(allScoreModels->get_config(),peak_model_str.c_str(),true))
		{
			cout << "Error: couldn't read peak model " << peak_model_str << endl;
			exit(1);
		}
	}

	if (! compAssigner.get_ind_was_initialized())
		compAssigner.read_and_init_from_tables(allScoreModels->get_config(),comp_assign_str.c_str());
	

	string score_model_name = allScoreModels->get_model_name();

	int num_models_read=0;
	init_tables(silent_ind);
	int c=0;
	for (c=0; c<dnv_part_models.size(); c++)
	{
		int size_idx;
		for (size_idx=0; size_idx<dnv_part_models[c].size(); size_idx++)
		{
			ostringstream oss;
			oss << allScoreModels->get_config()->get_resource_dir() << "/" <<score_model_name << "_" << 
				type_string << "/" << dnv_model_name << "_" << c << "_" << size_idx << "_model.txt";

		//	cout << oss.str() << endl;

			ifstream ifs(oss.str().c_str());
			if (! ifs.is_open())
				continue;
			ifs.close();
			dnv_part_models[c][size_idx]=new DeNovoPartitionModel;
			if (dnv_part_models[c][size_idx]->read_denovo_part_model(oss.str().c_str(),allScoreModels->get_config()))
			{
				if (! silent_ind)
					cout << "Read de novo rank model " << c << " " << size_idx << " " << oss.str() << endl;
				num_models_read++;

			//	cout << dnv_part_models[c][size_idx]->ind_was_initialized <<  dnv_part_models[c][size_idx] << endl;
			}
		}
	}

	if (! silent_ind)
		cout << "Read " << num_models_read << " de novo rank models..." << endl;
}



void PeptideRankScorer::write_denovo_rank_scorer_model(char *name)
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);

	if (! allScoreModels || ! allScoreModels->get_ind_pmcsqs_was_intialized())
	{
		cout << "Error: model not initialized!" << endl;
		exit(1);
	}

	PeakRankModel*& peakRankModel = allScoreModels->get_peak_prediction_model_ptr(model_type);

	if (! peakRankModel )
	{
		cout << "Error: peak model not initialized!" << endl;
		exit(1);
	}

	PeptideCompAssigner& compAssigner = allScoreModels->getPeptideCompositionAssigner();


	if (! compAssigner.get_ind_was_initialized())
	{
		cout << "Error: comp assigner not initialized!" << endl;
		exit(1);
	}

	string dir = allScoreModels->get_config()->get_resource_dir();
	string base_path = dir + "/" + name;
	string model_file_name = base_path + "_model.txt";

	ofstream ofs(model_file_name.c_str());
	if (! ofs.is_open() || ! ofs.good())
	{
		cout << "Error: couldn't open file for reading: " << model_file_name << endl;
		exit(1);
	}

	ofs << name << " " << model_type << endl;
	ofs << allScoreModels->get_max_score_model_charge() << endl;
	ofs << allScoreModels->get_model_name() << endl;
	ofs << peakRankModel->get_peak_rank_model_name() << endl;
	ofs << compAssigner.get_model_name() << endl;

	int c;
	for (c=0; c<this->dnv_part_models.size(); c++)
	{
		int size_idx;
		for (size_idx=0; size_idx<dnv_part_models[c].size(); size_idx++)
		{
			if (dnv_part_models[c][size_idx] && 
				dnv_part_models[c][size_idx]->ind_was_initialized)
			{
				ostringstream oss;
				oss << base_path << "_" << c << "_" << size_idx << ".txt";
				string path = oss.str();
				dnv_part_models[c][size_idx]->write_denovo_partition_model(model_type, path.c_str());
			}
		}
	}
}



/*********************************************************************
Adds the counts for peaks around a breakage to their respective bins

**********************************************************************/
void add_offset_counts_for_unannotated_peaks_arround_mass(vector<int>& counts, 
									AnnotatedSpectrum *spec,
									mass_t min_mass, 
									mass_t max_mass, 
									mass_t bin_coef, 
									mass_t break_mass,
									int charge)
{
	const vector< vector<PeakAnnotation> >& peak_anns = spec->get_peak_annotations();
	const mass_t low_range = (break_mass + min_mass+1)/(mass_t)charge;
	const mass_t high_range = (break_mass + max_mass)/(mass_t)charge;
	const PeakRange pr = spec->findPeaksInRange(low_range, high_range);

//	cout << spec->get_peptide().as_string(spec->get_config()) << " " << break_mass << endl;
	
	if (pr.num_peaks<=0)
		return;

	int skipped=0;

	// add counts
	int p_idx;
	for (p_idx = pr.low_idx; p_idx<=pr.high_idx; p_idx++)
	{
		if (spec->get_peak_iso_level(p_idx)>0)
			continue;

		if (peak_anns[p_idx].size()>0)
		{
			skipped++;
			continue;
		}
	
		const mass_t peak_mass = spec->getPeakMass(p_idx);
		const mass_t b_mass = break_mass/(mass_t)charge;
		const int bin = (int)((peak_mass - b_mass - min_mass)*bin_coef);
		
		counts[bin]+=10;
		counts[bin-1]+=9;
		counts[bin+1]+=9;
		counts[bin-2]+=6;
		counts[bin+2]+=6;
		counts[bin-3]+=4;
		counts[bin+3]+=4;
		counts[bin-4]+=2;
		counts[bin+4]+=2;
	}
}


/***********************************************************************
Uses the offset count method to determine if there are any fragments
that are special for a given PTM (function only prints a list of
the most interesting cases).
************************************************************************/
void find_special_PTM_frags_using_offset_counts(
										   const string& PTM_label,
										   FileManager& fm,
										   const vector<SingleSpectrumFile *>& all_ssfs,
										   AllScoreModels *model,
										   int max_charge)
{
	Config *config = model->get_config();
	const mass_t min_offset_mass = -120;
	const mass_t max_offset_mass = 120;
	const mass_t tolerance = config->getTolerance();
	const mass_t bin_size = tolerance * 0.1;
	const mass_t bin_coef = 1.0 / bin_size;
	const int count_size = (int)((max_offset_mass - min_offset_mass + 1) / bin_size);
	vector< vector< vector<int> > > prefix_counts, suffix_counts; // charge, distance from cut, bin_idx
	vector< vector< int > > pre_instances, suf_instances;
	int i,c,d;

	float min_frag_prob = 0.02;

	const int ptm_aa = config->get_aa_from_label(PTM_label);
	if (ptm_aa <0)
	{
		cout << "Error: PTM not supported in this model: " << PTM_label << endl;
		exit(1);
	}

	const int max_distance = 1;
	
	prefix_counts.resize(max_charge+1);
	suffix_counts.resize(max_charge+1);
	pre_instances.resize(max_charge+1);
	suf_instances.resize(max_charge+1);

	for (c=1; c<=max_charge; c++)
	{
		prefix_counts[c].resize(max_distance+1);
		suffix_counts[c].resize(max_distance+1);
		pre_instances[c].resize(max_distance+1,0);
		suf_instances[c].resize(max_distance+1,0);
		for (d=0; d<=max_distance; d++)
		{
			prefix_counts[c][d].resize(count_size,0);
			suffix_counts[c][d].resize(count_size,0);
		}
	}

	vector<QCPeak> peaks;
	BasicSpecReader bsr;

	peaks.resize(5000);

	int spectra_used=0;
	for (i=0; i<all_ssfs.size(); i++)
	{
		AnnotatedSpectrum as;
		vector<mass_t> break_masses;

		const Peptide& pep = all_ssfs[i]->peptide;
		const vector<int>& amino_acids = pep.get_amino_acids();
		int aa_idx;
		for (aa_idx = 0; aa_idx<amino_acids.size(); aa_idx++)
			if (amino_acids[aa_idx] == ptm_aa)
				break;
		if (aa_idx == amino_acids.size())
			continue;
		
		int num_peaks = bsr.read_basic_spec(config,fm,all_ssfs[i],&peaks[0]);
		as.init_from_QCPeaks(config,&peaks[0],num_peaks,all_ssfs[i]);
		as.set_peptide(pep);
		as.annotate_spectrum(pep.get_mass_with_19(),0,true);
		spectra_used++;

		
		pep.calc_expected_breakage_masses(config,break_masses);
		const mass_t true_mass = as.getPeptide().get_mass();
		for (aa_idx=0; aa_idx<amino_acids.size(); aa_idx++)
		{
			if (amino_acids[aa_idx] != ptm_aa)
				continue;

			int d;
			for (d=0; d<=max_distance; d++)
			{
				int cut_idx = aa_idx + 1 -d;
				if (cut_idx<1)
					break;

				int c;
				const mass_t break_mass = break_masses[cut_idx];
				for (c=1; c<=max_charge; c++)
				{
				//	cout << "P: " << c << " " << d << " ";
					add_offset_counts_for_unannotated_peaks_arround_mass(prefix_counts[c][d], &as,
						min_offset_mass, max_offset_mass, bin_coef, break_mass,c);
					pre_instances[c][d]++;
				}
			}

			for (d=0; d<=max_distance; d++)
			{
				int cut_idx = aa_idx + d;
				if (cut_idx>amino_acids.size())
					break;
				
				int c;
				const mass_t break_mass = break_masses[cut_idx];
				for (c=1; c<=max_charge; c++)
				{
				//	cout << "S: " << c << " " << d << " ";
					add_offset_counts_for_unannotated_peaks_arround_mass(suffix_counts[c][d], &as,
						min_offset_mass, max_offset_mass, bin_coef, (true_mass - break_mass),c);
					suf_instances[c][d]++;
				}
			}
		}
	}

	cout << "Using: " << spectra_used << " spectra for offset counts..." << endl;

	// select 30 top fragments becasue many are likely to be caused by previous/next
	// amino acids and will be later removed
	vector<FragmentTypeSet> pre_fts , suf_fts;

	pre_fts.resize(max_distance+1);
	suf_fts.resize(max_distance+1);
	for (c=1; c<=max_charge; c++)
	{
		int d;
		for (d=0; d<=max_distance; d++)
		{
			select_fragments_from_bins(prefix_counts[c][d],pre_fts[d],30,c,PREFIX,min_offset_mass,bin_coef,tolerance);
			select_fragments_from_bins(suffix_counts[c][d],suf_fts[d],30,c,SUFFIX,min_offset_mass,bin_coef,tolerance);

		/*	int f;
			for (f=0; f<pre_fts[d].get_num_fragments(); f++)
			{
				FragmentType& frag = pre_fts[d].get_non_const_fragment(f);
				frag.prob = frag.prob / (float)pre_instances[frag.charge][d];
			}
			pre_fts[d].sort_fragments_according_to_probs();

			for (f=0; f<suf_fts[d].get_num_fragments(); f++)
			{
				FragmentType& frag = suf_fts[d].get_non_const_fragment(f);
				frag.prob = frag.prob / (float)suf_instances[frag.charge][d];
			}
			suf_fts[d].sort_fragments_according_to_probs();*/
		}
	}


	for (d=0; d<=max_distance; d++)
	{
	//	calculate_true_fragment_probabilities(fm,config,fts[d], min_frag_prob);

		cout << "Fragments selected from spectra [distance " << d << "]:" << endl;
		for (i=0; i<pre_fts[d].get_num_fragments() && i<5; i++)
		{
			const FragmentType& frag = pre_fts[d].get_fragment(i);
			cout << left << setw(3) << i << right << setw(5) << frag.spec_count << " ";
			cout << setw(6) << setprecision(3) << right << frag.prob << " ";
			frag.write_fragment(cout);
		}
		cout << endl;

		for (i=0; i<suf_fts[d].get_num_fragments() && i<5; i++)
		{
			const FragmentType& frag = suf_fts[d].get_fragment(i);
			cout << left << setw(3) << i << right << setw(5) << frag.spec_count << " ";
			cout << setw(6) << setprecision(3) << right << frag.prob << " ";
			frag.write_fragment(cout);
		}
		cout << endl;
	}
}





/****************************************************************************
Runs benchmark on spectra with a complete de novo solution in there
*****************************************************************************/
void benchmark_ranking_on_full_denovo(AllScoreModels *model,
									  char *mgf_test_file, 
									  int max_num_spectra,
									  int num_solutions,
									  char *report_path,
									  int min_length,
									  int max_length)
{
	const int charge = 2;

	Config *config;

	char report_file_path[512];
	char *report_name=NULL;

	ofstream out_file_stream;
	
	if (report_path)
	{
		sprintf(report_file_path,"%s_bench_test.txt",report_path);
		report_name = report_file_path;
		out_file_stream.open(report_file_path);
		cout << "Writing results to " << report_file_path << endl;
	}
	else
		cout << "Writing results to cout..." << endl;


	
	ostream& sum_stream	= (report_path ? out_file_stream : cout);

	FileManager fm;
	FileSet		fs;
	vector<PrmGraph *> prm_ptrs;
	PrmGraph opt_prm;

	config = model->get_config();
	config->set_use_spectrum_charge(2);
	const mass_t tolerance = config->getTolerance();


	fm.init_from_file(config,mgf_test_file);
	fs.select_all_files(fm);
	fs.randomly_reduce_ssfs((int)(max_num_spectra * 1.5));
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	
	// 0 all seqs , 1 complete seqs, 2 all cuts ,3 complete cuts
	vector< vector<int> > top_rerank, top_denovo;
	vector<double> total_corr_denovo, total_pred_denovo, total_corr_rerank, total_pred_rerank;

	top_rerank.resize(4);
	top_denovo.resize(4);
	total_corr_denovo.resize(4,0);
	total_pred_denovo.resize(4,0);
	total_corr_rerank.resize(4,0);
	total_pred_rerank.resize(4,0);

	BasicSpecReader bsr;
	QCPeak peaks[4000];

	int num_spec_tested=0;
	int num_top_denovo_correct=0;
	int num_top_rerank_correct=0;

	sum_stream << "Test file     : " << mgf_test_file << endl;
	sum_stream << "Num cases     : " << max_num_spectra << endl;
	sum_stream << "Num solutions : " << num_solutions << endl << endl;
	sum_stream << "START TESTING..." << endl;
	double start_t = time(NULL);

	double num_pred_aa=0;

	int spec_idx;
	for (spec_idx=0; spec_idx<all_ssf.size(); spec_idx++)
	{
	
		MGF_single *ssf = (MGF_single *)all_ssf[spec_idx];
		const string corr_pep_str = ssf->peptide.as_string(config);
		const Peptide& full_pep = ssf->peptide;
		const mass_t true_mass_with_19 = full_pep.get_mass_with_19();
		const int correct_pep_length = ssf->peptide.get_num_aas();
		BasicSpectrum bs;
		Spectrum s;

		// exclude pepitdes with M+16/Q-17
		const vector<int>& aas = ssf->peptide.get_amino_acids();
		int a;
		for (a=0; a<aas.size(); a++)
			if (aas[a]>Val)
				break;
		if (a<aas.size())
			continue;


		bs.ssf = all_ssf[spec_idx];
		bs.peaks = peaks;
		bs.num_peaks = bsr.read_basic_spec(config,fm,bs.ssf,bs.peaks);
		s.init_from_QCPeaks(config,bs.peaks,bs.num_peaks,ssf);


		// create the correct pm_graph, test for presence of optimal solutions
		opt_prm.clear();
		opt_prm.create_graph_from_spectrum(model,&s,s.get_true_mass_with_19(),s.getCharge());
		SeqPath longest = opt_prm.get_longest_subpath(ssf->peptide,0);

	//	if (longest.get_num_aa()<6)
	//		continue;

		bool has_complete_path = false;
		if (longest.get_num_aa() == correct_pep_length)
			has_complete_path = true;

		// generate de novo solutions
		vector<PmcSqsChargeRes> pmc_sqs_res;
		vector<mass_t> pms_with_19;
		vector<int>    charges;
		model->select_pms_and_charges(config,bs,pms_with_19,charges);
		

		if (prm_ptrs.size()<pms_with_19.size())
			prm_ptrs.resize(pms_with_19.size(),NULL);
		
		vector<SeqPath> solutions;
		generate_denovo_solutions_from_several_pms(
				prm_ptrs,
				model,
				&s,
				true, 
				num_solutions,
				6,
				14,
				pms_with_19,
				charges,
				solutions,
				false);

		vector<mass_t> exp_cut_masses;
		full_pep.calc_expected_breakage_masses(config, exp_cut_masses);


//		if (solutions.size()>0)
//			num_pred_aa += solutions[0].get_num_aa();

		int min_denovo_correct_rank=NEG_INF;
		int min_denovo_cut_correct_rank=POS_INF;

		int s_idx;
		bool has_correct_path_in_results=false;
		for (s_idx=0; s_idx<solutions.size(); s_idx++)
		{
			if (solutions[s_idx].check_if_correct(corr_pep_str,config))
			{
				if (min_denovo_correct_rank<0)
					min_denovo_correct_rank=s_idx;

				solutions[s_idx].is_correct = true;
			}
			else
				solutions[s_idx].is_correct = false;

			if (solutions[s_idx].check_if_cut_correct(exp_cut_masses, tolerance))
			{
				if (min_denovo_cut_correct_rank == POS_INF)
					min_denovo_cut_correct_rank=s_idx;

				solutions[s_idx].is_cut_correct = true;
			}
			else
				solutions[s_idx].is_correct = false;

			solutions[s_idx].org_rank = s_idx;
		}

		// only include spectra that have a good path in the top xxx
	//	if (min_denovo_correct_rank==NEG_INF)
	//		continue;

		// if the correct path is not in the solutions, look for one in the prm graphs and add it
		if (min_denovo_correct_rank==NEG_INF)
		{
			cout << endl << spec_idx << " skipped..." << endl;
			if (solutions.size()>0)
			{
				total_pred_rerank[0]+=solutions[0].get_num_aa();
				total_pred_denovo[0]+=solutions[0].get_num_aa();
				int corr=solutions[0].get_num_correct_aas(ssf->peptide,config);
				total_corr_rerank[0]+=corr;
				total_corr_denovo[0]+=corr;
			}
			continue;
		}

		cout << endl << spec_idx << "\tDNV rank: " << min_denovo_correct_rank << endl;



		PeptideRankScorer *drs = (PeptideRankScorer *)model->get_rank_model_ptr(1);
		vector<score_pair> scores;
		drs->score_denovo_sequences(solutions,ssf,peaks,bs.num_peaks,scores,-1);

		sort(scores.begin(),scores.end());

		int min_rerank_correct=NEG_INF;
		int min_rerank_cut_correct=POS_INF;

		int j;
		for (j=0; j<scores.size(); j++)
		{
			int org_rank = scores[j].idx;
			if (min_rerank_correct<0 && solutions[org_rank].is_correct)
			{
				min_rerank_correct = j;
			}

			if (min_rerank_cut_correct==POS_INF && solutions[org_rank].is_cut_correct)
			{
				min_rerank_cut_correct = j;
			}

			if (j<10 || min_rerank_correct == j)
			{
				cout << j << "\t" << org_rank << "\t" << scores[j].score << endl;
			}
		}
		cout << "Best ranks:  Dnv " << min_denovo_correct_rank << "  Rnk " << min_rerank_correct << endl;
		cout << endl;

		if (min_rerank_correct == NEG_INF)
			min_rerank_correct = POS_INF;

		if (min_rerank_cut_correct > min_rerank_correct)
			min_rerank_cut_correct = min_rerank_correct;

		if (min_denovo_correct_rank == NEG_INF)
			min_denovo_correct_rank = POS_INF;

		if (min_denovo_cut_correct_rank > min_denovo_correct_rank)
			min_denovo_cut_correct_rank = min_denovo_correct_rank;

		top_denovo[0].push_back(min_denovo_correct_rank);
		top_rerank[0].push_back(min_rerank_correct);

		top_denovo[1].push_back(min_denovo_cut_correct_rank);
		top_rerank[1].push_back(min_rerank_cut_correct);


		if (min_denovo_correct_rank==0)
			num_top_denovo_correct++;
		if (min_rerank_correct==0)
			num_top_rerank_correct++;


		int num_corr_denovo=solutions[0].get_num_correct_aas(ssf->peptide,config);
		int num_pred_denovo=solutions[0].get_num_aa();
		int num_corr_rerank = solutions[scores[0].idx].get_num_correct_aas(ssf->peptide,config);
		int num_pred_rerank = solutions[scores[0].idx].get_num_aa();

		total_corr_denovo[0] += num_corr_denovo;
		total_pred_denovo[0] += num_pred_denovo;
		total_corr_rerank[0] += num_corr_rerank;
		total_pred_rerank[0] += num_pred_rerank;

		num_spec_tested++;

		if (num_spec_tested % 10 == 0)
		{
			double curr_t = time(NULL);
			sum_stream << num_spec_tested << "\t" << curr_t - start_t << "\t" <<  fixed << 
				setprecision(3) << (float)num_top_denovo_correct/spec_idx <<
				"\t" << num_top_rerank_correct/(float)spec_idx << "\t(" << spec_idx << ")" <<  endl;
		}

		if (num_spec_tested == max_num_spectra)
			break;
	}

	vector<int> b_vals;
	b_vals.push_back(1);
	b_vals.push_back(2);
	b_vals.push_back(5);
	b_vals.push_back(10);
	b_vals.push_back(20);
	b_vals.push_back(50);
	b_vals.push_back(100);
	b_vals.push_back(200);
	b_vals.push_back(500);
	b_vals.push_back(1000);
	b_vals.push_back(2000);
	b_vals.push_back(5000);
	b_vals.push_back(10000);

	int b=b_vals.size()-1;
	while (b>0 && b_vals[b]>=num_solutions)
	{
		b_vals.pop_back();
		b--;
	}
	b_vals.push_back(num_solutions);

	int rep;
	for (rep=0; rep<1; rep++)
	{
		if (rep==0)
		{
			sum_stream << endl << "De novo correctness results" ;
		}
		else
		{
			sum_stream << "Cut correctness results";
		}

		sum_stream << " for " << fixed << setprecision(0) << mgf_test_file << endl << endl;

		vector<int> d_counts,r_counts;
		r_counts.resize(b_vals.size()+1,0);
		d_counts.resize(b_vals.size()+1,0);

		int i;
		for (i=0; i<top_denovo[rep].size(); i++)
		{
			int b;
			int d_rank = top_denovo[rep][i];
			for (b=0; b<b_vals.size(); b++)
				if (d_rank<b_vals[b])
					break;
			int k;
			for (k=d_counts.size()-1; k>=b; k--)
				d_counts[k]++;

			int r_rank = top_rerank[rep][i];
			for (b=0; b<b_vals.size(); b++)
				if (r_rank<b_vals[b])
					break;
			for (k=d_counts.size()-1; k>=b; k--)
				r_counts[k]++;
		}

	//	const double total = top_denovo[rep].size();
		const double total = all_ssf.size();

		sum_stream << endl << "results for " << setprecision(0) << total << " test spectra" << endl;
		for (i=0; i<b_vals.size(); i++)
		{
			sum_stream << setprecision(0) << b_vals[i] << "\t" << setprecision(3) << d_counts[i] << "\t" << (double)d_counts[i]/total << "\t" <<
				r_counts[i] << "\t" << (double)r_counts[i]/total << endl;
		}
		sum_stream << endl;

		if (rep == 0)
		{
			sum_stream << "reg  denovo: aa correct " << 100*total_corr_denovo[rep]/total_pred_denovo[rep] << "% ( avg pred " <<
				total_pred_denovo[rep] / total << " aa)" << endl;
			sum_stream << "rank denovo: aa correct " << 100*total_corr_rerank[rep]/total_pred_rerank[rep]<< "% ( avg pred " <<
				total_pred_rerank[rep]/ total << " aa)" << endl;
		}

		for (i=0; i<b_vals.size(); i++)
			sum_stream << b_vals[i] <<"\t";
		sum_stream << endl;

		for (i=0; i<b_vals.size(); i++)
			sum_stream <<(double)d_counts[i]/total <<"\t";
		sum_stream << endl;

		for (i=0; i<b_vals.size(); i++)
			sum_stream <<(double)r_counts[i]/total <<"\t";
		sum_stream << endl;
		
	}
	if (out_file_stream.is_open())
		out_file_stream.close();
}



/****************************************************************************
Runs benchmark on spectra with a complete de novo solution in there
*****************************************************************************/
void make_ranking_examples(AllScoreModels *model,
						   char *mgf_test_file)
{
	const int charge = 2;

	Config *config;

	FileManager fm;
	FileSet		fs;
	vector<PrmGraph *> prm_ptrs;
	PrmGraph opt_prm;

	config = model->get_config();
	config->set_use_spectrum_charge(2);
	const mass_t tolerance = config->getTolerance();


	fm.init_from_file(config,mgf_test_file);
	fs.select_all_files(fm);
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	

	BasicSpecReader bsr;
	QCPeak peaks[4000];

	cout << "Test file     : " << mgf_test_file << endl;
	
	double start_t = time(NULL);

	double num_pred_aa=0;

	int spec_idx;
	for (spec_idx=0; spec_idx<all_ssf.size(); spec_idx++)
	{
	
		MGF_single *ssf = (MGF_single *)all_ssf[spec_idx];
		const string corr_pep_str = ssf->peptide.as_string(config);
		const Peptide& full_pep = ssf->peptide;
		const mass_t true_mass_with_19 = full_pep.get_mass_with_19();
		const int correct_pep_length = ssf->peptide.get_num_aas();
		BasicSpectrum bs;
		Spectrum s;

		// exclude pepitdes with M+16/Q-17
		const vector<int>& aas = ssf->peptide.get_amino_acids();
		int a;
		for (a=0; a<aas.size(); a++)
			if (aas[a]>Val)
				break;
		if (a<aas.size())
			continue;


		bs.ssf = all_ssf[spec_idx];
		bs.peaks = peaks;
		bs.num_peaks = bsr.read_basic_spec(config,fm,bs.ssf,bs.peaks);
		s.init_from_QCPeaks(config,bs.peaks,bs.num_peaks,ssf);


		// create the correct pm_graph, test for presence of optimal solutions
		opt_prm.clear();
		opt_prm.create_graph_from_spectrum(model,&s,s.get_true_mass_with_19(),s.getCharge());
		SeqPath longest = opt_prm.get_longest_subpath(ssf->peptide,0);

	//	if (longest.get_num_aa()<6)
	//		continue;

		bool has_complete_path = false;
		if (longest.get_num_aa() == correct_pep_length)
			has_complete_path = true;

		// generate de novo solutions
		vector<PmcSqsChargeRes> pmc_sqs_res;
		vector<mass_t> pms_with_19;
		vector<int>    charges;
		model->select_pms_and_charges(config,bs,pms_with_19,charges);
		

		if (prm_ptrs.size()<pms_with_19.size())
			prm_ptrs.resize(pms_with_19.size(),NULL);
		
		vector<SeqPath> solutions;
		generate_denovo_solutions_from_several_pms(
				prm_ptrs,
				model,
				&s,
				true, 
				400,
				6,
				14,
				pms_with_19,
				charges,
				solutions,
				true);

		vector<mass_t> exp_cut_masses;
		full_pep.calc_expected_breakage_masses(config, exp_cut_masses);


//		if (solutions.size()>0)
//			num_pred_aa += solutions[0].get_num_aa();

		int min_denovo_correct_rank=NEG_INF;
		int min_denovo_cut_correct_rank=POS_INF;

		int s_idx;
		bool has_correct_path_in_results=false;
		for (s_idx=0; s_idx<solutions.size(); s_idx++)
		{
			if (solutions[s_idx].check_if_correct(corr_pep_str,config))
			{
				if (min_denovo_correct_rank<0)
					min_denovo_correct_rank=s_idx;

				solutions[s_idx].is_correct = true;
			}
			else
				solutions[s_idx].is_correct = false;

			if (solutions[s_idx].check_if_cut_correct(exp_cut_masses, tolerance))
			{
				if (min_denovo_cut_correct_rank == POS_INF)
					min_denovo_cut_correct_rank=s_idx;

				solutions[s_idx].is_cut_correct = true;
			}
			else
				solutions[s_idx].is_correct = false;

			solutions[s_idx].org_rank = s_idx;
		}

		// only include spectra that have a good path in the top xxx
	//	if (min_denovo_correct_rank==NEG_INF)
	//		continue;

		// if the correct path is not in the solutions, look for one in the prm graphs and add it
		if (min_denovo_correct_rank==NEG_INF)
		{
			cout << endl << spec_idx << " skipped..." << endl;
			continue;
		}

		PeptideRankScorer *drs = (PeptideRankScorer *)model->get_rank_model_ptr(1);
		vector<score_pair> scores;
		drs->score_denovo_sequences(solutions,ssf,peaks,bs.num_peaks,scores,-1);

		sort(scores.begin(),scores.end());

		int min_rerank_correct=NEG_INF;
		int min_rerank_cut_correct=POS_INF;
		int sol_rank=-1;

		int j;
		for (j=0; j<scores.size(); j++)
		{
			int org_rank = scores[j].idx;
			if (min_rerank_correct<0 && solutions[org_rank].is_correct)
			{
				min_rerank_correct = j;
				sol_rank = org_rank;
			}

		}
		


		if (min_rerank_correct == 0 && min_denovo_correct_rank>0)
		{
			vector<SeqPath> example_paths;

			example_paths.push_back(solutions[0]);
			example_paths.push_back(solutions[sol_rank]);

			cout << endl;
			bs.ssf->print_ssf_stats(config);
			cout << "Denovo rank 0 vs. rerank 0 ( original rank " << sol_rank << ")" << endl << endl;
			
			vector<int> aas;
			solutions[0].get_amino_acids(aas);
			Peptide p;
			p.set_peptide_aas(aas);

			cout << "0\tdnv : " << p.as_string(config) << endl;
			cout << sol_rank << "\trnk : " << bs.ssf->peptide.as_string(config) << endl;
			

			drs->list_feature_differences(example_paths,bs.ssf,peaks,bs.num_peaks);
			cout << endl;
		}
	}
}



void PeptideRankScorer::give_de_novo_and_peak_match_examples(
				const string& db_dir,
				const string& correct_dir,
				const string& denovo_dir,
				const string& mgf_list,
				const int charge,
				const int size_idx)
{
	vector<bool>   file_indicators;
	PeptideSetMap	 psm;

	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel *&peak_model = allScoreModels->get_peak_prediction_model_ptr(model_type);
	const vector< vector<mass_t> >& size_thresholds = peak_model->get_size_thresholds();
	

	Config *config = allScoreModels->get_config();

	if (dnv_part_models.size()<=charge)
		init_tables();

	if (! dnv_part_models[charge][size_idx])
		dnv_part_models[charge][size_idx] = new DeNovoPartitionModel;

	vector<int> ppp_frag_type_idxs;
	ppp_frag_type_idxs.clear();
	ppp_frag_type_idxs.push_back(0);
	ppp_frag_type_idxs.push_back(1);
	ppp_frag_type_idxs.push_back(2);
	ppp_frag_type_idxs.push_back(3);

	DeNovoPartitionModel *part_model = dnv_part_models[charge][size_idx];
	part_model->init_features(model_type,charge,size_idx,ppp_frag_type_idxs,config);
	

	// read sample peptide sequences
	create_complete_denovo_set_map(config,mgf_list,db_dir,correct_dir,
		denovo_dir,charge,size_idx,psm, file_indicators);
	

	mass_t min_mz=0;
	mass_t max_mz = 10000;
	if (size_idx>0)
		min_mz = (size_thresholds[charge][size_idx-1]+charge-1)/(mass_t)charge;
	if (size_idx< size_thresholds[charge].size())
		max_mz = (size_thresholds[charge][size_idx]+charge-1)/(mass_t)charge;

//	max_mz = min_mz + 30;
	
	cout << "Charge " << charge << " size " << size_idx << endl;
	cout << "Min m/z " << setprecision(2) << fixed << min_mz  << "  Max m/z " << max_mz << endl;

		// read spectra
	FileManager fm;
	FileSet	    fs;
	BasicSpecReader bsr;
	QCPeak		peaks[4000];

	fm.init_from_list_file(allScoreModels->get_config(),mgf_list.c_str());
	fs.select_files_in_mz_range(fm,min_mz,max_mz,charge);
	
	static vector<PrmGraph *> prm_ptrs;
	static vector<SeqPath> solutions;
	const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();
	int num_spectra_read=0;

	int num_sams=0;

	// Generate various types of samples from spectra
	int i;
	for (i=0; i<all_ssfs.size(); i++)
	{
		PeptideSetMap::const_iterator it;
		MGF_single *ssf = (MGF_single *)all_ssfs[i];
		scan_pair key(ssf->file_idx,ssf->idx_in_file);

		if (ssf->peptide.get_num_aas() != 12)
			continue;

		it = psm.find(key);
		if (it == psm.end())
			continue;

		const PeptideSet& set = (*it).second;
		BasicSpectrum     bs;
		AnnotatedSpectrum as;
		vector<PmcSqsChargeRes> pmc_sqs_res;

		int charge1=0,charge2=0;
		mass_t mz1=0,mz2=0;
		float prob1=0,prob2=0;

		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.num_peaks = num_peaks;
		bs.ssf = ssf;

		allScoreModels->get_best_mz_charge(config,bs, &mz1, &charge1, &prob1,
								  &mz2, &charge2, &prob2, &pmc_sqs_res);

		const mass_t mass_with_19 = mz1*charge1 - (charge1-1.0);
		as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
		as.set_corrected_pm_with_19(mass_with_19);
	
	
	
	
		vector<int> idxs;
		idxs.push_back(part_model->ann_peak_start_idx+6);
		idxs.push_back(part_model->ann_peak_start_idx+7);

		RankBoostSample corr_rbs;
		fill_complete_peptide_rbs(set.correct_sol, peaks, num_peaks, as, pmc_sqs_res, corr_rbs, size_idx);
	
		cout << ++num_sams << 
				"\t" << part_model->feature_names[idxs[0]] << 
				"\t" << part_model->feature_names[idxs[1]] << endl;

		cout << set.correct_sol.type;
		int k;
		for (k=0; k<idxs.size(); k++)
		{
			float val=NEG_INF;
			corr_rbs.get_feature_val(idxs[k],&val);
			cout << "\t" << val;
		}
		cout << "\t" << set.correct_sol.pep.as_string(config) << endl;


		// add db samples
		int num_db=0;
		int num_dnv=0;
		int j;
		for (j=0; j<set.incorrect_sols.size(); j++)
		{
			
			if (set.incorrect_sols[j].type == 1)
			{ 
				num_db++;
				if (num_db>5)
					continue;
			}

			if (set.incorrect_sols[j].type != 1)
			{ 
				num_dnv++;
				if (num_dnv==1)
					cout << endl;
			}
			if (num_dnv>5)
				break;

			RankBoostSample bad_rbs;
			fill_complete_peptide_rbs(set.incorrect_sols[j], peaks, num_peaks, as, pmc_sqs_res, 
				bad_rbs, size_idx);

			cout << set.incorrect_sols[j].type;
			int k;
			for (k=0; k<idxs.size(); k++)
			{
				float val=NEG_INF;
				bad_rbs.get_feature_val(idxs[k],&val);
				cout << "\t" << val + 0.05 * (set.incorrect_sols[j].type == 1 ? num_db : num_dnv);
			}
			cout << "\t" << set.incorrect_sols[j].pep.as_string(config) << endl;

			
		}

		if (num_sams==50)
			break;
	
	}
}



void create_bench_mgf(Config *config, 
					  char *file_list, 
					  char *out_name, 
					  int num_spectra,
					  bool use_exact_pm)
{	
	mass_t min_mz=0;
	mass_t max_mz = 1300;
	int charge =2;

	ofstream ofs(out_name);

	// read spectra
	FileManager fm;
	FileSet	    fs;
	BasicSpecReader bsr;
	QCPeak		peaks[4000];

	fm.init_from_list_file(config,file_list);
	fs.select_files_in_mz_range(fm,min_mz,max_mz,charge);
	fs.randomly_reduce_ssfs(int(num_spectra*1.07));

	const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();
	int spec_idx;
	int num_written=0;
	vector<int> length_counts;
	length_counts.resize(100,0);
	for (spec_idx=0; spec_idx<all_ssfs.size() && num_written<num_spectra; spec_idx++)
	{
		MGF_single *ssf = (MGF_single *)all_ssfs[spec_idx];
		BasicSpectrum     bs;
		AnnotatedSpectrum as;
		PrmGraph prm;

		const vector<int>& aas = ssf->peptide.get_amino_acids();
		int a;
		for (a=0; a<aas.size(); a++)
			if (aas[a]>Val)
				break;
		if (a<aas.size())
			continue;

		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks,false,true);
		bs.peaks = peaks;
		bs.num_peaks = num_peaks;
		bs.ssf = ssf;
		as.init_from_QCPeaks(config,peaks,num_peaks,ssf,false);

		ostringstream oss;
		oss << "spec_" << num_written;
		as.setTitle(oss.str());

		if (use_exact_pm)
		{
			as.set_corrected_pm_with_19(as.get_true_mass_with_19());
			as.set_org_pm_with_19(as.get_true_mass_with_19());
			as.set_m_over_z((as.get_true_mass_with_19()+MASS_PROTON*(as.getCharge()-1))/as.getCharge());
		}

		as.output_as_MGF(ofs);

	//	cout << as.get_org_pm_with_19() << "\t" << as.get_peptide().get_mass_with_19() << "\t" <<
	//		as.get_peptide().as_string(config) << endl;

		num_written++;

		length_counts[as.getPeptide().get_num_aas()]++;
	}

	ofs.close();

	cout << "Wrote " << num_written << " to " << out_name << endl;
	cout << "Histogram: " << endl;
	int i;
	for (i=0; i<length_counts.size(); i++)
	{
		if (length_counts[i]>0)
		{
			cout << i << "\t" << length_counts[i] << "\t" << setprecision(1) << fixed << 100.0*length_counts[i]/(float)num_written << endl;
		}
	}
}


void run_peak_benchmark(AllScoreModels *model, char *benchmark_file)
{
	const int max_rank = 7;
	const int num_ranks_to_consider = 50;

	Config *config = model->get_config();
	PeptideRankScorer *dnv_rank = (PeptideRankScorer *)model->get_rank_model_ptr(1);
	PeakRankModel *prm = model->get_peak_prediction_model_ptr(3);
	FileManager fm;
	FileSet		fs;
	int num_exact_first=0;
	vector<int> correct_rank_counts;
	correct_rank_counts.resize(10,0);

	fm.init_from_file(config,benchmark_file);
	fs.select_all_files(fm);
	fs.randomly_reduce_ssfs(2000);

	
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	vector<QCPeak> peaks;
	peaks.resize(5000);

	BasicSpecReader bsr;
	int i;
	for (i=0; i<all_ssf.size(); i++)
	{
		vector< vector<intensity_t> > ann_intens;
		vector< vector<mass_t> >	  ann_masses;
		AnnotatedSpectrum as;
		Peptide pep = all_ssf[i]->peptide;
		PeptideSolution sol;
		sol.pep = pep;
		sol.reaches_n_terminal=true;
		sol.reaches_c_terminal=true;
		sol.charge = all_ssf[i]->charge;
		sol.pm_with_19 = pep.get_mass_with_19();

		const int num_peaks = bsr.read_basic_spec(config,fm,all_ssf[i],&peaks[0]);
		as.init_from_QCPeaks(config,&peaks[0],num_peaks,all_ssf[i]);
		as.set_peptide(sol.pep);
		as.annotate_spectrum(sol.pm_with_19, sol.charge, true);
		as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

		PeptidePeakPrediction ppp;
		prm->calc_peptide_predicted_scores(sol, ppp);


		// reduce intensities to the same dimensionality
		const int num_frags = ppp.frag_idxs.size();
		vector< vector< float> > observed_intens;
		observed_intens.resize(num_frags);

		int f;
		for (f=0; f<num_frags; f++)
		{
			const int frag_idx = ppp.frag_idxs[f];
			observed_intens[f]=ann_intens[frag_idx]; 
		}

		// calculate the ranks and mapping between predicted and observed
		vector< vector<int> > observed_ranks, predicted_ranks;
		calc_combined_peak_ranks(observed_intens, observed_ranks);
		calc_combined_peak_ranks(ppp.rank_scores, predicted_ranks);

		vector<int> pred2obs, obs2pred;
		vector<int> num_obs_for_frag, num_pred_for_frag;
		vector<float> ordered_scores,  // scores sorted according to their value
					  obs_ordered_scores;  // scores sorted according to the observed intensity rank
		pred2obs.resize(num_ranks_to_consider,999); // look at top 50 peaks
		obs2pred.resize(num_ranks_to_consider,999);
		ordered_scores.resize(num_ranks_to_consider,NEG_INF);
		obs_ordered_scores.resize(num_ranks_to_consider,NEG_INF);
		num_obs_for_frag.resize(num_frags,0);
		num_pred_for_frag.resize(num_frags,0);


		for (f=0; f<num_frags; f++)
		{
			if (observed_ranks[f].size() != predicted_ranks[f].size())
			{
				cout << "#obs  frags: " << observed_ranks.size() << endl;
				cout << "#pred frags: " << predicted_ranks.size() << endl;
				cout << "Error: mismatch in rank dimensionalities!" << endl;
				cout << f << "\tobs : " << observed_ranks[f].size() << "   pred " << predicted_ranks[f].size() << endl;
				exit(1);
			}
			const int num_ranks = predicted_ranks[f].size();
			const vector<float>& frag_rank_scores = ppp.rank_scores[f];
			const vector<float>& frag_intens = observed_intens[f];
			int j;
			for (j=0; j<num_ranks; j++)
			{
				const int obs_rank  = observed_ranks[f][j];
				const int pred_rank = predicted_ranks[f][j];
				const float pred_score = frag_rank_scores[j];

				if (pred_rank<num_ranks_to_consider)
				{
					pred2obs[pred_rank]=obs_rank;
					ordered_scores[pred_rank]=pred_score;
				}

				if (obs_rank<num_ranks_to_consider)
				{
					obs2pred[obs_rank]=pred_rank;
					obs_ordered_scores[obs_rank]=pred_score;
				}

				if (frag_intens[j]>0)
					num_obs_for_frag[f]++;

				if (frag_rank_scores[j]>NEG_INF)
					num_pred_for_frag[f]++;
			}
		}


		int j;
	/*	cout << i << ":\t";
		for (j=0; j<7; j++)
			cout << obs2pred[j] << "\t";
		cout << endl;*/

		if (obs2pred[0] == 0)
		{
			num_exact_first++;
		}

		
		for (j=0; j<max_rank; j++)
		{
			int min=j-3;
			int max=j+3;
		//	if (min<0)
		//		max-=min;

			if (pred2obs[j]>=min && pred2obs[j]<=max)
				correct_rank_counts[j]++;
				
		}
	}

	cout << i;
	cout << setprecision(3) << fixed;
//	cout << num_exact_first/(float)all_ssf.size() << "\t";
	int j;
	for (j=0; j<max_rank; j++)
	{
		cout << " & " << (correct_rank_counts[j]/(float)all_ssf.size());
	}
	cout << "\\\\" << endl;
}


void make_peak_hist_for_obs_rank(AllScoreModels *model, 
								 char *benchmark_file, 
								 int obs_rank)
{
	const int num_ranks_to_consider = 50;
	const int max_hist_rank = 20;

	Config *config = model->get_config();
	PeptideRankScorer *dnv_rank = (PeptideRankScorer *)model->get_rank_model_ptr(1);
	PeakRankModel *prm = model->get_peak_prediction_model_ptr(3);
	FileManager fm;
	FileSet		fs;
	int num_exact_first=0;
	vector<int> rank_hist;
	rank_hist.resize(max_hist_rank+1,0);

	fm.init_from_file(config,benchmark_file);
	fs.select_all_files(fm);
	fs.randomly_reduce_ssfs(10000);

	
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	vector<QCPeak> peaks;
	peaks.resize(5000);

	BasicSpecReader bsr;
	int count =0;
	int i;
	for (i=0; i<all_ssf.size(); i++)
	{
		vector< vector<intensity_t> > ann_intens;
		vector< vector<mass_t> >	  ann_masses;
		AnnotatedSpectrum as;
		Peptide pep = all_ssf[i]->peptide;
		PeptideSolution sol;
		sol.pep = pep;
		sol.reaches_n_terminal=true;
		sol.reaches_c_terminal=true;
		sol.charge = all_ssf[i]->charge;
		sol.pm_with_19 = pep.get_mass_with_19();

		const int num_peaks = bsr.read_basic_spec(config,fm,all_ssf[i],&peaks[0]);
		as.init_from_QCPeaks(config,&peaks[0],num_peaks,all_ssf[i]);
		as.set_peptide(sol.pep);
		as.annotate_spectrum(sol.pm_with_19, sol.charge, true);
		as.extract_annotated_intens_and_masses(ann_intens,ann_masses);

		PeptidePeakPrediction ppp;
		prm->calc_peptide_predicted_scores(sol, ppp);


		// reduce intensities to the same dimensionality
		const int num_frags = ppp.frag_idxs.size();
		vector< vector< float> > observed_intens;
		observed_intens.resize(num_frags);

		int f;
		for (f=0; f<num_frags; f++)
		{
			const int frag_idx = ppp.frag_idxs[f];
			observed_intens[f]=ann_intens[frag_idx]; 
		}

		// calculate the ranks and mapping between predicted and observed
		vector< vector<int> > observed_ranks, predicted_ranks;
		calc_combined_peak_ranks(observed_intens, observed_ranks);
		calc_combined_peak_ranks(ppp.rank_scores, predicted_ranks);

		vector<int> pred2obs, obs2pred;
		vector<int> num_obs_for_frag, num_pred_for_frag;
		vector<float> ordered_scores,  // scores sorted according to their value
					  obs_ordered_scores;  // scores sorted according to the observed intensity rank
		pred2obs.resize(num_ranks_to_consider,999); // look at top 50 peaks
		obs2pred.resize(num_ranks_to_consider,999);
		ordered_scores.resize(num_ranks_to_consider,NEG_INF);
		obs_ordered_scores.resize(num_ranks_to_consider,NEG_INF);
		num_obs_for_frag.resize(num_frags,0);
		num_pred_for_frag.resize(num_frags,0);


		for (f=0; f<num_frags; f++)
		{
			if (observed_ranks[f].size() != predicted_ranks[f].size())
			{
				cout << "#obs  frags: " << observed_ranks.size() << endl;
				cout << "#pred frags: " << predicted_ranks.size() << endl;
				cout << "Error: mismatch in rank dimensionalities!" << endl;
				cout << f << "\tobs : " << observed_ranks[f].size() << "   pred " << predicted_ranks[f].size() << endl;
				exit(1);
			}
			const int num_ranks = predicted_ranks[f].size();
			const vector<float>& frag_rank_scores = ppp.rank_scores[f];
			const vector<float>& frag_intens = observed_intens[f];
			int j;
			for (j=0; j<num_ranks; j++)
			{
				const int obs_rank  = observed_ranks[f][j];
				const int pred_rank = predicted_ranks[f][j];
				const float pred_score = frag_rank_scores[j];

				if (pred_rank<num_ranks_to_consider)
				{
					pred2obs[pred_rank]=obs_rank;
					ordered_scores[pred_rank]=pred_score;
				}

				if (obs_rank<num_ranks_to_consider)
				{
					obs2pred[obs_rank]=pred_rank;
					obs_ordered_scores[obs_rank]=pred_score;
				}

				if (frag_intens[j]>0)
					num_obs_for_frag[f]++;

				if (frag_rank_scores[j]>NEG_INF)
					num_pred_for_frag[f]++;
			}
		}

		int pred_rank = obs2pred[obs_rank];
		if (pred_rank>= max_hist_rank)
		{
			rank_hist[max_hist_rank]++;
		}
		else
			rank_hist[pred_rank]++;

		count++;

	}

	cout << i << endl;
	cout << setprecision(3) << fixed;
//	cout << num_exact_first/(float)all_ssf.size() << "\t";
	int j;
	for (j=0; j<max_hist_rank; j++)
	{
		cout << j+1 << "\t" << rank_hist[j]/(float)count << endl;
	}
	cout << ">" << max_hist_rank << "\t" << rank_hist[max_hist_rank]/(float)count << endl;
}



/***************************************************************************************
This function touches up inspect search results by rescoring the sequences returned by
inspect. The function produces a new inspect results file with the scores (and delta scores)
replaced.
****************************************************************************************/
void PeptideRankScorer::make_peak_table_examples(char *spectra_file) const
{
	AllScoreModels* allScoreModels = static_cast<AllScoreModels*>(allScoreModelsPtr_);
	PeakRankModel*& peakRankModel = allScoreModels->get_peak_prediction_model_ptr(model_type);
	Config *config = allScoreModels->get_config();
	int max_num_examples=100,num_examples=0;

	FileManager fm;
	FileSet     fs;
	fm.init_from_file(config,spectra_file);
	fs.select_all_files(fm);
	const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();

	cout << "Read " <<  all_ssfs.size() << " spectra headers..." << endl;

	BasicSpecReader bsr;
	QCPeak *peaks = new QCPeak[5000];


	int i;
	for (i=0; i<all_ssfs.size(); i++)
	{
		SingleSpectrumFile *ssf = all_ssfs[i];
		
		const int num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		AnnotatedSpectrum as;
		vector< vector< float > > intens;
		vector< vector< mass_t > > masses;

		as.init_from_QCPeaks(config,peaks,num_peaks,ssf);
		as.annotate_spectrum(ssf->peptide.get_mass_with_19());
		as.extract_annotated_intens_and_masses(intens,masses);

		PeptideSolution sol;

		sol.pep = ssf->peptide;
		sol.pm_with_19 = sol.pep.get_mass_with_19();
		sol.charge = ssf->charge;
		sol.reaches_n_terminal = true;
		sol.reaches_c_terminal = true;

		
		if (peakRankModel->make_peak_prediction_table(sol,intens,3))
		{
			cout << i << "\t" << ssf->get_scan() << "\t" <<ssf->peptide.as_string(config) << endl;
			num_examples++;
			cout << endl << endl << endl;
		}


		if (num_examples>max_num_examples)
			break;
	}
}

