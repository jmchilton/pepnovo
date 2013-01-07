#include "FragmentSelection.h"
#include "SpectraList.h"

/*********************************************************************
Adds the counts for peaks around a breakage to their respective bins

**********************************************************************/
void add_offset_counts_arround_mass(vector<int>& counts, 
									Spectrum *spec,
									mass_t min_mass, 
									mass_t max_mass, 
									mass_t bin_coef, 
									mass_t breakageMass,
									int charge)
{
	const mass_t low_range = (breakageMass + min_mass+1)/(mass_t)charge;
	const mass_t high_range = (breakageMass + max_mass-1)/(mass_t)charge;
	const PeakRange pr = spec->findPeaksInRange(low_range, high_range);
	
	if (pr.num_peaks<=0)
		return;

	// add counts
	int peakIndex;
	for (peakIndex = pr.low_idx; peakIndex<=pr.high_idx; peakIndex++)
	{
		if (spec->get_peak_iso_level(peakIndex)>0)
		{
			continue;
		}

		mass_t peak_mass = spec->getPeakMass(peakIndex);
		mass_t b_mass = breakageMass/(mass_t)charge;
		int bin = (int)((peak_mass - b_mass - min_mass)*bin_coef);

		if (bin>4 && bin<counts.size()-5)
		{
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
}


/**************************************************************************
Selects the bins that have the highest counts. If a bin is chosen, then
the bins near it are removed from further selection.
***************************************************************************/
void select_fragments_from_bins(vector<int>& counts, 
								FragmentTypeSet& fts, 
								int maxNumFragments, 
								int charge, 
								int orientation, 
								mass_t minimalMassOffset, 
								mass_t bin_coef, 
								mass_t tolerance)
{
	// select top max_numFragments fragment bins
	int erase_bins_size = (int)((bin_coef * 0.75)/charge);
	cout << "EBS: " << erase_bins_size << endl;
	cout << "bin_coef: " <<bin_coef << "  charge " << charge << endl;
	int i;
	for (i=0; i<maxNumFragments; i++)
	{
		int j;
		int max_count =0;
		int max_idx = -1;
		for (j=0; j<counts.size(); j++)
			if (counts[j]>max_count)
			{
				max_count = counts[j];
				max_idx = j;
			}

		if (max_idx<0)
			break;

		FragmentType ft;
		ft.charge=charge;
		ft.orientation = orientation;
		ft.offset      = minimalMassOffset + (static_cast<mass_t>(max_idx) / bin_coef);
		ft.prob = max_count; // for now use the probability
		ft.make_frag_label(tolerance);

		cout << i << "\t" << max_idx << "\t" << counts[max_idx] << "\t" << ft.label << endl;

		fts.add_fragment_type(ft);
	
		for (j=0; j<=erase_bins_size; j++)
		{
			counts[max_idx+j] = -1;
			counts[max_idx-j] = -1;
		}
	}
}


/***********************************************************************
Checks for each selected fragment what is the true count of fragments.
This is done by annotaing the highest count fragments first (and thus
removing fragments that are due to previous/next amino acid peaks).
This function also corrects the offset from the rounded-off value of
the bin to the average of the peak offsets.

If the selected fragment looks like an isotopic peak of a previously 
chosen fragment, it is removed from the list.
************************************************************************/
void calculate_true_fragment_probabilities(
							 const FileManager& fm, 
							 Config *config,
							 FragmentTypeSet& fts,  
							 float min_prob)
{
	const mass_t halfTolerance = config->getTolerance() * 0.5;
	const int numFragments = fts.get_num_fragments();
	FileSet fs;
	vector< vector<double> > trueOffsets;
	vector< vector<double> > totalBreakageCounts; // counts how many breakages are relevant for a charge and the fragment
	vector< vector<double> > perChargeCounts;
	vector< vector<int> > spectraCounts;
	int maxCharge = fm.get_max_charge();
	int c;

	totalBreakageCounts.resize(maxCharge+1);
	spectraCounts.resize(maxCharge+1);
	perChargeCounts.resize(maxCharge+1);
	trueOffsets.resize(maxCharge+1);
	for (c=0; c<=maxCharge; c++)
	{
		perChargeCounts[c].resize(numFragments,0);
		totalBreakageCounts[c].resize(numFragments,0);
		spectraCounts[c].resize(numFragments,0);
		trueOffsets[c].resize(numFragments,0);
	}


	fts.sort_fragments_according_to_probs();
	fs.select_all_files(fm);

//	int f;
//	for (f=0; f<numFragments; f++)
//		fts.get_non_const_fragment(f).spec_count =0;

	while(1)
	{
		Spectrum s;
		vector<mass_t> breakageMasses;
		vector<bool> indUsedPeaks;     // flag to indicate if peaks were already used by another fragment
		mass_t truePeptideMassWith19,truePeptideMass;
		vector<int> indFragmentsFound;

		if (! fs.get_next_spectrum(fm,config,&s))
			break;
		
		s.init_spectrum();
		s.getPeptide().calc_expected_breakage_masses(config,breakageMasses);
		const int charge = s.getCharge();

		truePeptideMass=s.getPeptide().get_mass();
		truePeptideMassWith19 =  truePeptideMass + MASS_OHHH;
		indUsedPeaks.resize(s.getNumPeaks(),false);
		indFragmentsFound.resize(numFragments,0);
	
		// loop on fragments first, so high count fragments get precedence over
		// low count fragments that are actually due to b/y ions of previous or
		// next amino acids
		int f;
		for (f=0; f<numFragments; f++)
		{
			const FragmentType& frag = fts.get_fragment(f);
			if (frag.charge>charge)
				continue;

			int nb=0;
			int b;
			for (b=1; b<breakageMasses.size()-1; b++)
			{
				const mass_t breakageMass = breakageMasses[b];
				const mass_t expectedPeakMass = frag.calc_expected_mass(breakageMass,truePeptideMassWith19);
				if (expectedPeakMass>2000)
					continue;

				if (expectedPeakMass<s.getPeak(0).mass)
					continue;

				nb++;
				const int peakIndex = s.findPeakWithMaxIntensity(expectedPeakMass,halfTolerance);
				if (peakIndex>=0 && ! indUsedPeaks[peakIndex] && s.get_peak_iso_level(peakIndex)==0)
				{
					perChargeCounts[charge][f]++;
					indUsedPeaks[peakIndex]=true;
					mass_t baseMass = ((frag.orientation == PREFIX) ? breakageMass : truePeptideMass - breakageMass);
					baseMass /= frag.charge;

					trueOffsets[charge][f] += (s.getPeakMass(peakIndex)-baseMass);

					indFragmentsFound[f]=1;
				}
			}
			totalBreakageCounts[charge][f]+= nb;
		}

		for (f=0; f<numFragments; f++)
			spectraCounts[charge][f] += indFragmentsFound[f];
	}

	int f;
	for (f=0; f<numFragments; f++)
	{
		FragmentType& frag = fts.get_non_const_fragment(f);

		// find the highest probability for this fragments (amongst all charges)
		frag.prob = -1;
		int c;
		for (c=0; c<=maxCharge; c++)
		{
			if (perChargeCounts[c][f]>0)
			{
				float prob = perChargeCounts[c][f] / totalBreakageCounts[c][f];
				if (prob>frag.prob)
				{
					frag.prob = prob;
					frag.spec_count = spectraCounts[c][f];
					frag.offset  = trueOffsets[c][f] / perChargeCounts[c][f];
				}

			}
		}
		frag.make_frag_label(config->getTolerance());
	}


	// remove fragments that are isotopic peaks of previously selected fragments
	// (e.g. p+2 which is basicly b+1), this doesn't affect -NH3 losses which
	// might appear as an isotope of -H2O

	fts.remove_isotopic_fragments(config->getTolerance());

}


void calculateTrueFragmentProbabilities(const SpectraAggregator& sa,
										FragmentTypeSet& fts,  
										float minProbability)
{
	const Config* const config = sa.getConfig();
	const mass_t halfTolerance = config->getTolerance() * 0.5;
	const int numFragments = fts.get_num_fragments();
	
	vector< vector<double> > trueOffsets;
	vector< vector<double> > totalBreakageCounts; // counts how many breakages are relevant for a charge and the fragment
	vector< vector<double> > perChargeCounts;
	vector< vector<int> >    spectraCounts;

	const int maxCharge = sa.getMaxSepctrumCharge();
	int c;

	totalBreakageCounts.resize(maxCharge+1);
	spectraCounts.resize(maxCharge+1);
	perChargeCounts.resize(maxCharge+1);
	trueOffsets.resize(maxCharge+1);

	for (c=0; c<=maxCharge; c++)
	{
		perChargeCounts[c].resize(numFragments,0);
		totalBreakageCounts[c].resize(numFragments,0);
		spectraCounts[c].resize(numFragments,0);
		trueOffsets[c].resize(numFragments,0);
	}

	fts.sort_fragments_according_to_probs();

	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();



	int i;
	for (i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);

		Spectrum s;
		if (! s.readSpectrum(sa, header))
			continue;

		vector<mass_t> breakageMasses;
		vector<bool> indUsedPeaks;     // flag to indicate if peaks were already used by another fragment
		mass_t truePeptideMassWith19, truePeptideMass;
		vector<int> indFragmentsFound;

		s.getPeptide().calc_expected_breakage_masses(config, breakageMasses);
		const int charge = s.getCharge();

		truePeptideMass=s.getPeptide().get_mass();
		truePeptideMassWith19 =  truePeptideMass + MASS_OHHH;
		indUsedPeaks.resize(s.getNumPeaks(),false);
		indFragmentsFound.resize(numFragments,0);
	
		// loop on fragments first, so high count fragments get precedence over
		// low count fragments that are actually due to b/y ions of previous or
		// next amino acids
		int f;
		for (f=0; f<numFragments; f++)
		{
			const FragmentType& frag = fts.get_fragment(f);

			if (frag.charge> s.getCharge())
				continue;

			int numBreakagesExamined=0;
			int b;
			for (b=1; b<breakageMasses.size()-1; b++)
			{
				const mass_t breakageMass = breakageMasses[b];
				const mass_t expectedPeakMass = frag.calc_expected_mass(breakageMass, truePeptideMassWith19);
				if (expectedPeakMass > s.getMaximalPeakMassToConsider() )
					continue;

				if (expectedPeakMass < s.get_min_peak_mass() )
					continue;

				numBreakagesExamined++;
				const int peakIndex = s.findPeakWithMaxIntensity(expectedPeakMass, halfTolerance);
				if (peakIndex>=0 && ! indUsedPeaks[peakIndex] && s.get_peak_iso_level(peakIndex)==0)
				{
					perChargeCounts[charge][f]++;
					indUsedPeaks[peakIndex]=true;
					mass_t baseMass = ((frag.orientation == PREFIX) ? breakageMass : truePeptideMass - breakageMass);
					baseMass /= frag.charge;

					trueOffsets[charge][f] += (s.getPeakMass(peakIndex)-baseMass);
					indFragmentsFound[f]=1;
				}
			}
			totalBreakageCounts[charge][f]+= numBreakagesExamined;
		}

		for (f=0; f<numFragments; f++)
			spectraCounts[charge][f] += indFragmentsFound[f];
	}

	int f;
	for (f=0; f<numFragments; f++)
	{
		FragmentType& frag = fts.get_non_const_fragment(f);

		// find the highest probability for this fragments (amongst all charges)
		frag.prob = -1;
		int c;
		for (c=1; c<=maxCharge; c++)
		{
			if (perChargeCounts[c][f]>0)
			{
				float prob = perChargeCounts[c][f] / totalBreakageCounts[c][f];
				if (prob>frag.prob)
				{
					frag.prob = prob;
					frag.spec_count = spectraCounts[c][f];
					frag.offset  = trueOffsets[c][f] / perChargeCounts[c][f];
				}

			}
		}
		frag.make_frag_label(config->getTolerance());
	}


	// remove fragments that are isotopic peaks of previously selected fragments
	// (e.g. p+2 which is basicly b+1), this doesn't affect -NH3 losses which
	// might appear as an isotope of -H2O

	fts.remove_isotopic_fragments(config->getTolerance(),false);
}



bool cmp( int a, int b ) { return a > b; }		

/***********************************************************************
Uses the offset count method to determine which fragments are most likely
to occur.
************************************************************************/
void select_frags_using_frag_offset_counts(const FileManager& fm, 
										   Config *config,
										   FragmentTypeSet& final_frags,
										   float min_frag_prob)
{
	FileSet fs;
	FragmentTypeSet fts;
	const mass_t minimalMassOffset = -50;
	const mass_t maximalMassOffset = 50;
	const mass_t tolerance = config->getTolerance();
	const mass_t binSize = tolerance * 0.1;
	const mass_t bin_coef = 1.0 / binSize;
	const int countVectorSize = (int)((maximalMassOffset - minimalMassOffset + 1) / binSize);
	vector< vector<int> > prefixCounts, suffixCounts; // charge, bin_idx
	int i,c;


	fs.select_files(fm,0,10000,-1,-1,0);

	int maxCharge = fm.get_max_charge();
	
	prefixCounts.resize(maxCharge+1);
	suffixCounts.resize(maxCharge+1);
	for (c=1; c<=maxCharge; c++)
	{
		prefixCounts[c].resize(countVectorSize,0);
		suffixCounts[c].resize(countVectorSize,0);
	}

	
	int numSpectraUsed=0;
	while(1)
	{
		Spectrum s;
		vector<mass_t> breakageMasses;
		mass_t truePeptideMass;

		if (! fs.get_next_spectrum(fm,config,&s))
			break;

		numSpectraUsed++;
		s.init_spectrum();
	//	cout << s.get_file_name() << endl;
	
		s.getPeptide().calc_expected_breakage_masses(config,breakageMasses);
		truePeptideMass = s.getPeptide().get_mass();

		int b;
		for (b=1; b<breakageMasses.size()-1; b++)
		{
			int c;
			mass_t breakageMass = breakageMasses[b];
			for (c=1; c<=maxCharge; c++)
			{
				add_offset_counts_arround_mass(prefixCounts[c], &s,
					minimalMassOffset, maximalMassOffset, bin_coef, breakageMass,c);

				add_offset_counts_arround_mass(suffixCounts[c], &s,
					minimalMassOffset, maximalMassOffset, bin_coef, (truePeptideMass - breakageMass),c);
			}
		}
	}



	cout << "Using " << numSpectraUsed << " spectra for offset counts..." << endl;

	// select 30 top fragments becasue many are likely to be caused by previous/next
	// amino acids and will be later removed
	for (c=0; c<=maxCharge; c++)
	{
		select_fragments_from_bins(prefixCounts[c],fts,20,c,PREFIX,minimalMassOffset,bin_coef,tolerance);
		select_fragments_from_bins(suffixCounts[c],fts,20,c,SUFFIX,minimalMassOffset,bin_coef,tolerance);
	}


	cout << "BEFORE: " << fts.get_fragments().size() << endl;
	calculate_true_fragment_probabilities(fm,config,fts, min_frag_prob);
	cout << "AFTER:  " << fts.get_fragments().size() << endl;

	final_frags.clear_set();

	cout << "Fragments selected from spectra:" << endl;
	for (i=0; i<fts.get_num_fragments(); i++)
	{
		const FragmentType& frag = fts.get_fragment(i);
		if (frag.prob >= min_frag_prob)
		{
			final_frags.add_fragment_type(frag);
			cout << left << setw(3) << i << right << setw(5) << frag.spec_count << " ";
			cout << setw(6) << setprecision(3) << right << frag.prob << " ";
			frag.write_fragment(cout);
		}
	}
	cout << endl;
}



void createFragmentsAccordingToOffsetCounts(const SpectraAggregator& sa, 
										    const Config* const config,
										    FragmentTypeSet& finalFragmentTypes, 
										    float mininmalFragmentProbability)
{
	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();

	
	

	const mass_t minimalMassOffset = -50;
	const mass_t maximalMassOffset = 50;
	const mass_t tolerance = config->getTolerance();
	const mass_t binSize = tolerance * 0.1;

	const mass_t bin_coef = 1.0 / binSize;

	const int countVectorSize = (int)((maximalMassOffset - minimalMassOffset + 1) / binSize);
	vector< vector<int> > prefixCounts, suffixCounts; // charge, bin_idx
	
	const int maxCharge = sa.getMaxSepctrumCharge();
	
	prefixCounts.resize(maxCharge+1);
	suffixCounts.resize(maxCharge+1);

	int c;
	for (c=1; c<=maxCharge; c++)
	{
		prefixCounts[c].resize(countVectorSize,0);
		suffixCounts[c].resize(countVectorSize,0);
	}

	
	int numSpectraUsed=0;
	int i;
	for (i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);

		Spectrum s;
		if (! s.readSpectrum(sa, header))
			continue;
		
		numSpectraUsed++;

		vector<mass_t> breakageMasses;
		mass_t truePeptideMass;

		s.getPeptide().calc_expected_breakage_masses(s.getConfig(), breakageMasses);

		truePeptideMass = s.getPeptide().get_mass();

		int b;
		for (b=1; b<breakageMasses.size()-1; b++)
		{
			int c;
			for (c=1; c<=maxCharge; c++)
			{
				add_offset_counts_arround_mass(prefixCounts[c], &s,
					minimalMassOffset, maximalMassOffset, bin_coef, breakageMasses[b],c);

				add_offset_counts_arround_mass(suffixCounts[c], &s,
					minimalMassOffset, maximalMassOffset, bin_coef, (truePeptideMass - breakageMasses[b]), c);
			}
		}
	}



	cout << "Using: " << numSpectraUsed << " spectra for offset counts..." << endl;

	// select 30 top fragments becasue many are likely to be caused by previous/next
	// amino acids and will be later removed
	FragmentTypeSet fts;
	for (c=1; c<=maxCharge; c++)
	{
		select_fragments_from_bins(prefixCounts[c], fts, 20, c, PREFIX, minimalMassOffset, bin_coef, tolerance);
		select_fragments_from_bins(suffixCounts[c], fts, 20, c, SUFFIX, minimalMassOffset, bin_coef, tolerance);
	}

	calculateTrueFragmentProbabilities(sa, fts, mininmalFragmentProbability);

	finalFragmentTypes.clear_set();

//	cout << "Fragments selected from spectra:" << endl;
	for (i=0; i<fts.get_num_fragments(); i++)
	{
		const FragmentType& frag = fts.get_fragment(i);
		if (frag.prob >= mininmalFragmentProbability)
		{
			finalFragmentTypes.add_fragment_type(frag);
		}
	}
//	cout << endl;
}


// returns the probablility of observing a the different types of
// fragments for different charges/sizes/regions
void collectProbabilitiesOfFragments(const SpectraAggregator& sa, 
									 Config *config,
									 vector< vector< vector<double> > >& fragmentProbabilities, // charge , size, region, frag_idx
									 vector< vector< vector<double> > >& inRangeCounts,
									 vector< vector< vector<int> > >& spectraCounts,
									 double& avgRand, 
									 int& numSpectraUsed, 
									 int charge, 
									 bool verbose)
{
	bool veryVerbose=false;
	SpectraList sl(sa);
	sl.selectHeaders(0,POS_INF,charge,charge);

	const vector < vector< vector< RegionalFragments> > >& regionalFragmentSts = 
														config->get_regional_fragment_sets();

	const vector<FragmentType>& all_fragments = config->get_all_fragments();

	vector<int> num_regions;

	int i;

	if (charge<1)
	{
		cout << "Error: charge must be > 1" << endl;
		exit(1);
	}

	sl.selectHeaders(0, POS_INF, charge, charge);

	numSpectraUsed = sl.getNumHeaders();
	

	// resize the fragmentProbabilities
	fragmentProbabilities.resize(regionalFragmentSts[charge].size());
	inRangeCounts.resize(regionalFragmentSts[charge].size());
	spectraCounts.resize(regionalFragmentSts[charge].size());
	num_regions.resize(regionalFragmentSts[charge].size(),0);

	int j;
	for (j=0; j<regionalFragmentSts[charge].size(); j++)
	{
		int k;
		fragmentProbabilities[j].resize(regionalFragmentSts[charge][j].size());
		inRangeCounts[j].resize(regionalFragmentSts[charge][j].size());
		spectraCounts[j].resize(regionalFragmentSts[charge][j].size());
		num_regions[j] = regionalFragmentSts[charge][j].size();

		for (k=0; k<regionalFragmentSts[charge][j].size(); k++)
		{
			fragmentProbabilities[j][k].resize(all_fragments.size(),0);
			inRangeCounts[j][k].resize(all_fragments.size(),0);
			spectraCounts[j][k].resize(all_fragments.size(),0);

//				cout << "SIZE " << fragmentProbabilities[j][k].size() << endl;
		}
	}
	
	cout << "Collecting fragment statistics from " << sl.getNumHeaders() << " spectra:" << endl;
	// fill in the stats
	avgRand=0;
	for (i=0; i<sl.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
		AnnotatedSpectrum as;
		if (! as.readSpectrum(sa, header))
			continue;

		if (veryVerbose)
		{
			cout << i << "  " << as.getTitle();
		}
		else if (verbose && i>0 && i % 5000 == 0)
		{
			cout << i << "/" << sl.getNumHeaders() << endl;
		}


		as.annotate_spectrum(as.get_true_mass_with_19());

	//	as.print_annotated();
	
		const int charge = as.getCharge();
		const int size_idx = as.get_size_idx();
		const vector<Breakage>& breakages = as.get_breakages();
		const mass_t min_mass = as.get_min_peak_mass() ;
		const mass_t max_mass = as.get_max_peak_mass() ;

		vector< vector<int> > frag_ind;
		frag_ind.clear();
		frag_ind.resize(num_regions[size_idx]);

		int j;
		for (j=0; j<num_regions[size_idx]; j++)
			frag_ind[j].resize(all_fragments.size(),0);

		for (j=0; j<breakages.size(); j++)
		{
			const Breakage& breakage = breakages[j];
		
			int f;
			
			if (veryVerbose)
				cout << " " << breakage.fragments.size();
		//	cout << j << " (" << breakage.breakage_region << ")   ";
			for (f=0; f<breakage.fragments.size(); f++)
			{
				if (breakage.fragments[f].mass>0)
				{
				//	cout << all_fragments[breakage.fragments[f].frag_type_idx].frag_label
				//		 << "    ";
					fragmentProbabilities[size_idx][breakage.region_idx]
						      [breakage.fragments[f].frag_type_idx]++;

					frag_ind[breakage.region_idx][breakage.fragments[f].frag_type_idx]=1;
				}
			}

			
			for (f=0; f<all_fragments.size(); f++)
			{
				mass_t expectedPeakMass = all_fragments[f].calc_expected_mass(breakage.mass,as.get_true_mass_with_19());
				if (expectedPeakMass >= min_mass  && expectedPeakMass <= max_mass)
					inRangeCounts[size_idx][breakage.region_idx][f]++;

			}
		}
		
		if (veryVerbose)
			cout << endl;


		int f;
		for (f=0; f<all_fragments.size(); f++)
		{
			int r;
			for (r=0; r<num_regions[size_idx]; r++)
				spectraCounts[size_idx][r][f]+= frag_ind[r][f];
		}
	
		avgRand += (as.getNumPeaks() * config->getTolerance() *2.0) / (max_mass - min_mass) ;
	}

	const double num_spectra = static_cast<double>(sl.getNumHeaders());
	avgRand /= num_spectra;

	int s;
	for (s=0; s<fragmentProbabilities.size(); s++)
	{
		int r;
		for (r=0; r<fragmentProbabilities[s].size(); r++)
		{
			int f;
			for (f=0; f<fragmentProbabilities[s][r].size(); f++)
				if (fragmentProbabilities[s][r][f]>0)
				{
				//	cout << fragmentProbabilities[c][s][r][f] << "   " << 
				//			inRangeCounts[c][s][r][f] << endl;
					fragmentProbabilities[s][r][f] /= inRangeCounts[s][r][f];
						
				}
		}
	}
}


// returns the probablility of observing a the different types of
// fragments for different charges/sizes/regions
void collect_probs_of_known_fragments(const FileManager& fm, Config *config,
				vector< vector< vector<double> > >& fragmentProbabilities, // charge , size, region, frag_idx
				vector< vector< vector<double> > >& inRangeCounts,
				vector< vector< vector<int> > >& spectraCounts,
				double &avgRand, int &numSpectraUsed, int charge, bool verbose)
{
	bool very_verbose=false;
	FileSet fs;
	const vector < vector< vector< RegionalFragments> > >& regional_fragment_sets = 
														config->get_regional_fragment_sets();

	const vector<FragmentType>& all_fragments = config->get_all_fragments();

	vector<int> num_regions;

	int i;

	if (charge<1)
	{
		cout << "Error: charge must be > 1" << endl;
		exit(1);
	}

	fs.select_files(fm,0,5000,-1,-1,charge);

	numSpectraUsed = fs.get_total_spectra();
	

	// resize the fragmentProbabilities
	fragmentProbabilities.resize(regional_fragment_sets[charge].size());
	inRangeCounts.resize(regional_fragment_sets[charge].size());
	spectraCounts.resize(regional_fragment_sets[charge].size());
	num_regions.resize(regional_fragment_sets[charge].size(),0);

	int j;
	for (j=0; j<regional_fragment_sets[charge].size(); j++)
	{
		int k;
		fragmentProbabilities[j].resize(regional_fragment_sets[charge][j].size());
		inRangeCounts[j].resize(regional_fragment_sets[charge][j].size());
		spectraCounts[j].resize(regional_fragment_sets[charge][j].size());
		num_regions[j] = regional_fragment_sets[charge][j].size();

		for (k=0; k<regional_fragment_sets[charge][j].size(); k++)
		{
			fragmentProbabilities[j][k].resize(all_fragments.size(),0);
			inRangeCounts[j][k].resize(all_fragments.size(),0);
			spectraCounts[j][k].resize(all_fragments.size(),0);

//				cout << "SIZE " << fragmentProbabilities[j][k].size() << endl;
		}
	}
	


	// fill in the stats
	avgRand=0;
	for (i=0; i<fs.get_total_spectra(); i++)
	{
		int j;
		AnnotatedSpectrum as;
		fs.get_next_spectrum(fm,config,&as);

		if (very_verbose)
		{
			cout << i << "  " << as.getTitle();
		}
		else if (verbose && i>0 && i % 1000 == 0)
		{
			cout << i << "/" << fs.get_total_spectra() << endl;
		}


		as.annotate_spectrum(as.get_true_mass_with_19());

	//	as.print_annotated();
	
		const int charge = as.getCharge();
		const int size_idx = as.get_size_idx();
		const vector<Breakage>& breakages = as.get_breakages();
		const mass_t min_mass = as.get_min_peak_mass() ;
		const mass_t max_mass = as.get_max_peak_mass() ;

		vector< vector<int> > frag_ind;
		frag_ind.clear();
		frag_ind.resize(num_regions[size_idx]);
		for (j=0; j<num_regions[size_idx]; j++)
			frag_ind[j].resize(all_fragments.size(),0);

		for (j=0; j<breakages.size(); j++)
		{
			const Breakage& breakage = breakages[j];
		
			int f;
			
			if (very_verbose)
				cout << " " << breakage.fragments.size();
		//	cout << j << " (" << breakage.breakage_region << ")   ";
			for (f=0; f<breakage.fragments.size(); f++)
			{
				if (breakage.fragments[f].mass>0)
				{
				//	cout << all_fragments[breakage.fragments[f].frag_type_idx].frag_label
				//		 << "    ";
					fragmentProbabilities[size_idx][breakage.region_idx]
						      [breakage.fragments[f].frag_type_idx]++;

					frag_ind[breakage.region_idx][breakage.fragments[f].frag_type_idx]=1;
				}
			}

			
			for (f=0; f<all_fragments.size(); f++)
			{
				mass_t expectedPeakMass = all_fragments[f].calc_expected_mass(breakage.mass,as.get_true_mass_with_19());
				if (expectedPeakMass >= min_mass  && expectedPeakMass <= max_mass)
					inRangeCounts[size_idx][breakage.region_idx][f]++;

			}
		}
		
		if (very_verbose)
			cout << endl;


		int f;
		for (f=0; f<all_fragments.size(); f++)
		{
			int r;
			for (r=0; r<num_regions[size_idx]; r++)
				spectraCounts[size_idx][r][f]+= frag_ind[r][f];
		}
	
		avgRand += (as.getNumPeaks() * config->getTolerance() *2.0) / (max_mass - min_mass) ;
	}

	const double num_spectra = fs.get_total_spectra();
	avgRand /= num_spectra;

	int s;
	for (s=0; s<fragmentProbabilities.size(); s++)
	{
		int r;
		for (r=0; r<fragmentProbabilities[s].size(); r++)
		{
			int f;
			for (f=0; f<fragmentProbabilities[s][r].size(); f++)
				if (fragmentProbabilities[s][r][f]>0)
				{
				//	cout << fragmentProbabilities[c][s][r][f] << "   " << 
				//			inRangeCounts[c][s][r][f] << endl;
					fragmentProbabilities[s][r][f] /= inRangeCounts[s][r][f];
						
				}
		}
	}
}



struct frag_per {
	bool operator< (const frag_per& other) const
	{
		return (percent > other.percent);
	}
	int frag_idx;
	double percent;
	int in_range_count;
};

/********************************************************************
Determines which fragments should belong in a regional fragment set.
These are the actual fragments that get scored when a breakage falls in
that regions.
*********************************************************************/
void selectRegionalFragmentSets(const SpectraAggregator& sa, 
								Config *config, 
								int charge, 
								bool verbose)
{
	vector< vector< vector<double> > > fragmentProbabilities, inRangeCounts;
	vector< vector< vector<int> > >    spectraCounts;
	double avgRand;
	int    numSpectraUsed;


	collectProbabilitiesOfFragments(sa, config, fragmentProbabilities, inRangeCounts, 
									spectraCounts,avgRand, numSpectraUsed, charge, verbose);

	const double minNumInRange = 0.05 * numSpectraUsed;

	int size_idx;
	for (size_idx=0; size_idx < fragmentProbabilities.size(); size_idx++)
	{
		int region_idx;
		for (region_idx=0; region_idx<fragmentProbabilities[size_idx].size(); region_idx++)
		{
			if (verbose)
			{
				cout << "Charge " << charge << ", Size " << size_idx <<", Region " << 
					region_idx << endl;
				cout << "Random " << fixed << setprecision(4) << avgRand << endl;
			}

			const vector<int>& frag_type_idxs = config->get_regional_fragment_type_idxs(charge,
				size_idx,region_idx);

			vector<score_t> probs;
			probs.resize(frag_type_idxs.size(),0);
			int f;

			for (f=0; f<frag_type_idxs.size(); f++)
			{
				const int frag_type_idx = frag_type_idxs[f];
				if (inRangeCounts[size_idx][region_idx][frag_type_idx]>=minNumInRange)
				{
					probs[f] = fragmentProbabilities[size_idx][region_idx][frag_type_idx];
				}
			}
			
			config->sort_accoriding_to_fragment_probs(probs, charge, size_idx, region_idx);
			config->set_regional_random_probability(charge, size_idx, region_idx, (float)avgRand);
		
			if (verbose)
			{
				const vector<score_t>& fragmentProbabilities = config->get_regional_fragments(charge,size_idx,region_idx).get_frag_probs();
				
				cout << "Fragments observed for charge " << charge << " size " << size_idx << 
						" region " << region_idx << endl;
				for (f=0; f<frag_type_idxs.size(); f++)
				{

					if (inRangeCounts[size_idx][region_idx][frag_type_idxs[f]]>0)
					{
						cout << f<< "\t" << setw(12) << left << config->get_fragment(frag_type_idxs[f]).label <<
							setprecision(5) << setw(10) << config->get_fragment(frag_type_idxs[f]).offset << 
							"   " << fixed << setprecision(5) << fragmentProbabilities[f] << 
							"   " << setprecision(0) << inRangeCounts[size_idx][region_idx][frag_type_idxs[f]] 
							<< " ( " << spectraCounts[size_idx][region_idx][frag_type_idxs[f]] <<  ")"<< endl;
					}
				}
				cout << endl << endl;
			}
		}
	}
}

/********************************************************************
Determines which fragments should belong in a regional fragment set.
These are the actual fragments that get scored when a breakage falls in
that regions.
*********************************************************************/
void select_regional_fragments(const FileManager& fm, Config *config, int charge,
						bool verbose)
{
	vector< vector< vector<double> > > fragmentProbabilities, inRangeCounts;
	vector< vector< vector<int> > >    spectraCounts;
	double avgRand;
	int    numSpectraUsed;

	collect_probs_of_known_fragments(fm, config, fragmentProbabilities, inRangeCounts, spectraCounts,
		avgRand, numSpectraUsed, charge, verbose);
	
	const double min_num_in_range = 0.05 * numSpectraUsed;

	int size_idx;
	for (size_idx=0; size_idx < fragmentProbabilities.size(); size_idx++)
	{
		int region_idx;
		for (region_idx=0; region_idx<fragmentProbabilities[size_idx].size(); region_idx++)
		{
			if (verbose)
			{
				cout << "Charge " << charge << ", Size " << size_idx <<", Region " << 
					region_idx << endl;
				cout << "Random " << fixed << setprecision(4) << avgRand << endl;
			}

			const vector<int>& frag_type_idxs = config->get_regional_fragment_type_idxs(charge,
				size_idx,region_idx);

			vector<score_t> probs;
			probs.resize(frag_type_idxs.size(),0);
			int f;

			cout << "FRAG_TYPES::: " << frag_type_idxs.size() << endl;
			for (f=0; f<frag_type_idxs.size(); f++)
			{
				const int frag_type_idx = frag_type_idxs[f];
				if (inRangeCounts[size_idx][region_idx][frag_type_idx]>=min_num_in_range)
				{
					probs[f] = fragmentProbabilities[size_idx][region_idx][frag_type_idx];
				}
			}
			
			config->sort_accoriding_to_fragment_probs(probs,charge,size_idx,region_idx);
			config->set_regional_random_probability(charge,size_idx,region_idx,(float)avgRand);
		
			if (verbose)
			{
				const vector<score_t>& fragmentProbabilities = config->get_regional_fragments(charge,size_idx,region_idx).get_frag_probs();
				
				for (f=0; f<frag_type_idxs.size(); f++)
				{

					cout << setw(12) << left << config->get_fragment(frag_type_idxs[f]).label <<
							setprecision(5) << setw(10) << config->get_fragment(frag_type_idxs[f]).offset << 
							"   " << fixed << setprecision(5) << fragmentProbabilities[f] << 
							"   " << setprecision(0) << inRangeCounts[size_idx][region_idx][frag_type_idxs[f]] 
							<< " ( " << spectraCounts[size_idx][region_idx][frag_type_idxs[f]] <<  ")"<< endl;
				}
				cout << endl;
			}
		}
	}
}


struct combo_pair {
	bool operator< (const combo_pair& other) const
	{
		return count>other.count;
	}

	int count, f_idx1, f_idx2;
};




void selectFragmentCombinations(const SpectraAggregator& sa, 
								Config* config,
								int charge, 
								int maxNumCombos)
{
	const vector< vector< RegionalFragments> > & regionalFragmentSets = 
												config->get_regional_fragment_sets()[charge];

	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	const int num_all_frags = all_fragments.size();

	if (charge<1)
	{
		cout << "Error: charge must be >= 1" << endl;
		exit(1);
	}

	int j;
	for (j=0; j<regionalFragmentSets.size(); j++)
	{
		int k;
		for (k=0; k<regionalFragmentSets[j].size(); k++)
			config->clear_combos(charge,j,k);
	}
	
	

	int c; // combo counter
	for (c=0; c<maxNumCombos && c<6; c++)
	{
		SpectraList sl(sa);
		sl.selectHeaders(0,POS_INF,charge,charge);

		vector< vector< vector< vector< int > > > > counts;
		counts.resize(regionalFragmentSets.size());
		int j;
		for (j=0; j<regionalFragmentSets.size(); j++)
		{
			int k;
			counts[j].resize(regionalFragmentSets[j].size());
			for (k=0; k<regionalFragmentSets[j].size(); k++)
			{
				int f;
				counts[j][k].resize(num_all_frags);
				for (f=0; f<num_all_frags; f++)
					counts[j][k][f].resize(num_all_frags,0);
			}
		}

//		cout << endl << "Round " << c+1 << endl << endl;

		// check breakage instances
		int i;
		for (i=0; i<sl.getNumHeaders(); i++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
			AnnotatedSpectrum as;
			if (! as.readSpectrum(sa, header))
				continue;

			as.annotate_spectrum(as.get_true_mass_with_19());

			const int size_idx = as.get_size_idx();
			const vector<Breakage>& breakages = as.get_breakages();
			const Peak* const peaks = as.getPeaks();
			const vector<int>& strong_peak_idxs = as.get_strong_peak_idxs();
			vector<int> strong_ind;
			strong_ind.resize(as.getNumPeaks(),0);
			for (j=0; j<strong_peak_idxs.size(); j++)
				strong_ind[strong_peak_idxs[j]]=1;

			for (j=0; j<breakages.size(); j++)
			{
				const Breakage& breakage = breakages[j];
				if (breakage.num_frags_detected<2)
					continue;

				// check if breakage is found by a strong frag

				const int region_idx = breakage.region_idx;
				const RegionalFragments& rf = regionalFragmentSets[size_idx][region_idx];
				int f;

				for (f=0; f<breakage.fragments.size(); f++)
					if (rf.is_a_strong_frag_type(breakage.fragments[f].frag_type_idx) &&
						strong_ind[breakage.fragments[f].peak_idx])
						break;

				if (f<breakage.fragments.size())
					continue;

				// check if breakage is found by an existing combo
				const vector<FragmentCombo>& combos = rf.get_frag_type_combos();
				bool found_by_combo = false;
				if (combos.size()>0)
				{
					for (f=0; f<combos.size(); f++)
					{
						int i;

						for (i=0; i<combos[f].frag_inten_idxs.size(); i++)
							if (breakage.get_position_of_frag_idx(combos[f].frag_inten_idxs[i])<0)
								break;

						if (i<combos[f].frag_inten_idxs.size())
							continue;

						if (combos[f].frag_no_inten_idxs.size() == 0)
						{
							found_by_combo = true;
							break;
						}

						for (i=0; i<combos[f].frag_no_inten_idxs.size(); i++)
							if (breakage.get_position_of_frag_idx(combos[f].frag_no_inten_idxs[i])>0)
								break;

						if (i<combos[f].frag_inten_idxs.size())
							continue;
						
						found_by_combo = true;
						break;
					}
				}
				if (found_by_combo)
					continue;

				// add to count for every pair
				for (f=0; f<breakage.fragments.size()-1; f++)
				{
					int g;
					for (g=f+1; g<breakage.fragments.size(); g++)
					{
						int f_idx = breakage.fragments[f].frag_type_idx;
						int g_idx = breakage.fragments[g].frag_type_idx;

						if (! rf.is_a_strong_frag_type(f_idx) &&
							! rf.is_a_strong_frag_type(g_idx) )
							continue;

						if (f_idx < g_idx)
						{
							counts[size_idx][region_idx][f_idx][g_idx]++;
						//	cout << size_idx << " " << region_idx << " " << f_idx << " " << g_idx <<
						//		" >> " << counts[size_idx][region_idx][f_idx][g_idx] << endl;
						}
						else
						{
							counts[size_idx][region_idx][g_idx][f_idx]++;
						//	cout << size_idx << " " << region_idx << " " << g_idx << " " << f_idx <<
						//		" >> " << counts[size_idx][region_idx][g_idx][f_idx] << endl;
						}
					}
				}
			}
		}

		for (j=0; j<regionalFragmentSets.size(); j++)
		{
			int k;
			for (k=0; k<regionalFragmentSets[j].size(); k++)
			{
				vector<combo_pair> pairs;
				pairs.clear();
				int f,g;
				for (f=0; f<num_all_frags-1; f++)
					for (g=f+1; g<num_all_frags; g++)
						if ( counts[j][k][f][g]>0)
						{
							combo_pair p;
							p.count = counts[j][k][f][g];
							p.f_idx1 = f;
							p.f_idx2 = g;
							pairs.push_back(p);
						}
				sort(pairs.begin(),pairs.end());

				if (pairs.size() == 0)
					continue;

				const RegionalFragments& rf = regionalFragmentSets[j][k];
				FragmentCombo combo;

				if (rf.is_a_strong_frag_type(pairs[0].f_idx1))
				{
					combo.frag_inten_idxs.push_back(pairs[0].f_idx1);
					combo.frag_inten_idxs.push_back(pairs[0].f_idx2);				
				}
				else
				{
					combo.frag_inten_idxs.push_back(pairs[0].f_idx2);
					combo.frag_inten_idxs.push_back(pairs[0].f_idx1);
				}
				config->get_non_const_regional_fragments(charge,j,k).add_combo(combo);

			/*	int i;
				for (i=0; i<pairs.size(); i++)
					cout << pairs[i].count << "  " << config->get_fragment(pairs[i].f_idx1).label
						<< " " << config->get_fragment(pairs[i].f_idx2).label << endl;
				cout << endl;*/
			}
		}
	}
}



// chooses for each regional model the pair of fragments that identify
// the largest number of cuts not already found using strong fragments
void select_frag_combos(const FileManager& fm, 
						Config *config, 
						int charge,
						int max_num_combos)
{
	const vector< vector< RegionalFragments> > & regional_fragment_sets = 
														config->get_regional_fragment_sets()[charge];

	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	const int num_all_frags = all_fragments.size();

	if (charge<1)
	{
		cout << "Error: charge must be > 1" << endl;
		exit(1);
	}

	int j;
	for (j=0; j<regional_fragment_sets.size(); j++)
	{
		int k;
		for (k=0; k<regional_fragment_sets[j].size(); k++)
			config->clear_combos(charge,j,k);
	}
	
	

	int c; // combo counter
	for (c=0; c<max_num_combos && c<6; c++)
	{
		FileSet fs;
		fs.select_files(fm,0,10000,-1,1,charge);

		vector< vector< vector< vector< int > > > > counts;
		counts.resize(regional_fragment_sets.size());
		int j;
		for (j=0; j<regional_fragment_sets.size(); j++)
		{
			int k;
			counts[j].resize(regional_fragment_sets[j].size());
			for (k=0; k<regional_fragment_sets[j].size(); k++)
			{
				int f;
				counts[j][k].resize(num_all_frags);
				for (f=0; f<num_all_frags; f++)
					counts[j][k][f].resize(num_all_frags,0);
			}
		}

//		cout << endl << "Round " << c+1 << endl << endl;

		// check breakage instances
		while(1)
		{
			int j;
			AnnotatedSpectrum as;
			if (! fs.get_next_spectrum(fm,config,&as))
				break;

			as.annotate_spectrum(as.get_true_mass_with_19());

			const int size_idx = as.get_size_idx();
			const vector<Breakage>& breakages = as.get_breakages();
			const Peak* const peaks = as.getPeaks();
			const vector<int>& strong_peak_idxs = as.get_strong_peak_idxs();
			vector<int> strong_ind;
			strong_ind.resize(as.getNumPeaks(),0);
			for (j=0; j<strong_peak_idxs.size(); j++)
				strong_ind[strong_peak_idxs[j]]=1;

			for (j=0; j<breakages.size(); j++)
			{
				const Breakage& breakage = breakages[j];
				if (breakage.num_frags_detected<2)
					continue;

				// check if breakage is found by a strong frag

				const int region_idx = breakage.region_idx;
				const RegionalFragments& rf = regional_fragment_sets[size_idx][region_idx];
				int f;

				for (f=0; f<breakage.fragments.size(); f++)
					if (rf.is_a_strong_frag_type(breakage.fragments[f].frag_type_idx) &&
						strong_ind[breakage.fragments[f].peak_idx])
						break;

				if (f<breakage.fragments.size())
					continue;

				// check if breakage is found by an existing combo
				const vector<FragmentCombo>& combos = rf.get_frag_type_combos();
				bool found_by_combo = false;
				if (combos.size()>0)
				{
					for (f=0; f<combos.size(); f++)
					{
						int i;

						for (i=0; i<combos[f].frag_inten_idxs.size(); i++)
							if (breakage.get_position_of_frag_idx(combos[f].frag_inten_idxs[i])<0)
								break;

						if (i<combos[f].frag_inten_idxs.size())
							continue;

						if (combos[f].frag_no_inten_idxs.size() == 0)
						{
							found_by_combo = true;
							break;
						}

						for (i=0; i<combos[f].frag_no_inten_idxs.size(); i++)
							if (breakage.get_position_of_frag_idx(combos[f].frag_no_inten_idxs[i])>0)
								break;

						if (i<combos[f].frag_inten_idxs.size())
							continue;
						
						found_by_combo = true;
						break;
					}
				}
				if (found_by_combo)
					continue;

				// add to count for every pair
				for (f=0; f<breakage.fragments.size()-1; f++)
				{
					int g;
					for (g=f+1; g<breakage.fragments.size(); g++)
					{
						int f_idx = breakage.fragments[f].frag_type_idx;
						int g_idx = breakage.fragments[g].frag_type_idx;

						if (! rf.is_a_strong_frag_type(f_idx) &&
							! rf.is_a_strong_frag_type(g_idx) )
							continue;

						if (f_idx < g_idx)
						{
							counts[size_idx][region_idx][f_idx][g_idx]++;
						//	cout << size_idx << " " << region_idx << " " << f_idx << " " << g_idx <<
						//		" >> " << counts[size_idx][region_idx][f_idx][g_idx] << endl;
						}
						else
						{
							counts[size_idx][region_idx][g_idx][f_idx]++;
						//	cout << size_idx << " " << region_idx << " " << g_idx << " " << f_idx <<
						//		" >> " << counts[size_idx][region_idx][g_idx][f_idx] << endl;
						}
					}
				}
			}
		}

		for (j=0; j<regional_fragment_sets.size(); j++)
		{
			int k;
			for (k=0; k<regional_fragment_sets[j].size(); k++)
			{
				vector<combo_pair> pairs;
				pairs.clear();
				int f,g;
				for (f=0; f<num_all_frags-1; f++)
					for (g=f+1; g<num_all_frags; g++)
						if ( counts[j][k][f][g]>0)
						{
							combo_pair p;
							p.count = counts[j][k][f][g];
							p.f_idx1 = f;
							p.f_idx2 = g;
							pairs.push_back(p);
						}
				sort(pairs.begin(),pairs.end());

				if (pairs.size() == 0)
					continue;

				const RegionalFragments& rf = regional_fragment_sets[j][k];
				FragmentCombo combo;

				if (rf.is_a_strong_frag_type(pairs[0].f_idx1))
				{
					combo.frag_inten_idxs.push_back(pairs[0].f_idx1);
					combo.frag_inten_idxs.push_back(pairs[0].f_idx2);				
				}
				else
				{
					combo.frag_inten_idxs.push_back(pairs[0].f_idx2);
					combo.frag_inten_idxs.push_back(pairs[0].f_idx1);
				}
				config->get_non_const_regional_fragments(charge,j,k).add_combo(combo);

			/*	int i;
				for (i=0; i<pairs.size(); i++)
					cout << pairs[i].count << "  " << config->get_fragment(pairs[i].f_idx1).label
						<< " " << config->get_fragment(pairs[i].f_idx2).label << endl;
				cout << endl;*/
			}
		}
	}
}





struct f_pair {
	bool operator< (const f_pair& other) const
	{
		return (prob>other.prob);
	}
	int idx;
	score_t prob;
};

void explore_fragment_set(FileManager& fm, Config *config)
{
	FileSet fs;
	vector<FragmentType> all_frags = config->get_all_fragments();
	const int numFragments = all_frags.size();
	const mass_t tolerance = 0.0075;
	vector<double> exp_counts, obs_counts, spectraCounts;

	exp_counts.resize(numFragments,0);
	obs_counts.resize(numFragments,0);
	spectraCounts.resize(numFragments,0);

	
	fs.select_all_files(fm);

	while (1)
	{
		Spectrum s;

		if (! fs.get_next_spectrum(fm,config,&s))
			break;

		vector<mass_t> breakageMasses;

		s.getPeptide().calc_expected_breakage_masses(config,breakageMasses);
		const int num_peaks = s.getNumPeaks();
		vector<int> used_peak_ind;
		vector<int> frag_ind;

		used_peak_ind.resize(num_peaks,0);
		frag_ind.resize(numFragments,0);
		const mass_t min_mass = s.get_min_peak_mass() - tolerance;
		const mass_t max_mass = s.get_max_peak_mass() + tolerance;
		const mass_t truePeptideMassWith19 = s.get_true_mass_with_19();
		// annotate peaks according to their order
		int i;
		for (i=0; i<all_frags.size(); i++)
		{
			const FragmentType& frag = all_frags[i];

			int j;

			for (j=0; j<breakageMasses.size(); j++)
			{
				mass_t expectedPeakMass = frag.calc_expected_mass(breakageMasses[j],truePeptideMassWith19);
				if (expectedPeakMass>min_mass && expectedPeakMass<max_mass)
				{
					exp_counts[i]++;
					int idx = s.findPeakWithMaxIntensity(expectedPeakMass,tolerance);
					if (idx>=0)
					{
						if (! used_peak_ind[idx])
						{
							obs_counts[i]++;
							used_peak_ind[idx]=1;
							frag_ind[i]=1;
						}
					}

				}
			}
		}

		for (i=0; i<all_frags.size(); i++)
			spectraCounts[i]+=frag_ind[i];
	}

	vector<f_pair> pairs;
	int f;
	for (f=0; f<numFragments; f++)
	{
		all_frags[f].prob = obs_counts[f] /exp_counts[f];
		if (all_frags[f].prob<0.001)
			continue;

		f_pair p;
		p.idx=f;
		p.prob = all_frags[f].prob;
		pairs.push_back(p);
	}

	sort(pairs.begin(),pairs.end());

	int i;
	for (i=0; i<pairs.size(); i++)
	{
		int f=pairs[i].idx;

		cout << all_frags[f].label << " & " ;
		cout << all_frags[f].charge<< " & " ;
		string term = "C";
		if (all_frags[f].orientation==PREFIX)
			term = "N";
		cout << term << " & ";
		cout << all_frags[f].offset << " &   & ";
		cout << setprecision(0) << obs_counts[f] << "/" << exp_counts[f] << " & " << spectraCounts[f] << " & ";
		cout << setprecision(3) << all_frags[f].prob << " \\\\" << endl;
		
	}
}


void show_occurences(FileManager& fm, Config *config, string label)
{
	FileSet fs;
	vector<FragmentType> all_frags = config->get_all_fragments();
	const vector<string>& aa2label = config->get_aa2label();
	const int numFragments = all_frags.size();
	const mass_t tolerance = 0.0075;
	int f_idx = config->get_frag_idx_from_label(label);

	fs.select_all_files(fm);

	while (1)
	{
		Spectrum s;

		if (! fs.get_next_spectrum(fm,config,&s))
			break;

		vector<mass_t> breakageMasses;

		s.getPeptide().calc_expected_breakage_masses(config,breakageMasses);
		const int num_peaks = s.getNumPeaks();
		vector<int> used_peak_ind;
		used_peak_ind.resize(num_peaks,0);

		const mass_t truePeptideMassWith19 = s.get_true_mass_with_19();
		// annotate peaks according to their order
		int i;
		for (i=0; i<all_frags.size(); i++)
		{
			const FragmentType& frag = all_frags[i];

			int j;

			for (j=0; j<breakageMasses.size(); j++)
			{
				mass_t expectedPeakMass = frag.calc_expected_mass(breakageMasses[j],truePeptideMassWith19);
				
				int idx = s.findPeakWithMaxIntensity(expectedPeakMass,tolerance);
				if (idx>=0)
				{
					if (! used_peak_ind[idx])
					{
						used_peak_ind[idx]=1;
						
						if (i == f_idx)
						{
							const vector<int>& amino_acids = s.getPeptide().get_amino_acids();
							int k;
							for (k=0; k<j; k++)
								cout << aa2label[amino_acids[k]];

							cout <<" ";
							for ( ; k<amino_acids.size(); k++)
								cout << aa2label[amino_acids[k]];
							cout << endl;
						}
					}
				}
			}
		}
	}
}


void make_frag_rank_histogram(FileManager& fm, Config *config)
{
	FileSet fs;
	vector<FragmentType> all_frags = config->get_all_fragments();
	const mass_t tolerance = 0.45;
	int i,f;

	vector<int> rank_levels;
	rank_levels.push_back(1);
	rank_levels.push_back(2);
	rank_levels.push_back(3);
	rank_levels.push_back(5);
	rank_levels.push_back(7);
	rank_levels.push_back(10);
	rank_levels.push_back(15);
	rank_levels.push_back(20);
	rank_levels.push_back(30);
	rank_levels.push_back(55);
	rank_levels.push_back(99999);

	const int num_levels = rank_levels.size();

	const int numFragments = 6;
	const string labels[numFragments]={"y","b","y-H2O","b-H2O","s2+10.2","b2"};
	//	"b-H2O","y-H2O","y2","b-NH3","a","y2-H2O",
	//	"y2-H2OH2O","b-H2OH2O","y-NH3","y2-H2ONH3",,"b-NH3H2O","a-NH3","a-H2O",
	//	"b-NH3NH3","b2","y-NH3H2O","y-H2OH2O"};

	vector<int> frag_idxs;
	for (f=0; f<numFragments; f++)
		frag_idxs.push_back(config->get_frag_idx_from_label(labels[f]));
	


	vector<double> level_counts;
	vector< vector<double> > frag_level_counts;

	level_counts.resize(num_levels,0);
	frag_level_counts.resize(num_levels);
	for (i=0; i<num_levels; i++)
		frag_level_counts[i].resize(numFragments,0);

	
	
	fs.select_all_files(fm);

	while (1)
	{
		Spectrum s;

		if (! fs.get_next_spectrum(fm,config,&s))
			break;

		vector<mass_t> breakageMasses;

		s.getPeptide().calc_expected_breakage_masses(config,breakageMasses);
		const int num_peaks = s.getNumPeaks();
		vector<int> used_peak_ind;

		used_peak_ind.resize(num_peaks,0);
	
		const mass_t min_mass = s.get_min_peak_mass() - tolerance;
		const mass_t max_mass = s.get_max_peak_mass() + tolerance;
		const mass_t truePeptideMassWith19 = s.get_true_mass_with_19();
		// annotate peaks according to their order
		int i;
		for (i=0; i<frag_idxs.size(); i++)
		{
			const int frag_idx = frag_idxs[i];
			const FragmentType& frag = all_frags[frag_idx];

			int j;

			for (j=0; j<breakageMasses.size(); j++)
			{
				mass_t expectedPeakMass = frag.calc_expected_mass(breakageMasses[j],truePeptideMassWith19);
				if (expectedPeakMass>min_mass && expectedPeakMass<max_mass)
				{
					int idx = s.findPeakWithMaxIntensity(expectedPeakMass,tolerance);
					if (idx>=0)
					{
						if (! used_peak_ind[idx])
						{
							int rank =s.get_peak_rank(idx);
							int l;
							for (l=0; l<num_levels; l++)
								if (rank<=rank_levels[l])
									break;

							frag_level_counts[l][i]++;		
						}
					}

				}
			}
		}

		// add counts for all peaks
		for (i=0; i<s.getNumPeaks(); i++)
		{
			int rank =s.get_peak_rank(i);
			int l;
			for (l=0; l<num_levels; l++)
				if (rank<=rank_levels[l])
					break;
			level_counts[l]++;
		}

	}


	for (i=0; i<numFragments; i++)
	{
		const int frag_idx = frag_idxs[i];
		const FragmentType& frag = all_frags[frag_idx];

		cout << setw(15) << left << frag.label ;
		int l;

		for (l=0; l<num_levels; l++)
		{
	//		cout << " & " << setw(7) << setprecision(3) << frag_level_counts[l][i] / level_counts[l];
			cout << "\t" << setw(7) << setprecision(3) << frag_level_counts[l][i] / level_counts[l];

		}
	
		cout << endl;		
	}

//	cout << setw(15) << left << "Unexplained";
	cout << setw(15) << left << "Other";
	int l;
	for (l=0; l<num_levels; l++)
	{
		int f;
		double exp_sum=0;
		for (f=0; f<numFragments; f++)
			exp_sum += frag_level_counts[l][f];

		double exp_prob = exp_sum/level_counts[l];
		//cout << " & "<< setw(7) << setprecision(3) << 1.0 - exp_prob ;
		cout << "\t"<< setw(7) << setprecision(3) << 1.0 - exp_prob ;
	}
	cout  << endl;
}



// calculates the average random probability according to the offset frequency function
double calc_avg_rand(FileManager& fm, Config *config)
{
	FileSet fs;
	vector<FragmentType> all_frags = config->get_all_fragments();
	const mass_t tolerance = 0.0075;
	int i;
	vector<double> bins;
	double count=0;

	bins.resize(201,0);
	
	fs.select_all_files(fm);

	while (1)
	{
		Spectrum s;

		if (! fs.get_next_spectrum(fm,config,&s))
			break;

		vector<mass_t> breakageMasses;

		s.getPeptide().calc_expected_breakage_masses(config,breakageMasses);
		const int num_peaks = s.getNumPeaks();
		vector<int> used_peak_ind;

		used_peak_ind.resize(num_peaks,0);
	
		const mass_t min_mass = s.get_min_peak_mass() - tolerance;
		const mass_t max_mass = s.get_max_peak_mass() + tolerance;
		const mass_t truePeptideMassWith19 = s.get_true_mass_with_19();
	
		// annotate peaks according to their order
		int i;
		for (i=0; i<all_frags.size() && i<10; i++)
		{
			const FragmentType& frag = all_frags[i];

			int j;

			for (j=0; j<breakageMasses.size(); j++)
			{
				mass_t expectedPeakMass = frag.calc_expected_mass(breakageMasses[j],truePeptideMassWith19);
				if (expectedPeakMass>min_mass && expectedPeakMass<max_mass)
				{
					int idx = s.findPeakWithMaxIntensity(expectedPeakMass,tolerance);
					if (idx>=0)
						used_peak_ind[idx]=1;
				}
			}
		}

		int j;
		for (j=0; j<breakageMasses.size(); j++)
		{
			mass_t expectedPeakMass = breakageMasses[j];

			if (expectedPeakMass>min_mass && expectedPeakMass<max_mass)
			{
				int k;
				for (k=0; k<s.getNumPeaks(); k++)
				{
					mass_t p_mass = s.getPeak(i).mass;
					int b_idx = (int)(expectedPeakMass + 100 - p_mass);
					if (b_idx>=0 && b_idx<=200)
						bins[b_idx]++;
				}
				count++;
			}
		}
	}

	double avg = 0;
	for (i=0; i<bins.size(); i++)
		avg+=bins[i];

	avg/=201;
	avg/=count;

	avg*= config->getTolerance()*2.0;

	cout << "AVG: " << setprecision(5) << avg << endl;

	return avg;
}


void make_bin_histograms(FileManager& fm, Config *config)
{
	FileSet fs;
	vector<FragmentType> all_frags = config->get_all_fragments();
	const mass_t mass_inc = 0.001;
	const mass_t max_off = 20*mass_inc;

	int i;
	vector<double> bins;

	bins.resize(41,0);
	
	fs.select_all_files(fm);

	while (1)
	{
		Spectrum s;

		if (! fs.get_next_spectrum(fm,config,&s))
			break;

		vector<mass_t> breakageMasses;

		s.getPeptide().calc_expected_breakage_masses(config,breakageMasses);
	
		const mass_t min_mass = s.get_min_peak_mass() - 0.1;
		const mass_t max_mass = s.get_max_peak_mass() + 0.1;
		const mass_t truePeptideMassWith19 = s.get_true_mass_with_19();
		// annotate peaks according to their order
		int i;
		for (i=0; i<all_frags.size() && i<2; i++)
		{
			const FragmentType& frag = all_frags[i];

			int j;

			for (j=0; j<breakageMasses.size(); j++)
			{
				mass_t expectedPeakMass = frag.calc_expected_mass(breakageMasses[j],truePeptideMassWith19);
				if (expectedPeakMass>min_mass && expectedPeakMass<max_mass)
				{
					int idx = s.findPeakWithMaxIntensity(expectedPeakMass,max_off);
					if (idx>=0)
					{
						int bin = (int)((expectedPeakMass - s.getPeak(idx).mass)/mass_inc + 20);
						if (bin>=0 && bin<=40)
							bins[bin]++;
					}
				}
			}
		}
	}

	double count = 0;
	for (i=0; i<bins.size(); i++)
		count+=bins[i];

	for (i=0; i<bins.size(); i++)
		cout << setprecision(3) << i*mass_inc - max_off << "  " << bins[i]/count << endl;


}



