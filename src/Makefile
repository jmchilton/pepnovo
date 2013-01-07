SHELL = /bin/sh
CC = g++
CFLAGS = -O3 
LDFLAGS = -lm

OBJS = PepNovo_main.o AllScoreModels.o AnnotatedSpecturm.o auxfun.o base64.o BasicDataStructs.o Config.o\
    ConversionTables.o DeNovoDp.o DeNovoSolutions.o EdgeModel.o FileManagement.o\
    Fragmentation.o FragmentSelection.o FragProbTable.o Isotopes.o MassRangeMap.o ME_REG.o\
    ME_Regression_DataSet.o ME_Regression_Model.o MultiPath.o MZXML_parsing.o\
    PrmGraph.o PTM_stuff.o RegionalPepNovoModel.o AminoAcidProbs.o CumulativeSeqProb.o\
    Spectrum.o PMCSQS.o PMCSQS_model.o MSBSequence.o MSBSequenceSet.o MSBlastCollector.o MSBlastKeywordFile.o MSBlastMapFile.o\
        QCoutputter.o PMC_rank.o RankBoost.o RankBoostSample.o AdvancedScoreModel.o RegularFragmentModel.o StrongFragmentModel.o\
	AdvancedScoreModel_regional.o PeptideComp.o PartitionModel.o PeakRankModel.o PrmNodeScoreModel.o RegionalPrmNodeScoreModel.o\
	PeakRankPartialDenovo.o PeakRank_combined.o PeakRank_advanced.o PeakRank_simple.o PeakRank_combined_dnv.o\
	DeNovoRankTrain.o DeNovoPartModel.o ReScoreDB.o SpectraAggregator.o SpectraFile.o SpectraList.o PeakList.o SingleSpectrumHeader.o\
	DatFile.o PeptideRankScorer.o QuickClusteringSpectra.o QuickClustering.o QCBasicSpecReader.o

HDRS = AllScoreModels.h AnnotatedSpectrum.h auxfun.h base64.h BasicDataStructs.h Config.h ConversionTables.h\
       DeNovoDp.h DeNovoSolutions.h DiscretePeakModel.h EdgeModel.h FileManagement.h Fragmentation.h RegularFragmentModel.h\
	StrongFragmentModel.h FragmentSelection.h FragProbTable.h includes.h Isotopes.h ME_REG.h Model.h PrmNodeScoreModel.h\
       PrmGraph.h RegionalPepNovoModel.h RegionalPrmNodeScoreModel.h  RegularRankModel.h AminoAcidProbs.h CumulativeSeqProb.h\
	SpectraAggregator.h SpectraList.h SpectraFile.h PeakList.h SingleSpectrumHeader.h DatFile.h PeptideRankScorer.h\
       Spectrum.h TagModel.h PMCSQS.h RankBoost.h AdvancedScoreModel.h PeptideComp.h PeakRankModel.h DeNovoRankScore.h\
	QuickClustering.h MSBlast.h 

EXE = PepNovo_bin

.SUFFIXES: .cpp 

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(EXE) $(OBJS)

.cpp.o: 
	$(CC) -c $(CFLAGS)-o $@ $<

clean-objs:
	@rm -f $(OBJS) 

clean-exe:
	@rm -f $(EXE) 

clean: clean-objs clean-exe
