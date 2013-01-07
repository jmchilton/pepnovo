#ifndef __MSBLAST_H__
#define __MSBLAST_H__

#include <set>
#include "PrmGraph.h"
#include "includes.h"

class AllScoreModels; // fwd dclr

const size_t MIN_MSB_DENOVO_SEQ_LENGTH = 7;
const size_t MAX_MSB_DENOVO_SEQ_LENGTH = 16;


const int B_SYM = 997; // this is R or K
const int X_SYM = 998; // this is a wild card 
const int Z_SYM = 999; // this is Q or K


//
struct MSBSequence {
	MSBSequence() : nGap(0.0), cGap(0.0), rankScore(NEG_INF), msbScore(NEG_INF),
		fileIndex(-1), mz(0.0), scanNumber(-1) {}
				
	bool operator== (const MSBSequence& rhs) const;

	bool operator< (const MSBSequence& rhs) const
	{
		return (msbScore > rhs.msbScore);
	}

	void print(const Config* config) const;

	vector<int>	  seq;
	vector<float> aaProbs;
	vector<bool>  markedAas;
	mass_t	nGap;
	mass_t	cGap;
	score_t rankScore;
	score_t msbScore;
	
	int    fileIndex;
	mass_t mz;
	int    scanNumber;

	void pushToNSide(int aa, float prob, bool mark)
	{
		seq.insert(seq.begin(), aa);
		aaProbs.insert(aaProbs.begin(), prob);
		markedAas.insert(markedAas.begin(), mark);
	}

	void pushToCSide(int aa, float prob, bool mark)
	{
		seq.push_back(aa);
		aaProbs.push_back(prob);
		markedAas.push_back(mark);
	}


	/*! \brief Computes the MS-Blast score for a sequence (approximately number of expected correct amino acids).
	*/
	float calcExpectedMSBScore(const Config* config);


	/*! \brief Computes the optimal alignment with the true peptide. */
	float calcMatchScore(const Config* config, const string& pep) const;


	/*! \brief creates a new sequence by replacing amino acids in the first*/
	void  cloneAndReplace(const MSBSequence& org, size_t orgPos, size_t orgSegmentLength, int *new_aas, size_t newSegmentLength);



	/*! \brief cerates a string of chars from the amino acid ints*/
	string makeSeqString(const Config* config) const;
};



// sets of sequences that should be outputted together to avoid penalty for 
// possible shifts in sequence alignment such as N->GG, Q->DA, etc.
class MSBSequenceSet {
public:

	MSBSequenceSet() :   config_(0) {}



	void convertSeqPathsToMSBSequences(const Config* config, vector<SeqPath>& path, size_t maxNumSequences = 7);

	/*! \brief Creates the output string that is written to the MS-Blast query files.*/
	void createMSBFullLine(string& msbString);

	const vector<MSBSequence>& getMSBSequences() const { return sequences_; }

	/*! \brief adds the sequence only if it is not present (if present and has a better score, replaces the score
		 and source.*/
	void addToSet(const MSBSequence& newSeq);

	void printSet() const;



private:
	vector<MSBSequence> sequences_;

	const Config* config_;

	/*! \brief adds the symbols for B_SYM (with possible gaps) */
	void addMSBlastSymbols();


	/*! \brief Generates all possible msb sequences that can arise from
			   double edges being replaced by single amino acids and vice versa.
	*/
	void generateDoubleAndSingleVariants();

	/*! \brief makes single aa substitutions like Q/K->Z */
	void makeSingleAaSubstitutions();

	/*! \brief replaces stretches of low probability amino acids with the likely number of amino acids XX */
	void replaceLowProbAasWithGaps();

};


struct MSBString {
	MSBString() : seqStr(std::string()), msbScore(-999.0), fileIndex(-1), mz(0.0), scanNumber(-1), indFoundInMap(0) {}

	bool operator< (const MSBString& rhs) const
	{
		return (seqStr.compare(rhs.seqStr)<0);
	}

	bool operator== (const MSBString& rhs) const
	{
		return (seqStr.compare(rhs.seqStr) == 0);
	}

	string seqStr;
	float  msbScore;
	int    fileIndex;
	mass_t mz;
	int    scanNumber;
	bool   indFoundInMap;
};


struct MSBStringSet {
	MSBStringSet() : indMarkedForExclusion(false), indMarkedForKeep(false), indFromMainFullFile(true) {}

	bool indMarkedForExclusion;
	bool indMarkedForKeep;
	bool indFromMainFullFile;

	vector<MSBString> strings;
};




class MSBlastMapFile {
public:

	MSBlastMapFile() : iterationIdx_(-1) {}

	void readFile(const char* file);

	void writeFile(const char* file) const;

	void parseAndAddResultsFile(const char* file);

	const map<string, int>& getProteinNames() const { return proteinNames_; }

	const map<string, vector<int> >& getPeptides() const { return peptides_; }

	int getIterationIdx() const { return iterationIdx_; }

	bool checkIfPeptideInMap(const string& pep) const { return (peptides_.find(pep) != peptides_.end()); }

	bool checkIfProteinInMap(const string& prot) const { return (proteinNames_.find(prot) != proteinNames_.end()); }

private:
	int iterationIdx_;

	map<string, int> proteinNames_;
	map<string, vector<int> > peptides_;
};


class MSBlastKeywordFile {
public:

	bool readKeywordFile(const char* file);

	bool checkForMatch(const string& str) const;

private:
	vector<string> keywords_;
};

class MSBlastIdFile {
public:
	bool readIdFile(const char* file);

	bool checkForIdMatch(const string& str) const;

private:

	set<string> ids_;
};


class MSBlastCollector {
public:

	MSBlastCollector() : config_(0) { msbStringSets.clear(); coveredStrings.clear(); }


	// generates the msb sequences (outputs dnv.txt, msb_full.txt files)
	void generateMsBlastSequences(const char* msb_name,
								  const vector<string>& list_vector, 
								  AllScoreModels *model, 
								  float min_filter_prob,
								  size_t maxSequencesPerSet = 7,
								  bool	indOutputCumulativeProbs = false);


	// generates the msb query (outputs msb_query.txt)
	void writeMsBlastQuery(const char* name,  size_t maxQuerySize, float minMsbScore );




	// generates the msb query (outputs msb_query.txt and updates the "pep_hist.txt" file)
	void generateMsBlastQueryWithHistory(const char* name,  map<string,int>& peptideHistory,
										 size_t maxQuerySize = 137000, float minMsbScore = 3.0);


	// Takes a bunch of MS-Blast "_full.txt" files and creates a single one
	// good for merging parallel jobs or adding to an existing query
	void mergeAndSplitQueries(const char* msb_name,
							  const vector<string>& list_vector,
							  size_t maxSequencesPerSet = 7,
							  size_t maxQuerySize = 500000, 
							  float minMsbScore = 3.0);

	void generateMsBlastFinalQuery(const char* msb_name,
								   const vector<string>& list_vector);


	


	void markAccordingToKeywords(const MSBlastMapFile& mapFile, const set<string>& goodPeptides, const set<string>& badPeptides);

private:

	// adds the unfiltered query file "_full.txt" into the collector
	void readFullFileIntoCollector(const char* file, size_t maxSequencesPerSet, float minMsbScore);

	const Config* config_;
	vector< MSBStringSet >  msbStringSets;
	map<string, size_t>	coveredStrings;  // if string is present holds the set index
	

	void addToExistingSets(MSBSequenceSet& set, size_t maxSequencesPerSet =20);
	void addToExistingSets(const char* line, size_t maxSequencesPerSet =22, float minMsbScore = 3.0);
};





void create_iterative_query(const char* msb_query_name,
							const char* msb_map_file,
							const char* full_file, 
							const char* seconday_file = NULL,
							const char* id_file = NULL,
							const char* keep_file = NULL,
							const char* exclude_file = NULL,
							size_t maxSequencesPerSet = 7,
							size_t maxQuerySize = 137000, 
							float minMsbScore = 3.0);

void examine_msb_results(const char* msb_results_file,
						 const char* msb_map_file,
						 const char* id_file = NULL,
						 const char* keep_file = NULL,
						 const char* exclude_file = NULL);

void create_msb_reprot( const char* list_file,
						const char* msb_query_name,
						const char* map_file,
						const char* msb_id_file,
						const char* msb_keep_file,
						const char* msb_exclude_file);


#endif





