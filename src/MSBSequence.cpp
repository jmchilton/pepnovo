#include "MSBlast.h"



bool MSBSequence::operator== (const MSBSequence& rhs) const
{
	if (seq.size() != rhs.seq.size())
		return false;

	for (size_t i=0; i<seq.size(); i++)
		if (seq[i] != rhs.seq[i])
			return false;

	return true;
}

/**********************************************************
Calcs the expected score for an MSB seq. Bases calculation
on the probabilities of the predicted amino acids, the 
presence of marked aa's and the presense of problematic
combinations of amino acids. Uses lots of empirically 
derived thresholds.
***********************************************************/
float MSBSequence::calcExpectedMSBScore(const Config *config)
{
	// Emprical probabilities 
	static const float prob_R = 0.80;
	static const float prob_W = 0.70;
	static const float prob_N = 0.96;
	static const float prob_Q = 0.96;
	static const float prob_GV = 0.85;
	static const float prob_VS = 0.95;
	static const float prob_SV = 0.95;
	static const float prob_DA = 0.90;
	static const float prob_AD = 0.95;
	static const float prob_AG = 0.96;
                                        
	static const float prob_F_vs_Met_ox = 0.9; // probs for single aa mixup
	static const float prob_K_first = 0.45;
	static const float prob_K_mid   = 0.25;
	static const float prob_K_last  = 1.0;  

	static const float minProbForAa = 0.275; // if the probability of the amino acid is below this, it will get
										   // a negative score

	static const float multProbForAa = 1.0/(1.0-minProbForAa);

	const int oxidizedMetAA = config->get_aa_from_label("M+16");
	const size_t seqLength = seq.size();

	if (seqLength<MIN_MSB_DENOVO_SEQ_LENGTH)
		return NEG_INF;

	vector<float> prefixScores(seqLength+1, 0.0);
	vector<float> suffixScores(seqLength+1, 0.0);
	
	if (seq[0] == X_SYM)
	{
		prefixScores[0]=0.0;
	}
	else if (seq[0] == Z_SYM)
	{
		prefixScores[0] = aaProbs[0]*0.5;
	}
	else
		prefixScores[0] = aaProbs[0];

	for (size_t i=1; i<seqLength; i++)
	{
		prefixScores[i] = prefixScores[i-1];
		if (seq[i] == X_SYM)
			continue;

		float modProb = (aaProbs[i]- minProbForAa)*multProbForAa;
		if (seq[i]==Z_SYM)
		{
			prefixScores[i]+= modProb*0.5;
		}
		else
			prefixScores[i] += modProb;
	}

	suffixScores[seqLength]=0.0;
	for (int i=seqLength-1; i>=0; i--)
	{
		suffixScores[i] = suffixScores[i+1];
		if (seq[i]==X_SYM)
			continue;

		float modProb = (aaProbs[i]- minProbForAa)*multProbForAa;
		if (seq[i]==Z_SYM)
		{
			suffixScores[i] += modProb * 0.5;
		}
		else
			suffixScores[i] += modProb;
	}

	float simpleMatchScore = prefixScores[seqLength-1];

	// discount occurences of substituteable single amino acids
	// F <=> M* ,  K <=>Q
	// the score adjustment works as follows:
	// with prob 1-k we have the wrong amino acid, we first need to substract
	// the positive score given (aaProbs[i]), and then add the mismatch
	// penalty (+1).
	if (seq[0] == B_SYM && seq[1] ==Lys)
	{
		simpleMatchScore -= ((1-prob_K_first)*(aaProbs[1]+1));
	}
	else if (seq[0] == B_SYM && seq[1] ==Gln)
	{
		simpleMatchScore -= (prob_K_first*(aaProbs[1]+1));
	}

	size_t start = (seq[0] == B_SYM) ? 2 : 0;
	for (size_t i=start; i<seqLength-1; i++)
	{
		if (seq[i] == Lys)
		{
			simpleMatchScore -= ((1-prob_K_mid)*(aaProbs[i]+1));
		}
		else if (seq[i] == Gln)
		{
			simpleMatchScore -= (prob_K_mid*(aaProbs[i]+1));
		}
	}

	if (seq[seqLength-1] == Lys)
	{
		simpleMatchScore -= ((1-prob_K_last)*(aaProbs[seqLength-1]+1));
	}
	else if (seq[seqLength-1] == Gln)
	{
		simpleMatchScore -= (prob_K_last*(aaProbs[seqLength-1]+1));
	}


	if (oxidizedMetAA >0)
	{
		for (size_t i=0; i<seqLength; i++)
		{
			if (seq[i] == oxidizedMetAA)
			{
				simpleMatchScore -= ((prob_F_vs_Met_ox)*(aaProbs[i]+1));
			}
			else if (seq[i] == Phe )
			{
				simpleMatchScore -= ((1-prob_F_vs_Met_ox)*(aaProbs[i]+1));
			}
		}
	}
	
	float min_score = simpleMatchScore;

	// first consider single amino acids that could be replaced by doubles
	// look only in the center, because if it appears on the edge, it can do
	// much damage.
	for (size_t i=1; i< seqLength-1; i++)
	{
		if (seq[i] == X_SYM)
			continue;
		
		float prob = -1.0;

		switch (seq[i])
		{
			case Arg : prob = prob_R; break;
			case Asn : prob = prob_N; break;
			case Gln : prob = prob_Q; break;
			case Trp : prob = prob_W; break;
		}
		
		if (prob <0)
			continue;

		float pre_match = prefixScores[i-1] - seqLength+ i;
		float suf_match = suffixScores[i+1] - i - 1;
		float mismatch_score = (pre_match>suf_match ) ? pre_match : suf_match;
		float exp_score = prob * simpleMatchScore + (1-prob)*mismatch_score;
		if (exp_score<min_score)
			min_score=exp_score;
	}

	// look for double combos that appear within the peptide
	for (size_t i=1; i<seqLength-2; i++)
	{
		if (seq[i]== X_SYM)
			continue;

		float prob=-1.0;

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

		const float pre_match = prefixScores[i-1] - seqLength+ i;
		const float suf_match = suffixScores[i+1] - i - 1;
		float mismatch_score = (pre_match>suf_match ) ? pre_match : suf_match;
		float exp_score = prob * simpleMatchScore + (1-prob)*mismatch_score;
		if (exp_score<min_score)
			min_score=exp_score;


	//	printf("%d : pre score %.2f  suf_score: %.2f\n",i,pre_match,suf_match);
	//	printf("mismatach score: %.2f  prob: %.2f   exp_score %.2f\n",
	//		mismatch_score,prob,exp_score);
	}

	msbScore = min_score;
	return min_score;
}



void MSBSequence::cloneAndReplace(const MSBSequence& org, 
								  size_t orgPos, 
								  size_t orgSegmentLength, 
								  int *new_aas, 
								  size_t newSegmentLength)
{
	assert(orgSegmentLength>0 && newSegmentLength>0);

	int deltaLength = newSegmentLength - orgSegmentLength;
	const int newSize = org.seq.size() + deltaLength;
	const int orgSize = org.seq.size();
	float orgProb = 1.0;

	for (size_t i=0; i<orgSegmentLength; i++)
		if (org.aaProbs[orgPos+i]>0.0)
			orgProb *= org.aaProbs[orgPos+i];

	
	float newProb=pow(static_cast<float>(1.0-orgProb),static_cast<float>(1.0/newSegmentLength));

	if (newProb<0.03)
		newProb = 0.03;

	*this=org;

	seq.resize(newSize);
	aaProbs.resize(newSize);
	markedAas.resize(newSize);
	
	// update swap area
	for (size_t i=0; i<newSegmentLength; i++)
	{
		const size_t pos = i+orgPos;
		seq[pos]=new_aas[i];
		aaProbs[pos]=newProb;
		markedAas[pos]=true; // don't mess with these aas any more
		
		if (new_aas[i] == X_SYM)        // give fixed low prob for X
			aaProbs[i+orgPos] = 0.05;

	}
	
	// add right portion 
	for (size_t i=orgPos+orgSegmentLength; i<orgSize; i++)
	{
		const size_t pos = i+deltaLength;
		seq[pos] = org.seq[i];
		aaProbs[pos]=org.aaProbs[i];
		markedAas[pos] = org.markedAas[i];
	}

	assert(seq.size()>0);
	assert(aaProbs.size() == seq.size());
	assert(markedAas.size() == seq.size());
}




float MSBSequence::calcMatchScore(const Config* config, const string& pep) const
{
	const string msb = makeSeqString(config);

	float maxScore = -999.0;
	
	for (size_t i=0; i<msb.length(); i++)
	{
		float score = 0.0;
		for (size_t j=0; j<pep.length(); j++)
		{
			if (i+j>=msb.length())
				break;
			char m=msb[i+j];
			char p=pep[j];

			if (m==p)
			{
				score+=1.0;
				continue;
			}

			if (m=='X')
				continue;

			if ((m=='I' || m=='L') && (p=='I' || p=='L'))
			{
				score+=1.0;
				continue;
			}
			if (m == 'Z' && (p=='Q' || p=='K'))
			{
				score+=0.5;
				continue;
			}
			score -= 0.5;
		}
		if (score > maxScore)
			maxScore=score;
	}

	for (size_t i=0; i<pep.length(); i++)
	{
		float score = 0.0;
		for (size_t j=0; j<msb.length(); j++)
		{
			if (i+j>=pep.length())
				break;
			char p=pep[i+j];
			char m=msb[j];

			if (m==p)
			{
				score+=1.0;
				continue;
			}
			if ((m=='I' || m=='L') && (p=='I' || p=='L'))
			{
				score+=1.0;
				continue;
			}
			if (m == 'Z' && (p=='Q' || p=='K'))
			{
				score+=0.5;
				continue;
			}
			score -= 0.5;
		}
		if (score > maxScore)
			maxScore=score;
	}
	return maxScore;
}

string MSBSequence::makeSeqString(const Config* config) const
{
	ostringstream oss;
	if (seq.size()<=0)
		return (std::string());

	const vector<char>& aa2Char = config->get_aa2char();
	const vector<int>& orgAAs = config->get_org_aa();
	for (size_t i=0; i<seq.size(); i++)
	{
		if (seq[i] == B_SYM)
		{
			oss << "B";
			continue;
		}

		if (seq[i] == X_SYM)
		{
			oss << "X";
			continue;
		}

		if (seq[i] == Z_SYM)
		{
			oss << "Z";
			continue;
		}

		if (seq[i]>Val)
		{
			oss << aa2Char[orgAAs[seq[i]]];
			oss << ".";
			continue;
		}

		assert(seq[i]>=Ala && seq[i]<=Val);
		oss << aa2Char[seq[i]];
	}

	return oss.str();
}


void MSBSequence::print(const Config* config) const
{
	cout << setprecision(3) << fixed;
	cout << msbScore << "\t" << makeSeqString(config);
	for (size_t j=0; j<aaProbs.size(); j++)
		cout << " " << aaProbs[j];
	cout << "\t";
	for (size_t k=0; k<markedAas.size(); k++)
		if (markedAas[k])
		{
			cout << "1";
		}
		else
			cout << "0";
}







