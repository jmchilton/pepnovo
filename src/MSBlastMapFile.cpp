#include "MSBlast.h"
#include "auxfun.h"

void MSBlastMapFile::writeFile(const char* file) const
{
	ofstream ofs(file);
	if (! ofs.good())
		error("Could not open file for writing: ",file);

	ofs << iterationIdx_ << "\t" << proteinNames_.size() << "\t" << peptides_.size() << endl;
	for (map<string, int>::const_iterator it = proteinNames_.begin(); it != proteinNames_.end(); it++)
		ofs << it->second << "\t" << (it->first).substr(0,128) << endl;

	size_t numRefs=0;
	for (map<string, vector<int> >::const_iterator it = peptides_.begin(); it != peptides_.end(); it++)
	{
		ofs << it->first << "\t" << it->second.size();
		for (size_t i=0; i<it->second.size(); i++)
			ofs << "\t" << it->second[i];
		ofs << endl;
		numRefs+= it->second.size();
	}
	ofs.close();

	cout << "Wrote summary of iteration " << iterationIdx_ << " to " << file << endl;
	cout << "Num proteins : " << proteinNames_.size() << endl;
	cout << "Num peptides : " << peptides_.size() << endl;
	cout << "Num indexes  : " << numRefs << endl;
}



void MSBlastMapFile::readFile(const char* file)
{
	ifstream ifs(file);
	if (! ifs.good())
		error("Could not open map file for reading: ",file);

	iterationIdx_=-1;
	size_t numProteins=0, numPeptides=0;
	ifs >> iterationIdx_ >> numProteins >> numPeptides;
	if (iterationIdx_<0 || numProteins<0 || numPeptides<0)
		error("Problem parsing map file: ",file);
	
	string line;
	for (size_t i=0; i<numProteins; i++)
	{
		int idx=-1;
		ifs >> idx;
		getline(ifs, line);
		assert(idx>=0 && line.length()>1);
		proteinNames_[line.substr(1,string::npos)]=idx;
	}

	size_t totalLocs=0;
	for (size_t i=0; i<numPeptides; i++)
	{
		string pep;
		int numLocs = 0;
		ifs >> pep >> numLocs;
		assert(pep.length() && numLocs>0);
		vector<int> locs(numLocs,-1);
		for (size_t j=0; j<numLocs; j++)
		{
			ifs >> locs[j];
			assert(locs[j]>=0);
		}
		totalLocs+=numLocs;
		peptides_[pep]=locs;
	}

	cout << "Read map file : " << file << endl;
	cout << proteinNames_.size() << "\tproteins" << endl;
	cout << peptides_.size()     << "\tpeptides" << endl;
	cout << totalLocs << "\tpeptide-protein connections" << endl;
}


void MSBlastMapFile::parseAndAddResultsFile(const char *file)
{
	const size_t fileSize = getFileSize(file);
	ifstream ifs(file);
	if (! ifs.good())
		error("Could not open MSBlast results file for reading: ",file);

	cout << "Parsing : " << file << " (" << fileSize << " bytes)" << endl;

	char* buffer = new char[fileSize+1024];
	ifs.read(buffer, fileSize);
	buffer[fileSize]='\0';

	// remove terminating strings if they are somewhere in the buffer
	for (size_t i=0; i<fileSize; i++)
		if (buffer[i]=='\0')
			buffer[i]=' ';

	// find the protein results
	const char* alignmentsPtr = strstr(buffer,"Alignments:");
	if (! alignmentsPtr)
		error("Could not find \"Alignments:\" in results file");

	size_t numProteinsExamined=0;
	size_t numProteinsAdded=0;
	size_t numPeptidesExamined=0;
	size_t numPeptidesAdded=0;
	size_t numRefsAdded =0;

	// loop until all proteins are done
	const char* ptr = alignmentsPtr;
	while (ptr)
	{
		const char* protPtr = strstr(ptr, "^ =");
		if (! protPtr)
			break;
		
		const char* lengthPtr = strstr(protPtr,"Length = ");
		if (! lengthPtr)
			error("Could not find length of protein!");

		int protLength = atoi(lengthPtr+9);
		assert(protLength>0 && protLength < 1E7);

		// parse protein name
		assert(lengthPtr - 5 > protPtr);
		string orgName = std::string(protPtr+4, lengthPtr - protPtr - 5);
		string name = std::string();
		// convert white space
		for (size_t i=0; i<orgName.length(); i++)
		{
			if (orgName[i] == '\t' || orgName[i] == '\r' || orgName[i] == '\n')
				orgName[i] = ' ';
			if (i>0 && orgName[i-1]==' ' && orgName[i] == ' ')
				continue;
			name.push_back(orgName[i]);
		}

		// remove trailing whitespace
		assert(name.length()>5);
		size_t last=name.length()-1;
		while (last>0 && name[last] == ' ')
			last--;
		if (last<name.length()-1)
			name.erase(last+1);

	//	cout << ">> " << name << endl;
		
		const char* nextProt = strstr(lengthPtr,"^ =");
		if (! nextProt)
		{
			nextProt = strstr(lengthPtr,"Parameters:");
			if (! nextProt)
				error("Could not find teminating \"Parameters:\"");
		}

		numProteinsExamined++;

		// parse the query hits
		vector<string> sequences;
		const char* p = lengthPtr;
		while (1)
		{
			p=strstr(p,"Query:");
			if (! p || p>=nextProt)
				break;
			
			istringstream iss(std::string(p,64));
			string dummy=std::string(), idx=std::string(), pep=std::string();
			iss >> dummy >> idx >> pep;
			if (pep.length()<3)
			{
				char* pp=const_cast<char*>(nextProt);
				*pp='\0';
				error("Problem parsing line: ",p);
			}
			p+=16;
			numPeptidesExamined++;
			sequences.push_back(pep);
		}
		ptr = lengthPtr;

		// add protein if needed
		name=name.substr(0,128);
		int protIdx=-1;
		map<string, int>::const_iterator it = proteinNames_.find(name);
		if (it != proteinNames_.end())
		{
			protIdx = it->second;
		}
		else
		{
			protIdx = proteinNames_.size();
			proteinNames_[name]=protIdx;
			numProteinsAdded++;
		//	cout << protIdx << "\t(" << sequences.size() << ")\t" << name << endl;
		}

		// add peptides pointers if needed
		for (size_t i=0; i<sequences.size(); i++)
		{
			map<string, vector<int> >::iterator it = peptides_.find(sequences[i]);
			if (it != peptides_.end())
			{
				size_t j;
				for (j=0; j<it->second.size(); j++)
					if (it->second[j] == protIdx)
						break;
				if (j==it->second.size())
				{
					it->second.push_back(protIdx);
					numRefsAdded++;
				}
			}
			else
			{
				peptides_[sequences[i]]=vector<int>(1,protIdx);
				numPeptidesAdded++;
				numRefsAdded++;
			}
		}
	}

	iterationIdx_++;
	cout << "Done parsing MS-Blast results - iteration " << iterationIdx_ << endl;
	cout << "Num proteins examined: " << numProteinsExamined << " (" << numProteinsAdded << " new)" << endl;
	cout << "Num peptides examined: " << numPeptidesExamined << " (" << numPeptidesAdded << " new)" << endl;
	cout << "Num peptide refreneces added: " << numRefsAdded << endl;

	delete [] buffer;
	//-model CID_IT_TRYP -msb_query_name MSB/SPBOTH -msb_create_iterative_query MSB/SPBOTH_L1_full.txt

}


