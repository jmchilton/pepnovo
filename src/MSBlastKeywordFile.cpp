#include "MSBlast.h"
#include "auxfun.h"

bool MSBlastKeywordFile::readKeywordFile(const char *file)
{
	ifstream ifs(file);
	if (! ifs.good())
		return false;

	keywords_.clear();
	char buffer[1024];
	while (! ifs.eof())
	{
		ifs.getline(buffer,1024);
		if (ifs.gcount()>3)
		{
			string str = std::string(buffer);
			for (size_t i=0; i<str.length(); i++)
				str[i] = tolower(str[i]);
			keywords_.push_back(str);
		}
	}
	cout << "Read " << keywords_.size() << " keywords from " << file << endl;
	ifs.close();
	return true;
}

bool MSBlastKeywordFile::checkForMatch(const string& str) const
{
	for (size_t i=0; i<keywords_.size(); i++)
		if (str.find(keywords_[i]) != string::npos)
			return true;
	return false;
}

bool MSBlastIdFile::readIdFile(const char *file)
{
	ifstream ifs(file);
	if (! ifs.good())
		return false;

	ids_.clear();
	char buffer[1024];
	while (! ifs.eof())
	{
		ifs.getline(buffer,1024);
		if (ifs.gcount()>3)
		{
			size_t startIdx=0;
			if (buffer[0] == '>')
				startIdx=1;

			if (buffer[startIdx]=='g' && buffer[startIdx+1] == 'i' && buffer[startIdx+2] == '|')
			{
				string str = std::string(buffer);
				size_t endIdx=str.find_first_of('|',startIdx+3);
				if (endIdx != string::npos)
				{
					ids_.insert(str.substr(startIdx, endIdx-startIdx));
					continue;
				}	
			}
			ids_.insert(std::string(buffer));
		}
	}
	cout << "Read " << ids_.size() << " ids from " << file << endl;
	size_t k=0;
	for (set<string>::const_iterator it=ids_.begin(); it != ids_.end() && k<10; it++,k++)
	{
		cout << *it << endl;
	}
	ifs.close();
	return true;
}

bool MSBlastIdFile::checkForIdMatch(const string& str) const
{
	size_t startIdx=0;
	if (str[0] == '>')
		startIdx=1;

	if (str.length()>3 && str[startIdx]=='g' && str[startIdx+1] == 'i' && str[startIdx+2] == '|')
	{
		size_t endIdx=str.find_first_of('|',startIdx+3);
		if (endIdx != string::npos)
		{
			return (ids_.find(str.substr(startIdx, endIdx-startIdx)) != ids_.end());
		}
	}
	return (ids_.find(str) != ids_.end());
}




