#ifndef __INCLUDES_H__
#define __INCLUDES_H__


#pragma warning (disable:4786)
#pragma warning (disable:4305)
#pragma warning (disable:4503)

#ifdef WIN32
#pragma warning (disable:4996)
#pragma warning (disable:4018)
#pragma warning (disable:4244)
#pragma warning (disable:4267)
#pragma warning (disable:4305)
#endif


#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <time.h>
#include <assert.h>
#include <cstring>



using namespace std;
using std::vector;

#ifdef WIN32
// Improved vector that checks access violations - has overhead use only in DEBUG
template<class T, class A = allocator<T> >
class my_vector : public vector<T,A>
{
public:
	explicit my_vector(const A& al = A()) : vector<T,A> (al) {}
	explicit my_vector(size_type n, const T& v = T(), const A& al = A()) : vector<T,A>(n, v, al) {}
	
	my_vector(const my_vector& x) : vector<T,A>(x) {}
	my_vector(const_iterator first, const_iterator last, const A& al = A()) : vector<T,A>(first, last, al) {}

	const_reference operator[](size_type pos) const
	{
		return at(pos);
	}
	
	reference operator[](size_type pos)
	{
		return at(pos);
	}
};

#ifdef _DEBUG
#define vector my_vector
#endif

#endif // WIN32


#

#define BYTEORDER_LITTLE_ENDIAN

#define NEG_INF -999999999
#define POS_INF  999999999

#define MASS_PROTON 1.00728
#define MASS_ELECTRON 0.00055
#define MASS_ISO 1.0033	 // mass (13C) - mass (12C)
#define MASS_H2O 18.01056 
#define MASS_NH3 17.02655
#define MASS_CO	 27.99492
#define MASS_H2ONH3 35.03711
#define MASS_H2OH2O 36.02112
#define MASS_NH3NH3 34.0531
#define MASS_OHHH   19.0184

#define NUM_SIG_DIGITS 3

#define ESI_MASS_SPEC 1

typedef enum AminoAcids {N_TERM, C_TERM, Gap,Xle,Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,
						 Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val} AminoAcids;

typedef enum InputFileTypes_IFT {
	IFT_DTA, IFT_MGF, IFT_MZXML, IFT_DAT, IFT_PKL, IFT_MS2, IFT_TXT, IFT_ZIP, IFT_NUM_FILE_TYPES
} InputFileTypes_IFT;

// Data types for common variables

typedef float  mass_t;
typedef float  intensity_t;
typedef float  score_t;
typedef float  value_t;
typedef double weight_t;

typedef unsigned int   clusterIdx_t; // used mostly by MsCluster but needed here because used in DatFile
typedef long int  longInt8_t;  // used mostly by MsCluster but needed here because used in DatFile


const unsigned long MAX_ULONG = numeric_limits<unsigned long>::max();
const unsigned int  MAX_UINT  = numeric_limits<unsigned int>::max();
const int			MAX_INT   = numeric_limits<int>::max();
const int			MIN_INT   = numeric_limits<int>::min();
const float			MAX_FLOAT = numeric_limits<float>::max();
const float			MIN_FLOAT = -numeric_limits<float>::max();
const double		MAX_DOUBLE = numeric_limits<double>::max();
const double		MIN_DOUBLE = -numeric_limits<double>::max();
const size_t		MAX_SIZE_T = numeric_limits<size_t>::max();
const value_t		MAX_VALUE_T = numeric_limits<value_t>::max();
const value_t		MIN_VALUE_T = numeric_limits<value_t>::min();
const weight_t		MAX_WEIGHT_T = numeric_limits<weight_t>::max();
const weight_t		MIN_WEIGHT_T = numeric_limits<weight_t>::min();

const float NON_FLOAT   (static_cast<float>(1.1231231e36) - numeric_limits<float>::max()); // my float value that is considered 


typedef map< string , int , less<string> > STRING2INT_MAP;
typedef map< mass_t , mass_t, less<mass_t> > MASS_T_MAP;
typedef map< int , int , less<int> > INT_MAP;
typedef map< mass_t , int   , less<mass_t> > MASS_T2INT_MAP;


#endif


