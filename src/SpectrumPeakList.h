#ifndef __SPECTRUMPEAKLIST_H__
#define __SPECTRUMPEAKLIST_H__

#include "BasicDataStructs.h"
#include "SingleSpectrumHeader.h"

////////////////////////////////////////////////////////////////////////////
/// This is the basic class for all spectra in PepNovo/MSCluster. This allows
/// for the use of common function in the different applications, despite the
/// big differences in the types of the inherted spectrum classes.
/// It holds a pointer to the spectrum's header, number of peaks
/// and a pointer to the peaks in the spectrum. There is an optional
/// array of peaks that can be used if the peaks are not allocated from a central
/// peak list.
/// A note about object distruction. There is no new type allocation performed  
/// by this object, so none of the pointers should be deleted upon destruction.
class SpectrumPeakList {
public:
	PeakListSpectrum() : header_(NULL), numPeaks_(0), peaks_(NULL) {};

	SingleSpectrumHeader const* getHeader()   const { return header_; }
	SingleSpectrumHeader* getHeader()				{ return header_; }
	Peak*				  getPeaks()    const { return peaks_; }
	int					  getNumPeaks() const { return numPeaks_; }

protected:
	SingleSpectrumHeader* header_;
	int numPeaks_;
	Peak *peaks_;
	vector<Peak> localPeakAllocation;
};




#endif