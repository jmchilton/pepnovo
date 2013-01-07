


/*
Copyright 2006, The Regents of the University of California
All Rights Reserved

Permission to use, copy, modify and distribute any part of this 
program for educational, research and non-profit purposes, without fee, 
and without a written agreement is hereby granted, provided that the 
above copyright notice, this paragraph and the following three paragraphs 
appear in all copies.

Those desiring to incorporate this work into commercial 
products or use for commercial purposes should contact the Technology 
Transfer & Intellectual Property Services, University of California, 
San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910, 
Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY 
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, 
INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN 
IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY 
OF SUCH DAMAGE.

THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY 
OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, 
ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO 
REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR 
EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF 
THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/

#include "includes.h"
#include "QuickClustering.h"
#include "base64.h"
#include "auxfun.h"





// Parse an int - skip characters until you see digits or -, then read until you 
// see something else.
int ParseIntFromXML(char* AttributeString)
{
    char Buffer[256];
    int CharCount;
    //
    if (!AttributeString || !*AttributeString)
    {
        return 0;
    }
    CharCount = 0;
    while ((*AttributeString < '0' || *AttributeString > '9') && *AttributeString != '-')
    {
        if (!*AttributeString || CharCount > 256)
        {
            return 0; // too much non-digit garbage!
        }
        AttributeString++;
    }
    CharCount = 0;
    while (*AttributeString >= '0' && *AttributeString <= '9')
    {
        Buffer[CharCount++] = *AttributeString;
        if (CharCount > 10)
        {
            break;
        }
        AttributeString++;
    }
    Buffer[CharCount] = '\0';
    return atoi(Buffer);
}

mass_t ParseMassFromXML(char* AttributeString)
{
    char Buffer[256];
    int CharCount;
    //
    if (!AttributeString || !*AttributeString)
    {
        return 0;
    }
    CharCount = 0;
    while ((*AttributeString < '0' || *AttributeString > '9') && *AttributeString != '-')
    {
        if (!*AttributeString || CharCount > 256)
        {
            return 0; // too much non-digit garbage!
        }
        AttributeString++;
    }
    CharCount = 0;
    while ((*AttributeString >= '0' && *AttributeString <= '9') || *AttributeString == '.')
    {
        Buffer[CharCount++] = *AttributeString;
        if (CharCount > 10)
        {
            break;
        }
        AttributeString++;
    }
    Buffer[CharCount] = '\0';
    return (mass_t)atof(Buffer);
}







/*
// Remove peaks that are not reasonably high for their mass window.
// If WindowWidth and KeepCount are <= 0, use reasonable defaults.
void WindowFilterPeaks(MSSpectrum* Spectrum, float WindowWidth, int KeepCount)
{
    int FilterPeakIndex;
    int NewIndex;
    int OtherPeakIndex;
    float* Intensities;
    int Neighbors;
    float WindowStart;
    float WindowEnd;
    int FilteredCount = 0;
    //
    if (Spectrum->UnfilteredPeaks)
    {
        // We've already performed window filtering; don't do it again!
        return;
    }
    if (WindowWidth <= 0)
    {
        WindowWidth = 50;
    }
    if (KeepCount <= 0)
    {
        KeepCount = 6;
    }

    //
    Intensities = (float*)calloc(Spectrum->PeakCount, sizeof(float));
    for (FilterPeakIndex = 0; FilterPeakIndex < Spectrum->PeakCount; FilterPeakIndex++)
    {
        WindowStart = Spectrum->Peaks[FilterPeakIndex].Mass - (WindowWidth / (float)2.0);
        WindowEnd = Spectrum->Peaks[FilterPeakIndex].Mass + (WindowWidth / (float)2.0);
        Neighbors = 0;
        for (OtherPeakIndex = 0; OtherPeakIndex < Spectrum->PeakCount; OtherPeakIndex++)
        {
            if (Spectrum->Peaks[OtherPeakIndex].Mass > WindowEnd)
            {
                break;
            }
            if (Spectrum->Peaks[OtherPeakIndex].Mass > WindowStart)
            {
                Intensities[Neighbors] = Spectrum->Peaks[OtherPeakIndex].Intensity;
                Neighbors++;
            }
        }
        qsort(Intensities, Neighbors, sizeof(float), (QSortCompare)CompareFloats);
        if (Neighbors < KeepCount || Spectrum->Peaks[FilterPeakIndex].Intensity >= Intensities[KeepCount - 1])
        {
            Spectrum->Peaks[FilterPeakIndex].FilterScore = 1;
            FilteredCount++;
        }
    }
    SafeFree(Intensities);
    // New array:
    Spectrum->UnfilteredPeakCount = Spectrum->PeakCount;
    Spectrum->UnfilteredPeaks = Spectrum->Peaks;
    Spectrum->PeakCount = FilteredCount;
    Spectrum->Peaks = (SpectralPeak*)calloc(FilteredCount, sizeof(SpectralPeak));
    NewIndex = 0;
    for (FilterPeakIndex = 0; FilterPeakIndex < Spectrum->UnfilteredPeakCount; FilterPeakIndex++)
    {
        if (Spectrum->UnfilteredPeaks[FilterPeakIndex].FilterScore)
        {
            memcpy(Spectrum->Peaks + NewIndex, Spectrum->UnfilteredPeaks + FilterPeakIndex, sizeof(SpectralPeak));
            Spectrum->Peaks[NewIndex].Index = NewIndex;
            NewIndex++;
        }
    }

}
*/

/****************************************************************************
// filters the peaks, returns the new number of peaks (that appear
// in the new_peaks buffer
*****************************************************************************/
int join_and_filter_peak_list(const Config *config, 
							  mass_t m_over_z, 
							  float *org_peaks, 
							  int num_org_peaks, 
							  float *new_peaks)
{
	static vector<MassInten> peak_list;
	static int peak_list_size =0;
	int i;

	const float WindowWidth = 50.0;
	const int   KeepCount   = 7;
	int qq = peak_list_size;

	if (num_org_peaks>peak_list_size)
	{
		peak_list_size = (int)(num_org_peaks * 1.5);
		if (peak_list_size<2000)
			peak_list_size = 2000;

		peak_list.resize(peak_list_size);
	}

	// copy org_peaks to the temporary peak_list
	int f_idx=0;
	for (i=0; i<num_org_peaks; i++)
	{
		peak_list[i].mass=org_peaks[f_idx++];
		peak_list[i].intensity=org_peaks[f_idx++];
	}

	vector<bool> keep_inds;
	keep_inds.resize(num_org_peaks,false);

	int start_idx = 0;
	int f;
	for (f=0; f<num_org_peaks; f++)
	{
		if (peak_list[f].intensity<=0)
			continue;

		const mass_t WindowStart = peak_list[f].mass - (WindowWidth / (float)2.0);
        const mass_t WindowEnd =   peak_list[f].mass + (WindowWidth / (float)2.0);
    	
		const float current_inten = peak_list[f].intensity;
		bool first_flag=true;
		int num_above=0;
		int OtherPeakIndex;
        for (OtherPeakIndex = start_idx; OtherPeakIndex < num_org_peaks; OtherPeakIndex++)
        {
            if (peak_list[OtherPeakIndex].mass> WindowEnd)
                break;
         
            if (peak_list[OtherPeakIndex].mass > WindowStart)
			{
                if (peak_list[OtherPeakIndex].intensity>current_inten)
					num_above++;

				if (first_flag)
				{
					start_idx = OtherPeakIndex;
					first_flag = false;
				}
			}
        }

		if (num_above < KeepCount)
			keep_inds[f]=true;
	}


	f_idx=0;
	for (i=0; i<num_org_peaks; i++)
	{
		if (keep_inds[i])
		{
			new_peaks[f_idx++] = peak_list[i].mass;
			new_peaks[f_idx++] = peak_list[i].intensity;
		}

	} 

	return (f_idx/2);
}


// reads the spectrum info from an MZXML file. Assumes the file pointer is 
// at the right position.
int BasicSpecReader::get_peak_list_from_MZXML(FILE *mzxml_stream)
{
    static char* Buffer = NULL;
    char* PeakCountStr;
    int PeakCount;
    char* PeakStr;
    static char* PeakBuffer = NULL;
    static char* DecodedPeakBuffer = NULL;
    static float* Peaks = NULL;
    static int PeakBufferSize = 0;

    int PeakIndex;
    static char* PrecursorStr;
    int FloatIndex;
    char* ByteOrderStr;
    int ByteOrderLittle = 1;

	long start_pos = ftell(mzxml_stream);
	
    //
    if (!Buffer)
    {
        Buffer = (char*)calloc(XML_BUFFER_SIZE + 1, sizeof(char));
    }
    
	char *bf=Buffer;
    fread(Buffer, XML_BUFFER_SIZE, sizeof(char), mzxml_stream);
    PeakCountStr = strstr(Buffer, "peaksCount=\"");
    if (!PeakCountStr)
    {
        cout << "Error parsing peaks from mzxml! " << endl;
        return 0;
    }
    PeakCount = ParseIntFromXML(PeakCountStr);
    if (!PeakCount)
    {
        // A spectrum with zero peaks!  This is aberrant, but it can happen,
        // so bail out politely.  
        cout << "Error: mzXML spectrum contains no peaks!" << endl;
		exit(1);
    }

  
    PeakStr = strstr(Buffer, "<peaks");
    if (PeakStr)
    {
        // Get byte order:
        ByteOrderStr = strstr(PeakStr, "byteOrder=\"");
        if (ByteOrderStr)
        {
            ByteOrderStr += 11;
            if (!strncmp(ByteOrderStr, "network", 7))
            {
                ByteOrderLittle = 0;
            }
            if (!strncmp(ByteOrderStr, "big", 3))
            {
                ByteOrderLittle = 0;
            }
            if (!strncmp(ByteOrderStr, "little", 6))
            {
                ByteOrderLittle = 1;
            }
        }
        PeakStr = strstr(PeakStr, ">");
    }
    if (!PeakStr)
    {
        cout << "Error parsing peaks from mzxml! " << endl;
        return 0;
    }

    PeakStr++;
    if (PeakBufferSize < PeakCount)
    {
        if (PeakBuffer)
        {
            free(PeakBuffer);
            PeakBuffer = NULL;
            free(DecodedPeakBuffer);
            DecodedPeakBuffer = NULL;
            free(Peaks);
            Peaks = NULL;
        }
        PeakBufferSize = PeakCount;
        PeakBuffer = (char*)calloc(PeakBufferSize * 22 + 100, 1);
        DecodedPeakBuffer = (char*)calloc(PeakCount * 8 + 8, 1);
        Peaks = (float*)calloc(PeakCount * 2, sizeof(float));
    }

	char *pb = PeakBuffer;
	int pbs= PeakBufferSize;
    fseek(mzxml_stream, start_pos + (PeakStr - Buffer) +15, 0);
    fread(PeakBuffer, PeakBufferSize * 22 + 100, sizeof(char), mzxml_stream);
    int Trail = (PeakCount % 3);
    if (!(PeakCount % 3))
    {
        PeakBuffer[PeakCount * 32/3] = '\0';
    }
    else
    {
        PeakBuffer[(PeakCount * 32/3) + Trail + 1] = '\0';
    }
   // b64_decode_mio(PeakBuffer, DecodedPeakBuffer);
	b64_decode_mio( DecodedPeakBuffer, PeakBuffer);
    for (FloatIndex = 0; FloatIndex < (2 * PeakCount); FloatIndex++)
    {
#ifdef BYTEORDER_LITTLE_ENDIAN
        if (!ByteOrderLittle)
        {
            char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
            DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
            DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
            ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
            DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
            DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
        }
        memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
#else
        if (ByteOrderLittle)
        {
            char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
            DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
            DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
            ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
            DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
            DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
        }
        memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
#endif
    }
   
	int i_pos=0;
    for (PeakIndex = 0; PeakIndex < PeakCount; PeakIndex++)
    {
		peak_list[PeakIndex].mass = Peaks[i_pos++];
		peak_list[PeakIndex].intensity =  Peaks[i_pos++];

//		cout << PeakIndex << " " << peak_list[PeakIndex].mass << " " << peak_list[PeakIndex].intensity << endl;
    }

	// sanity check
	int i;
	for (i=1; i<PeakCount; i++)
		if (peak_list[i].intensity<0 || peak_list[i].mass < peak_list[i-1].mass)
			break;

	if (i<PeakCount || peak_list[0].mass<0 || peak_list[0].intensity<0)
	{
		cout << "Error parsing peaks in mzXML! i:" <<i << " peakCount:" << PeakCount << endl;
		exit(1);
	}

    return PeakCount;
}




// Iterate over scans from a file in .mzXML format.
// Ignore everything except scans at MS level 2 or higher.
// Make a note of where these scans are in the file, so we can go back and read the peaks later.
// (All we need to do for now is find <scan tags, and check the scan number and msLevel attribute!)


void MZXML_file::initial_read(const Config *config, int file_idx)
//void ParseSpectraFromMZXML(char* FileName, InputFileNode* InputFile, int FirstScan, int LastScan)
{
    int BytesToRead;
    char* Buffer;
    int BufferStartPos = 0;
    int Pos;
    int BytesRead;
    int BufferEnd = 0;
    FILE* MZXMLFile;
    int ParseState = 0;
    int FilePos = 0;
    char* ScanStr;
    char* ScanNumberStr;
    int ScanNumber;
    char* MSLevelStr;
    int MSLevel;
	char *retentionTimeStr;
	float retentionTime;
	char *precursorIntensityStr;
	float precursorIntensity;
	char *PrecursorStr;
	mass_t precursorMZ;

    //
    Buffer = (char*)calloc(XML_BUFFER_SIZE + 1, sizeof(char));
    MZXMLFile = fopen(mzxml_name.c_str(), "rb");
    if (!MZXMLFile)
    {
        cout << "Error: Can't open MZXML file " <<  mzxml_name << endl;
        exit(1);
    }

	bool first_round=true;
    while (1)
    {
		MZXML_single spec_header;
		spec_header.file_idx    = file_idx;

		// Now move the buffer forward a bit
		// put it here to allow continue from any point in the loop (if there is a problem with 
		// parsing the spectrum)
		if (! first_round)
		{
			memmove(Buffer, Buffer + Pos + 10, BufferEnd - (Pos + 10));
			BufferEnd -= Pos + 10;
			FilePos += Pos + 10;
		}
		first_round = false;
		
        // Read more data, to fill up the buffer:
        BytesToRead = XML_BUFFER_SIZE - BufferEnd;
		if (BytesToRead>0)
		{
			BytesRead =   fread(Buffer + BufferEnd, sizeof(char), BytesToRead, MZXMLFile);
			BufferEnd += BytesRead;
		}
        Buffer[BufferEnd] = '\0';

        // Look for a new <scan tag opening:
        ScanStr = strstr(Buffer, "<scan");
        if (ScanStr)
        {
            Pos = ScanStr - Buffer;
        }
        else
        {
            Pos = 0;
        }
        if (!ScanStr || Pos > XML_BUFFER_HALF_SIZE)
        {
            // There's not a <scan tag in the first half of the buffer.  
            // If we're at EOF, then stop now:
            if (BufferEnd < XML_BUFFER_HALF_SIZE)
            {
                break;
            }
            // Shunt the tail of the buffer to the front, and carry on:
            memmove(Buffer, Buffer + XML_BUFFER_HALF_SIZE, BufferEnd - XML_BUFFER_HALF_SIZE);
            BufferEnd -= XML_BUFFER_HALF_SIZE;
            FilePos += XML_BUFFER_HALF_SIZE;
            continue;
        }

        ScanNumberStr = strstr(ScanStr, "num=");
        if (!ScanNumberStr)
        {
            printf("** Warning: mzXML parser encountered a scan with no scan number!  File %s Pos %d\n", 
				   mzxml_name.c_str(), FilePos + Pos);

            ScanNumber = -1;
        }
        else
        {
            ScanNumber = ParseIntFromXML(ScanNumberStr);
        }

		retentionTimeStr = strstr(ScanStr,"retentionTime=\"PT");
		if (! retentionTimeStr)
		{
		//	printf("Error: mzXML parser encountered a scan with no retnetion time: File %s Pos %d\n", 
		//		mzxml_name.c_str(), FilePos + Pos);
         //   exit(1);
			retentionTime = -1;
		}
		else
		{
			retentionTime = ParseMassFromXML(retentionTimeStr);
		}


		char *PeakCountStr = strstr(ScanStr, "peaksCount=\"");
		if (!PeakCountStr)
		{
			cout << "Warning: couldn't parse peaks from mzxml! " << endl;
			cout << "Scan: " << ScanNumber << " Pos: " << Pos << endl;
			continue;
		}
		int PeakCount = ParseIntFromXML(PeakCountStr);
		
        MSLevelStr = strstr(ScanStr, "msLevel=");
        if (!MSLevelStr)
        {
            printf("** Warning: mzXML parser encountered a scan with no MS level! File %s Scan %d Pos %d\n", 
				mzxml_name.c_str(), ScanNumber, FilePos + Pos);
            MSLevel = -1;
			continue;
        }
        else
        {
            MSLevel = ParseIntFromXML(MSLevelStr);
        }

		precursorIntensityStr = strstr(ScanStr,"precursorIntensity=");
		if (! precursorIntensityStr)
		{
			//if (MSLevel>1)
			//{
			//	printf("Error: mzXML parser encountered a scan with no precursor intenisty: File %s Pos %d\n", 
			//		mzxml_name.c_str(), FilePos + Pos);
			//	exit(1);
			//}
			precursorIntensity = 0;
		}
		else
		{
			precursorIntensity = ParseMassFromXML(precursorIntensityStr);
		}


		PrecursorStr = strstr(ScanStr, "<precursorMz");

		if (PrecursorStr)
		{
			PrecursorStr = strstr(PrecursorStr, ">");
			precursorMZ = ParseMassFromXML(PrecursorStr);
		}

		if (!PrecursorStr && MSLevel > 1)
		{
			printf("Warning: mzXML parser encountered a scan with no m/z: File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + Pos);
			continue;
		}
			
        if (MSLevel > 1 && ScanNumber >= 0 && PeakCount > 10)
        {
			spec_header.scan_number = ScanNumber;
			spec_header.MS_level = MSLevel;
			spec_header.precursor_intensity = precursorIntensity;
			spec_header.retention_time = retentionTime;
			spec_header.m_over_z = precursorMZ;
			spec_header.num_peaks = PeakCount;
			spec_header.charge = 0;
			spec_header.type = MZXML;
			spec_header.file_pos = FilePos + Pos;

			single_spectra.push_back(spec_header);
		}
    }


    free(Buffer);
	fclose(MZXMLFile);
//	cout << single_spectra.size() << " spectra..." << endl;
}


/******************************************************************************
	This is a special function designed to overcome parsing problems I have
	with mzXML in Linux enviornments. The function serially extracts spectra from
	an mzXML file and stores the peak lists (floats of pairs (mass,intensity)
*******************************************************************************/
int MZXML_file::extract_peak_lists_from_mzXML(const Config *config, 
								  string& mzxml_name, 
								  int file_idx,
								  mass_t min_m_over_z, 
								  mass_t max_m_over_z)
{
	static char* PeakBuffer = NULL;
    static char* DecodedPeakBuffer = NULL;
	static float* Peaks = NULL;
	static float* FilteredPeaks = NULL;
    static int PeakBufferSize = 0;

	int BytesToRead;
    char* Buffer;
    int BufferStartPos = 0;
    int Pos;
    int BytesRead;
    int BufferEnd = 0;
    FILE* MZXMLFile;
    int ParseState = 0;
    int FilePos = 0;
	int ByteOrderLittle = 1;

    char* ScanStr;
    char* ScanNumberStr;
	char* PeakStr;
    char* MSLevelStr;
	char *retentionTimeStr;
	char *precursorIntensityStr;
	char *PrecursorStr;
	char* ByteOrderStr;
	
	int ScanNumber;
    int MSLevel;
	float retentionTime;
	float precursorIntensity;
	mass_t precursorMZ;


    //
    Buffer = (char*)calloc(XML_BUFFER_SIZE + 1, sizeof(char));
    MZXMLFile = fopen(mzxml_name.c_str(), "rb");
    if (!MZXMLFile)
    {
        cout << "Error: Can't open MZXML file " <<  mzxml_name << endl;
        exit(1);
    }

//	cout << "Extracting peaks from: " << mzxml_name << endl;

	// initialize the file_peak_buff. Initialize to file_size*0.08 floats
	int file_size = getFileSize(mzxml_name.c_str());
	int num_floats = (int)(file_size*0.08);
	if (num_floats<20000)
		num_floats=20000;

	if (file_peak_buff.size()<(int)(num_floats*1.5))
		file_peak_buff.resize((int)(num_floats*1.5));

	
	int spec_counter=0;
	char *scan_start_ptr = NULL;

    while (1)
    {
		MZXML_single spec_header;

		spec_header.file_idx    = file_idx;

		// check if the peak buff should be extended
		if (file_peak_buff.size()-file_peak_buff_pos<5000)
		{
			file_peak_buff.resize((int)(file_peak_buff.size()*1.5));
		}
		
        // Read more data, to fill up the buffer:
     	if ( ! scan_start_ptr || 
			( (Buffer + BufferEnd - scan_start_ptr) < XML_BUFFER_HALF_SIZE) )
		{
			// try shunt half of the buffer
			if (scan_start_ptr)
			{
				if (BufferEnd - XML_BUFFER_HALF_SIZE>0)
				{
					memmove(Buffer, Buffer + XML_BUFFER_HALF_SIZE, BufferEnd - XML_BUFFER_HALF_SIZE);
					BufferEnd -= XML_BUFFER_HALF_SIZE;
					scan_start_ptr -= XML_BUFFER_HALF_SIZE;

//					cout << "MOVED!" << endl;
				}
			}
			else
				scan_start_ptr = Buffer;

			BytesToRead = XML_BUFFER_SIZE - BufferEnd;
			BytesRead = fread(Buffer + BufferEnd, sizeof(char), BytesToRead, MZXMLFile);

			if (BytesRead<5)
				break;

			BufferEnd += BytesRead;
			Buffer[BufferEnd] = '\0';

			FilePos += BytesRead;
		}
        // Look for a new <scan tag opening:
        const char *last_pos = Buffer + BufferEnd - 5;
		char *pos = scan_start_ptr;

		while (++pos<last_pos)
		{
			if (*pos != '<')
				continue;

			if (*(pos+1)=='s' && *(pos+2)=='c' && *(pos+3)=='a' && *(pos+4)=='n')
				break;
		}
		ScanStr =  (pos<last_pos) ? pos : NULL;

        if (ScanStr)
        {
            Pos = ScanStr - Buffer;
        }
        else
        {
            Pos = 0;
        }

        if (!ScanStr )
        {
			scan_start_ptr = Buffer + BufferEnd-5;
            continue;
        }

        ScanNumberStr = strstr(ScanStr, "num=");
        if (!ScanNumberStr)
        {
       //     printf("** Warning: mzXML parser encountered a scan with no scan number!  File %s Pos %d\n", 
	   //			   mzxml_name.c_str(), FilePos + Pos);

            ScanNumber = -1;
        }
        else
        {
            ScanNumber = ParseIntFromXML(ScanNumberStr);
        }

		retentionTimeStr = strstr(ScanStr,"retentionTime=\"PT");
		if (! retentionTimeStr)
		{
		//	printf("Error: mzXML parser encountered a scan with no retnetion time: File %s Pos %d\n", 
		//		mzxml_name.c_str(), FilePos + Pos);
         //   exit(1);
			retentionTime = -1;
		}
		else
		{
			retentionTime = ParseMassFromXML(retentionTimeStr);
		}


		char *PeakCountStr = strstr(ScanStr, "peaksCount=\"");
		if (!PeakCountStr)
		{
			cout << "Warninig: bad parsing of peaks from mzxml! " << endl;
			cout << "Scan: " << ScanNumber << "  Pos: " << FilePos << endl;
			scan_start_ptr += 50;
			continue;
			
		}
		int PeakCount = ParseIntFromXML(PeakCountStr);
		
        MSLevelStr = strstr(ScanStr, "msLevel=");
        if (!MSLevelStr)
        {
            printf("** Warning: mzXML parser encountered a scan with no MS level!  File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + Pos);
            scan_start_ptr += 50;
			continue;
        }
        else
        {
            MSLevel = ParseIntFromXML(MSLevelStr);
        }

		precursorIntensityStr = strstr(ScanStr,"precursorIntensity=");
		if (! precursorIntensityStr)
		{
			scan_start_ptr += 50;
			continue;
		//	if (MSLevel>1)
		//	{
		//		printf("Warning: mzXML parser encountered a scan with no precursor intenisty: File %s Pos %d\n", 
		//			mzxml_name.c_str(), FilePos + Pos);
		//		exit(1);
		//	}
		}
		else
		{
			precursorIntensity = ParseMassFromXML(precursorIntensityStr);
		}


		PrecursorStr = strstr(ScanStr, "<precursorMz");
		if (PrecursorStr)
		{
			PrecursorStr = strstr(PrecursorStr, ">");
			precursorMZ = ParseMassFromXML(PrecursorStr);
		}
		if (!PrecursorStr && MSLevel > 1)
		{
		
			printf("Warning: mzXML parser encountered a scan with no m/z: File %s Pos %d\n", 
				mzxml_name.c_str(), FilePos + Pos);
			scan_start_ptr += 50;
			continue;
		}

			
        if (MSLevel > 1 && ScanNumber >= 0 && PeakCount > 7 &&
			precursorMZ>=min_m_over_z && precursorMZ<=max_m_over_z)
        {
			// read peaks

			PeakStr = strstr(PrecursorStr, "<peaks");
			if (PeakStr)
			{
				// Get byte order:
				ByteOrderStr = strstr(PeakStr, "byteOrder=\"");
				if (ByteOrderStr)
				{
					ByteOrderStr += 11;
					if (!strncmp(ByteOrderStr, "network", 7))
					{
						ByteOrderLittle = 0;
					}
					if (!strncmp(ByteOrderStr, "big", 3))
					{
						ByteOrderLittle = 0;
					}
					if (!strncmp(ByteOrderStr, "little", 6))
					{
						ByteOrderLittle = 1;
					}
				}
				PeakStr = strstr(PeakStr, ">");
			}
			if (!PeakStr)
			{
				cout << "Warning: bad parsing of peaks from mzxml (scan " << ScanNumber << ") skipping..." << endl;
				scan_start_ptr += 50;
				continue;
			}

			PeakStr++;
			PeakBuffer = PeakStr;

			if (PeakBufferSize < PeakCount)
			{
				if (DecodedPeakBuffer)
				{
					char *dbf = DecodedPeakBuffer;
					free(DecodedPeakBuffer);
					DecodedPeakBuffer = NULL;
					free(Peaks);
					Peaks = NULL;
					free(FilteredPeaks);
					FilteredPeaks=NULL;
				}
				PeakBufferSize = (int)(PeakCount*1.5);
				DecodedPeakBuffer = (char*)calloc(PeakBufferSize * 8 + 8, 1);
				Peaks = (float*)calloc(PeakBufferSize * 2, sizeof(float));
				FilteredPeaks = (float*)calloc(PeakBufferSize * 2, sizeof(float));
			}
			
			int Trail = (PeakCount % 3);
			if (!(PeakCount % 3))
			{
				PeakBuffer[PeakCount * 32/3] = '\0';
			}
			else
			{
				PeakBuffer[(PeakCount * 32/3) + Trail + 1] = '\0';
			}
	
			b64_decode_mio( DecodedPeakBuffer, PeakBuffer);
			int FloatIndex;
			for (FloatIndex = 0; FloatIndex < (2 * PeakCount); FloatIndex++)
			{
		#ifdef BYTEORDER_LITTLE_ENDIAN
				if (!ByteOrderLittle)
				{
					char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
					DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
					DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
					ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
					DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
					DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
				}
				memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
		#else
				if (ByteOrderLittle)
				{
					char ByteSwap = DecodedPeakBuffer[FloatIndex*4];
					DecodedPeakBuffer[FloatIndex*4] = DecodedPeakBuffer[FloatIndex*4 + 3];
					DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
					ByteSwap = DecodedPeakBuffer[FloatIndex*4 + 1];
					DecodedPeakBuffer[FloatIndex*4 + 1] = DecodedPeakBuffer[FloatIndex*4 + 2];
					DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
				}
				memcpy(Peaks + FloatIndex, DecodedPeakBuffer + FloatIndex * 4, 4);
		#endif
			}


			
			int num_new_peaks = join_and_filter_peak_list(config,precursorMZ,Peaks,
				PeakCount, &file_peak_buff[file_peak_buff_pos]);

			spec_header.scan_number			= ScanNumber;
			spec_header.MS_level			= MSLevel;
			spec_header.precursor_intensity = precursorIntensity;
			spec_header.retention_time		= retentionTime;
			spec_header.m_over_z			= precursorMZ;
			spec_header.num_peaks			= num_new_peaks;
			spec_header.charge = 0;
			spec_header.type = MZXML;
			spec_header.file_pos = FilePos + Pos;

			spec_header.peak_buff_start_idx = file_peak_buff_pos;

			single_spectra.push_back(spec_header);

			spec_counter++;

			file_peak_buff_pos += num_new_peaks * 2;

			scan_start_ptr = PeakStr + 8 * PeakCount;
		}
		else
			scan_start_ptr = ScanStr +50;
    }
    free(Buffer);
	fclose(MZXMLFile);

//	cout << spec_counter << " spectra, " << file_peak_buff_pos << "/" << file_peak_buff.size() << endl;



	return spec_counter;
}









