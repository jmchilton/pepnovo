#include "QuickClustering.h"
#include "PMCSQS.h"
#include "auxfun.h"

// the sim matrix stores the similarity distances computed between clusters
int     num_sim_matrix_spectra = 0;
unsigned char  * sim_matrix = NULL;
unsigned char  * max_sim_addr = NULL;



void print_byte(unsigned char byte)
{
	int i;
	unsigned char mask=0X1;

	for (i=0; i<8; i++)
	{
		cout << ( (byte & mask) ? '1' : '0');
		mask <<= 1;
	}
}


void mark_bit_zero(unsigned char *addr, int position)
{
	const unsigned char ANDmasks[]={254,253,251,247,239,223,191,127};
	int bit_offset = (position & 0X7);
	int byte_offset = (position >> 3);
	
	*(addr+byte_offset) &= ANDmasks[bit_offset];
}


void mark_bit_one(unsigned char *addr, int position)
{
	const unsigned char ORmasks[] ={1,2,4,8,16,32,64,128};
	int bit_offset = (position & 0X7);
	int byte_offset = (position >> 3);

	*(addr+byte_offset) |= ORmasks[bit_offset];
}


// returns full int assumes we are at position 31 in the 4 bytes
int get_matrix_32_bits(unsigned char *row_start, int position)
{
	int cell_off = (position >> 5);
	return (*((int *)row_start+cell_off)); 
//	return 1;
}

int get_matrix_val(unsigned char *row_start, int position)
{
	const unsigned char masks[] ={1,2,4,8,16,32,64,128};
	int bit_offset = (position & 0X7);
	int byte_offset = (position >> 3);
	
	return ((*(row_start+byte_offset) & masks[bit_offset]));
//	return 1;
}


// global debugging variable
int wrongly_filtered_spectra;









