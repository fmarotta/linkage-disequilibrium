/* Interface definition
 *
 * Calculate linkage disequilibrium from a VCF file.
 */

// TODO: internally, use a single linked list to host all alleles, be them ref
// or alt. To the user, provide functions to access either the ref or all the
// alts.

#ifndef _LD_VCF_H_
#define _LD_VCF_H_
#include <stdbool.h>

//#define MAXFILTLEN 10
#define MAXVTLEN 5 // enough to accomodate `INDEL'.
#define MAXIDLEN 100
#define SEQLEN 150 // max length of allele seq

#define WINLEN 100000 // length of the window, in bases.

typedef struct vcf_allele {
	char allele_seq[SEQLEN]; // allocate space with malloc
	int allele_num; // number assigned to the allele (found in the genotype field)
	int ac; // number of ref/alt alleles in called genotypes
	float af; // ref/alt allele frequency in the range (0,1)
	char vt[MAXVTLEN]; // what type of variant the line represents; for ref alleles this is always `REF'.
	struct vcf_allele *next;
} VCF_ALLELE;

// This complicates a bit the structure. Should we use a linked list? I doubt
// it. Maybe some codes. In our vcf, however, all records are "PASS".
typedef struct vcf_filter {
	bool pass;
} VCF_FILTER;

typedef struct vcf_info {
	int ns; // total number of samples with data
	int an; // total number of alleles in called genotypes
	int _an; // number of different alleles (1 ref + N alt) (*)
} VCF_INFO;

typedef struct vcf_format_gt {
	int m; // maternal allele
	int p; // paternal allele
} VCF_FORMAT_GT;

// whether this sample is phased or not (*)
/* XXX The user shoud build this structure with as many fields as the format
 * field in the vcf, or she can limit to the relevant ones. Also, here we
 * consider diploid organisms, but this can be edited as well. We also add a
 * member for the phasing of the sample.
 */
typedef struct vcf_sample {
	VCF_FORMAT_GT gt;
	bool phased;
	struct vcf_sample *next;
} VCF_SAMPLE;

typedef struct vcf_locus {
	int chrom; // what about X and Y? -1 and -2? 23 and 24? enum??
	unsigned long pos;
	char id[MAXIDLEN]; // the complete ID field.
	VCF_ALLELE *alleles; // the first allele in this linked list is the ref.
	int qual;
	VCF_FILTER filter;
	VCF_INFO info;
	VCF_SAMPLE *samples;
	struct vcf_locus *next;
} VCF_LOCUS;

typedef struct vcf_window {
	VCF_LOCUS *head;
	VCF_LOCUS *tail;
	int nloci; // number of loci currently in the queue
	int winlen; // length of the sliding window
} VCF_WINDOW;


/* operation:		reads from the vcf file all the loci that fall within
 * 					<winlen> bases from the first locus.
 * precondition:	vcf_file is fopen'd, a pointer to window is defined;
 * 					winlen is the length of the window.
 * postcondition:	adds the first loci to the queue and initializes it. */
void Initialize_window(FILE *vcf_file, VCF_WINDOW *pwindow, int winlen);

/* operation:		moves the window forward one locus.
 * precondition:	pwindow is initialized.
 * postcondition:	removes the first locus from the queue; if appropriate
 * 					adds more locus from the file to the queue. returns true
 * 					until the window is full. */
bool Slide_window(FILE *vcf_file, VCF_WINDOW *pwindow);

/* operation:		removes all the loci from the window.
 * precondition:	pwindow is initialized.
 * postcondition:	all memory is freed. */
void Close_window(VCF_WINDOW *pwindow);

/* operation:		adds one locus to the window; do not bother using it,
 * 					please: use Initialize_window() and Slide_window().
 * precondition:	a locus must have been read and digested.
 * poscondition:	a new locus is added to the window. */
//bool Enqueue_locus(VCF_LOCUS locus, VCF_WINDOW *pwindow);

/* operation:		removes one locus from window; do not bother using it,
 * 					please: use Initialize_window() and Slide_window().
 * precondition:	
 * poscondition:	the memory occupied by the deleted locus is freed */
//bool Dequeue_locus(VCF_LOCUS *plocus, VCF_WINDOW *pwindow);

/* operation:		gets number of loci currently in the window.
 * precondition:	pwindow points to an initialized window.
 * poscondition:	returns nloci. */
unsigned int Nloci_in_window(const VCF_WINDOW *pwindow);

/* operation:		gets total number of alleles in a locus. The number
 * 					of alternate alleles is therefore (total - 1).
 * precondition:	plocus points to a locus (vcf row) of an initialized
 * 					window.
 * postcondition:	returns the number of alleles at that locus. */
int Nalleles_in_locus(const VCF_LOCUS *plocus);

/* operation:		finds all the alleles in a locus.
 * precondition:	plocus points to a locus (vcf row) of an initialized
 * 					window.
 * postcondition:	returns a linked list of alleles, the first being
 * 					the reference and with the alternates following. */
//VCF_ALLELE *alleles_in_locus(const VCF_LOCUS *plocus);

/* operation:		computes a unique id for an allele of a locus.
 * precondition:	plocus points to a locus.
 * postcondition:	returns a pointer to an ID string made out of the rs id
 * 					and the allele number. For triallelic loci, the ref
 * 					allele has both rs ids, the alt alleles have each the
 * 					appropriate one. */
char *Get_allele_id(int alnum, const VCF_LOCUS *plocus);

/* operation:		calculates the frequency of a genotype.
 * precondition:	pwindow points to an initialized window.
 * poscondition:	returns the frequency of the genotype. */
//float Genotype_freq(const VCF_FORMAT_GT *gt, const VCF_WINDOW *pwindow);

float Allele_freq(int alnum, const VCF_LOCUS *plocus);

// XXX this would be a great occasion to write a variable-argument-number
// function, if we knew how to calculate LD for more than two alleles!

/* operation:		calculates the frequency (p) of a pair of alleles at
 * 					different loci occurring together.
 * precondition:	alnum1 and alnum2 are two alleles of locus1 and locus2,
 * 					respectively.
 * postcondition:	returns the frequency of such event. */
float Linked_alleles_freq(int alnum1, const VCF_LOCUS *plocus1,
						  int alnum2, const VCF_LOCUS *plocus2);

float Calculate_D(float p_A, float p_B, float p_AB);

float Calculate_D_lewontin(float p_A, float p_B, float p_AB);

float Calculate_r_squared(float p_A, float p_B, float p_AB);

#endif
