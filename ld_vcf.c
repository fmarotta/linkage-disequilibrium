/* Interface implementation */

// TODO close fds and free malloc'd memory.

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ld_vcf.h"
#include "io_utils.h"
#include "type_utils.h"

#define FILTLEN 50 // in our case we can only have PASS
#define GTLEN 3 // max length of genotype
#define INFOLEN 200

static VCF_LOCUS buflocus;

static bool enqueue_locus(VCF_LOCUS locus, VCF_WINDOW *pwindow);

/* digest_line() returns 0 on success, -1 at EOF, 1 if memory failure,
 * and 2 when the line is malformed. */
static int digest_line(VCF_LOCUS *plocus, FILE *vcf_file);
static void vomit_line(const VCF_LOCUS *plocus);

static VCF_ALLELE *make_allele(const char *seq, int alnum);
static VCF_SAMPLE *make_sample(int m, int p, bool phased);
static void free_alleles(VCF_LOCUS *plocus);
static void free_samples(VCF_LOCUS *plocus);

static bool locus_is_in_window(const VCF_LOCUS *plocus, const VCF_WINDOW *pwindow);
static bool locus_is_valid(const VCF_LOCUS *plocus);

static void foreach_subfield(void (*fn)(char *subfield, VCF_LOCUS *plocus), const char *field, char sep, VCF_LOCUS *plocus);
static void parse_ac(char *subfield, VCF_LOCUS *plocus);
static void parse_af(char *subfield, VCF_LOCUS *plocus);
static void parse_vt(char *subfield, VCF_LOCUS *plocus);
static void parse_alt_seq(char *subfield, VCF_LOCUS *plocus);
static void parse_info(char *subfield, VCF_LOCUS *plocus);

static int compare_loci(VCF_LOCUS *plocus1, VCF_LOCUS *plocus2);


// Initialize_window {{{

/* We want to read a line from the file and decide whether it should be
 * added to the window according to the position of the locus and
 * perhaps the FILTER value. If a line is discarded because the locus
 * being read is too far from the first one, however, we should not
 * loose it forever, indeed that line can be added to the window upon
 * sliding. Therefore, we shall implement a one-locus buffer to keep the
 * locus, so that first the line is read into the buffer, then the
 * buffer is read and the decision whether to add the locus to the
 * window is taken, and if the buffer is emptied we read a new line from
 * the file.
 */
void Initialize_window(FILE *vcf_file, VCF_WINDOW *pwindow, int winlen) {
	char line[3];
	int status;

	// Initialize the window
	pwindow->head = pwindow->tail = NULL;
	pwindow->nloci = 0;
	pwindow->winlen = winlen;

	// Initialize the buffer pointers
	buflocus.alleles = NULL;
	buflocus.samples = NULL;

	// Discard header lines
	//
	// There is no need to read the whole line, the first two characters
	// are enough. Detect the change from ## to #.
	while (fgets(line, 3, vcf_file) != NULL && strcmp(line, "##") == 0)
		EATLINE(vcf_file);
	EATLINE(vcf_file);

	// Digest the first data line into the one-locus buffer
	if ((status = digest_line(&buflocus, vcf_file)) != 0)
	{
		if (status == -1)
		{
			fputs("ERROR: no data found in the vcf file.", stderr);
			exit(EXIT_FAILURE);
		}
		else
		{
			fputs("ERROR: we ran out of memory.", stderr);
			exit(EXIT_FAILURE);
		}
	}

	// Add loci to the window.
	while (locus_is_in_window(&buflocus, pwindow))
	{
		// filters for quality, number of alleles...
		if (locus_is_valid(&buflocus))
			// Add valid locus
			if (!enqueue_locus(buflocus, pwindow))
			{
				fputs("ERROR: could not add locus to the window.", stderr);
				exit(EXIT_FAILURE);
			}

		// Read the next line into the buffer
		if ((status = digest_line(&buflocus, vcf_file)) != 0)
		{
			if (status == -1)
			{
				return;
			}
			else
			{
				fputs("ERROR: we ran out of memory.", stderr);
				exit(EXIT_FAILURE);
			}
		}
	}
}
// }}}

// Slide_window {{{
bool Slide_window(FILE *vcf_file, VCF_WINDOW *pwindow)
{
	int status;
	VCF_LOCUS *tmp;

	if (pwindow->nloci == 0)
	{
		// do not remove first element, add the very next line from the
		// file. TODO
	}

	// Remove the first locus 
	tmp = pwindow->head;
	pwindow->head = pwindow->head->next;
	free_alleles(tmp);
	free_samples(tmp);
	free(tmp);
	pwindow->nloci--;
	if (pwindow->nloci == 0)
		pwindow->tail = NULL;

	// Add new locus if appropriate
	while (locus_is_in_window(&buflocus, pwindow))
	{
		// filters for quality, number of alleles...
		if (locus_is_valid(&buflocus))
			//Add valid locus
			if (!enqueue_locus(buflocus, pwindow))
			{
				fputs("ERROR: could not add locus to the window.", stderr);
				exit(EXIT_FAILURE);
			}

		// Read the next line into the buffer
		if ((status = digest_line(&buflocus, vcf_file)) != 0)
		{
			if (status == -1)
			{
				return false;
			}
			else
			{
				fputs("ERROR: we ran out of memory.", stderr);
				exit(EXIT_FAILURE);
			}
		}
	}

	return true;
}
// }}}

// Nalleles_in_locus {{{
int Nalleles_in_locus(const VCF_LOCUS *plocus)
{
	return plocus->info._an;
}
// }}}

// Linked_alleles_freq {{{
float Linked_alleles_freq(int alnum1, const VCF_LOCUS *plocus1,
						  int alnum2, const VCF_LOCUS *plocus2)
{
	VCF_SAMPLE *psample1, *psample2;
	int c_AB = 0; // count
	float p_AB; // frequency
	int ns;

	psample1 = plocus1->samples;
	psample2 = plocus2->samples;
	ns = (plocus1->info.ns <= plocus2->info.ns) ? plocus1->info.ns : plocus2->info.ns;
	for (int i = 0; i < ns; i++)
	{
		if ((psample1->gt.m == alnum1 && psample2->gt.m == alnum2) ||
			(psample1->gt.p == alnum1 && psample2->gt.p == alnum2))
			c_AB++;
		psample1 = psample1->next;
		psample2 = psample2->next;
	}

	p_AB = (float) c_AB / (2 * ns);

	return p_AB;
}
// }}}

// Allele_freq {{{
float Allele_freq(int alnum, const VCF_LOCUS *plocus)
{
	VCF_ALLELE *pallele;

	pallele = plocus->alleles;
	for (int i = 0; i < alnum; i++)
	{
		pallele = pallele->next;
		if (pallele == NULL)
			return -1;
	}

	return pallele->af;
}
// }}}

// Calculate_D {{{
float Calculate_D(float p_A, float p_B, float p_AB)
{
	return (p_AB - (p_A * p_B));
}
// }}}

// Calculate_D_lewontin {{{
float Calculate_D_lewontin(float p_A, float p_B, float p_AB)
{
	float D, Dmax;

	D = p_AB - (p_A*p_B);
	if (D < 0)
		Dmax = (-p_A*p_B > -(1-p_A)*(1-p_B)) ? -p_A*p_B : -(1-p_A)*(1-p_B);
	else
		Dmax = (p_A*(1-p_B) < p_B*(1-p_A)) ? p_A*(1-p_B) : p_B*(1-p_A);

	return D/Dmax;
}
// }}}

// Calculate_r_squared {{{
float Calculate_r_squared(float p_A, float p_B, float p_AB)
{
	float D, d;

	D = p_AB - (p_A*p_B);
	d = p_A*(1-p_A) * p_B*(1-p_B);

	return (D*D)/d;
}
// }}}

// compare_loci {{{

/* returns 0 if it is the same locus, -1 if locus1 is upstream, +1 if
 * locus2 is upstream. */
static int compare_loci(VCF_LOCUS *plocus1, VCF_LOCUS *plocus2)
{
	if (plocus1->pos == plocus2->pos)
		return 0;
	else if (plocus1->pos < plocus2->pos)
		return -1;
	else
		return 1;
}
// }}}

// enqueue_locus {{{
static bool enqueue_locus(VCF_LOCUS locus, VCF_WINDOW *pwindow)
{
	VCF_LOCUS *pnew;

	pnew = (VCF_LOCUS *) malloc(sizeof(VCF_LOCUS));
	if (pnew == NULL)
		return false;
	*pnew = locus;

	if (pwindow->head == NULL)
		pwindow->head = pnew;
	else
		pwindow->tail->next = pnew;
	pwindow->tail = pnew;
	pwindow->nloci++;

	return true;
}
// }}}

// digest_line {{{
static int digest_line(VCF_LOCUS *plocus, FILE *vcf_file)
{
	VCF_ALLELE *newallele, *lastallele;
	VCF_SAMPLE *newsample, *lastsample;

	char tmp_ref_seq[SEQLEN], tmp_alt_seq[SEQLEN];
	char tmp_filter[FILTLEN], tmp_info[INFOLEN];
	char tmp_gt[GTLEN];

	// XXX what about chr X and Y? are they integer?
	fscanf(vcf_file, "%d%ld%s%s%s%d%s%s", &plocus->chrom, &plocus->pos,
			&plocus->id, tmp_ref_seq, tmp_alt_seq,
			&plocus->qual, tmp_filter, tmp_info);

	// Filter
	if (strcmp(tmp_filter, "PASS") == 0)
		plocus->filter.pass = true;
	else
		plocus->filter.pass = false;

	// Free the previous allocations XXX DO NOT DO THAT!!!! otherwise
	// makeallele() is free to reallocate old locations, overwriting
	// alleles of previous loci!!!
	//free_alleles(plocus); free_samples(plocus);

	// Allocate space for the alleles
	plocus->info._an = 0;
	newallele = make_allele(tmp_ref_seq, plocus->info._an++);
	if (newallele == NULL)
		return 1;
	plocus->alleles = newallele;

	// alt seq
	foreach_subfield(parse_alt_seq, tmp_alt_seq, ',', plocus);

	foreach_subfield(parse_info, tmp_info, ';', plocus);

	// ref info
	VCF_ALLELE *ref;
	lastallele = plocus->alleles->next; // start from the first alt allele
	ref = plocus->alleles;
	ref->ac = plocus->info.an;
	ref->af = 1;
	strcpy(plocus->alleles->vt, "REF");
	while (lastallele != NULL)
	{
		ref->ac -= lastallele->ac;
		ref->af -= lastallele->af;
		lastallele = lastallele->next;
	}

	// skip one field (the format)
	fscanf(vcf_file, "%s", tmp_gt);

	// Read in the first sample
	fscanf(vcf_file, "%s", tmp_gt);
	newsample = make_sample(tmp_gt[0]-48, tmp_gt[2]-48, (tmp_gt[1] == '|'));
	if (newsample == NULL)
		return 1;
	plocus->samples = newsample;
	lastsample = plocus->samples;

	// Read the other samples
	if (plocus->info.ns != 2504)
		printf("%d\n",plocus->info.ns);
	for (int i = 1; i < plocus->info.ns; i++)
	{
		// XXX we assume that no locus has more than 10 alleles...
		fscanf(vcf_file, "%s", tmp_gt);
		newsample = make_sample(tmp_gt[0]-48, tmp_gt[2]-48, (tmp_gt[1] == '|'));
		if (newsample == NULL)
			return 1;
		lastsample->next = newsample;
		lastsample = lastsample->next;
	}

	EATLINE(vcf_file);
	//vomit_line(plocus);

	// after getting samples...
	if (feof(vcf_file))
		return -1;
	
	return 0;
}
// }}}

// vomit_line {{{
static void vomit_line(const VCF_LOCUS *plocus)
{
	VCF_ALLELE *pallele = plocus->alleles;
	VCF_SAMPLE *psample = plocus->samples;

	printf("LOCUS\n");
	printf("chrom %d\tpos %ld\tid %s\tn_samples %d\tn_haplotypes %d\tn_alleles %d\n",
			plocus->chrom, plocus->pos, plocus->id,
			plocus->info.ns, plocus->info.an, plocus->info._an);
	printf("ALLELES\n");
	while (pallele != NULL)
	{
		printf("allele %d: seq %s\tcount %d\tfreq %f\ttype %s\n",
				pallele->allele_num, pallele->allele_seq,
				pallele->ac, pallele->af, pallele->vt);
		pallele = pallele->next;
	}
	/*printf("SAMPLES\n");
	while (psample != NULL)
	{
		printf("%d%c%d\t", psample->gt.m, psample->phased ? '|' : '/', psample->gt.p);
		psample = psample->next;
	}
	putchar('\n');
	putchar('\n');
	*/
}
// }}}

// make_allele {{{
static VCF_ALLELE *make_allele(const char *seq, int alnum)
{
	VCF_ALLELE *newallele;

	newallele = (VCF_ALLELE *) malloc(sizeof(VCF_ALLELE));
	//printf("newly created allele: %p\n", newallele);
	if (newallele != NULL)
	{
		strcpy(newallele->allele_seq, seq);
		newallele->allele_num = alnum;
		newallele->ac = -1;
		newallele->af = -1;
		strcpy(newallele->vt, "");
		newallele->next = NULL;
	}
	
	return newallele;
}
// }}}

// make_sample {{{
static VCF_SAMPLE *make_sample(int m, int p, bool phased)
{
	VCF_SAMPLE *newsample;

	newsample = (VCF_SAMPLE *) malloc(sizeof(VCF_SAMPLE));
	if (newsample != NULL)
	{
		newsample->gt.m = m;
		newsample->gt.p = p;
		newsample->phased = phased;
		newsample->next = NULL;
	}

	return newsample;
}
// }}}

// free_alleles {{{
static void free_alleles(VCF_LOCUS *plocus)
{
	VCF_ALLELE *tmp;

	while (plocus->alleles != NULL)
	{
		tmp = plocus->alleles->next;
		free(plocus->alleles);
		plocus->alleles = tmp;
	}
}
// }}}

// free_samples {{{
static void free_samples(VCF_LOCUS *plocus)
{
	VCF_SAMPLE *tmp;

	while (plocus->samples != NULL)
	{
		tmp = plocus->samples->next;
		free(plocus->samples);
		plocus->samples = tmp;
	}
}
// }}}

// locus_is_in_window {{{
static bool locus_is_in_window(const VCF_LOCUS *plocus, const VCF_WINDOW *pwindow)
{
	if (pwindow->head == NULL)
		return true;
	else
		return (plocus->pos - pwindow->head->pos <= pwindow->winlen);
}
// }}}

// locus_is_valid {{{
static bool locus_is_valid(const VCF_LOCUS *plocus)
{
	// No conditions.
	return true;
}
// }}}

// foreach_subfield {{{
static void foreach_subfield(void (*fn)(char *subfield, VCF_LOCUS *plocus), const char *field, char sep, VCF_LOCUS *plocus)
{
	char *subfield, *find;

	while ((find = strchr(field, sep)) != NULL)
	{
		subfield = (char *) malloc(sizeof(char) * (find - field + 1));
		strncpy(subfield, field, find - field + 1);
		subfield[find-field] = '\0';
		(*fn)(subfield, plocus);
		field = find + 1;
		free(subfield);
	}
	subfield = (char *) malloc(sizeof(char) * (strlen(field) + 1));
	strncpy(subfield, field, strlen(field) + 1);
	subfield[strlen(field)] = '\0';
	(*fn)(subfield, plocus);
	free(subfield);
}
// }}}

// parse_ac {{{
static void parse_ac(char *subfield, VCF_LOCUS *plocus)
{
	VCF_ALLELE *tmp;

	//printf("AC subdfield: %s.\n", subfield);
	tmp = plocus->alleles->next; // start from the first alt allele

	//printf("expecting alt info: %d.\n", tmp->allele_info.ac);

	while (tmp->ac >= 0) 
	{
	//printf("expecting alt info: %d.\n", tmp->allele_info.ac);
	//printf("expecting next alt info: %d.\n",
	//tmp->next->allele_info.ac); printf("Whiling...\n"); if
	//(tmp->allele_info.ac < 0) break;
		tmp = tmp->next;
	}
	//printf("expercing tmp seq: %s.\n", tmp->allele_seq);
	tmp->ac = atoi(subfield);
}
// }}}

// parse_af {{{
static void parse_af(char *subfield, VCF_LOCUS *plocus)
{
	VCF_ALLELE *tmp;

	//printf("AF subdfield: %s.\n", subfield);
	tmp = plocus->alleles->next; // start from the first alt allele
	while (tmp->next != NULL && tmp->af >= 0)
		tmp = tmp->next;
	tmp->af = atof(subfield);
}
// }}}

// parse_vt {{{
static void parse_vt(char *subfield, VCF_LOCUS *plocus)
{
	VCF_ALLELE *tmp;

	tmp = plocus->alleles->next; // start from the first alt allele
	while (tmp->next != NULL && tmp->vt[0] != '\0') // stop at the first non-initialized member
		tmp = tmp->next;
	strcpy(tmp->vt, subfield);
}
// }}}

// parse_alt_seq {{{
static void parse_alt_seq(char *subfield, VCF_LOCUS *plocus)
{
	VCF_ALLELE *newallele, *tmp;

	newallele = make_allele(subfield, plocus->info._an++);
	// XXX handle this error
	if (newallele == NULL)
		exit(EXIT_FAILURE);

	tmp = plocus->alleles;
	//printf("first allele of locus %d: %p\n",plocus->pos, tmp);
	while (tmp->next != NULL)
	{
		//printf("seq: %s.\n", tmp->allele_seq);
		tmp = tmp->next;
	}
	//printf("last allele of locus %d: %p\n", plocus->pos, tmp);
	tmp->next = newallele;
	//printf("new last allele %p\n", tmp->next);
}
// }}}

// parse_info {{{
static void parse_info(char *subfield, VCF_LOCUS *plocus)
{
	char *datum;

	//printf("subfield %s.\n", subfield);

	if (strncmp(subfield, "NS=", 3) == 0)
	{
		datum = strchr(subfield, '=') + 1;
		plocus->info.ns = atoi(datum);
	}
	else if (strncmp(subfield, "AN=", 3) == 0)
	{
		datum = strchr(subfield, '=') + 1;
		plocus->info.an = atoi(datum);
	}
	else if (strncmp(subfield, "AC=", 3) == 0)
	{
		datum = strchr(subfield, '=') + 1;
		foreach_subfield(parse_ac, datum, ',', plocus);
	}
	else if (strncmp(subfield, "AF=", 3) == 0)
	{
		datum = strchr(subfield, '=') + 1;
		foreach_subfield(parse_af, datum, ',', plocus);
	}
	else if (strncmp(subfield, "VT=", 3) == 0)
	{
		datum = strchr(subfield, '=') + 1;
		foreach_subfield(parse_vt, datum, ',', plocus);
	}
}
// }}}

