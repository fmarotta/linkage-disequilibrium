/* Interface implementation */

// TODO close fds and free malloc'd memory.

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ld_vcf.h"
#include "io_utils.h"
#include "type_utils.h"

#define SEQLEN 150 // max length of allele seq
#define GTLEN 3 // max length of genotype
#define FILTLEN 50 // in our case we can only have PASS
#define INFOLEN 200

static VCF_LOCUS buflocus;

/* digest_line() returns 0 on success, -1 at EOF, 1 if memory failure,
 * and 2 when the line is malformed. */
static int digest_line(VCF_LOCUS *plocus, FILE *vcf_file);
static void vomit_line(const VCF_LOCUS *plocus);
static bool locus_is_in_window(const VCF_LOCUS *plocus, const VCF_WINDOW *pwindow);
static bool locus_is_valid(const VCF_LOCUS *plocus);

static void free_alleles(VCF_ALLELE **palleles);
static void free_samples(VCF_SAMPLE **psamples);
static VCF_ALLELE *make_allele(const char *seq, int alnum);
static VCF_SAMPLE *make_sample(int m, int p, bool phased);

static void foreach_subfield(void (*fn)(char *subfield, VCF_LOCUS *plocus), const char *field, char sep, VCF_LOCUS *plocus);
static void parse_info(char *subfield, VCF_LOCUS *plocus);
static void parse_ac(char *subfield, VCF_LOCUS *plocus);
static void parse_af(char *subfield, VCF_LOCUS *plocus);
static void parse_vt(char *subfield, VCF_LOCUS *plocus);
static void parse_alt_seq(char *subfield, VCF_LOCUS *plocus);
static void parse_gt(char *subfield, VCF_LOCUS *plocus);

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
			if (!Enqueue_locus(buflocus, pwindow))
			{
				fputs("ERROR: could not add locus to the window.", stderr);
				exit(EXIT_FAILURE);
			}

		// Read the next line into the buffer
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
	}
}

static bool locus_is_in_window(const VCF_LOCUS *plocus, const VCF_WINDOW *pwindow)
{
	if (pwindow->head == NULL)
		return true;
	else
		return (plocus->pos - pwindow->head->pos <= pwindow->winlen);
}

bool Enqueue_locus(VCF_LOCUS locus, VCF_WINDOW *pwindow)
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

	// Free the previous allocations
	if (plocus->alleles != NULL) free_alleles(&plocus->alleles);
	if (plocus->samples != NULL) free_samples(&plocus->samples);

	// Allocate space for the alleles
	plocus->info._an = 0;
	newallele = make_allele(tmp_ref_seq, plocus->info._an++);
	if (newallele == NULL)
		return 1;
	plocus->alleles = newallele;

	foreach_subfield(parse_alt_seq, tmp_alt_seq, ',', plocus);

	// alt info & general info
	foreach_subfield(parse_info, tmp_info, ';', plocus);

	// ref info
	VCF_ALLELE *ref;
	lastallele = plocus->alleles->next; // start from the first alt allele
	ref = plocus->alleles;
	ref->allele_info.ac = plocus->info.an;
	ref->allele_info.af = 1;
	strcpy(plocus->alleles->allele_info.vt, "REF");
	while (lastallele != NULL)
	{
		ref->allele_info.ac -= lastallele->allele_info.ac;
		ref->allele_info.af -= lastallele->allele_info.af;
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
	vomit_line(plocus);

	// after getting samples...
	if (feof(vcf_file))
		return -1;
	
	return 0;
}

static void foreach_subfield(void (*fn)(char *subfield, VCF_LOCUS *plocus), const char *field, char sep, VCF_LOCUS *plocus)
{
	char *subfield, *find;

	while ((find = strchr(field, sep)) != NULL)
	{
		subfield = (char *) malloc(sizeof(char) * (find - field + 1));
		strncpy(subfield, field, find - field);
		(*fn)(subfield, plocus);
		field = find + 1;
		free(subfield);
	}
	subfield = (char *) malloc(sizeof(char) * (strlen(field) + 1));
	strncpy(subfield, field, strlen(field) + 1);
	(*fn)(subfield, plocus);
	free(subfield);
}

static void parse_ac(char *subfield, VCF_LOCUS *plocus)
{
	VCF_ALLELE *tmp;

	tmp = plocus->alleles->next; // start from the first alt allele
	while (tmp->next != NULL && tmp->allele_info.ac >= 0)
		tmp = tmp->next;
	tmp->allele_info.ac = atoi(subfield);
}

static void parse_af(char *subfield, VCF_LOCUS *plocus)
{
	VCF_ALLELE *tmp;

	tmp = plocus->alleles->next; // start from the first alt allele
	while (tmp->next != NULL && tmp->allele_info.af >= 0)
		tmp = tmp->next;
	tmp->allele_info.af = atof(subfield);
}

static void parse_vt(char *subfield, VCF_LOCUS *plocus)
{
	VCF_ALLELE *tmp;

	tmp = plocus->alleles->next; // start from the first alt allele
	while (tmp->next != NULL && tmp->allele_info.vt[0] != '\0') // stop at the first non-initialized member
		tmp = tmp->next;
	strcpy(tmp->allele_info.vt, subfield);
}

static void parse_alt_seq(char *subfield, VCF_LOCUS *plocus)
{
	VCF_ALLELE *newallele, *tmp;

	newallele = make_allele(subfield, plocus->info._an++);
	// XXX handle this error
	if (newallele == NULL)
		exit(EXIT_FAILURE);

	tmp = plocus->alleles;
	while (tmp->next != NULL)
		tmp = tmp->next;
	tmp->next = newallele;
}

static void parse_info(char *subfield, VCF_LOCUS *plocus)
{
	char *datum;

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

static void free_alleles(VCF_ALLELE **palleles)
{
	VCF_ALLELE *tmp, *next;

	tmp = *palleles; // tmp points to the first node
	do {
		next = tmp->next;
		free(tmp);
		tmp = next;
	} while (tmp != NULL);

	*palleles = NULL;
}

static void free_samples(VCF_SAMPLE **psamples)
{
	VCF_SAMPLE *tmp, *next;

	tmp = *psamples; // tmp points to the first node
	do {
		next = tmp->next;
		free(tmp);
		tmp = next;
	} while (tmp != NULL);

	*psamples = NULL;
}

static VCF_ALLELE *make_allele(const char *seq, int alnum)
{
	VCF_ALLELE *newallele;

	newallele = (VCF_ALLELE *) malloc(sizeof(VCF_ALLELE));
	if (newallele != NULL)
	{
		newallele->allele_seq = (char *) malloc(sizeof(char) * (strlen(seq)+1));
		strcpy(newallele->allele_seq, seq);
		newallele->allele_num = alnum;
		newallele->allele_info.ac = -1;
		newallele->allele_info.af = -1;
		strcpy(newallele->allele_info.vt, "");
		newallele->next = NULL;
	}
	
	return newallele;
}

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
				pallele->allele_info.ac, pallele->allele_info.af, pallele->allele_info.vt);
		pallele = pallele->next;
	}
	printf("SAMPLES\n");
	while (psample != NULL)
	{
		printf("%d%c%d\t", psample->gt.m, psample->phased ? '|' : '/', psample->gt.p);
		psample = psample->next;
	}
	putchar('\n');
	putchar('\n');
}

static bool locus_is_valid(const VCF_LOCUS *plocus)
{
	// No conditions.
	return true;
}

