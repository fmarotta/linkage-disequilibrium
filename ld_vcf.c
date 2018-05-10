/* Interface implementation */

// TODO close fds and free malloc'd memory.

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ld_vcf.h"
#include "utils.h"

#define LLEN 30000 // max line length
#define FLEN 100 // max field length

#define DIGEST_NEXT_FIELD() \
	t = 0; \
	while (!isspace(line[l]) && t < FLEN-1) \
		tmpfield[t++] = line[l++]; \
	tmpfield[t] = '\0'; \
	while (isspace(line[l])) \
		l++;

static VCF_LOCUS buflocus;

static bool locus_is_in_window(const VCF_LOCUS *plocus, const VCF_WINDOW *pwindow);
static bool make_allele(VCF_ALLELE *tmpallele, const char *seq, int alnum);
static bool locus_is_valid(const VCF_LOCUS *plocus);
static bool digest_line(VCF_LOCUS *plocus, const char *line);
static void vomit_line(const VCF_LOCUS *plocus);
static char *digest_next_subfield(char *subfield, char *field);

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
	char line[LLEN];

	// Initialize the window
	pwindow->head = pwindow->tail = NULL;
	pwindow->nloci = 0;
	pwindow->winlen = winlen;

	// Discard comment and header lines
	//
	// FIXME there is no need to read the whole line, the first two
	// characters are enough. Detect the change from ## to #.
	while (fgets(line, LLEN, vcf_file) != NULL && line[0] == '#')
		continue;
	if (line[0] == '#')
	{
		fputs("ERROR: cannot find data in vcf file.", stderr);
		exit(EXIT_FAILURE);
	}

	// Digest the first data line into the one-locus buffer
	//
	// FIXME rewrite digest_line() using fscanf. Make fields like
	// tmpinfo. Use symbolic constants to "hard code" the strlengths.
	if (!digest_line(&buflocus, line))
	{
		fputs("ERROR: we ran out of memory.", stderr);
		exit(EXIT_FAILURE);
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
		//
		// FIXME we could be run out of file as well... Maybe this
		// function should return a buffer struct. At the end of each
		// iteration of the while below, we try to digest and inside the
		// body we check if the eof member is set to true and decide
		// whether to break the loop. OR, we check if we have run out of
		// file during the fgets...
		if (fgets(line, LLEN, vcf_file) == NULL)
			break;
		if (!digest_line(&buflocus, line))
		{
			// if buf.eof == true then break;
			fputs("ERROR: we ran out of memory.", stderr);
			exit(EXIT_FAILURE);
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
	if (pnew == NULL) return false;
	*pnew = locus;

	if (pwindow->head == NULL)
		pwindow->head = pnew;
	else
		pwindow->tail->next = pnew;
	pwindow->tail = pnew;
	pwindow->nloci++;

	return true;
}

static bool digest_line(VCF_LOCUS *plocus, const char *line)
{
	VCF_ALLELE *tmpallele, *last_allele;
	VCF_INFO *tmpinfo;
	VCF_SAMPLE *tmpsample;
	VCF_LOCUS *last_locus;
	char tmpfield[FLEN], tmpsubfield[FLEN];
	int l = 0, t;

	// Read the chromosome
	DIGEST_NEXT_FIELD();
	plocus->chrom = atoi(tmpfield);

	// Read the position
	DIGEST_NEXT_FIELD();
	plocus->pos = atoul(tmpfield);

	// Read the id
	DIGEST_NEXT_FIELD();
	plocus->id = (char *) malloc(strlen(tmpfield)+1);
	strcpy(plocus->id, tmpfield);

	// Read the REF seq
	//
	// Remember to add allele info when you read the INFO field
	DIGEST_NEXT_FIELD();
	int tmpalnum = 0;
	tmpallele = (VCF_ALLELE *) malloc(sizeof(VCF_ALLELE));
	if (tmpallele == NULL)
		return false;
	if (!make_allele(tmpallele, tmpfield, tmpalnum++))
		return false;
	plocus->alleles = tmpallele;
	last_allele = plocus->alleles;

	// Read the ALT seq
	DIGEST_NEXT_FIELD();
	char *ptmpfield = tmpfield;
	while ((ptmpfield = digest_next_subfield(tmpsubfield, ptmpfield)) != NULL)
	{
		tmpallele = (VCF_ALLELE *) malloc(sizeof(VCF_ALLELE));
		if (tmpallele == NULL)
			return false;
		if (!make_allele(tmpallele, tmpsubfield, tmpalnum++))
			return false;
		last_allele = last_allele->next = tmpallele;
	}

	// Read the quality
	DIGEST_NEXT_FIELD();
	plocus->qual = atoi(tmpfield);

	// TODO Read the info

	vomit_line(plocus);

	return true;
}

static bool make_allele(VCF_ALLELE *tmpallele, const char *seq, int alnum)
{
	tmpallele->allele_seq = (char *) malloc(sizeof(char) * (strlen(seq)+1));
	if (tmpallele->allele_seq == NULL)
		return false;
	strcpy(tmpallele->allele_seq, seq);
	tmpallele->allele_num = alnum;
	tmpallele->allele_info.ac = tmpallele->allele_info.af = 0;
	strcpy(tmpallele->allele_info.vt, "");
	tmpallele->next = NULL;

	return true;
}

static char *digest_next_subfield(char *subfield, char *field)
{
	int i = 0;

	if (field[i] == '\0')
		return NULL;

	while (field[i] != ',' && field[i] != '\0')
	{
		subfield[i] = field[i];
		i++;
	}
	subfield[i] = '\0';
	field += i;

	return field;
}

static void vomit_line(const VCF_LOCUS *plocus)
{
	printf("%d\t%ld\t%s\t%s\t%s\t%d\n", plocus->chrom, plocus->pos, plocus->id,
			plocus->alleles->allele_seq, plocus->alleles->next->allele_seq,
			plocus->qual);
}

static bool locus_is_valid(const VCF_LOCUS *plocus)
{
	// No conditions.
	return true;
}

