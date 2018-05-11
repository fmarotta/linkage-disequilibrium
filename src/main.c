#include <stdlib.h>
#include <stdio.h>
#include "ld_vcf.h"

#define R2_CUTOFF 0 // value under which we shall not print anything
#define WINLEN 10000 // length of the window, in bases.

int main(int argc, char *argv[])
{
	// TODO soft-code this parameters, maybe through command-line
	// options
	const int winlen = WINLEN;
	const float r2_cutoff = R2_CUTOFF;

	FILE *vcf_file;
	VCF_WINDOW window;
	VCF_LOCUS *plocus1, *plocus2;
	float p_AB, p_A, p_B;
	float D, D_lewontin, r_squared;

	if (argc != 2)
	{
		fprintf(stderr, "USAGE: %s <vcf_file>", argv[0]);
		exit(EXIT_FAILURE);
	}
	if ((vcf_file = fopen(argv[1], "r")) == NULL)
	{
		fprintf(stderr, "ERROR: could not read VCF: %s", argv[1]);
		exit(EXIT_FAILURE);
	}

	// Ensure that at least two loci are present in the window.
	Initialize_window(&window, vcf_file, winlen);
	while (window.nloci < 2 && !window.eow)
		Slide_window(&window);

	while (window.nloci >= 2)
	{
		plocus2 = plocus1 = window.head;

		// the loci must be biallelic in order for our formulae to work
		if (Nalleles_in_locus(plocus1) <= 2)
		while (plocus2 != window.tail)
		{
			plocus2 = plocus2->next;

			if (Nalleles_in_locus(plocus2) <= 2)
			for (int i = 0; i < Nalleles_in_locus(plocus1); i++)
				for (int j = 0; j < Nalleles_in_locus(plocus2); j++)
				{
					p_A = Allele_freq(i, plocus1);
					p_B = Allele_freq(j, plocus2);
					p_AB = Linked_alleles_freq(i, plocus1, j, plocus2);
					D = Calculate_D(p_A, p_B, p_AB);
					D_lewontin = Calculate_D_lewontin(p_A, p_B, p_AB); // a.k.a. D'
					r_squared = Calculate_r_squared(p_A, p_B, p_AB);
					if (r_squared >= r2_cutoff)
						printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\tD=%f\tD'=%f\tr^2=%f\n",
								i, plocus1->pos, j, plocus2->pos, p_A, p_B, p_AB,
								D, D_lewontin, r_squared);
				}
		}

		Slide_window(&window);
		// If we opened a window with less than two loci, we try again.
		while (window.nloci < 2 && !window.eow)
			Slide_window(&window);
	}

	Close_window(&window);
	fclose(vcf_file);

	return 0;
}
