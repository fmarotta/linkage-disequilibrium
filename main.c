#include <stdlib.h>
#include <stdio.h>
#include "ld_vcf.h"
#include "io_utils.h"

/* I am quite pleased with the result, though I still do not understand
 * why for some alleles pair I get D' or r^3 greater than 1, and by a
 * lot... However, it seems that of the four paris, two always get valid
 * D' and r^2, namely alleles 0-1 and 1-1. All the same, its two in the
 * morning and I haven't been doing anything other than this all day, so
 * instead of thinking about it, I am going to bed.
 */

int main(void)
{
	FILE *vcf_file;
	VCF_WINDOW window;
	VCF_LOCUS *plocus1, *plocus2;
	float p_AB, p_A, p_B;
	float D, D_lewontin, r_squared;

	if ((vcf_file = fopen("data/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", "r")) == NULL)
	{
		fputs("error reading vcf file", stderr);
		exit(EXIT_FAILURE);
	}

	Initialize_window(vcf_file, &window, WINLEN);
	//printf("head: %d, tail: %d, nloci: %d\n", window.head->pos,
	//window.tail->pos, window.nloci);
	while (window.nloci > 0)
	{
		plocus1 = window.head;
		plocus2 = plocus1->next;

		// the loci must be biallelic for our formulae to work
		if (Nalleles_in_locus(plocus1) <= 2)
		while (plocus2 != window.tail)
		{
			if (Nalleles_in_locus(plocus2) <= 2)
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
				{
					p_AB = Linked_alleles_freq(i, plocus1, j, plocus2);
					p_A = Allele_freq(i, plocus1);
					p_B = Allele_freq(j, plocus2);
					D = Calculate_D(p_A, p_B, p_AB);
					D_lewontin = Calculate_D_lewontin(p_A, p_B, p_AB);
					r_squared = Calculate_r_squared(p_A, p_B, p_AB);
					if (r_squared >= 0)
						printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\tD=%f\tD'=%f\tr^2=%f\n",
								i, plocus1->pos, j, plocus2->pos, p_A, p_B, p_AB,
								D, D_lewontin, r_squared);
				}
			plocus2 = plocus2->next;
		}

		Slide_window(vcf_file, &window);
		//printf("head: %d, tail: %d, nloci: %d\n", window.head->pos,
		//window.tail->pos, window.nloci);
	}

	fclose(vcf_file);
	puts("Done");
	return 0;
}
