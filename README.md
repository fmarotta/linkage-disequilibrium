# Linkage Disequilibrium

I wished to implement a simple program to calculate linkage 
disequilibrium from a VCF file. It is a very basic one as its purpose is 
merely didactical. Far be it from me to say this is a finished product 
-- it is not even a product, but only an exercise.

-----

Linkage disequilibrium can be measured in many ways, the most used of 
which are the following three:

* D = p\_AB - p\_A p\_B;
* D' = D / D\_max;
* r^2 = D^2 / (p\_A(1 - p\_A) p\_B(1 - p\_B)).

We just have to take every pair of alleles A and B, and compute p\_AB, 
p\_A and p\_B. Doing so literaly for every pair of alleles can be very 
slow, besides being useless (indeed LD decreases rather sharply with 
distance \[see e.g. 
[this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3587626)\]). 
Therefore we only consider a window, which by the way acts as the core 
of our ADT. Actually, the window is a queue, that is a FIFO data 
structure. 

We read the VCF first into a sort of buffer, then we decide whether to 
add the locus to the window or not, according to configurable filters 
(one being the distance from the first locus in the window). Should the 
locus pass the filters, it is added to the window. The buffer is needed 
as an intermediate storage not to lose records while we inch through the 
file. We ensure that at least two loci are present in the window so that 
we can parse the entire window computing the linkage, then slide over 
the loci, and repeat.

Each line is tokenized and digested into a structure; notably, the 
alleles and the samples are implemented as linked lists, since their 
number can vary. However, to compute linkage disequilibrium we use 
formulae which work only for biallelic loci.

There is still one thing that bothers me. If we have biallelic loci, 
then there are four possible pairs of alleles for each two loci. 
Nevertheless, most other programs I have seen provide a single number 
for each pair of loci. Moreover, many times our measures seem wrong, for 
the r^2 is greater than 1. For every pair of locus, though, there are 
two combination which are always right: alleles 0-1 and 1-1. Please, do 
contact me if you know the reason why. It may well be that I made a 
mistake, since it is almost three in the morning and I have been working 
non-stop since this morning.

As of the next day, this doesn't bother me anymore, for I realised my 
mistake... However, I still do not understand which of the four pairs is 
the relevant one, if any.

All in all, I think the purpose of this project has been accomplished; 
in doing it, I have (a) revised the awesome C language and some of its 
advanced features; (b) learned how to parse a file, in particular a VCF; 
(c) understood linkage disequilibrium at last, which now is not just a 
formula any more. I hope this experience will enable me to tackle 
similar and more difficult problems with more self-confidence in the 
future.
