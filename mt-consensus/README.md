# mt-consensus.py

This is a program for examining mtDNA sequencing data to try to identify positions with discordant reads - i.e., DNA from different haplotypes.
The program takes a bam file with reads aligned to the mtDNA genome, and outputs read counts for each base at each position.

It also attempts to identify discordant bases that could be the result of aDNA damage, and reports updated base counts with those bases removed.
Specifically, if the majority allele is a C or a G, then all T or A bases observed on the forward or reverse strands (respectively) are removed.

Example command:

`python mt-consensus/bin/mt-consensus.py -b data/Cap.K.4297.Hominidae.Homo_sapiens_deduped.bam -c NC_012920.1`
