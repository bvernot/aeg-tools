# pmd-filter-bam.py

`[super-alpha]`

A tool to reduce the effects of ancient DNA damage in reads that overlap SNPs of interest.
This program takes a bam file and a SNP file, and outputs a filtered bam file.
It identifies reads that overlap transition SNPs, and if the overlapping base could the result of ancient DNA damage,
then it is either masked or the read is removed entirely [`--mode mask | --mask filter`].
There is also the option to only filter/mask bases in the first and last X bases [`--X 3`].
Currently, the program assumes single-stranded libraries, and performs strand-aware detection of possible DNA damage.

For example, consider a C/T SNP (marked with an S below), and several overlapping reads
mapped to the forward strand (caps) or reverse strand (lower case):

```
   S
  ATCGGAGTACAGT...
GCACCGGAGTACAGT...

  atcggagtacagt...
gcaccggagtacagt...
   S
```

In this example, _only_ the first read would be masked or filtered, as only this read could be the result of aDNA damage.

### Notes:
 - Every read is processed, so it is much faster to pre-filter your bam file to only reads overlapping SNPs in the SNP file.
 

### Example command:

This command reads a set of snps (-s), filters reads from the bam file (-b), and saves them to a new bam file (-o).

```
python bin/pmd-filter-bam.py \
    --X 3 \
    -b ../data/Lib.F.9795.uniq.L35MQ25.bam \
    -s ../data/twist.1240k.burden.anno.txt \
    -o hey.bam \
    --mode filter
```

You can use `[ -h ]` to get a full description of options.
