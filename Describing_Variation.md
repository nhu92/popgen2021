# Describing Variation
Nan Hu / Feb. 4, 2020

---

In molecular population genetics, we tend to focus on those statistics that are estimators of the theoretical parameter θ. Chapter 3 offers us several common estimators and mathematical formula to calculate them. However, in practical researches, we usually face huge molecular population data in genomic level, which require us to use some efficient approaches to reach our goals.

This document provides some methods to calculate variation parameters in real population genetics data sets.

---
## Aims:
1. using vcftools to calculate per site pairwise differences (π) and windowed π
2. using angsd to calculate other θ estimators
3. using angsd to generate site frequency spectrum (SFS)
4. plotting results
---
## Calculating π using 'vcftools'
As described in manual of vcftools software, this package offers two mode to calculate π for `.vcf` file - per site π and windowed π.
### Per site π calculation
For per site π, vcftools will calculate pairwise differences based on polymorphic sites marked in `.vcf` files. It will skip sites that do not polymorphic in our population data sets.
```bash
#!/bin/bash
#SBATCH -J persitepi_vts
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 32

vcftools --vcf <vcf file> --site-pi --out <output file name>
```
Output file has the suffix of `.sites.pi`. It is a tablized text file with 3 columns - chromosome ID, position, and π, respectively.

### Windowed π calculation
However, per site π often fluctuates a lot so that it is very hard to figure out a diversity pattern in large scales. In order to see diversity level across a region, a chromosome or even the entire genome, we need stastistics to represent the diversity in a region. Moving average (will be described in 'plotting diversity') and diversity in windows are two ways to approach it.

Vcftools offers `--window-pi` function to generate π in pre-set windows. The size of windows is the only parameter need to be considered. It depends on genome size of target species, SNP density, resequencing methods (whole genome sequencing will cover most regions in genome; sequence capture is less 'random' covered than RAD-seq), etc. As a example here, I chose 25kb as the size of window but you can play around it (increase it into 50kb, 100kb ..., or decrease it into 20kb, 15kb, 10kb ...) to see different effects on estimating diversity.
```bash
#!/bin/bash
#SBATCH -J windowpi_vts
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 32

vcftools --vcf <vcf file> --window-pi 25000 --out <output file name>
```
Output file has the suffix of `.windowed.pi`. It is also a tablized text file with 6 columns - chromosome ID, position, bin start, bin end, number of segregating sites in this bin, windowed π. 

Notice here, the 'N_VARIANTS' column is the number of segregating sites in this window. It looks like θw·a - part of the estimator of diversity. However, it is not. The reason is that our data set is from sequence capture array. We only sampled genome regions that we have probes on. Thus, there will be some region contain SNPs but we just do not have reads covered on. Then, the number of segregating sites in this table is not the real segregating site number for each window. If your data is generated through whole genome resequencing and you have got a decent reference genome, then the 'N_VARIANTS' column might be perfect to calculate θw (by divided by a = 1/1 + 1/2 + 1/3 + ... + 1/Ne; assuming no ascertainment bias).

## Generate SFS and calculating diversity estimators using 'angsd'
> angsd requires `.bam` files at this stage. I am going back to generate these files. It might takes one extra day. (still waiting the BAM files, should be done in the evening of Saturday, Feb. 6, 2021 - Nan) It is done now - Nan 5:35PM, Feb. 6.

'angsd' takes `.bam` files to calculate site frequency spectrum (SFS) first. From the result of SFS (an `.idx` file), it estimate several estimators of diversity (θ). Since it takes `.bam` files, it may run longer than vcftools does, but normally less than 1 hour.
> This step requires the latest version of 'angsd'. I already installed under our course directory `/lustre/scratch/nhu/popgen2021/software/latest/angsd/`. You can directly use it without installing your own.
> bam file list are stored in `/lustre/scratch/nhu/popgen2021/scripts/diversity/`. There are `snivalis.bamlist`, `sreticulata.bamlist`, and `sphlebophylla.bamlist`. Choose the one that assigned to your group as the first input file for command below.
> We did a simple filter here to screen genotypes with minimum mapping quality and Phred quality. For vcftools we use `.vcf` file, they are already filtered so we do not need to do this.
```bash
#!/bin/bash
#SBATCH -J idx_angsd
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128

/lustre/scratch/nhu/popgen2021/software/latest/angsd/angsd -bam <bam.filelist> -doSaf 1 -anc <reference genome> -GL 1 -P 128 -minMapQ 30 -minQ 20 -out <output file name>
```
This command generate an `.idx` file record genotype likelihood and estimate of per site statistics of diversity. It is an intermediate file we will use for following analysis.

'angsd' contains a powerful package 'realSFS' to use `.idx` file to generate folded SFS and unfolded SFS. We suggest to generate folded SFS here because our reference is not the ancient states but share the same common ancestor.
```bash
#!/bin/bash
#SBATCH -J realsfs_angsd
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 64

# folded SFS
/lustre/scratch/nhu/popgen2021/software/latest/angsd/misc/realSFS [.idx file] -P 64 -fold 1 > [output file, .sfs]
/lustre/scratch/nhu/popgen2021/software/latest/angsd/misc/realSFS saf2theta [.idx file] -outname [output file name] -sfs [.sfs file] -fold 1

# unfolded SFS
/lustre/scratch/nhu/popgen2021/software/latest/angsd/misc/realSFS [.idx file] -P 64 > [output file, .sfs]
/lustre/scratch/nhu/popgen2021/software/latest/angsd/misc/realSFS saf2theta [.idx file] -outname [output file name] -sfs [.sfs file]
```
Output will be a file with suffix of `.thetas.idx`. Then, 'angsd' has another package 'thetaStat' to extract diversity estimation from SFS:
```bash
#!/bin/bash
#SBATCH -J theta_angsd
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 64

/lustre/scratch/nhu/popgen2021/software/latest/angsd/misc/thetaStat do_stat <.thetas.idx file>
```
'angsd' can also do sliding window estimation in diversity values. Here we still use 25k as the window size as an example. Steps are equal to window size to make un-overlapping windows.
```bash
#!/bin/bash
#SBATCH -J windowtheta_angsd
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 64

/lustre/scratch/nhu/popgen2021/software/latest/angsd/misc/thetaStat do_stat <.thetas.idx file> -win 25000 -step 25000 -outnames <output file name>
```

## plotting θπ, θw, SFS, and other estimators
See [Plotting Diversity](https://github.com/gudusanjiao/popgen2021/blob/main/plotting_diversity.R).






