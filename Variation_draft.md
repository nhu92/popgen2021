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

## Calculating diversity estimators using 'angsd'
> angsd requires `.bam` files at this stage. I am going back to generate these files. It might takes one extra day.





