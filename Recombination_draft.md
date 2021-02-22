# Recombination
Nan Hu / Feb. 7, 2021

---

In this exercise, we want to practice how to calculate recombination coefficients. In order to make results interpretable, we would like to compare them across different genomic regions, different chromosomes, or different genders. To begin with, there are some basic introductions about gender information in our data sets. Also, we will learn how to subset our data sets to work on specific regions and individuals.

---
## Aims:
1. Introducing males and females in datasets
2. Using vcftools to calculate r<sup>2</sup> statistics for specific region
3. Using ngsLD to calculate LD decay for specific region
4. Plotting r<sup>2</sup> metrics, LD decay curve
---
## Introducing males and females in datasets
### Dioecy in *Salicaceae*
Dioecy means having males and females individuals separately. In angiosperms, only 7% of species are dioecious. Specifically, we could find individuals with only male flowers (only have androecium include stamen, filament, anthur, etc.) and individuals with only female flowers (only have gynoecium include style, pistil, carpel, etc.) in a population. In our datasets, we nearly evenly sampled males and females for population study. For all of our three species, their sex determination system are ZW on chromosome 15 (we will learn how to find it in later chapters), which means for sex chromosome, males are homozygotes (ZZ) while females are heterozygotes (ZW). 

Knowing the gender of our data could bring some special perspectives in analyzing and interpretting results. For example, our previous diversity calculation on π and θ<sub>w</sub> could be done separately for males and females only. Major expected results will be high diversity level on female ZW due to high heterozygosities while low diversity level on male ZZ due to low effective population size ( Ne(Z) = 3/4 Ne(Autosome) ).

### How to do analysis based on part of our data
In practical research, we often want to subset our data to compare or study the part we are interested. In vcftools, it is very convenient to subset our data by using `--positions`, `--from-bp`, `--to-bp`, `-chr`, `--keep`, `--remove`, `--recode` parameters.
```bash
# Here are some examples for using vcftools to subset our data:
## calculate site pi for Chr05
vcftools --vcf snivalis.vcf --site-pi --chr Chr05 --out snivalis

## calculate windowed pi for Chr07 from 4Mbp to 7Mbp, window size is 50k
vcftools --vcf snivalis.vcf --window-pi 50000 --chr Chr07 --from-bp 4000000 --to-bp 7000000 --out snivalis

## subset vcf files to include only Chr15W
vcftools --vcf snivalis.vcf --chr Chr15W --recode --out snivalis.Chr15W
## the output will be snivalis.Chr15W.recode.vcf, which only has Chr15W.

## calculating hwe statistics only for Chr15W in males
vcftools --vcf snivalis.vcf --chr Chr15W --keep male.list --hardy --out snivalis
## male.list file is a pre-made text file. It contains all the males' names with each name in a line.

## subset vcf files that only include individuals with clear gender identifications
vcftools --vcf snivalis.vcf --remove nogender.list --recode --out snivalis.clean
## nogender.list is also a pre-made text file with all the names that do not have a gender identification.
```

For this exercise and later ones in this semester, we will use these command very often. For more usages and more fancy filters, you could visit [VCF manual](http://vcftools.sourceforge.net/man_latest.html).
## Calculating r<sup>2</sup> using 'vcftools'
In vcftools, we have `--geno-r2`, `--hap-r2`, `--geno-r2-positions`, `--hap-r2-positions` options to calculate r<sup>2</sup> in different situations. If our data are phased (means we know which alleles are in the same haplotype), we can calculate haplotype r<sup>2</sup>. For our data, since we sequenced with sequence capture technique, the genotype calling is less phased. Currently, we could only use `--geno-r2` or `--geno-r2-positions` to calculate r<sup>2</sup>.


## Calculate LD decay using 'ngsLD'

