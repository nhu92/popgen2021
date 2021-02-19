# Recombination
Nan Hu / Feb. 7, 2021

---

(description)

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


## Calculating r<sup>2</sup> using 'vcftools'


## Claculate LD decay using 'ngsLD'

