# Population Structure
Nan Hu / Feb. 25, 2021

---

(description)

---
## Aims:
1. Calculating Fst between males and females using 'vcftools'
2. Calculating Fst between *S. nivalis* and *S. reticulata* using 'vcftools'
3. Identifying Fst outliers
4. Calculating d<sub>XY</sub> between males and females
5. Plotting Fst manhattan plots, d<sub>XY</sub> plots
---
## Calculating Fst between males and females using 'vcftools'
Population structure focuses on comparisons among subpopulations with different isolation and migration histories. When inferring population structure, Fst is an estimator as the degree of measuring Wahlund effect. In our population data, we do not really have samples of subpopulations from the same species. Instead, luckily, all of our species are dioecious which having males and females. We could treat males and females as two different subpopulations and they certainly having some different patterns in selection and lifestyle. 
```bash
#!/bin/bash
#SBATCH -J fst
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 16

vcftools --vcf snivalis.vcf --weir-fst-pop male.list --weir-fst-pop female.list --out snivalis

```
The output will have a suffix of `.weir.fst`. It contains columns of positions and Fst estimation.

## Calculating Fst between *S. nivalis* and *S. reticulata* using 'vcftools'



