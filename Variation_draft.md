# Describing Variation
Nan Hu / Feb. 4, 2020

---

In molecular population genetics, we tend to focus on those statistics that are estimators of the theoretical parameter θ. Chapter 3 offers us several common estimators and mathematical formula to calculate them. However, in practical researches, we usually face huge molecular population data in genomic level, which require us to use some efficient approaches to reach our goals.

This document provides some methods to calculate variation parameters in real population genetics data sets.

---
## Aims:
1. using vcftools to calculate per site pairwise differences (π) and windowed π
2. using angsd to calculate per site π and windowed π
3. using angsd to generate site frequency spectrum (SFS)
4. plot results
---
## Calculating π using 'vcftools'

