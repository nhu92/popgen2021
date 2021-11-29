# HKA test
# adapted & edited from: https://github.com/eacooper400/sfsr/blob/master/Using_SFSR_w_VCF.md
# Original Creator of sfsr: Andrew Parker Morgan, 2016 https://github.com/andrewparkermorgan/sfsr
# Nan Hu
# Apr. 02, 2021

# install packages and load them
library(devtools)
install_github("eacooper400/sfsr")
library(sfsr)

# read input vcf
my.data=read.vcf("snivalis.vcf", header=TRUE, stringsAsFactors=FALSE)

#### Example run ####

# choose loci under neutral (random gene from chr06 that I assume it is under neutral evolution)
loc_neut = subset(my.data, (my.data$CHROM == 'Chr06' & my.data$POS>= 4794350 & my.data$POS<= 4798181))

# choose loci under selection (random genen from Chr15W inside SDR that I assume it is under selection)
loc_sel = subset(my.data, (my.data$CHROM == 'Chr15W' & my.data$POS>=7499612 & my.data$POS<=7505735))

# calculate SFS
sfs_neut=sfs_neut=sfs.fromVcf(list(loc_neut))
sfs_sel=sfs.fromVcf(list(loc_sel))

# performing HKA test
hka_test(sfs_neut, sfs_sel)
