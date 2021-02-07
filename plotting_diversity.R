library('gridExtra')
library('ggplot2')
library('cowplot')
theme_set(theme_cowplot())
library('magrittr')
library('dplyr')
require("zoo")

###################################################
# Plotting windowed Pi using output from vcftools #
###################################################
Pi <- read.csv("snivalis.windowed.pi", header=T, na.strings=c("", "NA"), sep="\t", stringsAsFactors=F)

## REMOVE NaN & SITES WITH < 5SNPs
Pi$PI[Pi$N_VARIANTS<5]<-NaN

head(Pi,40)
## remove scaffolds
Pi <- Pi %>% filter(grepl("^Chr", CHROM))

## calculate quantiles for rolling means for the entire genome
Pi$roll22<-rollmean(Pi$PI, 22, fill=NA, na.rm=TRUE)
rollQuantile.9975<-quantile(Pi$roll22,0.9975,fill=NA,na.rm=TRUE)
rollQuantile.0025<-quantile(Pi$roll22,0.0025,fill=NA,na.rm=TRUE)
rollQuantile.99<-quantile(Pi$roll22,0.99,fill=NA,na.rm=TRUE)
rollQuantile.01<-quantile(Pi$roll22,0.01,fill=NA,na.rm=TRUE)

head(Pi,40)

## plot
ggplot(Pi, aes(x=BIN_START, y=PI)) +
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=0.5) +
  facet_wrap(~ CHROM, nrow = 4) +
  scale_color_manual(values = rep("grey", 19)) +
  xlab("Genome Pos") +
  ylab("Pairwise Differences") +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_rect(fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(size=12), 
    axis.title.y = element_text(size=12), 
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(size=8) 
  ) +
  geom_line(aes(y=rollmean(Pi$PI, 22, fill=NA, na.rm=TRUE))) +
  geom_hline(yintercept=rollQuantile.01, color="red", size=0.25) +
  geom_hline(yintercept=rollQuantile.99, color="red", size=0.25)

###################################################
# Plotting windowed theta using output from angsd #
###################################################
theta <- read.csv("snivalis.25kwindow.pestPG", header=T, na.strings=c("", "NA"), sep="\t", stringsAsFactors=F)

## remove scaffolds
theta <- theta %>% filter(grepl("^Chr", Chr))

## create function related to theta_w (not used)
a_stats <- function(nsites) {
  sum <- 0
  for (i in seq(1, nsites)){
    sum <- sum + 1/i
  }
  return(sum)
}

## calculate quantiles for rolling means for the entire genome
theta$roll22<-rollmean(theta$tW, 22, fill=NA, na.rm=TRUE)
rollQuantile.9975<-quantile(theta$roll22,0.9975,fill=NA,na.rm=TRUE)
rollQuantile.0025<-quantile(theta$roll22,0.0025,fill=NA,na.rm=TRUE)
rollQuantile.99<-quantile(theta$roll22,0.99,fill=NA,na.rm=TRUE)
rollQuantile.01<-quantile(theta$roll22,0.01,fill=NA,na.rm=TRUE)

## plot
ggplot(theta, aes(x=WinCenter, y=tW)) +
  geom_point(aes(color=as.factor(Chr)), alpha=0.8, size=0.5) +
  facet_wrap(~ Chr, nrow = 4) +
  scale_color_manual(values = rep("grey", 19)) +
  xlab("Genome Pos") +
  ylab("Theta") +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_rect(fill=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(size=12), 
    axis.title.y = element_text(size=12), 
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(size=8) 
  ) + 
  geom_line(aes(y=rollmean(theta$tW, 22, fill=NA, na.rm=TRUE))) +
  geom_hline(yintercept=rollQuantile.01, color="red", size=0.25) +
  geom_hline(yintercept=rollQuantile.99, color="red", size=0.25)

########################################
# Plotting SFS using output from angsd #
########################################

#read data
sfs<-scan("snivalis.fold.sfs")

norm <- function(x) x/sum(x)

#the variable categories of the sfs
sfs <- norm(sfs[-c(1, length(sfs))]) 

# convert to dataframe
a_freq <- (0:(length(sfs)-1)/(length(sfs)-1))
sfs <- data.frame(a_freq, sfs)
sfs_df <- sfs %>% filter(sfs != 0)

ggplot(sfs_df, aes(a_freq, sfs)) +
  geom_bar(stat = "identity", fill = "orange") +
  xlab("Allele Frequency") +
  ylab("Proportions of SNPs")

