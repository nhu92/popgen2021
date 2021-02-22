library(dplyr)
library(ggplot2)

# -------------------------
# Function: plotPairwiseLD 
# take LD output from vcftools to draw a pairwise LD plot
# -------------------------
plotPairwiseLD <- function(dfr,chr,xlim=c(NA,NA),ylim=c(NA,NA),minr2) {
  if(missing(dfr)) stop("Input data.frame 'dfr' missing.")
  
  if(!missing(chr)) {
    ld <- filter(ld,CHROM_A==get("chr") & CHROM_B==get("chr"))
  }
  ld <- filter(ld,ld$POS_A<ld$POS_B)
  
  if(!missing(minr2)) {
    ld <- filter(ld,ld$R2>get("minr2"))
  }
  
  ld <- ld %>% arrange(ld$R2)
  
  ld$x <- ld$POS_A+((ld$POS_B-ld$POS_A)/2)
  ld$y <- ld$POS_B-ld$POS_A
  ld$r2c <- cut(ld$R2,breaks=seq(0,1,0.2),labels=c("0-0 - 0.2","0.2 - 0.4",
                                                   "0.4 - 0.6","0.6 - 0.8",
                                                   "0.8 - 1.0"))
  ld$r2c <- factor(ld$r2c,levels=rev(c("0-0 - 0.2","0.2 - 0.4",
                                       "0.4 - 0.6","0.6 - 0.8",
                                       "0.8 - 1.0")))
  
  ggplot(ld,aes(x=x,y=y,col=r2c))+
    geom_point(shape=20,size=0.1,alpha=0.8)+
# I have three examples of color system you may chose to display the result,
# you could have your own to make it nicer.
    # scale_color_manual(values=c("#0000ff","#0066ee","#ffcccc","#ffdddd","#efef00"))+
    # scale_color_manual(values=c("black","grey70","grey90","white","white"))+
    scale_color_manual(values=c("#ff0000","#efef00","#ffff1a","#ffff1a","#ffff1a"))+
    scale_x_continuous(limits=xlim)+
    scale_y_continuous(limits=ylim)+
    guides(colour=guide_legend(title="R2",override.aes=list(shape=20,size=8)))+
    labs(x="Chromosome (Bases)",y="")+
    theme_bw(base_size=14)+
    theme(panel.border=element_blank(),
          axis.ticks=element_blank()) %>%
    return()
}

# SDR plot
## load r2 output from vcftools
rinitial <- read.delim("snivalis.male.chr15.sdr.geno.ld")
colnames(rinitial)=c("chr","POS_A","POS_B","N_INDV","R2")
ld <- rinitial
g <-plotPairwiseLD(ld)

# #you could change the dimension of the plot
ggsave(filename ="R2plot_sdr.png", plot = g, height=45, width= 45)

# PAR1 plot
rinitial <- read.delim("snivalis.male.chr15.par1.geno.ld")
colnames(rinitial)=c("chr","POS_A","POS_B","N_INDV","R2")
ld <- rinitial
g <-plotPairwiseLD(ld)
ggsave(filename ="R2plot_par1.png", plot = g, height=45, width= 45)

# PAR2 plot
rinitial <- read.delim("snivalis.male.chr15.par2.geno.ld")
colnames(rinitial)=c("chr","POS_A","POS_B","N_INDV","R2")
ld <- rinitial
g <-plotPairwiseLD(ld)
ggsave(filename ="R2plot_par2.png", plot = g, height=45, width= 45)

