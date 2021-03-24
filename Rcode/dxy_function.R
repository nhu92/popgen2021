# dxy functions
rm(list = ls())
library(tidyverse)
library(PopGenome)
library(zoo)

dxy_dataframe <- function(vcf_path, population_list, chrlength, window_size, window_jump) {
  species <- readData(vcf_path,
                       format = "VCF",
                       include.unknown = TRUE,
                       FAST = TRUE)
  species_info <- read_delim(population_list, delim = "\t")
  
  # get the data for the popultions
  populations <- split(species_info$ind, species_info$pop)
  species <- set.populations(species, populations, diploid = T)
  # Setting up sliding windows
  # use seq to find the start points of each window
  window_start <- seq(from = 1, to = chrlength, by = window_jump)
  # add the size of the window to each start point 
  window_stop <- window_start + window_size
  
  # # no windows start before the end of chromosome 15
  # sum(window_start > chrlength)
  # # but some window stop positions do occur past the final point
  # sum(window_stop > chrlength)
  # 
  # remove windows from the start and stop vectors
  window_start <- window_start[which(window_stop < chrlength)]
  window_stop <- window_stop[which(window_stop < chrlength)]

  # save as a data.frame
  windows <- data.frame(start = window_start, stop = window_stop, 
                        mid = window_start + (window_stop-window_start)/2)
  
  # make a sliding window dataset
  species_sw <- sliding.window.transform(species,
                                          width = window_size,
                                          jump = window_jump,
                                          type = 2)
  
  # Calculating sliding window estimates of nucleotide diversity and differentiation
  
  # calculate diversity statistics
  species_sw <- diversity.stats(species_sw , pi = TRUE)
  
  species_sw <- F_ST.stats(species_sw , mode = "nucleotide")
  
  # extract nucleotide diversity and correct for window size
  nd <- species_sw@nuc.diversity.within/window_size
  
  # make population name vector
  pops <- names(populations)
  # set population names
  colnames(nd) <- paste0(pops, "_pi")
  
  # extract fst values
  fst <- t(species_sw @nuc.F_ST.pairwise)
  
  # extract dxy - pairwise absolute nucleotide diversity
  dxy <- get.diversity(species_sw , between = T)[[2]]/window_size
  
  # get column names 
  x <- colnames(fst)
  for (i in seq(1,length(pops))){
    x <- sub(paste0("pop",i), pops[i], x)
  }
  x <- sub("/", "_", x)
  
  colnames(fst) <- paste0(x, "_fst")
  colnames(dxy) <- paste0(x, "_dxy")
  
  windows <- head(windows, length(nd[,1]))
  return(as_tibble(data.frame(windows, nd, fst, dxy)))
}

# test run:
# vcf_path is a path to all the vcfs we want to include in
# population_list is a two-column tab-spaced list. 
# First column with header 'ind' refers to individuals name;
# Second column with header 'pop' refers to population it belongs to.

sr_table <- dxy_dataframe("sret", "sret_pops.txt", 15651726, 100000, 25000)



# rolling averages
sr_table$roll22<-rollmean(sr_table$Sret_F_Sret_M_dxy, 22, fill=NA, na.rm=TRUE)
rollQuantile.9975<-quantile(sr_table$roll22,0.9975,fill=NA,na.rm=TRUE)
rollQuantile.0025<-quantile(sr_table$roll22,0.0025,fill=NA,na.rm=TRUE)
rollQuantile.99<-quantile(sr_table$roll22,0.99,fill=NA,na.rm=TRUE)
rollQuantile.01<-quantile(sr_table$roll22,0.01,fill=NA,na.rm=TRUE)

head(sr_table,40)

## plot
ggplot(sr_table, aes(x=mid/10^6, y=Sret_F_Sret_M_dxy)) +
  geom_point(alpha=0.8, size=0.5) +
  scale_color_manual(values = rep("grey", 19)) +
  xlab("Genome Pos on Chr15W") +
  ylab(expression(italic(d)[xy])) +
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
  geom_line(aes(y=rollmean(sr_table$Sret_F_Sret_M_dxy, 22, fill=NA, na.rm=TRUE))) +
  geom_hline(yintercept=rollQuantile.01, color="red", size=0.25) +
  geom_hline(yintercept=rollQuantile.99, color="red", size=0.25)

# extract pairwise dxy

dxy_matrix <- sr_data %>% select(contains("dxy")) %>% summarise_all(mean)
write.table(dxy_matrix,"dxytable.txt")
