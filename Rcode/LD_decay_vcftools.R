# Libraries used in this notebook (R 3.6.2)
#from https://github.com/BrianSanderson/salix-nigra-slr/blob/master/notebooks/Salix-nigra-SLR-04-PopGen.ipynb
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(corrplot))
suppressMessages(library(reshape2))
suppressMessages(library(viridis))
suppressMessages(library(zoo))
suppressMessages(library(nlstools))
theme_set(theme_cowplot()); theme_update(plot.title = element_text(hjust = 0.5))


# read in LD values
# assuming SDR from 4M to 8M
r2_frame <- bind_rows(read_delim('snivalis.male.chr15.par1.geno.ld', delim="\t",  na =c("", "NA", "-nan")) %>%
                        mutate(., Sex = "Male", Region = "PAR1"),
                      read_delim('snivalis.male.chr15.4m_8m.geno.ld', delim="\t",  na =c("", "NA", "-nan")) %>%
                        mutate(., Sex = "Male", Region = "SDR"),
                      read_delim('snivalis.male.chr15.par2.geno.ld', delim="\t",  na =c("", "NA", "-nan")) %>%
                        mutate(., Sex = "Male", Region = "PAR2"))

# calculate distance between loci
r2_frame <- group_by(r2_frame, Sex, Region) %>%
  mutate(., dist = POS2 - POS1) %>%
  ungroup()

head(r2_frame,20)
tail(r2_frame,20)
#The following function implements the code described by Fabio Marroni in this blog post(https://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/), and described in greater detail in their 
#publication from 2011. It takes as input a vector of distances in base pairs, values of $r^2$ from vcftools, the number of chromosomes 
#(in our case with a minimum of 18 individuals I'm using 36), and a starting value for C (to be estimated by the nls model), which is the 
#product of the recombination fraction between sites and distance in bp.

estimate_decay <- function(data_frame, sex, region, n, HW.st) {
  require(nlstools)
  require(dplyr)
  distance = filter(data_frame, Sex == sex, Region == region) %>% 
    select(., dist) %>% as_vector() %>% unname()
  LD.data = filter(data_frame, Sex == sex, Region == region) %>% 
    select(., `R^2`) %>% as_vector() %>% unname()
  HW.nonlinear<-nls(LD.data~((10+C*distance)/
                               ((2+C*distance)*(11+C*distance)))*
                      (1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/
                         (n*(2+C*distance)*(11+C*distance))),
                    start=HW.st,control=nls.control(maxiter=100))
  tt<-summary(HW.nonlinear)
  new.rho<-tt$parameters[1]
  rho.ci <- confint2(HW.nonlinear, level = 0.95)
  fpoints<-((10+new.rho*distance)/
              ((2+new.rho*distance)*(11+new.rho*distance)))*
    (1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/
       (n*(2+new.rho*distance)*(11+new.rho*distance)))
  fpoints.lower <- ((10+rho.ci[1]*distance)/
                      ((2+rho.ci[1]*distance)*(11+rho.ci[1]*distance)))*
    (1+((3+rho.ci[1]*distance)*(12+12*rho.ci[1]*distance+(rho.ci[1]*distance)^2))/
       (n*(2+rho.ci[1]*distance)*(11+rho.ci[1]*distance)))
  fpoints.upper <- ((10+rho.ci[2]*distance)/
                      ((2+rho.ci[2]*distance)*(11+rho.ci[2]*distance)))*
    (1+((3+rho.ci[2]*distance)*(12+12*rho.ci[2]*distance+(rho.ci[2]*distance)^2))/
       (n*(2+rho.ci[2]*distance)*(11+rho.ci[2]*distance)))
  LD.decay<-data.frame(distance, LD.data, fpoints, fpoints.lower, fpoints.upper)
  return(LD.decay)
}

#Use the above function to estimate the decay of LD across the regions
estimate_frame <- bind_rows(estimate_decay(r2_frame, "Male", "PAR1", 36, c(C=0.1)) %>% 
                              mutate(., Sex = "Male", Region = "PAR1"),
                            estimate_decay(r2_frame, "Male", "SDR", 36, c(C=0.1)) %>% 
                              mutate(., Sex = "Male", Region = "SDR"),
                            estimate_decay(r2_frame, "Male", "PAR2", 36, c(C=0.1)) %>% 
                              mutate(., Sex = "Male", Region = "PAR2"))
#The following function also implements the code described by Fabio Marroni in this blog post, and calculates 
#the distance at which the estimates of LD have decayed by half.

half_distance <- function(data_frame, sex, region) {
  require(dplyr)
  sub_frame <- filter(data_frame, Sex == sex, Region == region)
  h.decay <- max(sub_frame$fpoints) / 2
  half.decay.distance <- sub_frame$distance[which.min(abs(sub_frame$fpoints - h.decay))]
  return(half.decay.distance)
}

half_distance(estimate_frame, "Male", "PAR1")
half_distance(estimate_frame, "Male", "SDR")
half_distance(estimate_frame, "Male", "PAR2")

#Calculate the maximum distance across the SLRs for males and females, to be used as the max X-axis in the plot of decay
filter(estimate_frame, Sex == "Male", Region == "SDR") %>% select(., distance) %>% max()


#Create a ggplot object to plot the decay of LD
#LINES <- c("SLR" = "solid", "Upstream" = "dashed", "Downstream" = "dotted" )

ld_plot <- estimate_frame %>% 
  mutate(.,Region = factor(Region, levels = c("SDR", "PAR1", "PAR2"))) %>%
  ggplot(., aes(x=distance, y=fpoints, linetype = Region)) +
  geom_ribbon(aes(ymin=fpoints.lower, ymax=fpoints.upper, size = 0.01), alpha = 0.3) +
  geom_line(aes(x=distance, y=fpoints, linetype = Region)) +
# replace limits with maximum distance calculated above
  scale_x_continuous(limits = c(0, 3995225),
                     labels = c(0, 1, 2, 3, 4),
                     breaks = c(0, 1000000, 2000000, 3000000, 4000000)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  xlab("Distance (MB)") + 
  ylab(expression(r^2)) +
  theme(legend.position="none",
        plot.margin = margin(10, 24, 10, 10),
        axis.title.x = element_text(size=36),
        axis.title.y = element_text(size=36),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24))
#head(ld_plot)
# options(repr.plot.width = 6.3, repr.plot.height = 3.15, repr.plot.res = 300)
#dev.off()

ggsave(filename="LD_Decay_color.png",plot = ld_plot, width = 36, height = 36)
