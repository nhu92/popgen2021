# library packages. If not installed, use install.packages() to install them
library(ggplot2)
library(tidyr)
library(dplyr)

# read output file from vcftools
hwe_tab <- read.delim("E:/OneDrive/OneDrive - Texas Tech University/2021Spring/BIOL6301/sexigua.hwe")

# remove NAs
hwe_tab_nohomozygous <- na.omit(hwe_tab)

# split genotype column into 3 genotypes
hwe_extract_geno_freq <- hwe_tab_nohomozygous %>% 
  separate(col = OBS.HOM1.HET.HOM2., into = c("a1a1", "a1a2", "a2a2"), sep = "/")
hwe_extract_geno_freq$a1a1 <- as.numeric(hwe_extract_geno_freq$a1a1)
hwe_extract_geno_freq$a1a2 <- as.numeric(hwe_extract_geno_freq$a1a2)
hwe_extract_geno_freq$a2a2 <- as.numeric(hwe_extract_geno_freq$a2a2)

# calculate allele frequencies
hwe_extract_allele_freq <- hwe_extract_geno_freq %>% 
  mutate(a1_freq = (a1a1 + a1a2 / 2) / (a1a1 + a1a2 + a2a2)) %>% 
  mutate(a2_freq = (a2a2 + a1a2 / 2) / (a1a1 + a1a2 + a2a2)) %>%
  mutate(hom_a1_freq = a1a1 / (a1a1 + a1a2 + a2a2)) %>%
  mutate(hom_a2_freq = a2a2 / (a1a1 + a1a2 + a2a2)) %>%
  mutate(het_freq = a1a2 / (a1a1 + a1a2 + a2a2))

# mark hwe outliers/equilibrium loci
hwe_mark_def <- hwe_extract_allele_freq %>%
  filter(P_HET_DEFICIT < 0.05 & P_HWE < 0.05) %>%
  mutate(het_status = "Heterozygous Deficiency")
hwe_mark_exc <- hwe_extract_allele_freq %>%
  filter(P_HET_EXCESS < 0.05 & P_HWE < 0.05) %>%
  mutate(het_status = "Heterozygous Excess")
hwe_mark_equi <- hwe_extract_allele_freq %>%
  filter(P_HWE >= 0.05) %>%
  mutate(het_status = "Equilibrium")

# combine and reshape data, keep useful column
hwe_clean <- rbind(hwe_mark_exc, hwe_mark_def, hwe_mark_equi) %>%
  select(CHR, POS, a1a1, a1a2, a2a2, a1_freq, a2_freq, hom_a1_freq, hom_a2_freq, het_freq, het_status)

# generate a standard equilibrium curve dataframe
hwe_expected_AF <- seq(0.0, 1.0, 0.01)
hwe_expected_Hom_AF <- hwe_expected_AF ^ 2
hwe_expected <- data.frame(hwe_expected_AF, hwe_expected_Hom_AF)

# plotting (you can change the title/color/size of the dot by yourself)
ggplot() + 
  geom_point(data = hwe_clean, aes(a1_freq, hom_a1_freq, col = het_status)) +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  scale_x_continuous(name = "Allele Frequency (AF)") +
  scale_y_continuous(name = "Homozygous Allele Frequency (Hom-AF)") +
  labs(color = "") +
  geom_line(data = hwe_expected, aes(hwe_expected_AF, hwe_expected_Hom_AF), size = 1) +
  theme_bw()

