
library(ggplot2)
library(zoo)
library(scales)
library(patchwork)
library(RColorBrewer)
library(tibble)
library(extrafont)
library(onewaytests)
library(boot)
library(car)
library(dplyr)


# check read depth and snp density across all chromosomes...

setwd("~/Desktop/Plots/FASTZ")

#--------------------------------------#
##### read in data and restructure #####
#--------------------------------------#

#snps...
swe_m_all_snps <- read.csv("P14502_103-filtered_200_snps_per_100kb_windows.txt", sep = " ")
swe_f_all_snps <- read.csv("P14502_104-filtered_200_snps_per_100kb_windows.txt", sep = " ")
mor_m_all_snps <- read.csv("P14458_108_S29-filtered_200_snps_per_100kb_windows.txt", sep = " ")
mor_f_all_snps <- read.csv("P14458_107_S28-filtered_200_snps_per_100kb_windows.txt", sep = " ")
dup_m_all_snps <- read.csv("P14458_104_S25-filtered_200_snps_per_100kb_windows.txt", sep = " ")
dup_f_all_snps <- read.csv("P14458_103_S24-filtered_200_snps_per_100kb_windows.txt", sep = " ")

colnames(swe_m_all_snps)[colnames(swe_m_all_snps) == "snp_count"] = "male_snp_count"
colnames(swe_f_all_snps)[colnames(swe_f_all_snps) == "snp_count"] = "female_snp_count"
colnames(mor_m_all_snps)[colnames(mor_m_all_snps) == "snp_count"] = "male_snp_count"
colnames(mor_f_all_snps)[colnames(mor_f_all_snps) == "snp_count"] = "female_snp_count"
colnames(dup_m_all_snps)[colnames(dup_m_all_snps) == "snp_count"] = "male_snp_count"
colnames(dup_f_all_snps)[colnames(dup_f_all_snps) == "snp_count"] = "female_snp_count"

# calculate normalization due to differens in snp count between samples...
swe_m_all_snps$male_snp_count_scaled <- swe_m_all_snps$male_snp_count/median(swe_m_all_snps$male_snp_count)
swe_f_all_snps$female_snp_count_scaled <- swe_f_all_snps$female_snp_count/median(swe_f_all_snps$female_snp_count)
mor_m_all_snps$male_snp_count_scaled <- mor_m_all_snps$male_snp_count/median(mor_m_all_snps$male_snp_count)
mor_f_all_snps$female_snp_count_scaled <- mor_f_all_snps$female_snp_count/median(mor_f_all_snps$female_snp_count)
dup_m_all_snps$male_snp_count_scaled <- dup_m_all_snps$male_snp_count/median(dup_m_all_snps$male_snp_count)
dup_f_all_snps$female_snp_count_scaled <- dup_f_all_snps$female_snp_count/median(dup_f_all_snps$female_snp_count)

# add 0.1 to allow division with 0...
swe_m_all_snps$male_snp_count_scaled <- swe_m_all_snps$male_snp_count_scaled+0.1
swe_f_all_snps$female_snp_count_scaled <- swe_f_all_snps$female_snp_count_scaled+0.1
mor_m_all_snps$male_snp_count_scaled <- mor_m_all_snps$male_snp_count_scaled+0.1
mor_f_all_snps$female_snp_count_scaled <- mor_f_all_snps$female_snp_count_scaled+0.1
dup_m_all_snps$male_snp_count_scaled <- dup_m_all_snps$male_snp_count_scaled+0.1
dup_f_all_snps$female_snp_count_scaled <- dup_f_all_snps$female_snp_count_scaled+0.1

# merge m/f data...
swe_all_snps <- merge(swe_m_all_snps, swe_f_all_snps, by=c("chr", "start", "end"))
mor_all_snps <- merge(mor_m_all_snps, mor_f_all_snps, by=c("chr", "start", "end"))
dup_all_snps <- merge(dup_m_all_snps, dup_f_all_snps, by=c("chr", "start", "end"))




#depth...
swe_m <- read.csv("LsinapisSweM-P14502_103-sorted.bam-unique.cov", sep = " ")
swe_f <- read.csv("LsinapisSweM-P14502_104-sorted.bam-unique.cov", sep = " ")
mor_m <- read.csv("LsinapisSweM-P14458_108_S29-sorted.bam-unique.cov", sep = " ")
mor_f <- read.csv("LsinapisSweM-P14458_107_S28-sorted.bam-unique.cov", sep = " ")
dup_m <- read.csv("LsinapisSweM-P14458_104_S25-sorted.bam-unique.cov", sep = " ")
dup_f <- read.csv("LsinapisSweM-P14458_103_S24-sorted.bam-unique.cov", sep = " ")


# plot variables...
#Z1_col="#89ad4e"
#Z2_col="#8e551e"
#Z3_col="#efbf5a"

Z1_col <- rgb(213/255,117/255,0)
Z2_col <- rgb(185/255,156/255,107/255)
Z3_col <- rgb(219/255,202/255,105/255)


auto_color="grey"

label_sizes=20
axis_text_sizes=12
x_max=0.8
x_min=-0.8
y_max=3
y_min=-3



#--------------------#
##### L. sinapis #####
#--------------------#

# snps...

# calculate ratio and mean...
swe_all_snps$snp_ratio <- log2(swe_all_snps$male_snp_count_scaled) - log2(swe_all_snps$female_snp_count_scaled)
swe_mf_snps_chr <- aggregate(snp_ratio ~ chr, data = swe_all_snps, FUN = mean)


# Bootstrap CI for snps...
swe_mf_snp_ratio_function <- function(swe_all_snps, index) {
  swe_mf_snp_ratio = mean(swe_all_snps[index, 8])
  return(swe_mf_snp_ratio)
}

swe_mf_snp_ratio_CI.res = data.frame()

for (i in 1:29) {
  
  swe_mf_snp_ratio.res <- boot(data = swe_all_snps[swe_all_snps$chr == paste0("HiC_scaffold_", i), ],
                               statistic = swe_mf_snp_ratio_function,
                               R = 10000)
  
  swe_mf_snp_ratio_CI <- boot.ci(boot.out = swe_mf_snp_ratio.res, type = "perc")
  print(swe_mf_snp_ratio_CI)
  
  output = c(paste0("HiC_scaffold_", i), swe_mf_snp_ratio_CI$percent[[4]], swe_mf_snp_ratio_CI$percent[[5]])
  
  swe_mf_snp_ratio_CI.res = rbind(swe_mf_snp_ratio_CI.res, output)
  
}

colnames(swe_mf_snp_ratio_CI.res) <- c("chr", "snp_CI_low", "snp_CI_high")

swe_mf_snp_ratio_CI.res$snp_CI_low <- as.numeric(swe_mf_snp_ratio_CI.res$snp_CI_low)
swe_mf_snp_ratio_CI.res$snp_CI_high <- as.numeric(swe_mf_snp_ratio_CI.res$snp_CI_high)



# depth ratio and mean...
swe_m_gene_median = aggregate(depth ~ gene_id + chr, data=swe_m, median)
swe_m_gene_median_auto <- swe_m_gene_median[swe_m_gene_median$chr != "HiC_scaffold_1" & swe_m_gene_median$chr != "HiC_scaffold_6" & swe_m_gene_median$chr != "HiC_scaffold_18", ]
swe_m_gene_mean = mean(swe_m_gene_median_auto$depth)
swe_m_gene_median$norm_median <- swe_m_gene_median$depth/swe_m_gene_mean
swe_m_gene_median$norm_median <- swe_m_gene_median$norm_median+0.1

swe_f_gene_median = aggregate(depth ~ gene_id + chr, data=swe_f, median)
swe_f_gene_median_auto <- swe_f_gene_median[swe_f_gene_median$chr != "HiC_scaffold_1" & swe_f_gene_median$chr != "HiC_scaffold_6" & swe_f_gene_median$chr != "HiC_scaffold_18", ]
swe_f_gene_mean = mean(swe_f_gene_median_auto$depth)
swe_f_gene_median$norm_median <- swe_f_gene_median$depth/swe_f_gene_mean
swe_f_gene_median$norm_median <- swe_f_gene_median$norm_median+0.1

swe_mf_gene_median <- merge(swe_m_gene_median, swe_f_gene_median, by = c("chr", "gene_id"))
swe_mf_gene_median$depth_ratio <- log2(swe_mf_gene_median$norm_median.x) - log2(swe_mf_gene_median$norm_median.y)
swe_mf_depth_chr <- aggregate(depth_ratio ~ chr, swe_mf_gene_median, FUN = mean)


# Bootstrap CI for depth...
swe_mf_depth_ratio_function <- function(swe_mf_gene_median, index) {
  swe_mf_depth_ratio = mean(swe_mf_gene_median[index, 7])
  return(swe_mf_depth_ratio)
}

swe_mf_depth_ratio_CI.res = data.frame()

for (i in 1:29) {
  
  swe_mf_depth_ratio.res <- boot(data = swe_mf_gene_median[swe_mf_gene_median$chr == paste0("HiC_scaffold_", i), ],
                                 statistic = swe_mf_depth_ratio_function,
                                 R = 10000)
  
  swe_mf_depth_ratio_CI <- boot.ci(boot.out = swe_mf_depth_ratio.res, type = "perc")
  print(swe_mf_depth_ratio_CI)
  
  output = c(paste0("HiC_scaffold_", i), swe_mf_depth_ratio_CI$percent[[4]], swe_mf_depth_ratio_CI$percent[[5]])
  
  swe_mf_depth_ratio_CI.res = rbind(swe_mf_depth_ratio_CI.res, output)
  
}

colnames(swe_mf_depth_ratio_CI.res) <- c("chr", "depth_CI_low", "depth_CI_high")

swe_mf_depth_ratio_CI.res$depth_CI_low <- as.numeric(swe_mf_depth_ratio_CI.res$depth_CI_low)
swe_mf_depth_ratio_CI.res$depth_CI_high <- as.numeric(swe_mf_depth_ratio_CI.res$depth_CI_high)


# merge data...
swe_mf_snp_depth_chr <- merge(swe_mf_snps_chr, swe_mf_depth_chr, by = c("chr"))
swe_mf_snp_depth_chr <- merge(swe_mf_snp_depth_chr, swe_mf_snp_ratio_CI.res, by = "chr")
swe_mf_snp_depth_chr <- merge(swe_mf_snp_depth_chr, swe_mf_depth_ratio_CI.res, by = "chr")

swe_mf_snp_depth_chr <- swe_mf_snp_depth_chr %>%  mutate(Z_vs_A = case_when(
                                       chr == 'HiC_scaffold_1' ~ 'Z1',
                                       chr == 'HiC_scaffold_2' ~ 'A',
                                       chr == 'HiC_scaffold_3' ~ 'A',
                                       chr == 'HiC_scaffold_4' ~ 'A',
                                       chr == 'HiC_scaffold_5' ~ 'A',
                                       chr == 'HiC_scaffold_6' ~ 'Z2',
                                       chr == 'HiC_scaffold_7' ~ 'A',
                                       chr == 'HiC_scaffold_8' ~ 'A',
                                       chr == 'HiC_scaffold_9' ~ 'A',
                                       chr == 'HiC_scaffold_10' ~ 'A',
                                       chr == 'HiC_scaffold_11' ~ 'A',
                                       chr == 'HiC_scaffold_12' ~ 'A',
                                       chr == 'HiC_scaffold_13' ~ 'A',
                                       chr == 'HiC_scaffold_14' ~ 'A',
                                       chr == 'HiC_scaffold_15' ~ 'A',
                                       chr == 'HiC_scaffold_16' ~ 'A',
                                       chr == 'HiC_scaffold_17' ~ 'A',
                                       chr == 'HiC_scaffold_18' ~ 'Z3',
                                       chr == 'HiC_scaffold_19' ~ 'A',
                                       chr == 'HiC_scaffold_20' ~ 'A',
                                       chr == 'HiC_scaffold_21' ~ 'A',
                                       chr == 'HiC_scaffold_22' ~ 'A',
                                       chr == 'HiC_scaffold_23' ~ 'A',
                                       chr == 'HiC_scaffold_24' ~ 'A',
                                       chr == 'HiC_scaffold_25' ~ 'A',
                                       chr == 'HiC_scaffold_26' ~ 'A',
                                       chr == 'HiC_scaffold_27' ~ 'A',
                                       chr == 'HiC_scaffold_28' ~ 'A',
                                       chr == 'HiC_scaffold_29' ~ 'A'))

# Z3 gets overplotted, split data in Z/A and plot Z on top of A to see Z3
swe_mf_snp_depth_chr_A <- swe_mf_snp_depth_chr[swe_mf_snp_depth_chr$Z_vs_A == "A", ]
swe_mf_snp_depth_chr_Z <- swe_mf_snp_depth_chr[swe_mf_snp_depth_chr$Z_vs_A != "A", ]

swe_plot <- ggplot(swe_mf_snp_depth_chr_A, aes(x=depth_ratio, y=snp_ratio, color=Z_vs_A)) +
  geom_point(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size=label_sizes),
        axis.text = element_text(size=label_sizes),
        axis.title = element_text(size=label_sizes),
        legend.title = element_blank(),
        legend.text = element_text(size=label_sizes)) +
  geom_linerange(aes(x=depth_ratio, ymin=snp_CI_low, ymax=snp_CI_high), linewidth=1) +
  geom_linerange(aes(y=snp_ratio, xmin=depth_CI_low, xmax=depth_CI_high), linewidth=1) +
  scale_color_manual(values=c(auto_color, Z1_col, Z2_col, Z3_col)) +
  coord_cartesian(ylim = c(y_min, y_max), xlim = c(x_min, x_max)) +
  xlab("\nlog2 M/F read depth") +
  ylab("log2 M/F snp density\n") +
  ggtitle("L. sinapis") +
  geom_point(data=swe_mf_snp_depth_chr_Z, aes(x=depth_ratio, y=snp_ratio), size=2) +
  geom_linerange(data=swe_mf_snp_depth_chr_Z, aes(x=depth_ratio, ymin=snp_CI_low, ymax=snp_CI_high), linewidth=1) +
  geom_linerange(data=swe_mf_snp_depth_chr_Z, aes(y=snp_ratio, xmin=depth_CI_low, xmax=depth_CI_high), linewidth=1)

#-------------------#
##### L. morsei #####
#-------------------#

# snps...

# calculate ratio and mean...
mor_all_snps$snp_ratio <- log2(mor_all_snps$male_snp_count_scaled) - log2(mor_all_snps$female_snp_count_scaled)
mor_mf_snps_chr <- aggregate(snp_ratio ~ chr, data = mor_all_snps, FUN = mean)


# Bootstrap CI for snps...
mor_mf_snp_ratio_function <- function(mor_all_snps, index) {
  mor_mf_snp_ratio = mean(mor_all_snps[index, 8])
  return(mor_mf_snp_ratio)
}

mor_mf_snp_ratio_CI.res = data.frame()

for (i in 1:29) {
  
  mor_mf_snp_ratio.res <- boot(data = mor_all_snps[mor_all_snps$chr == paste0("HiC_scaffold_", i), ],
                               statistic = mor_mf_snp_ratio_function,
                               R = 10000)
  
  mor_mf_snp_ratio_CI <- boot.ci(boot.out = mor_mf_snp_ratio.res, type = "perc")
  print(mor_mf_snp_ratio_CI)
  
  output = c(paste0("HiC_scaffold_", i), mor_mf_snp_ratio_CI$percent[[4]], mor_mf_snp_ratio_CI$percent[[5]])
  
  mor_mf_snp_ratio_CI.res = rbind(mor_mf_snp_ratio_CI.res, output)
  
}

colnames(mor_mf_snp_ratio_CI.res) <- c("chr", "snp_CI_low", "snp_CI_high")

mor_mf_snp_ratio_CI.res$snp_CI_low <- as.numeric(mor_mf_snp_ratio_CI.res$snp_CI_low)
mor_mf_snp_ratio_CI.res$snp_CI_high <- as.numeric(mor_mf_snp_ratio_CI.res$snp_CI_high)


# depth ratio and mean...
mor_m_gene_median = aggregate(depth ~ gene_id + chr, data=mor_m, median)
mor_m_gene_median_auto <- mor_m_gene_median[mor_m_gene_median$chr != "HiC_scaffold_1" & mor_m_gene_median$chr != "HiC_scaffold_6" & mor_m_gene_median$chr != "HiC_scaffold_18", ]
mor_m_gene_mean = mean(mor_m_gene_median$depth)
mor_m_gene_median$norm_median <- mor_m_gene_median$depth/mor_m_gene_mean
mor_m_gene_median$norm_median <- mor_m_gene_median$norm_median+0.1

mor_f_gene_median = aggregate(depth ~ gene_id + chr, data=mor_f, median)
mor_f_gene_median_auto <- mor_f_gene_median[mor_f_gene_median$chr != "HiC_scaffold_1" & mor_f_gene_median$chr != "HiC_scaffold_6" & mor_f_gene_median$chr != "HiC_scaffold_18", ]
mor_f_gene_mean = mean(mor_f_gene_median$depth)
mor_f_gene_median$norm_median <- mor_f_gene_median$depth/mor_f_gene_mean
mor_f_gene_median$norm_median <- mor_f_gene_median$norm_median+0.1

mor_mf_gene_median <- merge(mor_m_gene_median, mor_f_gene_median, by = c("chr", "gene_id"))
mor_mf_gene_median$depth_ratio <- log2(mor_mf_gene_median$norm_median.x) - log2(mor_mf_gene_median$norm_median.y)
mor_mf_depth_chr <- aggregate(depth_ratio ~ chr, mor_mf_gene_median, FUN = mean)


# Bootstrap CI for depth...
mor_mf_depth_ratio_function <- function(mor_mf_gene_median, index) {
  mor_mf_depth_ratio = mean(mor_mf_gene_median[index, 7])
  return(mor_mf_depth_ratio)
}

mor_mf_depth_ratio_CI.res = data.frame()

for (i in 1:29) {
  
  mor_mf_depth_ratio.res <- boot(data = mor_mf_gene_median[mor_mf_gene_median$chr == paste0("HiC_scaffold_", i), ],
                                 statistic = mor_mf_depth_ratio_function,
                                 R = 10000)
  
  mor_mf_depth_ratio_CI <- boot.ci(boot.out = mor_mf_depth_ratio.res, type = "perc")
  print(mor_mf_depth_ratio_CI)
  
  output = c(paste0("HiC_scaffold_", i), mor_mf_depth_ratio_CI$percent[[4]], mor_mf_depth_ratio_CI$percent[[5]])
  
  mor_mf_depth_ratio_CI.res = rbind(mor_mf_depth_ratio_CI.res, output)
  
}

colnames(mor_mf_depth_ratio_CI.res) <- c("chr", "depth_CI_low", "depth_CI_high")

mor_mf_depth_ratio_CI.res$depth_CI_low <- as.numeric(mor_mf_depth_ratio_CI.res$depth_CI_low)
mor_mf_depth_ratio_CI.res$depth_CI_high <- as.numeric(mor_mf_depth_ratio_CI.res$depth_CI_high)


# merge data...
mor_mf_snp_depth_chr <- merge(mor_mf_snps_chr, mor_mf_depth_chr, by = c("chr"))
mor_mf_snp_depth_chr <- merge(mor_mf_snp_depth_chr, mor_mf_snp_ratio_CI.res, by = "chr")
mor_mf_snp_depth_chr <- merge(mor_mf_snp_depth_chr, mor_mf_depth_ratio_CI.res, by = "chr")

mor_mf_snp_depth_chr <- mor_mf_snp_depth_chr %>% mutate(Z_vs_A = case_when(
      chr == 'HiC_scaffold_1' ~ 'Z1',
      chr == 'HiC_scaffold_2' ~ 'A',
      chr == 'HiC_scaffold_3' ~ 'A',
      chr == 'HiC_scaffold_4' ~ 'A',
      chr == 'HiC_scaffold_5' ~ 'A',
      chr == 'HiC_scaffold_6' ~ 'Z2',
      chr == 'HiC_scaffold_7' ~ 'A',
      chr == 'HiC_scaffold_8' ~ 'A',
      chr == 'HiC_scaffold_9' ~ 'A',
      chr == 'HiC_scaffold_10' ~ 'A',
      chr == 'HiC_scaffold_11' ~ 'A',
      chr == 'HiC_scaffold_12' ~ 'A',
      chr == 'HiC_scaffold_13' ~ 'A',
      chr == 'HiC_scaffold_14' ~ 'A',
      chr == 'HiC_scaffold_15' ~ 'A',
      chr == 'HiC_scaffold_16' ~ 'A',
      chr == 'HiC_scaffold_17' ~ 'A',
      chr == 'HiC_scaffold_18' ~ 'Z3',
      chr == 'HiC_scaffold_19' ~ 'A',
      chr == 'HiC_scaffold_20' ~ 'A',
      chr == 'HiC_scaffold_21' ~ 'A',
      chr == 'HiC_scaffold_22' ~ 'A',
      chr == 'HiC_scaffold_23' ~ 'A',
      chr == 'HiC_scaffold_24' ~ 'A',
      chr == 'HiC_scaffold_25' ~ 'A',
      chr == 'HiC_scaffold_26' ~ 'A',
      chr == 'HiC_scaffold_27' ~ 'A',
      chr == 'HiC_scaffold_28' ~ 'A',
      chr == 'HiC_scaffold_29' ~ 'A'))

# Z3 gets overplotted, split data in Z/A and plot Z on top of A to see Z3...
mor_mf_snp_depth_chr_A <- mor_mf_snp_depth_chr[mor_mf_snp_depth_chr$Z_vs_A == "A", ]
mor_mf_snp_depth_chr_Z <- mor_mf_snp_depth_chr[mor_mf_snp_depth_chr$Z_vs_A != "A", ]

mor_plot <- ggplot(mor_mf_snp_depth_chr_A, aes(x=depth_ratio, y=snp_ratio, color=Z_vs_A)) +
  geom_point(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size=label_sizes),
        axis.text = element_text(size=label_sizes),
        axis.title = element_text(size=label_sizes),
        legend.title = element_blank(),
        legend.text = element_text(size=label_sizes)) +
  geom_linerange(aes(x=depth_ratio, ymin=snp_CI_low, ymax=snp_CI_high), linewidth=1) +
  geom_linerange(aes(y=snp_ratio, xmin=depth_CI_low, xmax=depth_CI_high), linewidth=1) +
  scale_color_manual(values=c(auto_color, Z1_col, Z2_col, Z3_col)) +
  scale_alpha_manual(guide='none', values = list(a = 0.2, point = 1)) +
  coord_cartesian(ylim = c(y_min, y_max), xlim = c(x_min, x_max)) +
  xlab("\nlog2 M/F read depth") +
  ylab("log2 M/F snp density\n") +
  ggtitle("L. morsei") +
  geom_point(data=mor_mf_snp_depth_chr_Z, aes(x=depth_ratio, y=snp_ratio), size=2) +
  geom_linerange(data=mor_mf_snp_depth_chr_Z, aes(x=depth_ratio, ymin=snp_CI_low, ymax=snp_CI_high), linewidth=1) +
  geom_linerange(data=mor_mf_snp_depth_chr_Z, aes(y=snp_ratio, xmin=depth_CI_low, xmax=depth_CI_high), linewidth=1) 


#-----------------------#
##### L. duponcheli #####
#-----------------------#

# snps...

# calculate ratios and means...
dup_all_snps$snp_ratio <- log2(dup_all_snps$male_snp_count_scaled) - log2(dup_all_snps$female_snp_count_scaled)
dup_mf_snps_chr <- aggregate(snp_ratio ~ chr, data = dup_all_snps, FUN = mean)


# Bootstrap CI for snps...
dup_mf_snp_ratio_function <- function(dup_all_snps, index) {
  dup_mf_snp_ratio = mean(dup_all_snps[index, 8])
  return(dup_mf_snp_ratio)
}

dup_mf_snp_ratio_CI.res = data.frame()

for (i in 1:29) {
  
  dup_mf_snp_ratio.res <- boot(data = dup_all_snps[dup_all_snps$chr == paste0("HiC_scaffold_", i), ],
                               statistic = dup_mf_snp_ratio_function,
                               R = 10000)
  
  dup_mf_snp_ratio_CI <- boot.ci(boot.out = dup_mf_snp_ratio.res, type = "perc")
  print(dup_mf_snp_ratio_CI)
  
  output = c(paste0("HiC_scaffold_", i), dup_mf_snp_ratio_CI$percent[[4]], dup_mf_snp_ratio_CI$percent[[5]])
  
  dup_mf_snp_ratio_CI.res = rbind(dup_mf_snp_ratio_CI.res, output)
  
}

colnames(dup_mf_snp_ratio_CI.res) <- c("chr", "snp_CI_low", "snp_CI_high")

dup_mf_snp_ratio_CI.res$snp_CI_low <- as.numeric(dup_mf_snp_ratio_CI.res$snp_CI_low)
dup_mf_snp_ratio_CI.res$snp_CI_high <- as.numeric(dup_mf_snp_ratio_CI.res$snp_CI_high)


# depth ratios and means...
dup_m_gene_median = aggregate(depth ~ gene_id + chr, data=dup_m, median)
dup_m_gene_median_auto <- dup_m_gene_median[dup_m_gene_median$chr != "HiC_scaffold_1" & dup_m_gene_median$chr != "HiC_scaffold_6" & dup_m_gene_median$chr != "HiC_scaffold_18", ]
dup_m_gene_mean = mean(dup_m_gene_median$depth)
dup_m_gene_median$norm_median <- dup_m_gene_median$depth/dup_m_gene_mean
dup_m_gene_median$norm_median <- dup_m_gene_median$norm_median+0.1

dup_f_gene_median = aggregate(depth ~ gene_id + chr, data=dup_f, median)
dup_f_gene_median_auto <- dup_f_gene_median[dup_f_gene_median$chr != "HiC_scaffold_1" & dup_f_gene_median$chr != "HiC_scaffold_6" & dup_f_gene_median$chr != "HiC_scaffold_18", ]
dup_f_gene_mean = mean(dup_f_gene_median$depth)
dup_f_gene_median$norm_median <- dup_f_gene_median$depth/dup_f_gene_mean
dup_f_gene_median$norm_median <- dup_f_gene_median$norm_median+0.1

dup_mf_gene_median <- merge(dup_m_gene_median, dup_f_gene_median, by = c("chr", "gene_id"))
dup_mf_gene_median$depth_ratio <- log2(dup_mf_gene_median$norm_median.x) - log2(dup_mf_gene_median$norm_median.y)
dup_mf_depth_chr <- aggregate(depth_ratio ~ chr, dup_mf_gene_median, FUN = mean)


# Bootstrap CI for depth...
dup_mf_depth_ratio_function <- function(dup_mf_gene_median, index) {
  dup_mf_depth_ratio = mean(dup_mf_gene_median[index, 7])
  return(dup_mf_depth_ratio)
}

dup_mf_depth_ratio_CI.res = data.frame()

for (i in 1:29) {
  
  dup_mf_depth_ratio.res <- boot(data = dup_mf_gene_median[dup_mf_gene_median$chr == paste0("HiC_scaffold_", i), ],
                                 statistic = dup_mf_depth_ratio_function,
                                 R = 10000)
  
  dup_mf_depth_ratio_CI <- boot.ci(boot.out = dup_mf_depth_ratio.res, type = "perc")
  print(dup_mf_depth_ratio_CI)
  
  output = c(paste0("HiC_scaffold_", i), dup_mf_depth_ratio_CI$percent[[4]], dup_mf_depth_ratio_CI$percent[[5]])
  
  dup_mf_depth_ratio_CI.res = rbind(dup_mf_depth_ratio_CI.res, output)
  
}

colnames(dup_mf_depth_ratio_CI.res) <- c("chr", "depth_CI_low", "depth_CI_high")

dup_mf_depth_ratio_CI.res$depth_CI_low <- as.numeric(dup_mf_depth_ratio_CI.res$depth_CI_low)
dup_mf_depth_ratio_CI.res$depth_CI_high <- as.numeric(dup_mf_depth_ratio_CI.res$depth_CI_high)


# merge data...
dup_mf_snp_depth_chr <- merge(dup_mf_snps_chr, dup_mf_depth_chr, by = c("chr"))
dup_mf_snp_depth_chr <- merge(dup_mf_snp_depth_chr, dup_mf_snp_ratio_CI.res, by = "chr")
dup_mf_snp_depth_chr <- merge(dup_mf_snp_depth_chr, dup_mf_depth_ratio_CI.res, by = "chr")

dup_mf_snp_depth_chr <- dup_mf_snp_depth_chr %>%  mutate(Z_vs_A = case_when(
  chr == 'HiC_scaffold_1' ~ 'Z1',
  chr == 'HiC_scaffold_2' ~ 'A',
  chr == 'HiC_scaffold_3' ~ 'A',
  chr == 'HiC_scaffold_4' ~ 'A',
  chr == 'HiC_scaffold_5' ~ 'A',
  chr == 'HiC_scaffold_6' ~ 'Z2',
  chr == 'HiC_scaffold_7' ~ 'A',
  chr == 'HiC_scaffold_8' ~ 'A',
  chr == 'HiC_scaffold_9' ~ 'A',
  chr == 'HiC_scaffold_10' ~ 'A',
  chr == 'HiC_scaffold_11' ~ 'A',
  chr == 'HiC_scaffold_12' ~ 'A',
  chr == 'HiC_scaffold_13' ~ 'A',
  chr == 'HiC_scaffold_14' ~ 'A',
  chr == 'HiC_scaffold_15' ~ 'A',
  chr == 'HiC_scaffold_16' ~ 'A',
  chr == 'HiC_scaffold_17' ~ 'A',
  chr == 'HiC_scaffold_18' ~ 'Z3',
  chr == 'HiC_scaffold_19' ~ 'A',
  chr == 'HiC_scaffold_20' ~ 'A',
  chr == 'HiC_scaffold_21' ~ 'A',
  chr == 'HiC_scaffold_22' ~ 'A',
  chr == 'HiC_scaffold_23' ~ 'A',
  chr == 'HiC_scaffold_24' ~ 'A',
  chr == 'HiC_scaffold_25' ~ 'A',
  chr == 'HiC_scaffold_26' ~ 'A',
  chr == 'HiC_scaffold_27' ~ 'A',
  chr == 'HiC_scaffold_28' ~ 'A',
  chr == 'HiC_scaffold_29' ~ 'A'))

# Z3 gets overplotted, split data in Z/A and plot Z on top of A to see Z3
dup_mf_snp_depth_chr_A <- dup_mf_snp_depth_chr[dup_mf_snp_depth_chr$Z_vs_A == "A", ]
dup_mf_snp_depth_chr_Z <- dup_mf_snp_depth_chr[dup_mf_snp_depth_chr$Z_vs_A != "A", ]

dup_plot <- ggplot(dup_mf_snp_depth_chr_A, aes(x=depth_ratio, y=snp_ratio, color=Z_vs_A)) +
  geom_point(size=2) +
  theme_minimal() +
  theme(plot.title = element_text(size=label_sizes),
        axis.text = element_text(size=label_sizes),
        axis.title = element_text(size=label_sizes),
        legend.title = element_blank(),
        legend.text = element_text(size=label_sizes)) +
  geom_linerange(aes(x=depth_ratio, ymin=snp_CI_low, ymax=snp_CI_high), linewidth=1) +
  geom_linerange(aes(y=snp_ratio, xmin=depth_CI_low, xmax=depth_CI_high), linewidth=1) +
  scale_color_manual(values=c(auto_color, Z1_col, Z2_col, Z3_col)) +
  scale_alpha_manual(guide='none', values = list(a = 0.2, point = 1)) +
  coord_cartesian(ylim = c(y_min, y_max), xlim = c(x_min, x_max)) +
  xlab("\nlog2 M/F read depth") +
  ylab("log2 M/F snp density\n") +
  ggtitle("L. duponcheli") +
  geom_point(data=dup_mf_snp_depth_chr_Z, aes(x=depth_ratio, y=snp_ratio), size=2) +
  geom_linerange(data=dup_mf_snp_depth_chr_Z, aes(x=depth_ratio, ymin=snp_CI_low, ymax=snp_CI_high), linewidth=1) +
  geom_linerange(data=dup_mf_snp_depth_chr_Z, aes(y=snp_ratio, xmin=depth_CI_low, xmax=depth_CI_high), linewidth=1) 



#---------------------------#
##### combine all plots #####
#---------------------------#

all_plots <- swe_plot / mor_plot / dup_plot /  plot_layout(guides = 'collect')
all_plots

ggsave(filename = ("supplementary_figure_1.png"), width = 7.5, height = 15, units = "in", dpi = 300, limitsize = F)

