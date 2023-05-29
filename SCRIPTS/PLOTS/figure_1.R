
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

setwd("~/Desktop/Plots/FASTZ")


#-----------------------------------------------------------------#
##### read in coverage data and calculate averages and ratios #####
#-----------------------------------------------------------------#


## sinapis ##

swe_m <- read.csv("LsinapisSweM-P14502_103-sorted.bam-unique.cov", sep = " ")

# per cds...
swe_m_cds_median = aggregate(depth ~ cds_midpoint + chr, data=swe_m, median)

# calculate scaling factor based on autosomes...
swe_m_cds_median_auto <- swe_m_cds_median[swe_m_cds_median$chr != "HiC_scaffold_1" & swe_m_cds_median$chr != "HiC_scaffold_6" & swe_m_cds_median$chr != "HiC_scaffold_18", ]
swe_m_cds_mean = mean(swe_m_cds_median_auto$depth)
swe_m_cds_median$norm_median <- swe_m_cds_median$depth/swe_m_cds_mean


swe_f <- read.csv("LsinapisSweM-P14502_104-sorted.bam-unique.cov", sep = " ")

# per cds...
swe_f_cds_median = aggregate(depth ~ cds_midpoint + chr, data=swe_f, median)

# calculate scaling factor based on autosomes...
swe_f_cds_median_auto <- swe_f_cds_median[swe_f_cds_median$chr != "HiC_scaffold_1" & swe_f_cds_median$chr != "HiC_scaffold_6" & swe_f_cds_median$chr != "HiC_scaffold_18", ]
swe_f_cds_mean = mean(swe_f_cds_median_auto$depth)
swe_f_cds_median$norm_median <- swe_f_cds_median$depth/swe_f_cds_mean

# join m+f...
swe_mf_cds_median <- merge(swe_m_cds_median, swe_f_cds_median, by = c("chr", "cds_midpoint"))

swe_mf_cds_median$mf_ratio <- log(swe_mf_cds_median$norm_median.x/swe_mf_cds_median$norm_median.y)
swe_mf_cds_median <- swe_mf_cds_median[is.finite(swe_mf_cds_median$mf_ratio), ]


## duponcheli ##

dup_m <- read.csv("LsinapisSweM-P14458_104_S25-sorted.bam-unique.cov", sep = " ")

# per cds...
dup_m_cds_median = aggregate(depth ~ cds_midpoint + chr, data=dup_m, median)

# calculate scaling factor based on autosomes...
dup_m_cds_median_auto <- dup_m_cds_median[dup_m_cds_median$chr != "HiC_scaffold_1" & dup_m_cds_median$chr != "HiC_scaffold_6" & dup_m_cds_median$chr != "HiC_scaffold_18", ]
dup_m_cds_mean = mean(dup_m_cds_median_auto$depth)
dup_m_cds_median$norm_median <- dup_m_cds_median$depth/dup_m_cds_mean


dup_f <- read.csv("LsinapisSweM-P14458_103_S24-sorted.bam-unique.cov", sep = " ")

# per cds...
dup_f_cds_median = aggregate(depth ~ cds_midpoint + chr, data=dup_f, median)

# calculate scaling factor based on autosomes...
dup_f_cds_median_auto <- dup_f_cds_median[dup_f_cds_median$chr != "HiC_scaffold_1" & dup_f_cds_median$chr != "HiC_scaffold_6" & dup_f_cds_median$chr != "HiC_scaffold_18", ]
dup_f_cds_mean = mean(dup_f_cds_median$depth)
dup_f_cds_median$norm_median <- dup_f_cds_median$depth/dup_f_cds_mean

# join m+f...
dup_mf_cds_median <- merge(dup_m_cds_median, dup_f_cds_median, by = c("chr", "cds_midpoint"))

# calculate per cds log m/f ratio...
dup_mf_cds_median$mf_ratio <- log(dup_mf_cds_median$norm_median.x/dup_mf_cds_median$norm_median.y)
dup_mf_cds_median <- dup_mf_cds_median[is.finite(dup_mf_cds_median$mf_ratio), ]


## morsei ##

mor_m <- read.csv("LsinapisSweM-P14458_108_S29-sorted.bam-unique.cov", sep = " ")

# per cds...
mor_m_cds_median = aggregate(depth ~ cds_midpoint + chr, data=mor_m, median)

# calculate scaling factor based on autosomes...
mor_m_cds_median_auto <- mor_m_cds_median[mor_m_cds_median$chr != "HiC_scaffold_1" & mor_m_cds_median$chr != "HiC_scaffold_6" & mor_m_cds_median$chr != "HiC_scaffold_18", ]
mor_m_cds_mean = mean(mor_m_cds_median$depth)
mor_m_cds_median$norm_median <- mor_m_cds_median$depth/mor_m_cds_mean


mor_f <- read.csv("LsinapisSweM-P14458_107_S28-sorted.bam-unique.cov", sep = " ")

# per cds...
mor_f_cds_median = aggregate(depth ~ cds_midpoint + chr, data=mor_f, median)

# calculate scaling factor based on autosomes...
mor_f_cds_median_auto <- mor_f_cds_median[mor_f_cds_median$chr != "HiC_scaffold_1" & mor_f_cds_median$chr != "HiC_scaffold_6" & mor_f_cds_median$chr != "HiC_scaffold_18", ]
mor_f_cds_mean = mean(mor_f_cds_median$depth)
mor_f_cds_median$norm_median <- mor_f_cds_median$depth/mor_f_cds_mean

# join m+f...
mor_mf_cds_median <- merge(mor_m_cds_median, mor_f_cds_median, by = c("chr", "cds_midpoint"))

# calculate per cds log m/f ratio...
mor_mf_cds_median$mf_ratio <- log(mor_mf_cds_median$norm_median.x/mor_mf_cds_median$norm_median.y)
mor_mf_cds_median <- mor_mf_cds_median[is.finite(mor_mf_cds_median$mf_ratio), ]






#------------------#
##### snp data #####
#------------------#


# snp plot variables...
heat_max_y=1.5
heat_min_y=0.5
trans_method=sqrt
female_color="#ea5545"
male_color="#27aeef"

# chr max coordinate...
Z1_max_coord=34455698
Z2_max_coord=26593402
Z3_max_coord=17773115


# male sinapis #

swe_m_all_snps <- read.csv("P14502_103-filtered_200_snps_per_100kb_windows.txt", sep = " ")
swe_m_auto_snps <- swe_m_all_snps[swe_m_all_snps$chr != "HiC_scaffold_1" & swe_m_all_snps$chr != "HiC_scaffold_6" & swe_m_all_snps$chr != "HiC_scaffold_18", ]
swe_m_Z_snps <- swe_m_all_snps[swe_m_all_snps$chr == "HiC_scaffold_1" | swe_m_all_snps$chr == "HiC_scaffold_6" | swe_m_all_snps$chr == "HiC_scaffold_18", ]
swe_m_Z_snps$mid <- (swe_m_Z_snps$start+swe_m_Z_snps$end)/2

# calculate normalization due to differens in snp count between samples...
swe_m_Z_snps$snp_count <- swe_m_Z_snps$snp_count/median(swe_m_all_snps$snp_count)

swe_m_Z1_snps <- swe_m_Z_snps[swe_m_Z_snps$chr == "HiC_scaffold_1", ]
swe_m_Z2_snps <- swe_m_Z_snps[swe_m_Z_snps$chr == "HiC_scaffold_6", ]
swe_m_Z3_snps <- swe_m_Z_snps[swe_m_Z_snps$chr == "HiC_scaffold_18", ]

# female sinapis #

swe_f_all_snps <- read.csv("P14502_104-filtered_200_snps_per_100kb_windows.txt", sep = " ")
swe_f_auto_snps <- swe_f_all_snps[swe_f_all_snps$chr != "HiC_scaffold_1" & swe_f_all_snps$chr != "HiC_scaffold_6" & swe_f_all_snps$chr != "HiC_scaffold_18", ]
swe_f_Z_snps <- swe_f_all_snps[swe_f_all_snps$chr == "HiC_scaffold_1" | swe_f_all_snps$chr == "HiC_scaffold_6" | swe_f_all_snps$chr == "HiC_scaffold_18", ]
swe_f_Z_snps$mid <- (swe_f_Z_snps$start+swe_f_Z_snps$end)/2

# calculate normalization due to differens in snp count between samples...
swe_f_Z_snps$snp_count <- swe_f_Z_snps$snp_count/median(swe_f_all_snps$snp_count)

swe_f_Z1_snps <- swe_f_Z_snps[swe_f_Z_snps$chr == "HiC_scaffold_1", ]
swe_f_Z2_snps <- swe_f_Z_snps[swe_f_Z_snps$chr == "HiC_scaffold_6", ]
swe_f_Z3_snps <- swe_f_Z_snps[swe_f_Z_snps$chr == "HiC_scaffold_18", ]


# scale heatmaps based on male Z max (values in female above max male Z will be plotted as saturation)...
swe_snp_plot_scaling <- max(swe_m_Z_snps$snp_count)


# make snp plots...

swe_m_Z1_snps_plot <- ggplot(swe_m_Z1_snps, aes(x=mid, y=1, fill=trans_method(snp_count/swe_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z1_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0)) 

swe_m_Z2_snps_plot <- ggplot(swe_m_Z2_snps, aes(x=mid, y=1, fill=trans_method(snp_count/swe_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z2_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0))

swe_m_Z3_snps_plot <- ggplot(swe_m_Z3_snps, aes(x=mid, y=1, fill=trans_method(snp_count/swe_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z3_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0)) 


swe_f_Z1_snps_plot <- ggplot(swe_f_Z1_snps, aes(x=mid, y=1, fill=trans_method(snp_count/swe_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z1_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))

swe_f_Z2_snps_plot <- ggplot(swe_f_Z2_snps, aes(x=mid, y=1, fill=trans_method(snp_count/swe_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z2_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))

swe_f_Z3_snps_plot <- ggplot(swe_f_Z3_snps, aes(x=mid, y=1, fill=trans_method(snp_count/swe_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z3_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))
  



# male morsei #

mor_m_all_snps <- read.csv("P14458_108_S29-filtered_200_snps_per_100kb_windows.txt", sep = " ")
mor_m_auto_snps <- mor_m_all_snps[mor_m_all_snps$chr != "HiC_scaffold_1" & mor_m_all_snps$chr != "HiC_scaffold_6" & mor_m_all_snps$chr != "HiC_scaffold_18", ]
mor_m_Z_snps <- mor_m_all_snps[mor_m_all_snps$chr == "HiC_scaffold_1" | mor_m_all_snps$chr == "HiC_scaffold_6" | mor_m_all_snps$chr == "HiC_scaffold_18", ]
mor_m_Z_snps$mid <- (mor_m_Z_snps$start+mor_m_Z_snps$end)/2

# calculate normalization due to differens in snp count between samples...
mor_m_Z_snps$snp_count <- mor_m_Z_snps$snp_count/median(mor_m_all_snps$snp_count)

mor_m_Z1_snps <- mor_m_Z_snps[mor_m_Z_snps$chr == "HiC_scaffold_1", ]
mor_m_Z2_snps <- mor_m_Z_snps[mor_m_Z_snps$chr == "HiC_scaffold_6", ]
mor_m_Z3_snps <- mor_m_Z_snps[mor_m_Z_snps$chr == "HiC_scaffold_18", ]


# female morsei #
 
mor_f_all_snps <- read.csv("P14458_107_S28-filtered_200_snps_per_100kb_windows.txt", sep = " ")
mor_f_auto_snps <- mor_f_all_snps[mor_f_all_snps$chr != "HiC_scaffold_1" & mor_f_all_snps$chr != "HiC_scaffold_6" & mor_f_all_snps$chr != "HiC_scaffold_18", ]
mor_f_Z_snps <- mor_f_all_snps[mor_f_all_snps$chr == "HiC_scaffold_1" | mor_f_all_snps$chr == "HiC_scaffold_6" | mor_f_all_snps$chr == "HiC_scaffold_18", ]
mor_f_Z_snps$mid <- (mor_f_Z_snps$start+mor_f_Z_snps$end)/2

# calculate normalization due to differens in snp count between samples...
mor_f_Z_snps$snp_count <- mor_f_Z_snps$snp_count/median(mor_f_all_snps$snp_count)

mor_f_Z1_snps <- mor_f_Z_snps[mor_f_Z_snps$chr == "HiC_scaffold_1", ]
mor_f_Z2_snps <- mor_f_Z_snps[mor_f_Z_snps$chr == "HiC_scaffold_6", ]
mor_f_Z3_snps <- mor_f_Z_snps[mor_f_Z_snps$chr == "HiC_scaffold_18", ]


# scale based on male Z max...
mor_snp_plot_scaling <- max(mor_m_Z_snps$snp_count)



# make snp plots...

mor_m_Z1_snps_plot <- ggplot(mor_m_Z1_snps, aes(x=mid, y=1, fill=sqrt(snp_count/mor_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z1_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0))

mor_m_Z2_snps_plot <- ggplot(mor_m_Z2_snps, aes(x=mid, y=1, fill=sqrt(snp_count/mor_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z2_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0))

mor_m_Z3_snps_plot <- ggplot(mor_m_Z3_snps, aes(x=mid, y=1, fill=sqrt(snp_count/mor_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z3_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0))


mor_f_Z1_snps_plot <- ggplot(mor_f_Z1_snps, aes(x=mid, y=1, fill=sqrt(snp_count/mor_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z1_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))

mor_f_Z2_snps_plot <- ggplot(mor_f_Z2_snps, aes(x=mid, y=1, fill=sqrt(snp_count/mor_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z2_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))

mor_f_Z3_snps_plot <- ggplot(mor_f_Z3_snps, aes(x=mid, y=1, fill=sqrt(snp_count/mor_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z3_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))



# male duponcheli #

dup_m_all_snps <- read.csv("P14458_104_S25-filtered_200_snps_per_100kb_windows.txt", sep = " ")
dup_m_auto_snps <- dup_m_all_snps[dup_m_all_snps$chr != "HiC_scaffold_1" & dup_m_all_snps$chr != "HiC_scaffold_6" & dup_m_all_snps$chr != "HiC_scaffold_18", ]
dup_m_Z_snps <- dup_m_all_snps[dup_m_all_snps$chr == "HiC_scaffold_1" | dup_m_all_snps$chr == "HiC_scaffold_6" | dup_m_all_snps$chr == "HiC_scaffold_18", ]
dup_m_Z_snps$mid <- (dup_m_Z_snps$start+dup_m_Z_snps$end)/2

# calculate normalization due to differens in snp count between samples...
dup_m_Z_snps$snp_count <- dup_m_Z_snps$snp_count/median(dup_m_all_snps$snp_count)

dup_m_Z1_snps <- dup_m_Z_snps[dup_m_Z_snps$chr == "HiC_scaffold_1", ]
dup_m_Z2_snps <- dup_m_Z_snps[dup_m_Z_snps$chr == "HiC_scaffold_6", ]
dup_m_Z3_snps <- dup_m_Z_snps[dup_m_Z_snps$chr == "HiC_scaffold_18", ]


# female duponcheli #

dup_f_all_snps <- read.csv("P14458_103_S24-filtered_200_snps_per_100kb_windows.txt", sep = " ")
dup_f_auto_snps <- dup_f_all_snps[dup_m_all_snps$chr != "HiC_scaffold_1" & dup_f_all_snps$chr != "HiC_scaffold_6" & dup_f_all_snps$chr != "HiC_scaffold_18", ]
dup_f_Z_snps <- dup_f_all_snps[dup_m_all_snps$chr == "HiC_scaffold_1" | dup_f_all_snps$chr == "HiC_scaffold_6" | dup_f_all_snps$chr == "HiC_scaffold_18", ]
dup_f_Z_snps$mid <- (dup_f_Z_snps$start+dup_f_Z_snps$end)/2

# calculate normalization due to differens in snp count between samples...
dup_f_Z_snps$snp_count <- dup_f_Z_snps$snp_count/median(dup_f_all_snps$snp_count)

dup_f_Z1_snps <- dup_f_Z_snps[dup_f_Z_snps$chr == "HiC_scaffold_1", ]
dup_f_Z2_snps <- dup_f_Z_snps[dup_f_Z_snps$chr == "HiC_scaffold_6", ]
dup_f_Z3_snps <- dup_f_Z_snps[dup_f_Z_snps$chr == "HiC_scaffold_18", ]


# scale based on male Z max...
dup_snp_plot_scaling <- max(dup_m_Z_snps$snp_count)


# make snp plots...

dup_m_Z1_snps_plot <- ggplot(dup_m_Z1_snps, aes(x=mid, y=1, fill=sqrt(snp_count/dup_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "red") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z1_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0))

dup_m_Z2_snps_plot <- ggplot(dup_m_Z2_snps, aes(x=mid, y=1, fill=sqrt(snp_count/dup_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "red") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z2_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0))

dup_m_Z3_snps_plot <- ggplot(dup_m_Z3_snps, aes(x=mid, y=1, fill=sqrt(snp_count/dup_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = male_color, limits=c(0, 1), na.value = "red") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z3_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = male_color)) +
  scale_x_continuous(expand = c(0.001, 0))


dup_f_Z1_snps_plot <- ggplot(dup_f_Z1_snps, aes(x=mid, y=1, fill=sqrt(snp_count/dup_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z1_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))


dup_f_Z2_snps_plot <- ggplot(dup_f_Z2_snps, aes(x=mid, y=1, fill=sqrt(snp_count/dup_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z2_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))


dup_f_Z3_snps_plot <- ggplot(dup_f_Z3_snps, aes(x=mid, y=1, fill=sqrt(snp_count/dup_snp_plot_scaling))) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = female_color, limits=c(0, 1), na.value = "black") +
  coord_cartesian(ylim = c(heat_min_y, heat_max_y), xlim = c(0, Z3_max_coord)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = female_color)) +
  scale_x_continuous(expand = c(0.001, 0))




#----------------------------#
##### sub plot variables #####
#----------------------------#


# synteny coordinates & colors

#Z1
anc_Z_start <- 10536855
anc_Z_end <- 32728189
neo_17_start <- 25979
neo_17_end <- 10453464

anz_Z_color <- rgb(213/255,117/255,0)
neo_17_color <- rgb(102/255,141/255,60/255)

#Z2
neo_11_start <- 75688
neo_11_end <- 10209824
neo_7_start <- 11321879
neo_7_end <- 20054014
neo_24_start <- 20280782
neo_24_end <- 26584900

neo_11_color <- rgb(185/255,156/255,107/255)
neo_7_color <- rgb(129/255,108/255,91/255)
neo_24_color <- rgb(228/255,153/255,105/255)

#Z3
neo_19_start <- 106210
neo_19_end <- 1196438
neo_15_start <- 1460622
neo_15_end <- 6655142
neo_8_start <- 6699355
neo_8_end <- 12082468
neo_26_start <- 12448676
neo_26_end <- 15033783
neo_28_start <- 15090709
neo_28_end <- 17625616

neo_19_color <- rgb(131/255,146/255,159/255)
neo_15_color <- rgb(219/255,202/255,105/255)
neo_8_color <- rgb(160/255,161/255,140/255)
neo_26_color <- rgb(157/255,151/255,84/255)
neo_28_color <- rgb(183/255,166/255,173/255)


# set smooth span to 1mb...
Z1_smooth=1000000/Z1_max_coord
Z2_smooth=1000000/Z2_max_coord
Z3_smooth=1000000/Z3_max_coord


# visuals...
chr_size=8
chr_pos=3
chr_text_size=4
chr_text_pos=3
title_size=18
y_lab_size=12
axis_size=12

depth_point_color="grey"
depth_line_color="black"


# make gender symbols...
male_symbol <- ggplot() +
  theme_void() +
  annotate("text", label="♂", y=1, x=1, size=3.5, color=male_color, fontface="bold", angle='-45')

female_symbol <- ggplot() +
  theme_void() +
  annotate("text", label="♀", y=1, x=1, size=3.5, color=female_color, fontface="bold") 



#------------------#
##### subplots #####
#------------------#


swe_mf_Z1 <- swe_mf_cds_median[swe_mf_cds_median$chr == "HiC_scaffold_1", ]

swe_mf_Z1_plot <- ggplot(swe_mf_Z1, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3)) +
  geom_smooth(span = Z1_smooth, colour=depth_line_color, method = "loess") +
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z1_max_coord, size=chr_size/2, color = "grey") +
  annotate("segment", y=chr_pos, yend=chr_pos, x=anc_Z_start, xend=anc_Z_end, size=chr_size, color = anz_Z_color) +
  annotate("text", label="Z (anc)", y=chr_text_pos, x=(anc_Z_start+anc_Z_end)/2, size=chr_text_size) +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_17_start, xend=neo_17_end, size=chr_size, color = neo_17_color) +
  annotate("text", label="17", y=chr_text_pos, x=(neo_17_start+neo_17_end)/2, size=chr_text_size) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size),
        plot.title = element_text(size=title_size, hjust = 0.5)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6), expand = c(0, 0),
                     minor_breaks = seq(2500000, 25000000, by=2500000)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ggtitle("Z1") +
  ylab("log m/f depth ratio") 

swe_mf_Z1_plot <- swe_mf_Z1_plot + 
  inset_element(swe_m_Z1_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(swe_f_Z1_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) + 
  inset_element(male_symbol, left = -0.04, bottom = 0.05, right = 0, top = 0.16) +
  inset_element(female_symbol, left = -0.04, bottom = 0, right = 0, top = 0.07) 


swe_mf_Z2 <- swe_mf_cds_median[swe_mf_cds_median$chr == "HiC_scaffold_6", ]

swe_mf_Z2_plot <- ggplot(swe_mf_Z2, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3)) +
  geom_smooth(span = Z2_smooth, colour=depth_line_color, method = "loess") + 
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z2_max_coord, size=chr_size/2, color = "grey") +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_11_start, xend=neo_11_end, size=chr_size, color = neo_11_color) +
  annotate("text", label="11", y=chr_text_pos, x=(neo_11_start+neo_11_end)/2, size=chr_text_size) +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_7_start, xend=neo_7_end, size=chr_size, color = neo_7_color) +
  annotate("text", label="7", y=chr_text_pos, x=(neo_7_start+neo_7_end)/2, size=chr_text_size) +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_24_start, xend=neo_24_end, size=chr_size, color = neo_24_color) +
  annotate("text", label="24", y=chr_text_pos, x=(neo_24_start+neo_24_end)/2, size=chr_text_size) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size),
        plot.title = element_text(size=title_size, hjust = 0.5)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ggtitle("Z2") +
  ylab("")

swe_mf_Z2_plot <- swe_mf_Z2_plot + 
  inset_element(swe_m_Z2_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(swe_f_Z2_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) 



swe_mf_Z3 <- swe_mf_cds_median[swe_mf_cds_median$chr == "HiC_scaffold_18", ]

swe_mf_Z3_plot <- ggplot(swe_mf_Z3, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3)) +
  geom_smooth(span = Z3_smooth, color=depth_line_color, method = "loess") +
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z3_max_coord, size=chr_size/2, color = "grey") +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_19_start, xend=neo_19_end, size=chr_size, color = neo_19_color) +
  annotate("text", label="19", y=chr_text_pos, x=(neo_19_start+neo_19_end)/2, size=chr_text_size) +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_15_start, xend=neo_15_end, size=chr_size, color = neo_15_color) +
  annotate("text", label="15", y=chr_text_pos, x=(neo_15_start+neo_15_end)/2, size=chr_text_size) +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_8_start, xend=neo_8_end, size=chr_size, color = neo_8_color) +
  annotate("text", label="8", y=chr_text_pos, x=(neo_8_start+neo_8_end)/2, size=chr_text_size) +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_26_start, xend=neo_26_end, size=chr_size, color = neo_26_color) +
  annotate("text", label="26", y=chr_text_pos, x=(neo_26_start+neo_26_end)/2, size=chr_text_size) +
  annotate("segment", y=chr_pos, yend=chr_pos, x=neo_28_start, xend=neo_28_end, size=chr_size, color = neo_28_color) +
  annotate("text", label="28", y=chr_text_pos, x=(neo_28_start+neo_28_end)/2, size=chr_text_size) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size),
        plot.title = element_text(size=title_size, hjust = 0.5)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ggtitle("Z3") +
  ylab("") 

swe_mf_Z3_plot <- swe_mf_Z3_plot + 
  inset_element(swe_m_Z3_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(swe_f_Z3_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) 


mor_mf_Z1 <- mor_mf_cds_median[mor_mf_cds_median$chr == "HiC_scaffold_1", ]

mor_mf_Z1_plot <- ggplot(mor_mf_Z1, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3), expand = FALSE) +
  geom_smooth(span = Z1_smooth, color=depth_line_color, method = "loess") + 
  theme_minimal() +
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z1_max_coord, size=chr_size/2, color = "white", alpha=0) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6), expand = c(0, 0),
                     minor_breaks = seq(2500000, 25000000, by=2500000)) +  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("log m/f depth ratio") 

mor_mf_Z1_plot <- mor_mf_Z1_plot + 
  inset_element(mor_m_Z1_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(mor_f_Z1_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) +
  inset_element(male_symbol, left = -0.04, bottom = 0.05, right = 0, top = 0.16) +
  inset_element(female_symbol, left = -0.04, bottom = 0, right = 0, top = 0.07) 


mor_mf_Z2 <- mor_mf_cds_median[mor_mf_cds_median$chr == "HiC_scaffold_6", ]

mor_mf_Z2_plot <- ggplot(mor_mf_Z2, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3), expand = FALSE) +
  geom_smooth(span = Z2_smooth, color=depth_line_color, method = "loess") + 
  theme_minimal() +
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z2_max_coord, size=chr_size/2, color = "white", alpha=0) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("")

mor_mf_Z2_plot <- mor_mf_Z2_plot + 
  inset_element(mor_m_Z2_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(mor_f_Z2_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) 


mor_mf_Z3 <- mor_mf_cds_median[mor_mf_cds_median$chr == "HiC_scaffold_18", ]

mor_mf_Z3_plot <- ggplot(mor_mf_Z3, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3), expand = FALSE) +
  geom_smooth(span = Z3_smooth, color=depth_line_color, method = "loess") + 
  theme_minimal() +
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z3_max_coord, size=chr_size/2, color = "white", alpha=0) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("") 

mor_mf_Z3_plot <- mor_mf_Z3_plot + 
  inset_element(mor_m_Z3_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(mor_f_Z3_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) 



dup_mf_Z1 <- dup_mf_cds_median[dup_mf_cds_median$chr == "HiC_scaffold_1", ]

dup_mf_Z1_plot <- ggplot(dup_mf_Z1, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3), expand = FALSE) +
  geom_smooth(span = Z1_smooth, color=depth_line_color, method = "loess") + 
  theme_minimal() +
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z1_max_coord, size=chr_size/2, color = "white", alpha=0) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6), expand = c(0, 0),
                     minor_breaks = seq(2500000, 25000000, by=2500000)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("log m/f depth ratio") 

dup_mf_Z1_plot <- dup_mf_Z1_plot + 
  inset_element(dup_m_Z1_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(dup_f_Z1_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) +
  inset_element(male_symbol, left = -0.04, bottom = 0.05, right = 0, top = 0.16) +
  inset_element(female_symbol, left = -0.04, bottom = 0, right = 0, top = 0.07) 


dup_mf_Z2 <- dup_mf_cds_median[dup_mf_cds_median$chr == "HiC_scaffold_6", ]

dup_mf_Z2_plot <- ggplot(dup_mf_Z2, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3), expand = FALSE) +
  geom_smooth(span = Z2_smooth, color=depth_line_color, method = "loess") + 
  theme_minimal() +
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z2_max_coord, size=chr_size/2, color = "white", alpha=0) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("")

dup_mf_Z2_plot <- dup_mf_Z2_plot + 
  inset_element(dup_m_Z2_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(dup_f_Z2_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) 


dup_mf_Z3 <- dup_mf_cds_median[dup_mf_cds_median$chr == "HiC_scaffold_18", ]

dup_mf_Z3_plot <- ggplot(dup_mf_Z3, aes(x=cds_midpoint, y=mf_ratio)) +
  geom_point(aes(color = mf_ratio)) +
  scale_color_gradientn(colors=c(female_color, depth_point_color, male_color), limits=c(-3, 3)) +
  coord_cartesian(ylim = c(-3, 3), expand = FALSE) +
  geom_smooth(span = Z3_smooth, color=depth_line_color, method = "loess") + 
  theme_minimal() +
  annotate("segment", y=chr_pos, yend=chr_pos, x=0, xend=Z3_max_coord, size=chr_size/2, color = "white", alpha=0) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("") 

dup_mf_Z3_plot <- dup_mf_Z3_plot + 
  inset_element(dup_m_Z3_snps_plot, left = 0, bottom = 0.07, right = 1, top = 0.14) +
  inset_element(dup_f_Z3_snps_plot, left = 0, bottom = 0, right = 1, top = 0.07) 



#-------------------#
##### cladogram #####
#-------------------#

branch_size=2
node_size=1
tree_color="black"
species_label_size=5

tree_plot <- ggplot() +
  theme_void() +
  annotate("segment", y=2.5, yend=2.5, x=1, xend=2, size=branch_size, color=tree_color) +
  annotate("segment", y=1, yend=1, x=2, xend=4, size=branch_size, color=tree_color) +
  annotate("segment", y=4, yend=4, x=2, xend=3, size=branch_size, color=tree_color) +
  annotate("segment", y=3, yend=3, x=3, xend=4, size=branch_size, color=tree_color) +
  annotate("segment", y=5, yend=5, x=3, xend=4, size=branch_size, color=tree_color) +
  annotate("segment", y=1, yend=4, x=2, xend=2, size=branch_size, color=tree_color) +
  annotate("segment", y=3, yend=5, x=3, xend=3, size=branch_size, color=tree_color) +
  annotate("point", y=2.5, x=1, color=tree_color, size=node_size) +
  annotate("point", y=1, x=2, color=tree_color, size=node_size) +
  annotate("point", y=4, x=2, color=tree_color, size=node_size) +
  annotate("point", y=3, x=3, color=tree_color, size=node_size) +
  annotate("point", y=5, x=3, color=tree_color, size=node_size) +
  annotate("text", x=5.2, y=5, label="L. sinapis", fontface= "italic", hjust=0.5, size=species_label_size) +
  annotate("text", x=5.2, y=3, label="L. morsei", fontface= "italic", hjust=0.5, size=species_label_size) +
  annotate("text", x=5.2, y=1, label="L. duponcheli", fontface= "italic", hjust=0.5, size=species_label_size) +
  annotate("text", x=1, y=0.4, label="") +
  annotate("text", x=1, y=5.6, label="") +
  annotate("text", x=6.2, y=5, label="") +
  
  annotate("text", label="Z2+ ", y=3.25, x=1.5, size=chr_text_size) +
  annotate("segment", y=3.25, yend=3.25, x=1.75, xend=2.25, size=chr_size, color=neo_11_color) +
  annotate("text", label="11", y=3.25, x=2, size=chr_text_size) +
  
  annotate("text", label="Z3: ", y=4.5, x=1.5, size=chr_text_size) +
  annotate("segment", y=4.5, yend=4.5, x=1.75, xend=2.25, size=chr_size, color=neo_19_color) +
  annotate("text", label="19", y=4.5, x=2, size=chr_text_size) +
  annotate("segment", y=4.5, yend=4.5, x=2.25, xend=2.75, size=chr_size, color=neo_15_color) +
  annotate("text", label="15", y=4.5, x=2.5, size=chr_text_size) +
  annotate("segment", y=4.5, yend=4.5, x=2.75, xend=3.25, size=chr_size, color=neo_8_color) +
  annotate("text", label="8", y=4.5, x=3, size=chr_text_size) +
  annotate("segment", y=4.5, yend=4.5, x=3.25, xend=3.75, size=chr_size, color=neo_26_color) +
  annotate("text", label="26", y=4.5, x=3.5, size=chr_text_size) +
  annotate("segment", y=4.5, yend=4.5, x=3.75, xend=4.25, size=chr_size, color=neo_28_color) +
  annotate("text", label="28", y=4.5, x=4, size=chr_text_size) 
  




#---------------------------#
##### combine all plots #####
#---------------------------#


# scale chromosomes by size for plotting...

scale_Z1=34455698
scale_Z2=26593402
scale_Z3=17773115

total_size <- (scale_Z1+scale_Z2+scale_Z3)

Z_plots <- tree_plot +
  (swe_mf_Z1_plot + swe_mf_Z2_plot + swe_mf_Z3_plot + plot_layout(widths = c((scale_Z1/total_size*12), (scale_Z2/total_size*12), (scale_Z3/total_size*12)))) /
  (mor_mf_Z1_plot + mor_mf_Z2_plot + mor_mf_Z3_plot + plot_layout(widths = c((scale_Z1/total_size*12), (scale_Z2/total_size*12), (scale_Z3/total_size*12)))) /
  (dup_mf_Z1_plot + dup_mf_Z2_plot + dup_mf_Z3_plot + plot_layout(widths = c((scale_Z1/total_size*12), (scale_Z2/total_size*12), (scale_Z3/total_size*12)))) +
  plot_layout(widths = c(1, 4))

Z_plots


ggsave(filename = ("figure_1.png"), width = 15, height = 7.5, units = "in", dpi = 300, limitsize = F)

