
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


# sinapis

swe_m <- read.csv("LsinapisSweM-P14502_103-sorted.bam-unique.cov", sep = " ")

swe_m$gene_midpoint <- (swe_m$start+swe_m$stop)/2

# per gene
swe_m_gene_median = aggregate(depth ~ gene_id + chr + gene_midpoint, data=swe_m, median)

# calculate scaling factor based on autosomes...
swe_m_gene_median_auto <- swe_m_gene_median[swe_m_gene_median$chr != "HiC_scaffold_1" & swe_m_gene_median$chr != "HiC_scaffold_6" & swe_m_gene_median$chr != "HiC_scaffold_18", ]
swe_m_gene_mean = mean(swe_m_gene_median_auto$depth)
swe_m_gene_median$norm_median <- swe_m_gene_median$depth/swe_m_gene_mean


swe_f <- read.csv("LsinapisSweM-P14502_104-sorted.bam-unique.cov", sep = " ")

swe_f$gene_midpoint <- (swe_f$start+swe_f$stop)/2

# per gene
swe_f_gene_median = aggregate(depth ~ gene_id + chr + gene_midpoint, data=swe_f, median)

# calculate scaling factor based on autosomes...
swe_f_gene_median_auto <- swe_f_gene_median[swe_f_gene_median$chr != "HiC_scaffold_1" & swe_f_gene_median$chr != "HiC_scaffold_6" & swe_f_gene_median$chr != "HiC_scaffold_18", ]
swe_f_gene_mean = mean(swe_f_gene_median_auto$depth)
swe_f_gene_median$norm_median <- swe_f_gene_median$depth/swe_f_gene_mean

# join
swe_mf_gene_median <- merge(swe_m_gene_median, swe_f_gene_median, by = c("chr", "gene_midpoint"))

# calculate per gene log m/f ratio
swe_mf_gene_median$mf_ratio <- log(swe_mf_gene_median$norm_median.x/swe_mf_gene_median$norm_median.y)
swe_mf_gene_median <- swe_mf_gene_median[is.finite(swe_mf_gene_median$mf_ratio), ]


# duponcheli

dup_m <- read.csv("LsinapisSweM-P14458_104_S25-sorted.bam-unique.cov", sep = " ")

dup_m$gene_midpoint <- (dup_m$start+dup_m$stop)/2

# per gene
dup_m_gene_median = aggregate(depth ~ gene_id + chr + gene_midpoint, data=dup_m, median)

# calculate scaling factor based on autosomes...
dup_m_gene_median_auto <- dup_m_gene_median[dup_m_gene_median$chr != "HiC_scaffold_1" & dup_m_gene_median$chr != "HiC_scaffold_6" & dup_m_gene_median$chr != "HiC_scaffold_18", ]
dup_m_gene_mean = mean(dup_m_gene_median_auto$depth)
dup_m_gene_median$norm_median <- dup_m_gene_median$depth/dup_m_gene_mean

dup_f <- read.csv("LsinapisSweM-P14458_103_S24-sorted.bam-unique.cov", sep = " ")

dup_f$gene_midpoint <- (dup_f$start+dup_f$stop)/2

# per gene...
dup_f_gene_median = aggregate(depth ~ gene_id + chr + gene_midpoint, data=dup_f, median)

# calculate scaling factor based on autosomes...
dup_f_gene_median_auto <- dup_f_gene_median[dup_f_gene_median$chr != "HiC_scaffold_1" & dup_f_gene_median$chr != "HiC_scaffold_6" & dup_f_gene_median$chr != "HiC_scaffold_18", ]
dup_f_gene_mean = mean(dup_f_gene_median$depth)
dup_f_gene_median$norm_median <- dup_f_gene_median$depth/dup_f_gene_mean

# join
dup_mf_gene_median <- merge(dup_m_gene_median, dup_f_gene_median, by = c("chr", "gene_midpoint"))

# calculate per gene log m/f ratio
dup_mf_gene_median$mf_ratio <- log(dup_mf_gene_median$norm_median.x/dup_mf_gene_median$norm_median.y)
dup_mf_gene_median <- dup_mf_gene_median[is.finite(dup_mf_gene_median$mf_ratio), ]


# morsei

mor_m <- read.csv("LsinapisSweM-P14458_108_S29-sorted.bam-unique.cov", sep = " ")

mor_m$gene_midpoint <- (mor_m$start+mor_m$stop)/2

# per gene...
mor_m_gene_median = aggregate(depth ~ gene_id + chr + gene_midpoint, data=mor_m, median)

# calculate scaling factor based on autosomes...
mor_m_gene_median_auto <- mor_m_gene_median[mor_m_gene_median$chr != "HiC_scaffold_1" & mor_m_gene_median$chr != "HiC_scaffold_6" & mor_m_gene_median$chr != "HiC_scaffold_18", ]
mor_m_gene_mean = mean(mor_m_gene_median$depth)
mor_m_gene_median$norm_median <- mor_m_gene_median$depth/mor_m_gene_mean


mor_f <- read.csv("LsinapisSweM-P14458_107_S28-sorted.bam-unique.cov", sep = " ")

mor_f$gene_midpoint <- (mor_f$start+mor_f$stop)/2

# per gene...
mor_f_gene_median = aggregate(depth ~ gene_id + chr + gene_midpoint, data=mor_f, median)

# calculate scaling factor based on autosomes...
mor_f_gene_median_auto <- mor_f_gene_median[mor_f_gene_median$chr != "HiC_scaffold_1" & mor_f_gene_median$chr != "HiC_scaffold_6" & mor_f_gene_median$chr != "HiC_scaffold_18", ]
mor_f_gene_mean = mean(mor_f_gene_median$depth)
mor_f_gene_median$norm_median <- mor_f_gene_median$depth/mor_f_gene_mean

# join
mor_mf_gene_median <- merge(mor_m_gene_median, mor_f_gene_median, by = c("chr", "gene_midpoint"))

# calculate per gene log m/f ratio
mor_mf_gene_median$mf_ratio <- log(mor_mf_gene_median$norm_median.x/mor_mf_gene_median$norm_median.y)
mor_mf_gene_median <- mor_mf_gene_median[is.finite(mor_mf_gene_median$mf_ratio), ]



# plot variables

chr_size=8
chr_pos=3
chr_text_size=4
chr_text_pos=3
title_size=18
y_lab_size=12
axis_size=12

depth_point_color="grey"
depth_line_color="black"


#--------------------#
##### make plots #####
#--------------------#


## sin ##

plot_list <- list()   

for (i in 1:29) {
  
  subset_df <- swe_mf_gene_median[swe_mf_gene_median$chr == paste0("HiC_scaffold_", i), ]
  
  header=paste("Chr", i)
  
  chr_plot <- ggplot(subset_df, aes(x=gene_midpoint, y=mf_ratio)) +
    geom_point(color=depth_point_color) +
    coord_cartesian(ylim = c(-2, 2)) +
    geom_smooth(span = 0.1, method = "loess", color=depth_line_color) + 
    theme_minimal() +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_text(size=y_lab_size),
          axis.text.x = element_text(size=axis_size),
          axis.text.y = element_text(size=axis_size)) +
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(-1, 0, 1)) +
    xlab(header)
  
  plot_list[[i]] <- chr_plot
  
}

# remove Z...
plot_list <- plot_list[-c(1, 6, 18)]

label_plot <- ggplot() +
  annotate("text", x=1, y=1, label="A", size=8) +
  theme_void()

swe_plots <- (label_plot | plot_spacer() | plot_spacer() | plot_spacer() | 
                plot_spacer() | plot_spacer() |plot_spacer() | plot_spacer() | plot_spacer() | plot_spacer()) / 
  patchwork::wrap_plots(plot_list) + plot_layout(heights = c(1, 10))

swe_plots

ggsave(filename = ("supplementary_figure_1A.png"), width = 15, height = 15, units = "in", dpi = 300, limitsize = F)


## mor ##

plot_list <- list()   

for (i in 1:29) {

subset_df <- mor_mf_gene_median[mor_mf_gene_median$chr == paste0("HiC_scaffold_", i), ]

header=paste("Chr", i)

chr_plot <- ggplot(subset_df, aes(x=gene_midpoint, y=mf_ratio)) +
  geom_point(color=depth_point_color) +
  coord_cartesian(ylim = c(-2, 2)) +
  geom_smooth(span = 0.1, method = "loess", color=depth_line_color) + 
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=y_lab_size),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size)) +
  scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  xlab(header)

plot_list[[i]] <- chr_plot

}

# remove Z...
plot_list <- plot_list[-c(1, 6, 18)]

label_plot <- ggplot() +
  annotate("text", x=1, y=1, label="B", size=8) +
  theme_void()

mor_plots <- (label_plot | plot_spacer() | plot_spacer() | plot_spacer() | 
                plot_spacer() | plot_spacer() |plot_spacer() | plot_spacer() | plot_spacer() | plot_spacer()) / 
  patchwork::wrap_plots(plot_list) + plot_layout(heights = c(1, 10))

mor_plots

ggsave(filename = ("supplementary_figure_1B.png"), width = 15, height = 15, units = "in", dpi = 300, limitsize = F)




## dup ##

plot_list <- list()   

for (i in 1:29) {
  
  subset_df <- dup_mf_gene_median[dup_mf_gene_median$chr == paste0("HiC_scaffold_", i), ]
  
  header=paste("Chr", i)
  
  chr_plot <- ggplot(subset_df, aes(x=gene_midpoint, y=mf_ratio)) +
    geom_point(color=depth_point_color) +
    coord_cartesian(ylim = c(-2, 2)) +
    geom_smooth(span = 0.1, method = "loess", color=depth_line_color) + 
    theme_minimal() +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_text(size=y_lab_size),
          axis.text.x = element_text(size=axis_size),
          axis.text.y = element_text(size=axis_size)) +
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(-1, 0, 1)) +
    xlab(header)
  
  plot_list[[i]] <- chr_plot
  
}

# remove Z...
plot_list <- plot_list[-c(1, 6, 18)]

label_plot <- ggplot() +
  annotate("text", x=1, y=1, label="C", size=8) +
  theme_void()

dup_plots <- (label_plot | plot_spacer() | plot_spacer() | plot_spacer() | 
                plot_spacer() | plot_spacer() |plot_spacer() | plot_spacer() | plot_spacer() | plot_spacer()) / 
  patchwork::wrap_plots(plot_list) + plot_layout(heights = c(1, 10))

dup_plots

ggsave(filename = ("supplementary_figure_1C.png"), width = 15, height = 15, units = "in", dpi = 300, limitsize = F)
