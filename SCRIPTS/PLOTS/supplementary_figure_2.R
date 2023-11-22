
library(ggplot2)
library(rstatix)
library(scales)


# analyse Z3 with putative W gametologs

# check gene wise read depth difference between male and female L. sinapis
# equal read depth suggests presence of W gametolog in females

# check snp counts per gene
# higher snp count in female suggests W copy aligned to Z3 gene

setwd("~/Desktop/Plots/FASTZ")

df <- read.csv("lep_sin_dnds_data.txt", sep = "\t")

df$DoS <- (df$Dn/(df$Dn+df$Ds)) - (df$Pn/(df$Pn+df$Ps))

df_filtered <- df[df$dnds <999, ]

Z3 <- df[df$chrN == "Z3", ]


# use normalization from whole genome data to account for sequencing depth and mapping bias... 

swe_m <- read.csv("LsinapisSweM-P14502_103-sorted.bam-unique.cov", sep = " ")
swe_m$gene_midpoint <- (swe_m$start+swe_m$stop)/2

# per gene...
swe_m_gene_median = aggregate(depth ~ gene_midpoint + chr + gene_id, data=swe_m, median)

# calculate scaling factor based on autosomes...
swe_m_gene_median_auto <- swe_m_gene_median[swe_m_gene_median$chr != "HiC_scaffold_1" & swe_m_gene_median$chr != "HiC_scaffold_6" & swe_m_gene_median$chr != "HiC_scaffold_18", ]
swe_m_gene_mean = mean(swe_m_gene_median_auto$depth)
swe_m_gene_median$norm_depth <- swe_m_gene_median$depth/swe_m_gene_mean

swe_m_gene_median_Z3 <- swe_m_gene_median[swe_m_gene_median$chr == "HiC_scaffold_18", ]
colnames(swe_m_gene_median_Z3) <- c("gene_midpoint", "chr", "gene_id", "depth", "swe_m_norm_depth")



swe_f <- read.csv("LsinapisSweM-P14502_104-sorted.bam-unique.cov", sep = " ")
swe_f$gene_midpoint <- (swe_f$start+swe_f$stop)/2

# per gene...
swe_f_gene_median = aggregate(depth ~ gene_midpoint + chr + gene_id, data=swe_f, median)

# calculate scaling factor based on autosomes...
swe_f_gene_median_auto <- swe_f_gene_median[swe_f_gene_median$chr != "HiC_scaffold_1" & swe_f_gene_median$chr != "HiC_scaffold_6" & swe_f_gene_median$chr != "HiC_scaffold_18", ]
swe_f_gene_mean = mean(swe_f_gene_median_auto$depth)
swe_f_gene_median$norm_depth <- swe_f_gene_median$depth/swe_f_gene_mean

swe_f_gene_median_Z3 <- swe_f_gene_median[swe_f_gene_median$chr == "HiC_scaffold_18", ]
colnames(swe_f_gene_median_Z3) <- c("gene_midpoint", "chr", "gene_id", "depth", "swe_f_norm_depth")


# join m+f...
swe_mf_gene_median_Z3 <- merge(swe_m_gene_median_Z3, swe_f_gene_median_Z3, by = c("chr", "gene_id", "gene_midpoint"))




# spanish samples...

spa_m <- read.csv("LsinapisSweM-P14502_105-sorted.bam-unique.cov", sep = " ")
spa_m$gene_midpoint <- (spa_m$start+spa_m$stop)/2

# per gene...
spa_m_gene_median = aggregate(depth ~ gene_midpoint + chr + gene_id, data=spa_m, median)

# calculate scaling factor based on autosomes...
spa_m_gene_median_auto <- spa_m_gene_median[spa_m_gene_median$chr != "HiC_scaffold_1" & spa_m_gene_median$chr != "HiC_scaffold_6" & spa_m_gene_median$chr != "HiC_scaffold_18", ]
spa_m_gene_mean = mean(spa_m_gene_median_auto$depth)
spa_m_gene_median$norm_depth <- spa_m_gene_median$depth/spa_m_gene_mean

spa_m_gene_median_Z3 <- spa_m_gene_median[spa_m_gene_median$chr == "HiC_scaffold_18", ]
colnames(spa_m_gene_median_Z3) <- c("gene_midpoint", "chr", "gene_id", "depth", "spa_m_norm_depth")


spa_f <- read.csv("LsinapisSweM-P14502_106-sorted.bam-unique.cov", sep = " ")
spa_f$gene_midpoint <- (spa_f$start+spa_f$stop)/2

# per gene...
spa_f_gene_median = aggregate(depth ~ gene_midpoint + chr + gene_id, data=spa_f, median)

# calculate scaling factor based on autosomes...
spa_f_gene_median_auto <- spa_f_gene_median[spa_f_gene_median$chr != "HiC_scaffold_1" & spa_f_gene_median$chr != "HiC_scaffold_6" & spa_f_gene_median$chr != "HiC_scaffold_18", ]
spa_f_gene_mean = mean(spa_f_gene_median_auto$depth)
spa_f_gene_median$norm_depth <- spa_f_gene_median$depth/spa_f_gene_mean

spa_f_gene_median_Z3 <- spa_f_gene_median[spa_f_gene_median$chr == "HiC_scaffold_18", ]
colnames(spa_f_gene_median_Z3) <- c("gene_midpoint", "chr", "gene_id", "depth", "spa_f_norm_depth")

# join m+f...
spa_mf_gene_median_Z3 <- merge(spa_m_gene_median_Z3, spa_f_gene_median_Z3, by = c("chr", "gene_id", "gene_midpoint"))
spa_mf_gene_median <- merge(spa_m_gene_median, spa_f_gene_median, by = c("chr", "gene_id", "gene_midpoint"))

#spa_mf_gene_median$mf_ratio <- log2(spa_mf_gene_median$norm_median.x/spa_mf_gene_median$norm_median.y)

# join swe+spa...
all_mf_gene_median_Z3 <- merge(swe_mf_gene_median_Z3, spa_mf_gene_median_Z3, by = c("chr", "gene_id", "gene_midpoint"))

# remove genes without depth >10 in all samples...

all_mf_gene_median_Z3 <- all_mf_gene_median_Z3[all_mf_gene_median_Z3$depth.x.x > 10 |
                                               all_mf_gene_median_Z3$depth.y.x > 10 |
                                               all_mf_gene_median_Z3$depth.x.y > 10 |
                                               all_mf_gene_median_Z3$depth.y.y > 10, ]



Z_threshold <- log2(1/0.80)
ZW_threshold <- log2(1/0.80)

# for Z unique genes, require more than 80% in spanish male (swe/spa ratio < ZW threshold)
# and less than 80% in females (m/f ratio > Z threshold)

Z_genes <- all_mf_gene_median_Z3[log2((all_mf_gene_median_Z3$swe_m_norm_depth)/(all_mf_gene_median_Z3$spa_m_norm_depth)) < ZW_threshold &
                                   log2((all_mf_gene_median_Z3$swe_m_norm_depth)/(all_mf_gene_median_Z3$swe_f_norm_depth)) > Z_threshold &
                                   log2((all_mf_gene_median_Z3$swe_m_norm_depth)/(all_mf_gene_median_Z3$spa_f_norm_depth)) > Z_threshold, ]

Z_genes <- na.omit(Z_genes)

# for ZWgenes, require more than 80% in all compared to swe male (swe/n ratio < ZW threshold)

ZW_genes <- all_mf_gene_median_Z3[log2((all_mf_gene_median_Z3$swe_m_norm_depth)/(all_mf_gene_median_Z3$spa_m_norm_depth)) < ZW_threshold &
                                    log2((all_mf_gene_median_Z3$swe_m_norm_depth)/(all_mf_gene_median_Z3$swe_f_norm_depth)) < Z_threshold &
                                    log2((all_mf_gene_median_Z3$swe_m_norm_depth)/(all_mf_gene_median_Z3$spa_f_norm_depth)) < Z_threshold, ]

ZW_genes <- na.omit(ZW_genes)

# list all Z unique and ZW genes based on depth...

ZW_depth_list <- list(ZW_genes$gene_id)
Z_depth_list <- list(Z_genes$gene_id)

Z_depth_genes_df <- Z3[Z3$gene_id %in% Z_depth_list[[1]], ]
ZW_depth_genes_df <- Z3[Z3$gene_id %in% ZW_depth_list[[1]], ]

nrow(ZW_depth_genes_df)
nrow(Z_depth_genes_df)

# Counts of ZW and Z unique genes based on read depth:
# ZW: 251 
# Z:   24


#Update lists after filtering...

ZW_depth_list <- list(ZW_depth_genes_df$gene_id)
Z_depth_list <- list(Z_depth_genes_df$gene_id)


##### check gene wise snp density #####

# read in snp data...

swe_male_snps <- read.csv("P14502_103-filtered_200_snps_per_gene.txt", sep="\t")
swe_female_snps <- read.csv("P14502_104-filtered_200_snps_per_gene.txt", sep="\t")

spa_male_snps <- read.csv("P14502_105-filtered_200_snps_per_gene.txt", sep="\t")
spa_female_snps <- read.csv("P14502_106-filtered_200_snps_per_gene.txt", sep="\t")

# add sex variable and sum counts per gene...
#swe_male_snps$sex <- "swe_male"
swe_male_snps_agg <- aggregate(count ~ gene_id, swe_male_snps[swe_male_snps$chr == "HiC_scaffold_18", ], FUN=sum)
colnames(swe_male_snps_agg) <- c("gene_id", "swe_m_snps")

#swe_female_snps$sex <- "swe_female"
swe_female_snps_agg <- aggregate(count ~ gene_id, swe_female_snps[swe_female_snps$chr == "HiC_scaffold_18", ], FUN=sum)
colnames(swe_female_snps_agg) <- c("gene_id", "swe_f_snps")

#spa_male_snps$sex <- "spa_male"
spa_male_snps_agg <- aggregate(count ~ gene_id, spa_male_snps[spa_male_snps$chr == "HiC_scaffold_18", ], FUN=sum)
colnames(spa_male_snps_agg) <- c("gene_id", "spa_m_snps")

#spa_female_snps$sex <- "spa_female"
spa_female_snps_agg <- aggregate(count ~ gene_id, spa_female_snps[spa_female_snps$chr == "HiC_scaffold_18", ], FUN=sum)
colnames(spa_female_snps_agg) <- c("gene_id", "spa_f_snps")

swe_snps_agg <- merge(swe_male_snps_agg, swe_female_snps_agg, by="gene_id")
spa_snps_agg <- merge(spa_male_snps_agg, spa_female_snps_agg, by="gene_id")

all_snps_agg <- merge(swe_snps_agg, spa_snps_agg, by="gene_id")

# keep if snp in both swe and spa female...
ZW_snps_agg <- all_snps_agg[all_snps_agg$swe_f_snps > 0 & all_snps_agg$spa_f_snps > 0, ]
Z_snps_agg <- all_snps_agg[all_snps_agg$swe_f_snps == 0 & all_snps_agg$spa_f_snps == 0, ]

# combine lists of ZW candidate genes...
ZW_snp_genes_list <- list(ZW_snps_agg$gene_id)
Z_snp_genes_list <- list(Z_snps_agg$gene_id)

ZW_snp_genes_df <- Z3[Z3$gene_id %in% ZW_snp_genes_list[[1]], ]
Z_snp_genes_df <- Z3[Z3$gene_id %in% Z_snp_genes_list[[1]], ]

nrow(ZW_snp_genes_df)
nrow(Z_snp_genes_df)

# Counts of ZW and Z unique genes based on SNP density:

# ZW: 279
# Z:   12



# all genes that have both female snps and female depth
all_ZW_intersect_list <- Reduce(intersect, list(ZW_snp_genes_list[[1]], ZW_depth_list[[1]]))

length(all_ZW_intersect_list)

# ZW: 244

# all genes that don't have both female snps and female depth
all_Z_intersect_list <- Reduce(intersect, list(Z_snp_genes_list[[1]], Z_depth_list[[1]]))

length(all_Z_intersect_list)

# Z: 7

#--------------#
##### Plot #####
#--------------#


swe_m <- swe_m_gene_median[swe_m_gene_median$chr == "HiC_scaffold_18", ]
swe_f <- swe_f_gene_median[swe_f_gene_median$chr == "HiC_scaffold_18", ]
spa_m <- spa_m_gene_median[spa_m_gene_median$chr == "HiC_scaffold_18", ]
spa_f <- spa_f_gene_median[spa_f_gene_median$chr == "HiC_scaffold_18", ]

swe_m$sample <- "Swe male"
swe_f$sample <- "Swe female"
spa_m$sample <- "Spa male"
spa_f$sample <- "Spa female"

all_Z3_depth <- rbind(swe_m, swe_f, spa_m, spa_f)

all_Z3_Z <- all_Z3_depth[all_Z3_depth$gene_id %in% all_Z_intersect_list, ]
all_Z3_ZW <- all_Z3_depth[all_Z3_depth$gene_id %in% all_ZW_intersect_list, ]

all_Z3_Z$ZW <- "no"
all_Z3_ZW$ZW <- "yes"

all_Z3_depth <- rbind(all_Z3_Z, all_Z3_ZW)


#Z3

Z3_max_coord=17773115

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
Z3_smooth=1000000/Z3_max_coord


# visuals...
chr_size=8
chr_pos=-0.25
chr_text_size=4
chr_text_pos=-0.25
title_size=20
y_lab_size=18
axis_size=16
point_size=2.5

depth_point_color="grey"
depth_line_color="black"

all_Z3_depth_plot <- ggplot(all_Z3_depth, aes(x=gene_midpoint, y=norm_depth)) +
  geom_point(aes(color = sample, shape=ZW), size=point_size, alpha=0.5) +
  scale_color_manual(values = c("#eec76b", "#a9be77", "#ac1917", "#156b8a"), labels=c("Cat female", "Cat male", "Swe female", "Swe male")) +
  scale_shape_manual(values = c(0, 16), guide = "none") +
  geom_smooth(data=all_Z3_depth, aes(y=norm_depth, color=sample), 
              span = Z3_smooth, method = "loess", se=FALSE, size=0.5) +
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
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_size),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_size),
        axis.text.y = element_text(size=axis_size),
        axis.title.y = element_text(size=y_lab_size),
        plot.title = element_text(size=title_size, hjust = 0.5)) +
  scale_x_continuous(labels = unit_format(unit = "mb", scale = 1e-6, accuracy = 1L), expand = c(0, 0)) +
  xlab("") +
  ylab("Normalized read depth\n")

all_Z3_depth_plot

ggsave(filename = ("supplementary_figure_2.png"), width = 15, height = 10, units = "in", dpi = 300, limitsize = F)


