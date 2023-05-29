
library(ggplot2)
library(ggpubr)
library(gtools)
library(RColorBrewer)
library(fmsb)
library(rstatix)
library(patchwork)
library(car)
library(cowplot)
library(dplyr)
library(zoo)
library(scales)
library(extrafont)
library(tidyverse)
library(ggbreak)
library(boot)
library(ggbeeswarm)
library(see)


#---------------------------------------------#
##### read in data and set plot variables #####
#---------------------------------------------#


setwd("~/Desktop/Plots/FASTZ")

df <- read.csv("lep_sin_dnds_data.txt", sep = "\t")

label_sizes=20
axis_text_sizes=12
sign_size=12
auto_color <- "grey"
fem_col="#ea5545"
male_col="#27aeef"

Z1_col="#89ad4e"
Z2_col="#8e551e"
Z3_col="#efbf5a"

anz_Z_color <- rgb(213/255,117/255,0)


#---------------------------#
##### Filter high dN/dS #####
#---------------------------#


df_filtered <- df[df$dnds < 999, ]


#--------------------------#
##### Restructure data #####
#--------------------------#


Z_dnds <- df_filtered[df_filtered$Z_vs_A == "Z", ]

Z1 <- df_filtered[df_filtered$chrN == "Z1", ]
Z2 <- df_filtered[df_filtered$chrN == "Z2", ]
Z3 <- df_filtered[df_filtered$chrN == "Z3", ]

A_dnds <- df_filtered[df_filtered$anc_vs_neo == "A", ]
median_dnds_A <- median(A_dnds$dnds)
Z_anc_dnds <- df_filtered[df_filtered$anc_vs_neo == "anc Z", ]
median_dnds_anc_Z <- median(Z_anc_dnds$dnds)
Z_neo_dnds <- df_filtered[df_filtered$anc_vs_neo == "neo Z", ]
median_dnds_neo_Z <- median(Z_neo_dnds$dnds)


FBG <- df_filtered[df_filtered$sex_bias == "FBG", ]
FBG_A <- FBG[FBG$Z_vs_A == "A", ]
FBG_Z <- FBG[FBG$Z_vs_A == "Z", ]
FBG_anc_Z <- FBG[FBG$anc_vs_neo == "anc Z", ]
FBG_neo_Z <- FBG[FBG$anc_vs_neo == "neo Z", ]

MBG <- df_filtered[df_filtered$sex_bias == "MBG", ]
MBG_A <- MBG[MBG$Z_vs_A == "A", ]
MBG_Z <- MBG[MBG$Z_vs_A == "Z", ]
MBG_anc_Z <- MBG[MBG$anc_vs_neo == "anc Z", ]
MBG_neo_Z <- MBG[MBG$anc_vs_neo == "neo Z", ]

UBG <- df_filtered[df_filtered$sex_bias == "UBG", ]
UBG_A <- UBG[UBG$Z_vs_A == "A", ]
UBG_Z <- UBG[UBG$Z_vs_A == "Z", ]
UBG_anc_Z <- UBG[UBG$anc_vs_neo == "anc Z", ]
UBG_neo_Z <- UBG[UBG$anc_vs_neo == "neo Z", ]


df_filtered_neo_Z <- df_filtered[df_filtered$anc_vs_neo == "neo Z", ]
df_filtered_anc_Z <- df_filtered[df_filtered$anc_vs_neo == "anc Z", ]
df_filtered_A <- df_filtered[df_filtered$anc_vs_neo == "A", ]



#-------------------#
##### Figure 3A #####
#-------------------#


# dN/dS for SBG #

## A ##

# remove non-called genes (non expressed etc.)
A_dnds_sbg <- A_dnds[!(A_dnds$sex_bias == ""), ]

kruskal.test(dnds ~ sex_bias, data = A_dnds_sbg)
# p-value = 0.008883

pairwise.wilcox.test(A_dnds_sbg$dnds, A_dnds_sbg$sex_bias,
                     p.adjust.method = "bonferroni")
#    FBG    MBG   
#MBG 0.2270 -     
#UBG 1.0000 0.0067

# MBG faster than UBG on A...


## anc Z ##

Z_anc_dnds_sbg <- Z_anc_dnds[!(Z_anc_dnds$sex_bias == ""), ]

kruskal.test(dnds ~ sex_bias, data = Z_anc_dnds_sbg)
# p-value = 0.008102

pairwise.wilcox.test(Z_anc_dnds_sbg$dnds, Z_anc_dnds_sbg$sex_bias,
                     p.adjust.method = "bonferroni")
#    FBG   MBG  
#MBG 0.017 -    
#UBG 0.213 0.067

# FBG faster than MBG on anc Z...


## neo Z ##

Z_neo_dnds_sbg <- Z_neo_dnds[!(Z_neo_dnds$sex_bias == ""), ]

kruskal.test(dnds ~ sex_bias, data = Z_neo_dnds_sbg)
# p-value = 0.1298

# Neo Z no difference...



# plot...

df_filtered_sbg <- df_filtered[!(df_filtered$sex_bias == ""), ]

dnds_A_SBG_plot <- ggplot(data=A_dnds_sbg, aes(y=dnds, x=sex_bias, fill=sex_bias)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(name = "dN/dS") +
  coord_cartesian(ylim=c(0, 0.85)) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.text.x=element_text(size=axis_text_sizes),
        axis.title = element_text(size=axis_text_sizes),
        axis.text.y = element_text(size=axis_text_sizes),
        panel.grid.minor.y = element_line(linewidth=.1, color="black" ),
        panel.grid.major.y = element_line(linewidth=.1, color="black" )) +
  xlab("\nAutosomes") +
  annotate("text", label = "*", x=2.5, y=0.625, size = sign_size) +
  annotate("segment", x=2, xend=3, y=0.625, yend = 0.625) 


dnds_Z_anc_SBG_plot <- ggplot(data=Z_anc_dnds_sbg, aes(y=dnds, x=sex_bias, fill=sex_bias)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(name = "dN/dS") +
  coord_cartesian(ylim=c(0, 0.85)) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
        theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title=element_text(size=axis_text_sizes),
        axis.text=element_text(size=axis_text_sizes),
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor.y = element_line(linewidth=.1, color="black" ),
        panel.grid.major.y = element_line(linewidth=.1, color="black" )) +
  xlab("\nAncestral Z") +
  ylab("") +
  annotate("text", label = "*", x=1.5, y=0.75, size = sign_size) +
  annotate("segment", x=1, xend=2, y=0.75, yend = 0.75) 


Z_neo_dnds_sbg$chrN <- factor(Z_neo_dnds_sbg$chrN, levels = c("Z1", "Z2", "Z3"))

dnds_Z_neo_SBG_plot <- ggplot(data=Z_neo_dnds_sbg, aes(y=dnds, x=sex_bias, fill=sex_bias, color=chrN)) +
  geom_boxplot(outlier.shape = NA, lwd=1) +
  scale_y_continuous(name = "dN/dS") +
  coord_cartesian(ylim=c(0, 0.85)) +
  scale_color_manual("Chr", values = c(Z1_col, Z2_col, Z3_col), labels = c("Z1", "Z2", "Z3")) +
  scale_fill_manual(guide = "none", "Sex bias", values = c(fem_col, male_col, "wheat")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_text(size=axis_text_sizes),
        axis.text.x=element_text(size=axis_text_sizes),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor.y = element_line(linewidth=.1, color="black" ),
        panel.grid.major.y = element_line(linewidth=.1, color="black" )) +
  xlab("\nNeo Z") +
  ylab("")
 


#-------------------#
##### Figure 3B #####
#-------------------#


all_Z <- df_filtered[df_filtered$Z_vs_A == "Z", ]
all_Z <- all_Z[all_Z$sex_bias != "", ]
Z_dnds <- all_Z[all_Z$Z_vs_A == "Z", ]
Z_dnds <- Z_dnds[Z_dnds$sex_bias != "", ]


Z_dnds_FBG <- Z_dnds[Z_dnds$sex_bias == "FBG", ]
Z_dnds_MBG <- Z_dnds[Z_dnds$sex_bias == "MBG", ]


A_proportions = aggregate(gene_id ~ sex_bias, data=A_dnds_sbg, length)
A_proportions$chrN <- "A"
SBG_proportions = aggregate(gene_id ~ sex_bias + chrN, data=all_Z, length)
SBG_proportions <- SBG_proportions[SBG_proportions$chrN != "Z1", ]

Z1_proportions = aggregate(gene_id ~ sex_bias + anc_vs_neo, data=Z1, length)
Z1_proportions <- Z1_proportions[Z1_proportions$sex_bias != "", ]
names(Z1_proportions)[names(Z1_proportions) == "anc_vs_neo"] <- "chrN"

all_proportions <- rbind(SBG_proportions, A_proportions, Z1_proportions)

all_proportions <- all_proportions %>% group_by(chrN) %>% mutate(percent = gene_id/sum(gene_id))

all_proportions$chrN <- factor(all_proportions$chrN, levels = c("A", "Z3", "Z2", "neo Z", "anc Z"))
all_proportions$sex_bias <- factor(all_proportions$sex_bias, levels = c("MBG", "UBG", "FBG"))



SBG_proportions_plot <- ggplot(all_proportions, aes(x=chrN, y=percent, fill=sex_bias, width = 0.95)) +
  coord_cartesian(ylim = c(0,1)) +
  theme(legend.margin=margin(0,0,0,0),
        legend.box.spacing = unit(0, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=axis_text_sizes, color = "white"),
        axis.title.x=element_text(size=axis_text_sizes),
        axis.text.x=element_text(size=axis_text_sizes)) +
  geom_col() +
  scale_fill_manual("", breaks=c("MBG", "UBG", "FBG"), values = c(male_col, "wheat", fem_col)) +
  scale_x_discrete(labels=c("A", "Z3", "Z2", "Neo Z1", "Anc Z1"), ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  annotate("text", x=1, y=-0.075, label="A", size=5, fontface="plain", color="grey30") +
  annotate("text", x=2, y=-0.075, label="Z3", size=5, fontface="plain", color="grey30") +
  annotate("text", x=3, y=-0.075, label="Z2", size=5, fontface="plain", color="grey30") +
  annotate("text", x=4, y=-0.075, label="Neo Z1", size=5, fontface="plain", color="grey30") +
  annotate("text", x=5, y=-0.075, label="Anc Z1", size=5, fontface="plain", color="grey30") +
  annotate("text", x=5, y=-0.1, label="", size=4, fontface="plain") +
  coord_flip() +
  labs(y = "Proportion of genes")



#-------------------#
##### Figure 3C #####
#-------------------#



##### ALL_Z_A #####

ALL_Z_A <- log2(median(Z_dnds$dnds)/median(A_dnds$dnds))

# BS #

ALL_Z_A_median_ratio_function <- function(df_filtered, index) {
  ALL_Z = median(subset(df_filtered[index, 16], df_filtered[index, 11] == "Z"))
  ALL_A = median(subset(df_filtered[index, 16], df_filtered[index, 11] == "A"))
  ALL_Z_A_ratio = log2(ALL_Z/ALL_A)
  return(ALL_Z_A_ratio)
}

ALL_Z_A_median_ratio.res <- boot(data = df_filtered,
                                 statistic = ALL_Z_A_median_ratio_function,
                                 R = 10000,
                                 strata = factor(df_filtered$Z_vs_A))

print(ALL_Z_A_median_ratio.res)
plot(ALL_Z_A_median_ratio.res)

ALL_Z_A_median_ratio_CI <- boot.ci(boot.out = ALL_Z_A_median_ratio.res, type = "perc")
print(ALL_Z_A_median_ratio_CI)

# get range of CI for plott...
ALL_Z_A_CI_low <- ALL_Z_A_median_ratio_CI$percent[[1,4]]
ALL_Z_A_CI_high <- ALL_Z_A_median_ratio_CI$percent[[1,5]]



##### FBG_Z_A #####

FBG_Z_A <- log2(median(FBG_Z$dnds)/median(FBG_A$dnds))

# BS #

FBG_Z_A_median_ratio_function <- function(FBG, index) {
  FBG_Z = median(subset(FBG[index, 16], FBG[index, 11] == "Z"))
  FBG_A = median(subset(FBG[index, 16], FBG[index, 11] == "A"))
  FBG_Z_A_ratio = log2(FBG_Z/FBG_A)
  return(FBG_Z_A_ratio)
}

FBG_Z_A_median_ratio.res <- boot(data = FBG,
                                 statistic = FBG_Z_A_median_ratio_function,
                                 R = 10000, 
                                 strata = factor(FBG$Z_vs_A))
print(FBG_Z_A_median_ratio.res)
plot(FBG_Z_A_median_ratio.res)

FBG_Z_A_median_ratio_CI <- boot.ci(boot.out = FBG_Z_A_median_ratio.res, type = "perc")
print(FBG_Z_A_median_ratio_CI)

# get range of CI for plott...
FBG_Z_A_CI_low <- FBG_Z_A_median_ratio_CI$percent[[1,4]]
FBG_Z_A_CI_high <- FBG_Z_A_median_ratio_CI$percent[[1,5]]


##### MBG_Z_A #####

MBG_Z_A <- log2(median(MBG_Z$dnds)/median(MBG_A$dnds))

# BS #

MBG_Z_A_median_ratio_function <- function(MBG, index) {
  MBG_Z = median(subset(MBG[index, 16], MBG[index, 11] == "Z"))
  MBG_A = median(subset(MBG[index, 16], MBG[index, 11] == "A"))
  MBG_Z_A_ratio = log2(MBG_Z/MBG_A)
  return(MBG_Z_A_ratio)
}

MBG_Z_A_median_ratio.res <- boot(data = MBG,
                                 statistic = MBG_Z_A_median_ratio_function,
                                 R = 10000, 
                                 strata = factor(MBG$Z_vs_A))
print(MBG_Z_A_median_ratio.res)
plot(MBG_Z_A_median_ratio.res)

MBG_Z_A_median_ratio_CI <- boot.ci(boot.out = MBG_Z_A_median_ratio.res, type = "perc")
print(MBG_Z_A_median_ratio_CI)

# get range of CI for plott...
MBG_Z_A_CI_low <- MBG_Z_A_median_ratio_CI$percent[[1,4]]
MBG_Z_A_CI_high <- MBG_Z_A_median_ratio_CI$percent[[1,5]]


##### UBG_Z_A #####

UBG_Z_A <- log2(median(UBG_Z$dnds)/median(UBG_A$dnds))

# BS #

UBG_Z_A_median_ratio_function <- function(UBG, index) {
  UBG_Z = median(subset(UBG[index, 16], UBG[index, 11] == "Z"))
  UBG_A = median(subset(UBG[index, 16], UBG[index, 11] == "A"))
  UBG_Z_A_ratio = log2(UBG_Z/UBG_A)
  return(UBG_Z_A_ratio)
}

UBG_Z_A_median_ratio.res <- boot(data = UBG,
                                 statistic = UBG_Z_A_median_ratio_function,
                                 R = 10000, 
                                 strata = factor(UBG$Z_vs_A))
print(UBG_Z_A_median_ratio.res)
plot(UBG_Z_A_median_ratio.res)

UBG_Z_A_median_ratio_CI <- boot.ci(boot.out = UBG_Z_A_median_ratio.res, type = "perc")
print(UBG_Z_A_median_ratio_CI)

# get range of CI for plott...
UBG_Z_A_CI_low <- UBG_Z_A_median_ratio_CI$percent[[1,4]]
UBG_Z_A_CI_high <- UBG_Z_A_median_ratio_CI$percent[[1,5]]



##### ALL_anc_Z_A #####

ALL_anc_Z_A <- log2(median(Z_anc_dnds$dnds)/median(A_dnds$dnds))

# BS #

ALL_anc_Z_A_median_ratio_function <- function(df_filtered, index) {
  ALL_anc_Z = median(subset(df_filtered[index, 16], df_filtered[index, 12] == "anc Z"))
  ALL_anc_A = median(subset(df_filtered[index, 16], df_filtered[index, 12] == "A"))
  ALL_anc_Z_A_ratio = log2(ALL_anc_Z/ALL_anc_A)
  return(ALL_anc_Z_A_ratio)
}

ALL_anc_Z_A_median_ratio.res <- boot(data = df_filtered,
                                     statistic = ALL_anc_Z_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(df_filtered$anc_vs_neo))
print(ALL_anc_Z_A_median_ratio.res)
plot(ALL_anc_Z_A_median_ratio.res)

ALL_anc_Z_A_median_ratio_CI <- boot.ci(boot.out = ALL_anc_Z_A_median_ratio.res, type = "perc")
print(ALL_anc_Z_A_median_ratio_CI)

# get range of CI for plott...
ALL_anc_Z_A_CI_low <- ALL_anc_Z_A_median_ratio_CI$percent[[1,4]]
ALL_anc_Z_A_CI_high <- ALL_anc_Z_A_median_ratio_CI$percent[[1,5]]


##### FBG_anc_Z_A #####

FBG_anc_Z_A <- log2(median(FBG_anc_Z$dnds)/median(FBG_A$dnds))

# BS #

FBG_anc_Z_A_median_ratio_function <- function(FBG, index) {
  FBG_anc_Z = median(subset(FBG[index, 16], FBG[index, 12] == "anc Z"))
  FBG_anc_A = median(subset(FBG[index, 16], FBG[index, 12] == "A"))
  FBG_anc_Z_A_ratio = log2(FBG_anc_Z/FBG_anc_A)
  return(FBG_anc_Z_A_ratio)
}

FBG_anc_Z_A_median_ratio.res <- boot(data = FBG,
                                     statistic = FBG_anc_Z_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(FBG$anc_vs_neo))
print(FBG_anc_Z_A_median_ratio.res)
plot(FBG_anc_Z_A_median_ratio.res)

FBG_anc_Z_A_median_ratio_CI <- boot.ci(boot.out = FBG_anc_Z_A_median_ratio.res, type = "perc")
print(FBG_anc_Z_A_median_ratio_CI)

# get range of CI for plott...
FBG_anc_Z_A_CI_low <- FBG_anc_Z_A_median_ratio_CI$percent[[1,4]]
FBG_anc_Z_A_CI_high <- FBG_anc_Z_A_median_ratio_CI$percent[[1,5]]


##### MBG_anc_Z_A #####

MBG_anc_Z_A <- log2(median(MBG_anc_Z$dnds)/median(MBG_A$dnds))

# BS #

MBG_anc_Z_A_median_ratio_function <- function(MBG, index) {
  MBG_anc_Z = median(subset(MBG[index, 16], MBG[index, 12] == "anc Z"))
  MBG_anc_A = median(subset(MBG[index, 16], MBG[index, 12] == "A"))
  MBG_anc_Z_A_ratio = log2(MBG_anc_Z/MBG_anc_A)
  return(MBG_anc_Z_A_ratio)
}

MBG_anc_Z_A_median_ratio.res <- boot(data = MBG,
                                     statistic = MBG_anc_Z_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(MBG$anc_vs_neo))
print(MBG_anc_Z_A_median_ratio.res)
plot(MBG_anc_Z_A_median_ratio.res)

MBG_anc_Z_A_median_ratio_CI <- boot.ci(boot.out = MBG_anc_Z_A_median_ratio.res, type = "perc")
print(MBG_anc_Z_A_median_ratio_CI)

# get range of CI for plott...
MBG_anc_Z_A_CI_low <- MBG_anc_Z_A_median_ratio_CI$percent[[1,4]]
MBG_anc_Z_A_CI_high <- MBG_anc_Z_A_median_ratio_CI$percent[[1,5]]


##### UBG_anc_Z_A #####

UBG_anc_Z_A <- log2(median(UBG_anc_Z$dnds)/median(UBG_A$dnds))

# BS #

UBG_anc_Z_A_median_ratio_function <- function(UBG, index) {
  UBG_anc_Z = median(subset(UBG[index, 16], UBG[index, 12] == "anc Z"))
  UBG_anc_A = median(subset(UBG[index, 16], UBG[index, 12] == "A"))
  UBG_anc_Z_A_ratio = log2(UBG_anc_Z/UBG_anc_A)
  return(UBG_anc_Z_A_ratio)
}

UBG_anc_Z_A_median_ratio.res <- boot(data = UBG,
                                     statistic = UBG_anc_Z_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(UBG$anc_vs_neo))
print(UBG_anc_Z_A_median_ratio.res)
plot(UBG_anc_Z_A_median_ratio.res)

UBG_anc_Z_A_median_ratio_CI <- boot.ci(boot.out = UBG_anc_Z_A_median_ratio.res, type = "perc")
print(UBG_anc_Z_A_median_ratio_CI)

# get range of CI for plott...
UBG_anc_Z_A_CI_low <- UBG_anc_Z_A_median_ratio_CI$percent[[1,4]]
UBG_anc_Z_A_CI_high <- UBG_anc_Z_A_median_ratio_CI$percent[[1,5]]



##### ALL_neo_Z_A #####

ALL_neo_Z_A <- log2(median(Z_neo_dnds$dnds)/median(A_dnds$dnds))

# BS #

ALL_neo_Z_A_median_ratio_function <- function(df_filtered, index) {
  ALL_neo_Z = median(subset(df_filtered[index, 16], df_filtered[index, 12] == "neo Z"))
  ALL_neo_A = median(subset(df_filtered[index, 16], df_filtered[index, 12] == "A"))
  ALL_neo_Z_A_ratio = log2(ALL_neo_Z/ALL_neo_A)
  return(ALL_neo_Z_A_ratio)
}

ALL_neo_Z_A_median_ratio.res <- boot(data = df_filtered,
                                     statistic = ALL_neo_Z_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(df_filtered$anc_vs_neo))
print(ALL_neo_Z_A_median_ratio.res)
plot(ALL_neo_Z_A_median_ratio.res)

ALL_neo_Z_A_median_ratio_CI <- boot.ci(boot.out = ALL_neo_Z_A_median_ratio.res, type = "perc")
print(ALL_neo_Z_A_median_ratio_CI)

# get range of CI for plott...
ALL_neo_Z_A_CI_low <- ALL_neo_Z_A_median_ratio_CI$percent[[1,4]]
ALL_neo_Z_A_CI_high <- ALL_neo_Z_A_median_ratio_CI$percent[[1,5]]


##### FBG_neo_Z_A #####

FBG_neo_Z_A <- log2(median(FBG_neo_Z$dnds)/median(FBG_A$dnds))

# BS #

FBG_neo_Z_A_median_ratio_function <- function(FBG, index) {
  FBG_neo_Z = median(subset(FBG[index, 16], FBG[index, 12] == "neo Z"))
  FBG_neo_A = median(subset(FBG[index, 16], FBG[index, 12] == "A"))
  FBG_neo_Z_A_ratio = log2(FBG_neo_Z/FBG_neo_A)
  return(FBG_neo_Z_A_ratio)
}

FBG_neo_Z_A_median_ratio.res <- boot(data = FBG,
                                     statistic = FBG_neo_Z_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(FBG$anc_vs_neo))
print(FBG_neo_Z_A_median_ratio.res)
plot(FBG_neo_Z_A_median_ratio.res)

FBG_neo_Z_A_median_ratio_CI <- boot.ci(boot.out = FBG_neo_Z_A_median_ratio.res, type = "perc")
print(FBG_neo_Z_A_median_ratio_CI)

# get range of CI for plott...
FBG_neo_Z_A_CI_low <- FBG_neo_Z_A_median_ratio_CI$percent[[1,4]]
FBG_neo_Z_A_CI_high <- FBG_neo_Z_A_median_ratio_CI$percent[[1,5]]


##### MBG_neo_Z_A #####

MBG_neo_Z_A <- log2(median(MBG_neo_Z$dnds)/median(MBG_A$dnds))

# BS #

MBG_neo_Z_A_median_ratio_function <- function(MBG, index) {
  MBG_neo_Z = median(subset(MBG[index, 16], MBG[index, 12] == "neo Z"))
  MBG_neo_A = median(subset(MBG[index, 16], MBG[index, 12] == "A"))
  MBG_neo_Z_A_ratio = log2(MBG_neo_Z/MBG_neo_A)
  return(MBG_neo_Z_A_ratio)
}

MBG_neo_Z_A_median_ratio.res <- boot(data = MBG,
                                     statistic = MBG_neo_Z_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(MBG$anc_vs_neo))
print(MBG_neo_Z_A_median_ratio.res)
plot(MBG_neo_Z_A_median_ratio.res)

MBG_neo_Z_A_median_ratio_CI <- boot.ci(boot.out = MBG_neo_Z_A_median_ratio.res, type = "perc")
print(MBG_neo_Z_A_median_ratio_CI)

# get range of CI for plott...
MBG_neo_Z_A_CI_low <- MBG_neo_Z_A_median_ratio_CI$percent[[1,4]]
MBG_neo_Z_A_CI_high <- MBG_neo_Z_A_median_ratio_CI$percent[[1,5]]

##### UBG_anc_Z_A #####

UBG_neo_Z_A <- log2(median(UBG_neo_Z$dnds)/median(UBG_A$dnds))

# BS #

UBG_neo_Z_A_median_ratio_function <- function(UBG, index) {
  UBG_neo_Z = median(subset(UBG[index, 16], UBG[index, 12] == "neo Z"))
  UBG_neo_A = median(subset(UBG[index, 16], UBG[index, 12] == "A"))
  UBG_neo_Z_A_ratio = log2(UBG_neo_Z/UBG_neo_A)
  return(UBG_neo_Z_A_ratio)
}

UBG_neo_Z_A_median_ratio.res <- boot(data = UBG,
                                     statistic = UBG_neo_Z_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(UBG$anc_vs_neo))
print(UBG_neo_Z_A_median_ratio.res)
plot(UBG_neo_Z_A_median_ratio.res)

UBG_neo_Z_A_median_ratio_CI <- boot.ci(boot.out = UBG_neo_Z_A_median_ratio.res, type = "perc")
print(UBG_neo_Z_A_median_ratio_CI)

# get range of CI for plott...
UBG_neo_Z_A_CI_low <- UBG_neo_Z_A_median_ratio_CI$percent[[1,4]]
UBG_neo_Z_A_CI_high <- UBG_neo_Z_A_median_ratio_CI$percent[[1,5]]




# plot...

nudge=0.1
ci_end_size=0.05

Z_A_ratios <- data.frame(sex_bias = c("\nFBG", "\nMBG", "\nUBG", "\nAll",
                                      "\nFBG", "\nMBG", "\nUBG", "\nAll",
                                      "\nFBG", "\nMBG", "\nUBG", "\nAll"),
                         ratio  = c(FBG_Z_A, MBG_Z_A, UBG_Z_A, ALL_Z_A,
                                    FBG_anc_Z_A, MBG_anc_Z_A, UBG_anc_Z_A, ALL_anc_Z_A,
                                    FBG_neo_Z_A, MBG_neo_Z_A, UBG_neo_Z_A, ALL_neo_Z_A),
                         comparison = c("Z vs A", "Z vs A", "Z vs A", "Z vs A", 
                                        "anc-Z vs A", "anc-Z vs A", "anc-Z vs A", "anc-Z vs A", 
                                        "neo-Z vs A", "neo-Z vs A", "neo-Z vs A", "neo-Z vs A"),
                         CI_low = c(FBG_Z_A_CI_low, MBG_Z_A_CI_low, UBG_Z_A_CI_low, ALL_Z_A_CI_low,
                                    FBG_anc_Z_A_CI_low, MBG_anc_Z_A_CI_low, UBG_anc_Z_A_CI_low, ALL_anc_Z_A_CI_low,
                                    FBG_neo_Z_A_CI_low, MBG_neo_Z_A_CI_low, UBG_neo_Z_A_CI_low, ALL_neo_Z_A_CI_low),
                         CI_high = c(FBG_Z_A_CI_high, MBG_Z_A_CI_high, UBG_Z_A_CI_high, ALL_Z_A_CI_high,
                                     FBG_anc_Z_A_CI_high, MBG_anc_Z_A_CI_high, UBG_anc_Z_A_CI_high, ALL_anc_Z_A_CI_high,
                                     FBG_neo_Z_A_CI_high, MBG_neo_Z_A_CI_high, UBG_neo_Z_A_CI_high, ALL_neo_Z_A_CI_high))

Z_A_ratios$comparison <- factor(Z_A_ratios$comparison, levels = c("Z vs A", "anc-Z vs A", "neo-Z vs A"))


Z_A_ratios_plot <- ggplot(Z_A_ratios, aes(y=ratio, x=sex_bias, color=comparison)) +
  geom_point(position = position_dodge(0.4), size=4) +
  coord_cartesian(ylim = c(-2.5, 2.5)) +
  scale_color_manual(values = c("grey", anz_Z_color, "firebrick4")) +
  ylab("log2(median dN/dS Z / median dN/dS A)") +
  xlab("") +
  theme_classic() +
  theme(legend.title = element_blank(),
    axis.text.x = element_text(size=axis_text_sizes),
    axis.title.y = element_text(size=axis_text_sizes),
    axis.title = element_text(size=axis_text_sizes),
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(size=axis_text_sizes)) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), position = position_dodge(0.4), width = 0.2) +
  annotate("segment", x=0, xend = 5, y=0, yend = 0, linetype="dashed")




#---------------------------#
##### combine all plots #####
#---------------------------#


all_SBG_plot <- (dnds_A_SBG_plot + dnds_Z_anc_SBG_plot + dnds_Z_neo_SBG_plot + plot_layout(widths = c(1, 1, 2))) / 
                  (SBG_proportions_plot + Z_A_ratios_plot )

all_SBG_plot + plot_annotation(tag_levels = list(c("A", "", "", "B", "C"))) & 
  theme(plot.tag = element_text(face = 'plain', size = 20))


ggsave(filename = ("figure_3.png"), width = 15, height = 12, units = "in", dpi = 300, limitsize = F)

