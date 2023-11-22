
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
library(latex2exp)

#---------------------------------------------#
##### read in data and set plot variables #####
#---------------------------------------------#


setwd("~/Desktop/Plots/FASTZ")

df <- read.csv("lep_sin_dnds_data.txt", sep = "\t")

label_sizes=15
axis_text_sizes=12
sign_size=12
auto_color <- "grey"
fem_col="#ea5545"
male_col="#27aeef"

#Z1_col="#89ad4e"
#Z2_col="#8e551e"
#Z3_col="#efbf5a"

anz_Z_color <- rgb(213/255,117/255,0)
Z1_col <- rgb(213/255,117/255,0)
Z2_col <- rgb(185/255,156/255,107/255)
Z3_col <- rgb(219/255,202/255,105/255)


#---------------------------#
##### Filter high dN/dS #####
#---------------------------#


df_filtered <- df[df$dnds < 999, ]


#--------------------------#
##### Restructure data #####
#--------------------------#

df_filtered_SBG <- df_filtered[df_filtered$sex_bias != "", ]

all_Z <- df_filtered_SBG[df_filtered_SBG$Z_vs_A == "Z", ]
Z1 <- df_filtered_SBG[df_filtered_SBG$chrN == "Z1", ]
Z2 <- df_filtered_SBG[df_filtered_SBG$chrN == "Z2", ]
Z3 <- df_filtered_SBG[df_filtered_SBG$chrN == "Z3", ]
A_dnds <- df_filtered_SBG[df_filtered_SBG$anc_vs_neo == "A", ]



FBG <- df_filtered[df_filtered$sex_bias == "FBG", ]
FBG_A <- df_filtered[df_filtered$sex_bias == "FBG" & df_filtered$Z_vs_A == "A", ]
FBG_Z <- df_filtered[df_filtered$sex_bias == "FBG" & df_filtered$Z_vs_A == "Z", ]
FBG_Z1 <- FBG_Z[FBG_Z$chrN == "Z1", ]
FBG_Z2 <- FBG_Z[FBG_Z$chrN == "Z2", ]
FBG_Z3 <- FBG_Z[FBG_Z$chrN == "Z3", ]

MBG <- df_filtered[df_filtered$sex_bias == "MBG", ]
MBG_A <- df_filtered[df_filtered$sex_bias == "MBG" & df_filtered$Z_vs_A == "A", ]
MBG_Z <- df_filtered[df_filtered$sex_bias == "MBG" & df_filtered$Z_vs_A == "Z", ]
MBG_Z1 <- MBG_Z[MBG_Z$chrN == "Z1", ]
MBG_Z2 <- MBG_Z[MBG_Z$chrN == "Z2", ]
MBG_Z3 <- MBG_Z[MBG_Z$chrN == "Z3", ]


UBG <- df_filtered[df_filtered$sex_bias == "UBG", ]
UBG_A <- df_filtered[df_filtered$sex_bias == "UBG" & df_filtered$Z_vs_A == "A", ]
UBG_Z <- df_filtered[df_filtered$sex_bias == "UBG" & df_filtered$Z_vs_A == "Z", ]
UBG_Z1 <- UBG_Z[UBG_Z$chrN == "Z1", ]
UBG_Z2 <- UBG_Z[UBG_Z$chrN == "Z2", ]
UBG_Z3 <- UBG_Z[UBG_Z$chrN == "Z3", ]




#-------------------#
##### Figure 3A #####
#-------------------#


# dN/dS for SBG #

## A ##

kruskal.test(dnds ~ sex_bias, data = A_dnds)
# p-value = 0.008883

pairwise.wilcox.test(A_dnds$dnds, A_dnds$sex_bias,
                     p.adjust.method = "bonferroni")
#    FBG    MBG   
#MBG 0.2270 -     
#UBG 1.0000 0.0067

# MBG faster than UBG on A...


## Z1 ##

kruskal.test(dnds ~ sex_bias, data = Z1)
# p-value = 0.007557

# Overall difference...

pairwise.wilcox.test(Z1$dnds, Z1$sex_bias,
                     p.adjust.method = "bonferroni")

#    FBG    MBG   
#MBG 0.0075 -     
#UBG 0.0316 0.4263

# female biased genes higher dN/dS than both MBG and UBG, while MBG and UBG not different...

## Z2 ##

kruskal.test(dnds ~ sex_bias, data = Z2)
# p-value = 0.09128

# No difference...

## Z3 ##

kruskal.test(dnds ~ sex_bias, data = Z3)
# p-value = 0.9936

# No difference...

# On Z1, female biased genes have higher dN/dS but for the other Z there is no difference depending on sex bias...



# plot...

dnds_A_SBG_plot <- ggplot(data=A_dnds, aes(y=dnds, x=sex_bias, fill=sex_bias)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(name = expression(omega)) +
  coord_cartesian(ylim=c(0, 1)) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.text.x=element_text(size=axis_text_sizes),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=axis_text_sizes),
        axis.text.y = element_text(size=axis_text_sizes),
        panel.grid.minor.y = element_line(linewidth=.1, color="grey" ),
        panel.grid.major.y = element_line(linewidth=.1, color="grey" )) +
  xlab("\nAutosomes") +
  annotate("text", label = "*", x=2.5, y=0.625, size = sign_size) +
  annotate("segment", x=2, xend=3, y=0.625, yend = 0.625) 


dnds_A_SBG_plot

anz_Z_color <- rgb(213/255,117/255,0)
neo_17_color <- rgb(102/255,141/255,60/255)


dnds_Z1_SBG_plot <- ggplot(data=Z1, aes(y=dnds, x=sex_bias, fill=sex_bias, color=anc_vs_neo)) +
  geom_boxplot(outlier.shape = NA, lwd=1) +
  scale_y_continuous(name = "dN/dS") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_color_manual("Chr", values = c(anz_Z_color, neo_17_color), labels = c("Anc Z", "Bm 17")) +
  scale_fill_manual(guide = "none", "Sex bias", values = c(fem_col, male_col, "wheat")) +
        theme(legend.title = element_blank(),
              legend.position = c(0.8, 0.79),
              legend.text = element_text(size=axis_text_sizes),
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
        panel.grid.minor.y = element_line(linewidth=.1, color="grey" ),
        panel.grid.major.y = element_line(linewidth=.1, color="grey" )) +
  xlab("\nZ1") +
  ylab("") +
  annotate("text", label = "*", x=1.5, y=0.9, size = sign_size) +
  annotate("segment", x=1, xend=2, y=0.9, yend = 0.9) +
  annotate("text", label = "*", x=2, y=1, size = sign_size) +
  annotate("segment", x=1, xend=3, y=1, yend = 1) 


# split Z2 in older and newer regions...

#Z2
neo_11_end <- 10209824
neo_7_start <- 11321879
midpoint <- (10209824+11321879)/2

neo_11_color <- rgb(185/255,156/255,107/255)
neo_7_color <- rgb(129/255,108/255,91/255)
neo_24_color <- rgb(228/255,153/255,105/255)


Z2 <- Z2 %>% 
  mutate(age = case_when(gene_position > midpoint ~ "Bm 24+7",
                         gene_position < midpoint ~ "Bm 11"))

Z2$age <- factor(Z2$age, levels = c("Bm 24+7", "Bm 11"))

dnds_Z2_SBG_plot <- ggplot(data=Z2, aes(y=dnds, x=sex_bias, fill=sex_bias, color=age)) +
  geom_boxplot(outlier.shape = NA, lwd=1) +
  scale_y_continuous(name = "dN/dS") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_fill_manual(guide = "none", "Sex bias", values = c(fem_col, male_col, "wheat")) +
  scale_color_manual("Chr", values = c(neo_7_color, neo_11_color)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.79),
        legend.text = element_text(size=axis_text_sizes),
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
        panel.grid.minor.y = element_line(linewidth=.1, color="grey" ),
        panel.grid.major.y = element_line(linewidth=.1, color="grey" )) +
  xlab("\nZ2") +
  ylab("") 

dnds_Z2_SBG_plot


dnds_Z3_SBG_plot <- ggplot(data=Z3, aes(y=dnds, x=sex_bias, fill=sex_bias)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(name = "dN/dS") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  theme(legend.text = element_text(size=axis_text_sizes),
        legend.title = element_blank(),
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
        panel.grid.minor.y = element_line(linewidth=.1, color="grey" ),
        panel.grid.major.y = element_line(linewidth=.1, color="grey" )) +
  xlab("\nZ3") +
  ylab("") 





#-------------------#
##### Figure 3B #####
#-------------------#

all_Z_SBG <- df_filtered[df_filtered$Z_vs_A == "Z", ]
Autosomes_SBG_genes <- df_filtered[df_filtered$Z_vs_A == "A", ]

## FBG ##

FBG_Z <- all_Z_SBG[all_Z_SBG$sex_bias == "FBG" , ]
FBG_A <- Autosomes_SBG_genes[Autosomes_SBG_genes$sex_bias == "FBG" , ]

# rename chr numbers to A
FBG_A$chrN <- "A"

all_FBG <- rbind(FBG_Z, FBG_A)

# test all FBG vs A...
wilcox.test(FBG_Z$dnds, FBG_A$dnds)
#p-value = 1.229e-06

kruskal.test(dnds ~ chrN, data = all_FBG)
# p-value = 3.54e-06

# Overall difference...

pairwise.wilcox.test(all_FBG$dnds, all_FBG$chrN,
                     p.adjust.method = "bonferroni")

#   A      Z1     Z2    
#Z1 0.0010 -      -     
#Z2 0.0018 1.0000 -     
#Z3 0.8088 0.4952 0.2740

# female biased genes have higher dN/dS on Z1 and Z2 (which are not different) compared to Z3 and A (which are not different)...


## MBG ##

MBG_Z <- all_Z_SBG[all_Z_SBG$sex_bias == "MBG" , ]
MBG_A <- Autosomes_SBG_genes[Autosomes_SBG_genes$sex_bias == "MBG" , ]

# rename chr numbers to A
MBG_A$chrN <- "A"

all_MBG <- rbind(MBG_Z, MBG_A)

# test all MBG vs A...
wilcox.test(MBG_Z$dnds, MBG_A$dnds)
#p-value = 0.0294

kruskal.test(dnds ~ chrN, data = all_MBG)
# p-value = 0.002178

# Overall difference...

pairwise.wilcox.test(all_MBG$dnds, all_MBG$chrN,
                     p.adjust.method = "bonferroni")

#   A      Z1     Z2    
#Z1 1.0000 -      -     
#Z2 0.0014 0.0250 -     
#Z3 1.0000 1.0000 0.9850

# male biased genes on Z2 have higher dN/dS than A and Z1, no other comparison is sign. different...


## UBG ##

UBG_Z <- all_Z_SBG[all_Z_SBG$sex_bias == "UBG" , ]
UBG_A <- Autosomes_SBG_genes[Autosomes_SBG_genes$sex_bias == "UBG" , ]

# rename chr numbers to A
UBG_A$chrN <- "A"

all_UBG <- rbind(UBG_Z, UBG_A)

# test all UBG vs A...
wilcox.test(UBG_Z$dnds, UBG_A$dnds)
#p-value = 1.727e-14

kruskal.test(dnds ~ chrN, data = all_UBG)
# p-value = 2.523e-13

# Overall difference...

pairwise.wilcox.test(all_UBG$dnds, all_UBG$chrN,
                     p.adjust.method = "bonferroni")

#   A       Z1     Z2    
#Z1 3.5e-05 -      -     
#Z2 4.3e-08 0.6796 -     
#Z3 0.0018  1.0000 1.0000

# unbiased genes have higher dN/dS on all Z (which are not different) compared to A...




##### ALL_Z_A #####

ALL_Z_A <- log2(median(all_Z$dnds)/median(A_dnds$dnds))

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



##### ALL_Z1_A #####

ALL_Z1_A <- log2(median(Z1$dnds)/median(A_dnds$dnds))

# BS #

ALL_Z1_A_median_ratio_function <- function(df_filtered, index) {
  ALL_Z1 = median(subset(df_filtered[index, 16], df_filtered[index, 10] == "Z1"))
  ALL_A = median(subset(df_filtered[index, 16], df_filtered[index, 12] == "A"))
  ALL_Z1_A_ratio = log2(ALL_Z1/ALL_A)
  return(ALL_Z1_A_ratio)
}

ALL_Z1_A_median_ratio.res <- boot(data = df_filtered,
                                     statistic = ALL_Z1_A_median_ratio_function,
                                     R = 10000, 
                                     strata = factor(df_filtered$Z_vs_A))
print(ALL_Z1_A_median_ratio.res)
plot(ALL_Z1_A_median_ratio.res)

ALL_Z1_A_median_ratio_CI <- boot.ci(boot.out = ALL_Z1_A_median_ratio.res, type = "perc")
print(ALL_Z1_A_median_ratio_CI)

# get range of CI for plott...
ALL_Z1_A_CI_low <- ALL_Z1_A_median_ratio_CI$percent[[1,4]]
ALL_Z1_A_CI_high <- ALL_Z1_A_median_ratio_CI$percent[[1,5]]


##### ALL_Z2_A #####

ALL_Z2_A <- log2(median(Z2$dnds)/median(A_dnds$dnds))

# BS #

ALL_Z2_A_median_ratio_function <- function(df_filtered, index) {
  ALL_Z2 = median(subset(df_filtered[index, 16], df_filtered[index, 10] == "Z2"))
  ALL_A = median(subset(df_filtered[index, 16], df_filtered[index, 12] == "A"))
  ALL_Z2_A_ratio = log2(ALL_Z2/ALL_A)
  return(ALL_Z2_A_ratio)
}

ALL_Z2_A_median_ratio.res <- boot(data = df_filtered,
                                  statistic = ALL_Z2_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(df_filtered$Z_vs_A))
print(ALL_Z2_A_median_ratio.res)
plot(ALL_Z2_A_median_ratio.res)

ALL_Z2_A_median_ratio_CI <- boot.ci(boot.out = ALL_Z2_A_median_ratio.res, type = "perc")
print(ALL_Z2_A_median_ratio_CI)

# get range of CI for plott...
ALL_Z2_A_CI_low <- ALL_Z2_A_median_ratio_CI$percent[[1,4]]
ALL_Z2_A_CI_high <- ALL_Z2_A_median_ratio_CI$percent[[1,5]]


##### ALL_Z3_A #####

ALL_Z3_A <- log2(median(Z3$dnds)/median(A_dnds$dnds))

# BS #

ALL_Z3_A_median_ratio_function <- function(df_filtered, index) {
  ALL_Z3 = median(subset(df_filtered[index, 16], df_filtered[index, 10] == "Z3"))
  ALL_A = median(subset(df_filtered[index, 16], df_filtered[index, 12] == "A"))
  ALL_Z3_A_ratio = log2(ALL_Z3/ALL_A)
  return(ALL_Z3_A_ratio)
}

ALL_Z3_A_median_ratio.res <- boot(data = df_filtered,
                                  statistic = ALL_Z3_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(df_filtered$Z_vs_A))
print(ALL_Z3_A_median_ratio.res)
plot(ALL_Z3_A_median_ratio.res)

ALL_Z3_A_median_ratio_CI <- boot.ci(boot.out = ALL_Z3_A_median_ratio.res, type = "perc")
print(ALL_Z3_A_median_ratio_CI)

# get range of CI for plott...
ALL_Z3_A_CI_low <- ALL_Z3_A_median_ratio_CI$percent[[1,4]]
ALL_Z3_A_CI_high <- ALL_Z3_A_median_ratio_CI$percent[[1,5]]



##### FBG_Z1_A #####

FBG_Z1_A <- log2(median(FBG_Z1$dnds)/median(FBG_A$dnds))

# BS #

FBG_Z1_A_median_ratio_function <- function(FBG, index) {
  FBG_Z1 = median(subset(FBG[index, 16], FBG[index, 10] == "Z1"))
  FBG_A = median(subset(FBG[index, 16], FBG[index, 12] == "A"))
  FBG_Z1_A_ratio = log2(FBG_Z1/FBG_A)
  return(FBG_Z1_A_ratio)
}

FBG_Z1_A_median_ratio.res <- boot(data = FBG,
                                  statistic = FBG_Z1_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(FBG$Z_vs_A))
print(FBG_Z1_A_median_ratio.res)
plot(FBG_Z1_A_median_ratio.res)

FBG_Z1_A_median_ratio_CI <- boot.ci(boot.out = FBG_Z1_A_median_ratio.res, type = "perc")
print(FBG_Z1_A_median_ratio_CI)

# get range of CI for plott...
FBG_Z1_A_CI_low <- FBG_Z1_A_median_ratio_CI$percent[[1,4]]
FBG_Z1_A_CI_high <- FBG_Z1_A_median_ratio_CI$percent[[1,5]]


##### FBG_Z2_A #####

FBG_Z2_A <- log2(median(FBG_Z2$dnds)/median(FBG_A$dnds))

# BS #

FBG_Z2_A_median_ratio_function <- function(FBG, index) {
  FBG_Z2 = median(subset(FBG[index, 16], FBG[index, 10] == "Z2"))
  FBG_A = median(subset(FBG[index, 16], FBG[index, 12] == "A"))
  FBG_Z2_A_ratio = log2(FBG_Z2/FBG_A)
  return(FBG_Z2_A_ratio)
}

FBG_Z2_A_median_ratio.res <- boot(data = FBG,
                                  statistic = FBG_Z2_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(FBG$Z_vs_A))
print(FBG_Z2_A_median_ratio.res)
plot(FBG_Z2_A_median_ratio.res)

FBG_Z2_A_median_ratio_CI <- boot.ci(boot.out = FBG_Z2_A_median_ratio.res, type = "perc")
print(FBG_Z2_A_median_ratio_CI)

# get range of CI for plott...
FBG_Z2_A_CI_low <- FBG_Z2_A_median_ratio_CI$percent[[1,4]]
FBG_Z2_A_CI_high <- FBG_Z2_A_median_ratio_CI$percent[[1,5]]


##### FBG_Z3_A #####

FBG_Z3_A <- log2(median(FBG_Z3$dnds)/median(FBG_A$dnds))

# BS #

FBG_Z3_A_median_ratio_function <- function(FBG, index) {
  FBG_Z3 = median(subset(FBG[index, 16], FBG[index, 10] == "Z3"))
  FBG_A = median(subset(FBG[index, 16], FBG[index, 12] == "A"))
  FBG_Z3_A_ratio = log2(FBG_Z3/FBG_A)
  return(FBG_Z3_A_ratio)
}

FBG_Z3_A_median_ratio.res <- boot(data = FBG,
                                  statistic = FBG_Z3_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(FBG$Z_vs_A))
print(FBG_Z3_A_median_ratio.res)
plot(FBG_Z3_A_median_ratio.res)

FBG_Z3_A_median_ratio_CI <- boot.ci(boot.out = FBG_Z3_A_median_ratio.res, type = "perc")
print(FBG_Z3_A_median_ratio_CI)

# get range of CI for plott...
FBG_Z3_A_CI_low <- FBG_Z3_A_median_ratio_CI$percent[[1,4]]
FBG_Z3_A_CI_high <- FBG_Z3_A_median_ratio_CI$percent[[1,5]]




##### MBG_Z1_A #####

MBG_Z1_A <- log2(median(MBG_Z1$dnds)/median(MBG_A$dnds))

# BS #

MBG_Z1_A_median_ratio_function <- function(MBG, index) {
  MBG_Z1 = median(subset(MBG[index, 16], MBG[index, 10] == "Z1"))
  MBG_A = median(subset(MBG[index, 16], MBG[index, 12] == "A"))
  MBG_Z1_A_ratio = log2(MBG_Z1/MBG_A)
  return(MBG_Z1_A_ratio)
}

MBG_Z1_A_median_ratio.res <- boot(data = MBG,
                                  statistic = MBG_Z1_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(MBG$Z_vs_A))
print(MBG_Z1_A_median_ratio.res)
plot(MBG_Z1_A_median_ratio.res)

MBG_Z1_A_median_ratio_CI <- boot.ci(boot.out = MBG_Z1_A_median_ratio.res, type = "perc")
print(MBG_Z1_A_median_ratio_CI)

# get range of CI for plott...
MBG_Z1_A_CI_low <- MBG_Z1_A_median_ratio_CI$percent[[1,4]]
MBG_Z1_A_CI_high <- MBG_Z1_A_median_ratio_CI$percent[[1,5]]


##### MBG_Z2_A #####

MBG_Z2_A <- log2(median(MBG_Z2$dnds)/median(MBG_A$dnds))

# BS #

MBG_Z2_A_median_ratio_function <- function(MBG, index) {
  MBG_Z2 = median(subset(MBG[index, 16], MBG[index, 10] == "Z2"))
  MBG_A = median(subset(MBG[index, 16], MBG[index, 12] == "A"))
  MBG_Z2_A_ratio = log2(MBG_Z2/MBG_A)
  return(MBG_Z2_A_ratio)
}

MBG_Z2_A_median_ratio.res <- boot(data = MBG,
                                  statistic = MBG_Z2_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(MBG$Z_vs_A))
print(MBG_Z2_A_median_ratio.res)
plot(MBG_Z2_A_median_ratio.res)

MBG_Z2_A_median_ratio_CI <- boot.ci(boot.out = MBG_Z2_A_median_ratio.res, type = "perc")
print(MBG_Z2_A_median_ratio_CI)

# get range of CI for plott...
MBG_Z2_A_CI_low <- MBG_Z2_A_median_ratio_CI$percent[[1,4]]
MBG_Z2_A_CI_high <- MBG_Z2_A_median_ratio_CI$percent[[1,5]]


##### MBG_Z3_A #####

MBG_Z3_A <- log2(median(MBG_Z3$dnds)/median(MBG_A$dnds))

# BS #

MBG_Z3_A_median_ratio_function <- function(MBG, index) {
  MBG_Z3 = median(subset(MBG[index, 16], MBG[index, 10] == "Z3"))
  MBG_A = median(subset(MBG[index, 16], MBG[index, 12] == "A"))
  MBG_Z3_A_ratio = log2(MBG_Z3/MBG_A)
  return(MBG_Z3_A_ratio)
}

MBG_Z3_A_median_ratio.res <- boot(data = MBG,
                                  statistic = MBG_Z3_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(MBG$Z_vs_A))
print(MBG_Z3_A_median_ratio.res)
plot(MBG_Z3_A_median_ratio.res)

MBG_Z3_A_median_ratio_CI <- boot.ci(boot.out = MBG_Z3_A_median_ratio.res, type = "perc")
print(MBG_Z3_A_median_ratio_CI)

# get range of CI for plott...
MBG_Z3_A_CI_low <- MBG_Z3_A_median_ratio_CI$percent[[1,4]]
MBG_Z3_A_CI_high <- MBG_Z3_A_median_ratio_CI$percent[[1,5]]


##### UBG_Z1_A #####

UBG_Z1_A <- log2(median(UBG_Z1$dnds)/median(UBG_A$dnds))

# BS #

UBG_Z1_A_median_ratio_function <- function(UBG, index) {
  UBG_Z1 = median(subset(UBG[index, 16], UBG[index, 10] == "Z1"))
  UBG_A = median(subset(UBG[index, 16], UBG[index, 12] == "A"))
  UBG_Z1_A_ratio = log2(UBG_Z1/UBG_A)
  return(UBG_Z1_A_ratio)
}

UBG_Z1_A_median_ratio.res <- boot(data = UBG,
                                  statistic = UBG_Z1_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(UBG$Z_vs_A))
print(UBG_Z1_A_median_ratio.res)
plot(UBG_Z1_A_median_ratio.res)

UBG_Z1_A_median_ratio_CI <- boot.ci(boot.out = UBG_Z1_A_median_ratio.res, type = "perc")
print(UBG_Z1_A_median_ratio_CI)

# get range of CI for plott...
UBG_Z1_A_CI_low <- UBG_Z1_A_median_ratio_CI$percent[[1,4]]
UBG_Z1_A_CI_high <- UBG_Z1_A_median_ratio_CI$percent[[1,5]]



##### UBG_Z2_A #####

UBG_Z2_A <- log2(median(UBG_Z2$dnds)/median(UBG_A$dnds))

# BS #

UBG_Z2_A_median_ratio_function <- function(UBG, index) {
  UBG_Z2 = median(subset(UBG[index, 16], UBG[index, 10] == "Z2"))
  UBG_A = median(subset(UBG[index, 16], UBG[index, 12] == "A"))
  UBG_Z2_A_ratio = log2(UBG_Z2/UBG_A)
  return(UBG_Z2_A_ratio)
}

UBG_Z2_A_median_ratio.res <- boot(data = UBG,
                                  statistic = UBG_Z2_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(UBG$Z_vs_A))
print(UBG_Z2_A_median_ratio.res)
plot(UBG_Z2_A_median_ratio.res)

UBG_Z2_A_median_ratio_CI <- boot.ci(boot.out = UBG_Z2_A_median_ratio.res, type = "perc")
print(UBG_Z2_A_median_ratio_CI)

# get range of CI for plott...
UBG_Z2_A_CI_low <- UBG_Z2_A_median_ratio_CI$percent[[1,4]]
UBG_Z2_A_CI_high <- UBG_Z2_A_median_ratio_CI$percent[[1,5]]


##### UBG_Z3_A #####

UBG_Z3_A <- log2(median(UBG_Z3$dnds)/median(UBG_A$dnds))

# BS #

UBG_Z3_A_median_ratio_function <- function(UBG, index) {
  UBG_Z3 = median(subset(UBG[index, 16], UBG[index, 10] == "Z3"))
  UBG_A = median(subset(UBG[index, 16], UBG[index, 12] == "A"))
  UBG_Z3_A_ratio = log2(UBG_Z3/UBG_A)
  return(UBG_Z3_A_ratio)
}

UBG_Z3_A_median_ratio.res <- boot(data = UBG,
                                  statistic = UBG_Z3_A_median_ratio_function,
                                  R = 10000, 
                                  strata = factor(UBG$Z_vs_A))
print(UBG_Z3_A_median_ratio.res)
plot(UBG_Z3_A_median_ratio.res)

UBG_Z3_A_median_ratio_CI <- boot.ci(boot.out = UBG_Z3_A_median_ratio.res, type = "perc")
print(UBG_Z3_A_median_ratio_CI)

# get range of CI for plott...
UBG_Z3_A_CI_low <- UBG_Z3_A_median_ratio_CI$percent[[1,4]]
UBG_Z3_A_CI_high <- UBG_Z3_A_median_ratio_CI$percent[[1,5]]



# plot...

nudge=0.1
ci_end_size=0.05

Z_A_ratios <- data.frame(sex_bias = c("\nFBG", "\nMBG", "\nUBG",
                                      "\nFBG", "\nMBG", "\nUBG",
                                      "\nFBG", "\nMBG", "\nUBG",
                                      "\nFBG", "\nMBG", "\nUBG"),
                         ratio  = c(FBG_Z_A, MBG_Z_A, UBG_Z_A,
                                    FBG_Z1_A, MBG_Z1_A, UBG_Z1_A,
                                    FBG_Z2_A, MBG_Z2_A, UBG_Z2_A,
                                    FBG_Z3_A, MBG_Z3_A, UBG_Z3_A),
                         sign = c("Sign.", "Sign.", "Sign.",
                                  "Sign.", "Non sign.", "Sign.",
                                  "Sign.", "Sign.", "Sign.",
                                  "Non sign.", "Non sign.", "Sign."),
                         comparison = c("Z vs A", "Z vs A", "Z vs A", 
                                        "Z1 vs A", "Z1 vs A", "Z1 vs A", 
                                        "Z2 vs A", "Z2 vs A", "Z2 vs A",
                                        "Z3 vs A", "Z3 vs A", "Z3 vs A"),
                         CI_low = c(FBG_Z_A_CI_low, MBG_Z_A_CI_low, UBG_Z_A_CI_low,
                                    FBG_Z1_A_CI_low, MBG_Z1_A_CI_low, UBG_Z1_A_CI_low,
                                    FBG_Z2_A_CI_low, MBG_Z2_A_CI_low, UBG_Z2_A_CI_low,
                                    FBG_Z3_A_CI_low, MBG_Z3_A_CI_low, UBG_Z3_A_CI_low),
                         CI_high = c(FBG_Z_A_CI_high, MBG_Z_A_CI_high, UBG_Z_A_CI_high,
                                     FBG_Z1_A_CI_high, MBG_Z1_A_CI_high, UBG_Z1_A_CI_high,
                                     FBG_Z2_A_CI_high, MBG_Z2_A_CI_high, UBG_Z2_A_CI_high,
                                     FBG_Z3_A_CI_high, MBG_Z3_A_CI_high, UBG_Z3_A_CI_high))

Z_A_ratios$comparison <- factor(Z_A_ratios$comparison, levels = c("Z vs A", "Z1 vs A", "Z2 vs A", "Z3 vs A"))
Z_A_ratios$sign <- factor(Z_A_ratios$sign, levels = c("Sign.", "Non sign."))

Z_A_ratios_plot <- ggplot(Z_A_ratios, aes(y=ratio, x=sex_bias, color=comparison, shape=sign)) +
  geom_point(position = position_dodge(0.4), size=4) +
  coord_cartesian(ylim = c(-2.5, 2.5)) +
  scale_color_manual(values = c("grey", Z1_col, Z2_col, Z3_col)) +
  scale_shape_manual(values = c(17, 15)) +
  ylab(bquote(log2(median~omega[Z] / median~omega[A]))) +
  xlab("") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text.x = element_text(size=axis_text_sizes),
        axis.title.y = element_text(size=label_sizes),
        axis.title = element_text(size=axis_text_sizes),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size=axis_text_sizes)) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), position = position_dodge(0.4), width = 0.2) +
  annotate("segment", x=0, xend = 4, y=0, yend = 0, linetype="dashed") 

Z_A_ratios_plot


#-------------------#
##### Figure 3C #####
#-------------------#


# count genes and aggregate data...

A_proportions = aggregate(gene_id ~ sex_bias, data=A_dnds, length)
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
#  coord_cartesian(ylim = c(0,1)) +
  theme(legend.text = element_text(size=axis_text_sizes),
        legend.margin=margin(0,0,0,0),
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
  scale_fill_manual("", breaks=c("FBG", "UBG", "MBG"), values = c(fem_col, "wheat", male_col)) +
  scale_x_discrete(labels=c("A", "Z3", "Z2", "Neo Z1", "Anc Z1"), ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  annotate("text", x=1, y=-0.075, label="A", size=5, fontface="plain", color="black") +
  annotate("text", x=2, y=-0.075, label="Z3", size=5, fontface="plain", color="black") +
  annotate("text", x=3, y=-0.075, label="Z2", size=5, fontface="plain", color="black") +
  annotate("text", x=4, y=-0.075, label="Neo Z1", size=5, fontface="plain", color="black") +
  annotate("text", x=5, y=-0.075, label="Anc Z1", size=5, fontface="plain", color="black") +
  coord_flip() +
  labs(y = "Proportion of genes") +
  annotate("text", label = "", x=6, y=-0.15, size = 1) 

SBG_proportions_plot





#---------------------------#
##### combine all plots #####
#---------------------------#


all_SBG_plot <- (dnds_A_SBG_plot + dnds_Z1_SBG_plot + dnds_Z2_SBG_plot + dnds_Z3_SBG_plot + plot_layout(widths = c(1, 1, 1, 1))) / 
                  (Z_A_ratios_plot + SBG_proportions_plot)

all_SBG_plot + plot_annotation(tag_levels = list(c("A", "", "", "", "B", "C"))) & 
  theme(plot.tag = element_text(face = 'plain', size = 20))


ggsave(filename = ("figure_3.jpeg"), width = 15, height = 12, units = "in", dpi = 600, limitsize = F)

