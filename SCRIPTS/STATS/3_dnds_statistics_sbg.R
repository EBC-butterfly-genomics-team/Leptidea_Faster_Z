

# Statistical tests of dN/dS for sex biased genes #

# this corresponds to figure 3, table 2, supplementary table 5 & 6

library(boot)
library(rstatix)
library(chisq.posthoc.test)
library(ggplot2)
library(tidyr)

setwd("~/Desktop/Plots/FASTZ")
df <- read.csv("lep_sin_dnds_data.txt", sep = "\t")

#---------------------------#
##### filter dnds < 999 #####
#---------------------------#


df_filtered <- df[df$dnds < 999, ]


#-----------------------------------#
##### check distribution of SBG #####
#-----------------------------------#

# this relates to figure 3c...

Z_dnds <- df_filtered[df_filtered$Z_vs_A == "Z", ]
A_dnds <- df_filtered[df_filtered$Z_vs_A == "A", ]

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


Z1 <- df_filtered[df_filtered$chrN == "Z1", ]
Z1_SBG_genes <- aggregate(gene_id ~ sex_bias, data = Z1, length)
Z1_SBG_genes <- Z1_SBG_genes[!(Z1_SBG_genes$sex_bias == ""), ]
Z1_SBG_genes$chr <- "Z1"

Z2 <- df_filtered[df_filtered$chrN == "Z2", ]
Z2_SBG_genes <- aggregate(gene_id ~ sex_bias, data = Z2, length)
Z2_SBG_genes <- Z2_SBG_genes[!(Z2_SBG_genes$sex_bias == ""), ]
Z2_SBG_genes$chr <- "Z2"

Z3 <- df_filtered[df_filtered$chrN == "Z3", ]
Z3_SBG_genes <- aggregate(gene_id ~ sex_bias, data = Z3, length)
Z3_SBG_genes <- Z3_SBG_genes[!(Z3_SBG_genes$sex_bias == ""), ]
Z3_SBG_genes$chr <- "Z3"

A_dnds_sbg <- A_dnds[!(A_dnds$sex_bias == ""), ]
A_SBG_genes <- aggregate(gene_id ~ sex_bias, data = A_dnds_sbg, length)

# check if skew between Z chromosomes...

all_Z_genes <- rbind(Z1_SBG_genes, Z2_SBG_genes, Z3_SBG_genes)
all_Z_genes_wide <- spread(all_Z_genes, sex_bias, gene_id)
rownames(all_Z_genes_wide) <- all_Z_genes_wide$chr
all_Z_genes_wide$chr <- NULL

mosaicplot(all_Z_genes_wide, color = TRUE, main="Z SBG proportions")

chisq.test(all_Z_genes_wide)$expected
chisq.test(all_Z_genes_wide)
#X-squared = 26.848, df = 4, p-value = 2.133e-05


#-----------------------------------------------#
##### contrast SBG on different chromosomes #####
#-----------------------------------------------#

# this relates to figure 3a...

# remove non-called genes (non expressed)
A_dnds_sbg <- A_dnds[!(A_dnds$sex_bias == ""), ]

kruskal.test(dnds ~ sex_bias, data = A_dnds_sbg)
# p-value = 0.008883
pairwise.wilcox.test(A_dnds_sbg$dnds, A_dnds_sbg$sex_bias,
                     p.adjust.method = "bonferroni")
#    FBG    MBG   
#MBG 0.2270 -     
#UBG 1.0000 0.0067

# MBG faster than UBG on A


# Difference between Z genes depending on sex bias...
Z_dnds <- Z_dnds[!(Z_dnds$sex_bias == ""), ]

kruskal.test(dnds ~ sex_bias, data = Z_dnds)
# p-value = 0.04728

# Overall difference...

pairwise.wilcox.test(Z_dnds$dnds, Z_dnds$sex_bias,
                     p.adjust.method = "bonferroni")

#    FBG   MBG  
#MBG 0.069 -    
#UBG 0.054 1.000

# But no difference between individual groups after correction....


# contrast SBG per Z chromosome

Z1 <- Z1[!(Z1$sex_bias == ""), ]

kruskal.test(dnds ~ sex_bias, data = Z1)
# p-value = 0.007557

pairwise.wilcox.test(Z1$dnds, Z1$sex_bias,
                     p.adjust.method = "bonferroni")
#    FBG    MBG   
#MBG 0.0075 -     
#UBG 0.0316 0.4263

# On Z1 FBG evolve faster!

Z2 <- Z2[!(Z2$sex_bias == ""), ]
kruskal.test(dnds ~ sex_bias, data = Z2)
# p-value = 0.09128

Z3 <- Z3[!(Z3$sex_bias == ""), ]
kruskal.test(dnds ~ sex_bias, data = Z3)
# p-value = 0.9936

# No difference between U, F or M for Z2 and Z3



#-----------------------------------------------------------------------------#
##### Contrast different chromosome classes for each category of sex bias #####
#-----------------------------------------------------------------------------#

# this relates to figure 3b...


## FBG ##

# rename chr numbers to A
FBG_A$chrN <- "A"

all_FBG <- rbind(FBG_Z, FBG_A)

# test all FBG Z vs A...
#wilcox.test(FBG_Z$dnds, FBG_A$dnds)
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

# female biased genes have higher dN/dS on Z1 and Z2 (which are not different) compared A...



## MBG ##

# rename chr numbers to A
MBG_A$chrN <- "A"

all_MBG <- rbind(MBG_Z, MBG_A)

# test all MBG vs A...
#wilcox.test(MBG_Z$dnds, MBG_A$dnds)
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

# rename chr numbers to A
UBG_A$chrN <- "A"

all_UBG <- rbind(UBG_Z, UBG_A)

# test all UBG vs A...
#wilcox.test(UBG_Z$dnds, UBG_A$dnds)
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


#-----------------------------------#
##### Calculate medians for SBG #####
#-----------------------------------#

# table 2...

all_FBG_medians <- aggregate(dnds ~ chrN, data=all_FBG, FUN=median)

## bootstrap confidence intervals of median... ##

boot_dnds_FBG <- boot(FBG$dnds[FBG$Z_vs_A == "A"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_FBG, type = c("perc"))
plot(boot_dnds_FBG)
# ( 0.0490,  0.0725 )  

boot_dnds_FBG_Z1 <- boot(FBG$dnds[FBG$chrN == "Z1"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_FBG_Z1, type = c("perc"))
plot(boot_dnds_FBG_Z1)
# ( 0.0712,  0.3891 )  

boot_dnds_FBG_Z2 <- boot(FBG$dnds[FBG$chrN == "Z2"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_FBG_Z2, type = c("perc"))
plot(boot_dnds_FBG_Z2)
# ( 0.1043,  0.3428 )  

boot_dnds_FBG_Z3 <- boot(FBG$dnds[FBG$chrN == "Z3"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_FBG_Z3, type = c("perc"))
plot(boot_dnds_FBG_Z3)
# ( 0.0628,  0.1362 )  


all_MBG_medians <- aggregate(dnds ~ chrN, data=all_MBG, FUN=median)

boot_dnds_MBG <- boot(MBG$dnds[MBG$Z_vs_A == "A"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_MBG, type = c("perc"))
plot(boot_dnds_MBG)
# ( 0.0597,  0.0800 )  

boot_dnds_MBG_Z1 <- boot(MBG$dnds[MBG$chrN == "Z1"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_MBG_Z1, type = c("perc"))
plot(boot_dnds_MBG_Z1)
# ( 0.0001,  0.0959 )  

boot_dnds_MBG_Z2 <- boot(MBG$dnds[MBG$chrN == "Z2"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_MBG_Z2, type = c("perc"))
plot(boot_dnds_MBG_Z2)
# ( 0.1173,  0.2360 )  

boot_dnds_MBG_Z3 <- boot(MBG$dnds[MBG$chrN == "Z3"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_MBG_Z3, type = c("perc"))
plot(boot_dnds_MBG_Z3)
# ( 0.0671,  0.1486 )  



all_UBG_medians <- aggregate(dnds ~ chrN, data=all_UBG, FUN=median)

boot_dnds_UBG <- boot(UBG$dnds[UBG$Z_vs_A == "A"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_UBG, type = c("perc"))
plot(boot_dnds_UBG)
# ( 0.0474,  0.0562 )  

boot_dnds_UBG_Z1 <- boot(UBG$dnds[UBG$chrN == "Z1"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_UBG_Z1, type = c("perc"))
plot(boot_dnds_UBG_Z1)
# ( 0.0689,  0.1247 )  

boot_dnds_UBG_Z2 <- boot(UBG$dnds[UBG$chrN == "Z2"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_UBG_Z2, type = c("perc"))
plot(boot_dnds_UBG_Z2)
# ( 0.0788,  0.1321 )  

boot_dnds_UBG_Z3 <- boot(UBG$dnds[UBG$chrN == "Z3"], function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_UBG_Z3, type = c("perc"))
plot(boot_dnds_UBG_Z3)
# ( 0.0821,  0.1210 )  





