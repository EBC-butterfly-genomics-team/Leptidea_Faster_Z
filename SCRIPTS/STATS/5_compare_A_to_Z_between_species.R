
##### compare evolutionary rates between species for genes that are Z-linked in one species and autosomal in another #####

library(ggplot2)
library(dplyr)
library(rstatix)


setwd("C:/Users/larsf/Desktop/manuscripts/EVOLUTION")

df <- read.csv("lep_sin_mb_dnds_data.txt", sep = "\t")
df_filtered <- df[df$Lsin_dnds < 999, ]
df_filtered <- df_filtered[df_filtered$Lmor_dnds < 999, ]
df_filtered <- df_filtered[df_filtered$Ldup_dnds < 999, ]


Lsin <- df_filtered[, 1:10]
Lsin$species <- "L. sinapis"
names(Lsin)[names(Lsin) == 'Lsin_dnds'] <- 'dnds'
names(Lsin)[names(Lsin) == 'Lsin_dn'] <- 'dn'
names(Lsin)[names(Lsin) == 'Lsin_ds'] <- 'ds'

Lmor <- df_filtered[, c(1:7, 11:13)]
Lmor$species <- "L. morsei"
names(Lmor)[names(Lmor) == 'Lmor_dnds'] <- 'dnds'
names(Lmor)[names(Lmor) == 'Lmor_dn'] <- 'dn'
names(Lmor)[names(Lmor) == 'Lmor_ds'] <- 'ds'

Ldup <- df_filtered[, c(1:7, 14:16)]
Ldup$species <- "L. duponcheli"
names(Ldup)[names(Ldup) == 'Ldup_dnds'] <- 'dnds'
names(Ldup)[names(Ldup) == 'Ldup_dn'] <- 'dn'
names(Ldup)[names(Ldup) == 'Ldup_ds'] <- 'ds'


df_filtered <- rbind(Lsin, Lmor, Ldup)


Z1 <- na.omit(df_filtered[df_filtered$chrN == "Z1", ])
Z2 <- na.omit(df_filtered[df_filtered$chrN == "Z2", ])
Z3 <- na.omit(df_filtered[df_filtered$chrN == "Z3", ])

Auto <- na.omit(df_filtered[df_filtered$Z_vs_A == "A", ])

ggplot(Auto, aes(y=dnds, fill=species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 0.5))

kruskal.test(dnds ~ species, data = Auto)

#Kruskal-Wallis chi-squared = 1.8856, df = 2, p-value = 0.3895
# no difference for autosomes...


ggplot(Z1, aes(y=dnds, fill=species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 0.5))

kruskal.test(dnds ~ species, data = Z1)
#Kruskal-Wallis chi-squared = 0.07288, df = 2, p-value = 0.9642
# no difference for Z1...


#ggplot(Z2, aes(y=dnds, fill=species)) +
#  geom_boxplot(outlier.shape = NA) +
#  coord_cartesian(ylim = c(0, 0.5))

#kruskal.test(dnds ~ species, data = Z2)
#Kruskal-Wallis chi-squared = 10.779, df = 2, p-value = 0.004565

#wilcox_test(Z2, dnds ~ species, p.adjust.method = "bonferroni")

#group1        group2        n1    n2 statistic     p p.adj p.adj.signif
#L. duponcheli L. morsei    313   313    43020. 0.007 0.021 *           
#L. duponcheli L. sinapis   313   313    43142  0.009 0.026 *           
#L. morsei     L. sinapis   313   313    50560. 0.482 1     ns          


ggplot(Z3, aes(y=dnds, x=species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 0.5))

kruskal.test(dnds ~ species, data = Z3)
#Kruskal-Wallis chi-squared = 1.3506, df = 2, p-value = 0.509
# no difference for Z3...



Z2_11_end <- 10209824
Z2_7_start <- 11321879

Z2 <- Z2 %>% mutate(age = case_when(gene_position < Z2_11_end ~ "new",
                              gene_position > Z2_7_start ~ "old"))

Z2 <- na.omit(Z2)

ggplot(Z2, aes(y=dnds, x=species, color=age)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 0.5))

Z2_new <- Z2[Z2$age == "new", ]
Z2_old <- Z2[Z2$age == "old", ]

ggplot(Z2_new, aes(y=dnds, fill=species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 0.5))

kruskal.test(dnds ~ species, data = Z2_new)
#Kruskal-Wallis chi-squared = 6.57, df = 2, p-value = 0.03744

wilcox_test(Z2_new, dnds ~ species, p.adjust.method = "bonferroni")

#.y.   group1        group2        n1    n2 statistic     p p.adj p.adj.signif
#* <chr> <chr>         <chr>      <int> <int>     <dbl> <dbl> <dbl> <chr>       
# 1 dnds  L. duponcheli L. morsei    151   151     9678. 0.021 0.064 ns          
# 2 dnds  L. duponcheli L. sinapis   151   151    10760  0.392 1     ns          
# 3 dnds  L. morsei     L. sinapis   151   151    12889  0.048 0.144 ns        


Z2_new_medians <- aggregate(dnds ~ species, data=Z2_new, FUN=median)
Z2_new_medians

#species   dnds
#1 L. duponcheli 0.0994
#2     L. morsei 0.2047
#3    L. sinapis 0.1332

species   dnds
1 L. duponcheli 0.0994
2     L. morsei 0.2047
3    L. sinapis 0.1332

ggplot(Z2_old, aes(y=dnds, fill=species)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 0.5))

kruskal.test(dnds ~ species, data = Z2_old)
#Kruskal-Wallis chi-squared = 7.5591, df = 2, p-value = 0.02283

wilcox_test(Z2_old, dnds ~ species, p.adjust.method = "bonferroni")







