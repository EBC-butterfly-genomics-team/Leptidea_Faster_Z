
#--------------------------------#
### Statistical tests of dN/dS ###
#--------------------------------#


# This corresponds to figure 2, table 1, supplementary figure 4, supplementary table 2, 3 & 4


library(boot)
library(rstatix)
library(chisq.posthoc.test)
library(ggplot2)
library(patchwork)
library(tidyr)
library(car)

setwd("~/Desktop/Plots/FASTZ")
df <- read.csv("lep_sin_dnds_data.txt", sep = "\t")

#------------------------------------------------------------------------#
##### calculate dN/dS based on sum across all genes before filtering #####
#------------------------------------------------------------------------#

# sum Non-synonymous...
Dn_sum <- aggregate(Dn ~ anc_vs_neo, FUN = sum, data = df)
Ns_sum <- aggregate(Ns ~ anc_vs_neo, FUN = sum, data = df)
N_sum <- merge(Dn_sum, Ns_sum, by = "anc_vs_neo")

# sum Synonymous...
Ds_sum <- aggregate(Ds ~ anc_vs_neo, FUN = sum, data = df)
Ss_sum <- aggregate(Ss ~ anc_vs_neo, FUN = sum, data = df)
S_sum <- merge(Ds_sum, Ss_sum, by = "anc_vs_neo")

# calculate dN/dS
NS_sum <- merge(N_sum, S_sum, by = "anc_vs_neo")
NS_sum$dN <- NS_sum$Dn/NS_sum$Ns
NS_sum$dS <- NS_sum$Ds/NS_sum$Ss
NS_sum$dNdS <- NS_sum$dN/NS_sum$dS

# check how many genes with Ds = 0 are in each category...
Ds_0 <- df[df$Ds == 0, ]
number_of_genes <- aggregate(gene_id ~ anc_vs_neo, data = df, length)
colnames(number_of_genes)[2] = "before"
number_of_genes_filtered <- aggregate(gene_id ~ anc_vs_neo, data = Ds_0, length)
colnames(number_of_genes_filtered)[2] = "after"

compare_filtered <- cbind(number_of_genes, number_of_genes_filtered)
compare_filtered$frac <- ((compare_filtered$after)/(compare_filtered$before))


#----------------------------------------------------#
##### Pair-wise Permutation tests between groups #####
#----------------------------------------------------#


# test for difference in point estimate...
A_vs_ancZ <- df[df$anc_vs_neo != "neo Z", ]
A_vs_neoZ <- df[df$anc_vs_neo != "anc Z", ]
ancZ_vs_neoZ <- df[df$anc_vs_neo != "A", ]

n <- 9999

# A vs anc Z...

A_ancZ_res <- numeric(n)
for (i in 1:n) {
  # reshuffle variable...
  perm <- sample(nrow(A_vs_ancZ))
  shuffled <- transform(A_vs_ancZ,anc_vs_neo=anc_vs_neo[perm])
  # compute & store difference in point estimate...
  A_ancZ_res[i] <- ((sum(shuffled$Dn[shuffled$anc_vs_neo=="A"])/sum(shuffled$Ns[shuffled$anc_vs_neo=="A"]))/(sum(shuffled$Ds[shuffled$anc_vs_neo=="A"])/sum(shuffled$Ss[shuffled$anc_vs_neo=="A"]))) -
    ((sum(shuffled$Dn[shuffled$anc_vs_neo=="anc Z"])/sum(shuffled$Ns[shuffled$anc_vs_neo=="anc Z"]))/(sum(shuffled$Ds[shuffled$anc_vs_neo=="anc Z"])/sum(shuffled$Ss[shuffled$anc_vs_neo=="anc Z"]))) 
}

# calculate observed...
obs <- ((sum(A_vs_ancZ$Dn[A_vs_ancZ$anc_vs_neo=="A"])/sum(A_vs_ancZ$Ns[A_vs_ancZ$anc_vs_neo=="A"]))/(sum(A_vs_ancZ$Ds[A_vs_ancZ$anc_vs_neo=="A"])/sum(A_vs_ancZ$Ss[A_vs_ancZ$anc_vs_neo=="A"]))) -
  ((sum(A_vs_ancZ$Dn[A_vs_ancZ$anc_vs_neo=="anc Z"])/sum(A_vs_ancZ$Ns[A_vs_ancZ$anc_vs_neo=="anc Z"]))/(sum(A_vs_ancZ$Ds[A_vs_ancZ$anc_vs_neo=="anc Z"])/sum(A_vs_ancZ$Ss[A_vs_ancZ$anc_vs_neo=="anc Z"]))) 

A_ancZ_res_df <- as.data.frame(A_ancZ_res)

A_vs_ancZ_perm_plot <- ggplot(A_ancZ_res_df, aes(x=A_ancZ_res)) +
  geom_histogram(fill="grey80", color="black") +
  geom_vline(xintercept = obs, color="red") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  xlab("Permutation test distribution")
A_vs_ancZ_perm_plot

# observed diff is lower than sampled distribution...


# combine with resample for p-value...
A_ancZ_res <- c(A_ancZ_res,obs)

# two-tailed test 
3*((2*mean(obs>=A_neoZ_res)))

# p < 2e-04
# Bonferroni p < 6e-04


# A vs neo Z

A_neoZ_res <- numeric(n)
for (i in 1:n) {
  # reshuffle variable...
  perm <- sample(nrow(A_vs_neoZ))
  shuffled <- transform(A_vs_neoZ,anc_vs_neo=anc_vs_neo[perm])
  # compute & store difference in point estimate...
  A_neoZ_res[i] <- ((sum(shuffled$Dn[shuffled$anc_vs_neo=="A"])/sum(shuffled$Ns[shuffled$anc_vs_neo=="A"]))/(sum(shuffled$Ds[shuffled$anc_vs_neo=="A"])/sum(shuffled$Ss[shuffled$anc_vs_neo=="A"]))) -
    ((sum(shuffled$Dn[shuffled$anc_vs_neo=="neo Z"])/sum(shuffled$Ns[shuffled$anc_vs_neo=="neo Z"]))/(sum(shuffled$Ds[shuffled$anc_vs_neo=="neo Z"])/sum(shuffled$Ss[shuffled$anc_vs_neo=="neo Z"]))) 
}

# calculate observed...
obs <- ((sum(A_vs_neoZ$Dn[A_vs_neoZ$anc_vs_neo=="A"])/sum(A_vs_neoZ$Ns[A_vs_neoZ$anc_vs_neo=="A"]))/(sum(A_vs_neoZ$Ds[A_vs_neoZ$anc_vs_neo=="A"])/sum(A_vs_neoZ$Ss[A_vs_neoZ$anc_vs_neo=="A"]))) -
  ((sum(A_vs_neoZ$Dn[A_vs_neoZ$anc_vs_neo=="neo Z"])/sum(A_vs_neoZ$Ns[A_vs_neoZ$anc_vs_neo=="neo Z"]))/(sum(A_vs_neoZ$Ds[A_vs_neoZ$anc_vs_neo=="neo Z"])/sum(A_vs_neoZ$Ss[A_vs_neoZ$anc_vs_neo=="neo Z"]))) 

A_neoZ_res_df <- as.data.frame(A_neoZ_res)

A_vs_neoZ_perm_plot <- ggplot(A_neoZ_res_df, aes(x=A_neoZ_res)) +
  geom_histogram(fill="grey80", color="black") +
  geom_vline(xintercept = obs, color="red") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  xlab("Permutation test distribution")
A_vs_neoZ_perm_plot

# observed diff is lower than sampled distribution...

# combine with resample for p-value...
A_neoZ_res <- c(A_neoZ_res,obs)

# two-tailed test 
3*((2*mean(obs>=A_neoZ_res)))


# p < 2e-04 (same as previous test because we never observe a more extreme value)
# Bonferroni p < 6e-04

# neo vs anc Z...

ancZ_vs_neoZ_res <- numeric(n)
for (i in 1:n) {
  # reshuffle variable...
  perm <- sample(nrow(ancZ_vs_neoZ))
  shuffled <- transform(ancZ_vs_neoZ,anc_vs_neo=anc_vs_neo[perm])
  # compute & store difference in point estimate...
  ancZ_vs_neoZ_res[i] <- ((sum(shuffled$Dn[shuffled$anc_vs_neo=="anc Z"])/sum(shuffled$Ns[shuffled$anc_vs_neo=="anc Z"]))/(sum(shuffled$Ds[shuffled$anc_vs_neo=="anc Z"])/sum(shuffled$Ss[shuffled$anc_vs_neo=="anc Z"]))) -
    ((sum(shuffled$Dn[shuffled$anc_vs_neo=="neo Z"])/sum(shuffled$Ns[shuffled$anc_vs_neo=="neo Z"]))/(sum(shuffled$Ds[shuffled$anc_vs_neo=="neo Z"])/sum(shuffled$Ss[shuffled$anc_vs_neo=="neo Z"]))) 
}

# calculate observed...
obs <- ((sum(ancZ_vs_neoZ$Dn[ancZ_vs_neoZ$anc_vs_neo=="anc Z"])/sum(ancZ_vs_neoZ$Ns[ancZ_vs_neoZ$anc_vs_neo=="anc Z"]))/(sum(ancZ_vs_neoZ$Ds[ancZ_vs_neoZ$anc_vs_neo=="anc Z"])/sum(ancZ_vs_neoZ$Ss[ancZ_vs_neoZ$anc_vs_neo=="anc Z"]))) -
  ((sum(ancZ_vs_neoZ$Dn[ancZ_vs_neoZ$anc_vs_neo=="neo Z"])/sum(ancZ_vs_neoZ$Ns[ancZ_vs_neoZ$anc_vs_neo=="neo Z"]))/(sum(ancZ_vs_neoZ$Ds[ancZ_vs_neoZ$anc_vs_neo=="neo Z"])/sum(ancZ_vs_neoZ$Ss[ancZ_vs_neoZ$anc_vs_neo=="neo Z"]))) 


ancZ_vs_neoZ_res_df <- as.data.frame(ancZ_vs_neoZ_res)

ancZ_vs_neoZ_perm_plot <- ggplot(ancZ_vs_neoZ_res_df, aes(x=ancZ_vs_neoZ_res)) +
  geom_histogram(fill="grey80", color="black") +
  geom_vline(xintercept = obs, color="red") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  xlab("Permutation test distribution")
ancZ_vs_neoZ_perm_plot


# observed diff is not different from sampled distribution...

# combine with resample...
ancZ_vs_neoZ_res <- c(ancZ_vs_neoZ_res,obs)

# two-tailed test 
(2*mean(ancZ_vs_neoZ_res>=obs))

# p = 0.5664
# Bonferroni p = 1



## plot perm test ##

all_perm_plots <- A_vs_ancZ_perm_plot + A_vs_neoZ_perm_plot + ancZ_vs_neoZ_perm_plot + plot_annotation(tag_levels = "A")

all_perm_plots

ggsave(filename = ("supplementary_figure_4.png"), width = 15, height = 5, units = "in", dpi = 300, limitsize = F)


#---------------------------------------------------#
##### filter dnds < 999 for gene wise estimates #####
#---------------------------------------------------#


df_filtered <- df[df$dnds < 999, ]


# count genes before and after filtering #
number_of_genes <- aggregate(gene_id ~ anc_vs_neo, data = df, length)
colnames(number_of_genes)[2] = "before"
number_of_genes_filtered <- aggregate(gene_id ~ anc_vs_neo, data = df_filtered, length)
colnames(number_of_genes_filtered)[2] = "after"

compare_filtered <- cbind(number_of_genes, number_of_genes_filtered)
compare_filtered$frac <- (1-((compare_filtered$after)/(compare_filtered$before)))*100

# removed by filtering dnds < 999
# A: 5.64% ancZ: 7.83% neoZ: 8.22% 


#---------------#
##### dN/dS #####
#---------------#


## dnds A vs Z ##

Z_dnds <- df_filtered[df_filtered$Z_vs_A == "Z", ]
A_dnds <- df_filtered[df_filtered$Z_vs_A == "A", ]

median(A_dnds$dnds)
# 0.0566
median(Z_dnds$dnds)
# 0.1023


## bootstrap confidence intervals of median... ##

boot_dnds_A <- boot(A_dnds$dnds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_A, type = c("perc"))
plot(boot_dnds_A)
# ( 0.0525,  0.0608 )

boot_dnds_Z <- boot(Z_dnds$dnds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_Z, type = c("perc"))
plot(boot_dnds_Z)
# ( 0.0932,  0.1149 )  

## test difference of distributions Z vs A... ##
wilcox.test(dnds ~ Z_vs_A, data = df_filtered, paired=FALSE, conf.int=TRUE)
# p-value < 2.2e-16


## dn A vs Z ##

median(A_dnds$dn)
# 0.0011
median(Z_dnds$dn)
# 0.0021

boot_dn_A <- boot(A_dnds$dn, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dn_A, type = c("perc"))
plot(boot_dn_A)
# ( 0.0010,  0.0012 ) 

boot_dn_Z <- boot(Z_dnds$dn, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dn_Z, type = c("perc"))
plot(boot_dn_Z)
# ( 0.0018,  0.0023 )  


## test difference of distributions Z vs A... ##

wilcox.test(dn ~ Z_vs_A, data = df_filtered, paired=FALSE)
# p-value < 2.2e-16


# ds A vs Z #

median(A_dnds$ds)
# 0.0178
median(Z_dnds$ds)
# 0.0186


## bootstrap confidence intervals of median... ##

boot_ds_A <- boot(A_dnds$ds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_ds_A, type = c("perc"))
plot(boot_ds_A)
# ( 0.0175,  0.0181 )  

boot_ds_Z <- boot(Z_dnds$ds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_ds_Z, type = c("perc"))
plot(boot_ds_Z)
# ( 0.0178,  0.0192 )


## test difference of distributions Z vs A... ##

wilcox.test(ds ~ Z_vs_A, data = df_filtered, paired=FALSE)
# p-value = 0.001892



#--------------------------------#
##### dnds anc vs neo Z vs A #####
#--------------------------------#


Z_anc_dnds <- df_filtered[df_filtered$anc_vs_neo == "anc Z", ]
Z_neo_dnds <- df_filtered[df_filtered$anc_vs_neo == "neo Z", ]

median(Z_anc_dnds$dnds)
# 0.0891
median(Z_neo_dnds$dnds)
# 0.1101


## bootstrap confidence intervals of median... ##

boot_dnds_anc_Z <- boot(Z_anc_dnds$dnds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_anc_Z, type = c("perc"))
plot(boot_dnds_anc_Z)
# ( 0.0681,  0.1111 )  

boot_dnds_neo_Z <- boot(Z_neo_dnds$dnds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_neo_Z, type = c("perc"))
plot(boot_dnds_neo_Z)
# ( 0.0986,  0.1203 )  


## test difference dN/dS anc vs neo vs A... ##

kruskal.test(dnds ~ anc_vs_neo, data = df_filtered)
# p-value < 2.2e-16

pairwise.wilcox.test(df_filtered$dnds, df_filtered$anc_vs_neo,
                     p.adjust.method = "bonferroni")
#      A       anc Z
#anc Z 2.1e-05 -    
#neo Z 2.9e-16 0.71 

# dn anc vs neo Z vs A #

median(Z_anc_dnds$dn)
# 0.0014
median(Z_neo_dnds$dn)
# 0.0024


## bootstrap confidence intervals of median... ##

boot_dn_anc_Z <- boot(Z_anc_dnds$dn, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dn_anc_Z, type = c("perc"))
plot(boot_dn_anc_Z)
# ( 0.0011,  0.0018 )

boot_dn_neo_Z <- boot(Z_neo_dnds$dn, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dn_neo_Z, type = c("perc"))
plot(boot_dn_neo_Z)
# ( 0.0022,  0.0026 ) 


## test difference dN anc vs neo vs A... ##

kruskal.test(dn ~ anc_vs_neo, data = df_filtered)
# p-value < 2.2e-16

pairwise.wilcox.test(df_filtered$dn, df_filtered$anc_vs_neo,
                     p.adjust.method = "bonferroni")
#      A       anc Z  
#anc Z 0.00484 -      
#neo Z < 2e-16 0.00063



## ds anc vs neo Z vs A ##

median(Z_anc_dnds$ds)
# 0.0156
median(Z_neo_dnds$ds)
# 0.0206


## bootstrap confidence intervals of median... ##

boot_ds_anc_Z <- boot(Z_anc_dnds$ds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_ds_anc_Z, type = c("perc"))
plot(boot_ds_anc_Z)
# ( 0.0148,  0.0168 )  

boot_ds_neo_Z <- boot(Z_neo_dnds$ds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_ds_neo_Z, type = c("perc"))
plot(boot_ds_neo_Z)
# ( 0.0193,  0.0213 )  

kruskal.test(ds ~ anc_vs_neo, data = df_filtered)
# p-value = 2.695e-10

pairwise.wilcox.test(df_filtered$ds, df_filtered$anc_vs_neo,
                     p.adjust.method = "bonferroni")

#      A       anc Z  
#anc Z 0.02    -      
#neo Z 1.8e-08 2.0e-09



#---------------------------#
##### neo vs anc on Z1 ######
#---------------------------#


Z1 <- df_filtered[df_filtered$chrN == "Z1", ]
Z1_anc <- Z1[Z1$anc_vs_neo == "anc Z", ]
Z1_neo <- Z1[Z1$anc_vs_neo == "neo Z", ]

median(Z1_neo$dnds)
# 0.0961
median(Z1_neo$dn)
# 0.00175
median(Z1_neo$ds)
# 0.0177


## bootstrap confidence intervals of median... ##

boot_dnds_neo_Z1 <- boot(Z1_neo$dnds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dnds_neo_Z1, type = c("perc"))
plot(boot_dnds_neo_Z1)
# ( 0.0566,  0.1251 ) 

boot_dn_neo_Z1 <- boot(Z1_neo$dn, function(x,i) median(x[i]), R=10000)
boot.ci(boot_dn_neo_Z1, type = c("perc"))
plot(boot_dn_neo_Z1)
# ( 0.0012,  0.0022 )  

boot_ds_neo_Z1 <- boot(Z1_neo$ds, function(x,i) median(x[i]), R=10000)
boot.ci(boot_ds_neo_Z1, type = c("perc"))
plot(boot_ds_neo_Z1)
# ( 0.0157,  0.0202 )  

wilcox.test(dnds ~ anc_vs_neo, data = Z1, paired=FALSE)
# p-value = 0.631

wilcox.test(dn ~ anc_vs_neo, data = Z1, paired=FALSE)
# p-value = 0.8213

wilcox.test(ds ~ anc_vs_neo, data = Z1, paired=FALSE)
# p-value = 0.1134



#-----------------------------#
##### dnds individual chr #####
#-----------------------------#


# test all-vs-all between chromosomes...


kruskal.test(df_filtered$dnds ~ df_filtered$chrN, data = df_filtered)

dunn_res <- dunn_test(df_filtered, dnds ~ chrN, p.adjust.method	= "fdr")

write.table(dunn_res,"~/Desktop/Plots/FASTZ/dunn_test_dnds.txt", sep = "\t")


#-----------------------------------------------------------#
##### Model log(dn/ds) ~ chr type + log(fpkm) + gc4fold #####
#-----------------------------------------------------------#


df_filtered <- df[df$dnds < 999, ]
df_filtered <- subset(df_filtered, dn!=0)

dnds_model <- lm(log(df_filtered$dnds) ~ 
                   df_filtered$anc_vs_neo + 
                   df_filtered$gc4fold +
                   log(df_filtered$FPKM),
                 data = df_filtered)

layout(matrix(c(1,2,3,4),2,2))
plot(dnds_model)
hist(residuals(object = dnds_model))

ncvTest(dnds_model)
#Non-constant Variance Score Test (heteroscedasticity)
#Variance formula: ~ fitted.values 
#Chisquare = 0.4881859, Df = 1, p = 0.48474

summary(dnds_model)

#Call:
#lm(formula = log(df_filtered$dnds) ~ df_filtered$anc_vs_neo + 
#     df_filtered$gc4fold + log(df_filtered$FPKM), data = df_filtered)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-4.0068 -0.6661 -0.0167  0.6228  8.3012 

#Coefficients:
#                               Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)                 -1.556479   0.040115 -38.800  < 2e-16 ***
#  df_filtered$anc_vs_neoanc Z  0.177458   0.061191   2.900  0.00375 ** 
#  df_filtered$anc_vs_neoneo Z  0.228779   0.042441   5.391  7.3e-08 ***
#  df_filtered$gc4fold         -0.806730   0.093928  -8.589  < 2e-16 ***
#  log(df_filtered$FPKM)       -0.002042   0.006324  -0.323  0.74681    
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.9826 on 5806 degrees of freedom
#(417 observations deleted due to missingness)
#Multiple R-squared:  0.01813,	Adjusted R-squared:  0.01745 
#F-statistic:  26.8 on 4 and 5806 DF,  p-value: < 2.2e-16


Anova(dnds_model)

# Anova Table (Type II tests)

# Response: log(df_filtered$dnds)
#                         Sum Sq   Df F value    Pr(>F)    
#  df_filtered$anc_vs_neo   34.0    2 17.6247 2.338e-08 ***
#  df_filtered$gc4fold      71.2    1 73.7673 < 2.2e-16 ***
#  log(df_filtered$FPKM)     0.1    1  0.1042    0.7468    
#  Residuals              5605.4 5806                      
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
