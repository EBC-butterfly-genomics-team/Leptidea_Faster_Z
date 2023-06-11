
## Direction of selection and Tajima's D statistics ##

# This corresponds to Figure 4, Table 3, Supplementary tables 7, 8 & 9


setwd("~/Desktop/Plots/FASTZ")
df <- read.csv("lep_sin_dnds_data.txt", sep = "\t")


#-----------------------#
##### calculate DoS #####
#-----------------------#

df$DoS <- (df$Dn/(df$Dn+df$Ds)) - (df$Pn/(df$Pn+df$Ps))

# median...

aggregate(DoS ~ anc_vs_neo, data = df, FUN = median)

#anc_vs_neo         DoS
#A              -0.10391726
#anc Z          -0.04618474
#neo Z          -0.05528042


# Test difference between chromosome categories...

kruskal.test(DoS ~ anc_vs_neo, data = df)
# p-value = 2.284e-07

pairwise.wilcox.test(df$DoS, df$anc_vs_neo,
                     p.adjust.method = "bonferroni")
#      A       anc Z 
#anc Z 0.0014  -     
#neo Z 1.9e-05 1.0000

# Both anc and neo Z higher than A



#--------------------#
##### Tajima's D #####
#--------------------#

aggregate(TajD ~ anc_vs_neo, data = df, FUN = median)

#anc_vs_neo       TajD
#         A -0.3270450
#     anc Z -0.6848615
#     neo Z -0.5277810


# Test difference between chromosome categories...

kruskal.test(TajD ~ anc_vs_neo, data = df)
# p-value < 2.2e-16

pairwise.wilcox.test(df$TajD, df$anc_vs_neo,
                     p.adjust.method = "bonferroni")
#      A       anc Z
#anc Z 2.5e-14 -    
#neo Z 2.5e-09 0.018




#--------------#
##### SBGs #####
#--------------#


Z1 <- df[df$chrN == "Z1" & df$sex_bias != "", ]
Z2 <- df[df$chrN == "Z2" & df$sex_bias != "", ]
Z3 <- df[df$chrN == "Z3" & df$sex_bias != "", ]
Autosomes <- df[df$Z_vs_A == "A" & df$sex_bias != "", ]
all_Z_SBG <- df[df$Z_vs_A == "Z" & df$sex_bias != "", ]



#------------------------------------------------#
##### DoS: contrast SBG on A and different Z #####
#------------------------------------------------#


## A ##

kruskal.test(DoS ~ sex_bias, data = Autosomes)
# p-value = 0.009562

# difference...

pairwise.wilcox.test(Autosomes$DoS, Autosomes$sex_bias,
                     p.adjust.method = "bonferroni")

#      FBG   MBG  
#MBG 0.020 -    
#UBG 1.000 0.016

Autosomes_DoS_median <- aggregate(DoS ~ sex_bias, data = Autosomes, FUN = median)

boot_DoS_FBG_A <- boot(Autosomes$DoS[Autosomes$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_FBG_A, type = c("perc"))
plot(boot_DoS_FBG_A)
# (-0.0992, -0.0546 )  

boot_DoS_MBG_A <- boot(Autosomes$DoS[Autosomes$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_MBG_A, type = c("perc"))
plot(boot_DoS_MBG_A)
# (-0.1467, -0.1154 )  

boot_DoS_UBG_A <- boot(Autosomes$DoS[Autosomes$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_UBG_A, type = c("perc"))
plot(boot_DoS_UBG_A)
# (-0.1111, -0.0924 )  


## All Z ##

kruskal.test(DoS ~ sex_bias, data = all_Z_SBG)
#p-value = 0.007844

pairwise.wilcox.test(all_Z_SBG$DoS, all_Z_SBG$sex_bias,
                     p.adjust.method = "bonferroni")

#    FBG    MBG   
#MBG 0.3196 -     
#UBG 1.0000 0.0056

all_Z_DoS_median <- aggregate(DoS ~ sex_bias, data = all_Z_SBG, FUN = median)

boot_DoS_FBG_Z <- boot(all_Z_SBG$DoS[all_Z_SBG$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_FBG_Z, type = c("perc"))
plot(boot_DoS_FBG_Z)
# (-0.1186,  0.0000 ) 

boot_DoS_MBG_Z <- boot(all_Z_SBG$DoS[all_Z_SBG$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_MBG_Z, type = c("perc"))
plot(boot_DoS_MBG_Z)
# (-0.1469, -0.0649 )  

boot_DoS_UBG_Z <- boot(all_Z_SBG$DoS[all_Z_SBG$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_UBG_Z, type = c("perc"))
plot(boot_DoS_UBG_Z)
# (-0.0704, -0.0112 )



## Z1 ##

kruskal.test(DoS ~ sex_bias, data = Z1)
# p-value = 0.008895

# Overall difference...

pairwise.wilcox.test(Z1$DoS, Z1$sex_bias,
                     p.adjust.method = "bonferroni")

#    FBG   MBG  
#MBG 0.088 -    
#UBG 0.892 0.018

Z1_DoS_median <- aggregate(DoS ~ sex_bias, data = Z1, FUN = median)

boot_DoS_FBG_Z1 <- boot(Z1$DoS[Z1$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_FBG_Z1, type = c("perc"))
plot(boot_DoS_FBG_Z1)
# (-0.1096,  0.1446 )  

boot_DoS_MBG_Z1 <- boot(Z1$DoS[Z1$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_MBG_Z1, type = c("perc"))
plot(boot_DoS_MBG_Z1)
# (-0.1769, -0.0311 )  

boot_DoS_UBG_Z1 <- boot(Z1$DoS[Z1$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_UBG_Z1, type = c("perc"))
plot(boot_DoS_UBG_Z1)
# (-0.0695,  0.0000 )  



## Z2 ##

kruskal.test(DoS ~ sex_bias, data = Z2)
# p-value = 0.03515

pairwise.wilcox.test(Z2$DoS, Z2$sex_bias,
                     p.adjust.method = "bonferroni")

#    FBG   MBG  
#MBG 0.080 -    
#UBG 1.000 0.078

# No difference...

Z2_DoS_median <- aggregate(DoS ~ sex_bias, data = Z2, FUN = median)

boot_DoS_FBG_Z2 <- boot(Z2$DoS[Z2$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_FBG_Z2, type = c("perc"))
plot(boot_DoS_FBG_Z2)
# (-0.0437,  0.0857 )  

boot_DoS_MBG_Z2 <- boot(Z2$DoS[Z2$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_MBG_Z2, type = c("perc"))
plot(boot_DoS_MBG_Z2)
# (-0.1723, -0.0208 )  

boot_DoS_UBG_Z2 <- boot(Z2$DoS[Z2$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_UBG_Z2, type = c("perc"))
plot(boot_DoS_UBG_Z2)
# (-0.0833,  0.0000 )  




## Z3 ##

kruskal.test(DoS ~ sex_bias, data = Z3)
# p-value = 0.2637

# No difference...

Z3_DoS_median <- aggregate(DoS ~ sex_bias, data = Z3, FUN = median)

boot_DoS_FBG_Z3 <- boot(Z3$DoS[Z3$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_FBG_Z3, type = c("perc"))
plot(boot_DoS_FBG_Z3)
# (-0.0437,  0.0857 )  

boot_DoS_MBG_Z3 <- boot(Z3$DoS[Z3$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_MBG_Z3, type = c("perc"))
plot(boot_DoS_MBG_Z3)
# (-0.15,  0.00 )  

boot_DoS_UBG_Z3 <- boot(Z3$DoS[Z3$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_DoS_UBG_Z3, type = c("perc"))
plot(boot_DoS_UBG_Z3)
# (-0.1429, -0.0238 )  



#-------------------------------------------------------#
##### Tajima's D: contrast SBG on A and different Z #####
#-------------------------------------------------------#


## A ##

kruskal.test(TajD ~ sex_bias, data = Autosomes)
# p-value = 0.6894

# No difference...

Autosomes_TajD_median <- aggregate(TajD ~ sex_bias, data = Autosomes, FUN = median)

boot_TajD_FBG_A <- boot(Autosomes$TajD[Autosomes$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_FBG_A, type = c("perc"))
plot(boot_TajD_FBG_A)
# (-0.3934, -0.2574 )  

boot_TajD_MBG_A <- boot(Autosomes$TajD[Autosomes$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_MBG_A, type = c("perc"))
plot(boot_TajD_MBG_A)
# (-0.4067, -0.2851 )  

boot_TajD_UBG_A <- boot(Autosomes$TajD[Autosomes$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_UBG_A, type = c("perc"))
plot(boot_TajD_UBG_A)
# (-0.1111, -0.0924 )  




## All Z ##

kruskal.test(TajD ~ sex_bias, data = all_Z_SBG)
#p-value = 0.006544

pairwise.wilcox.test(all_Z_SBG$TajD, all_Z_SBG$sex_bias,
                     p.adjust.method = "bonferroni")

#    FBG   MBG  
#MBG 0.029 -    
#UBG 0.992 0.015

all_Z_TajD_median <- aggregate(TajD ~ sex_bias, data = all_Z_SBG, FUN = median)

boot_TajD_FBG_Z <- boot(all_Z_SBG$TajD[all_Z_SBG$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_FBG_Z, type = c("perc"))
plot(boot_TajD_FBG_Z)
# (-0.9666, -0.5276 )  

boot_TajD_MBG_Z <- boot(all_Z_SBG$TajD[all_Z_SBG$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_MBG_Z, type = c("perc"))
plot(boot_TajD_MBG_Z)
# (-0.5541, -0.2870 )  

boot_TajD_UBG_Z <- boot(all_Z_SBG$TajD[all_Z_SBG$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_UBG_Z, type = c("perc"))
plot(boot_TajD_UBG_Z)
# (-0.7480, -0.5462 )  


## Z1 ##

kruskal.test(TajD ~ sex_bias, data = Z1)
# p-value = 0.1276

# No difference...

Z1_TajD_median <- aggregate(TajD ~ sex_bias, data = Z1, FUN = median)

boot_TajD_FBG_Z1 <- boot(Z1$TajD[Z1$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_FBG_Z1, type = c("perc"))
plot(boot_TajD_FBG_Z1)
# (-1.2692, -0.5905 )  

boot_TajD_MBG_Z1 <- boot(Z1$TajD[Z1$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_MBG_Z1, type = c("perc"))
plot(boot_TajD_MBG_Z1)
# (-0.7747, -0.3775 )  

boot_TajD_UBG_Z1 <- boot(Z1$TajD[Z1$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_UBG_Z1, type = c("perc"))
plot(boot_TajD_UBG_Z1)
# (-0.8133, -0.5916 )  


## Z2 ##

kruskal.test(TajD ~ sex_bias, data = Z2)
# p-value = 0.05849

Z2_TajD_median <- aggregate(TajD ~ sex_bias, data = Z2, FUN = median)

boot_TajD_FBG_Z2 <- boot(Z2$TajD[Z2$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_FBG_Z2, type = c("perc"))
plot(boot_TajD_FBG_Z2)
# (-1.1407, -0.3372 )  

boot_TajD_MBG_Z2 <- boot(Z2$TajD[Z2$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_MBG_Z2, type = c("perc"))
plot(boot_TajD_MBG_Z2)
# (-0.4850,  0.0702 )  

boot_TajD_UBG_Z2 <- boot(Z2$TajD[Z2$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_UBG_Z2, type = c("perc"))
plot(boot_TajD_UBG_Z2)
# (-0.6674, -0.3565 )  


## Z3 ##

kruskal.test(TajD ~ sex_bias, data = Z3)
# p-value = 0.0569

# No difference...

Z3_TajD_median <- aggregate(TajD ~ sex_bias, data = Z3, FUN = median)

boot_TajD_FBG_Z3 <- boot(Z3$TajD[Z3$sex_bias == "FBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_FBG_Z3, type = c("perc"))
plot(boot_TajD_FBG_Z3)
# (-0.9687, -0.3248 )  

boot_TajD_MBG_Z3 <- boot(Z3$TajD[Z3$sex_bias == "MBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_MBG_Z3, type = c("perc"))
plot(boot_TajD_MBG_Z3)
# (-0.5656, -0.0297 )  

boot_TajD_UBG_Z3 <- boot(Z3$TajD[Z3$sex_bias == "UBG"], function(x,i) median(x[i], na.rm = TRUE), R=10000)
boot.ci(boot_TajD_UBG_Z3, type = c("perc"))
plot(boot_TajD_UBG_Z3)
# (-0.9317, -0.5117 )  




