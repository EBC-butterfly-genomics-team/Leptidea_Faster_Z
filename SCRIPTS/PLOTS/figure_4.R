
library(ggplot2)
library(patchwork)
library(introdataviz)

setwd("~/Desktop/Plots/FASTZ")
df <- read.csv("lep_sin_dnds_data.txt", sep = "\t")

label_sizes=20
axis_text_sizes=15
auto_color <- "grey"
fem_col="#ea5545"
male_col="#27aeef"
sign_size=12
line_size=0.75

#-----------------------#
##### calculate DoS #####
#-----------------------#

df$DoS <- (df$Dn/(df$Dn+df$Ds)) - (df$Pn/(df$Pn+df$Ps))

#---------------------------#
##### filter dnds < 999 #####
#---------------------------#

#df_filtered <- df[df$dnds < 999, ]


#-----------------------------#
##### plot Z, anc and neo #####
#-----------------------------#

auto_color <- "grey"
anz_Z_color <- rgb(213/255,117/255,0)

#-------------#
##### DoS #####
#-------------#

all_DoS_plot <- ggplot(df, aes(y=DoS, x=anc_vs_neo, fill=anc_vs_neo)) +
  geom_violin() +
  geom_boxplot(width=0.1, outlier.shape = NA, color="white") +
  coord_cartesian(ylim=c(-1, 1.3)) +
  scale_fill_manual(values = c(auto_color, anz_Z_color, "firebrick4")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_text_sizes),
        axis.title.y = element_text(size=axis_text_sizes),
        axis.text.y = element_text(size=axis_text_sizes)) +
  annotate("text", label = "*", x=1.5, y=1.1, size = sign_size) +
  annotate("text", label = "*", x=2, y=1.25, size = sign_size) +
  annotate("segment", x=1, xend=2, y=1.1, yend=1.1) +
  annotate("segment", x=1, xend=3, y=1.25, yend=1.25) 


#--------------------#
##### Tajima's D #####
#--------------------#

kruskal.test(TajD ~ anc_vs_neo, data = df)
# p-value < 2.2e-16

pairwise.wilcox.test(df$TajD, df$anc_vs_neo,
                     p.adjust.method = "bonferroni")
#      A       anc Z
#anc Z 2.5e-14 -    
#neo Z 2.5e-09 0.018

aggregate(TajD ~ anc_vs_neo, data = df, FUN = median)

#anc_vs_neo       TajD
#         A -0.3270450
#     anc Z -0.6848615
#     neo Z -0.5277810

all_TajD_plot <- ggplot(df, aes(y=TajD, x=anc_vs_neo, fill=anc_vs_neo)) +
  geom_violin() +
  geom_boxplot(width=0.1, outlier.shape = NA, color="white") +
  coord_cartesian(ylim=c(-2.5, 4.5)) +
  scale_fill_manual(values = c(auto_color, anz_Z_color, "firebrick4")) +
  scale_y_continuous(name = "Tajima's D") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_text_sizes),
        axis.title.y = element_text(size=axis_text_sizes),
        axis.text.y = element_text(size=axis_text_sizes)) +
  annotate("text", label = "*", x=1.5, y=3.5, size = sign_size) +
  annotate("text", label = "*", x=2, y=4, size = sign_size) +
  annotate("text", label = "*", x=2.5, y=4.5, size = sign_size) +
  annotate("segment", x=1, xend=2, y=3.5, yend=3.5) +
  annotate("segment", x=1, xend=3, y=4, yend=4) +
  annotate("segment", x=2, xend=3, y=4.5, yend=4.5) 

all_TajD_plot

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

## DoS A: MBG < FBG & UBG




## All Z ##

kruskal.test(DoS ~ sex_bias, data = all_Z_SBG)
#p-value = 0.007844

pairwise.wilcox.test(all_Z_SBG$DoS, all_Z_SBG$sex_bias,
                     p.adjust.method = "bonferroni")

#    FBG    MBG   
#MBG 0.3196 -     
#UBG 1.0000 0.0056

all_Z_DoS_median <- aggregate(DoS ~ sex_bias, data = all_Z_SBG, FUN = median)

## DoS all Z: MBG < UBG



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

## DoS Z1: MBG < UBG




## Z2 ##

kruskal.test(DoS ~ sex_bias, data = Z2)
# p-value = 0.03515

pairwise.wilcox.test(Z2$DoS, Z2$sex_bias,
                     p.adjust.method = "bonferroni")

# No difference...

Z2_DoS_median <- aggregate(DoS ~ sex_bias, data = Z2, FUN = median)




## Z3 ##

kruskal.test(DoS ~ sex_bias, data = Z3)
# p-value = 0.2637

# No difference...

Z3_DoS_median <- aggregate(DoS ~ sex_bias, data = Z3, FUN = median)




Autosomes_DoS <- ggplot(Autosomes, aes(fill=sex_bias, x=DoS)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-1,1), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = Autosomes_DoS_median, aes(xintercept=DoS, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  scale_y_continuous(limits = c(0, 0.006), breaks = c(0.000, 0.003, 0.006)) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  annotate("text", label = "*", x=0.5, y=0.004, size = sign_size, color=male_col) +
  annotate("text", label = "*", x=0.60, y=0.004, size = sign_size, color="wheat") +
  
  annotate("text", label = "*", x=0.5, y=0.0025, size = sign_size, color=male_col) +
  annotate("text", label = "*", x=0.60, y=0.0025, size = sign_size, color=fem_col) +
  xlab("")

all_Z_DoS <- ggplot(all_Z_SBG, aes(fill=sex_bias, x=DoS)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-1,1), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = all_Z_DoS_median, aes(xintercept=DoS, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  scale_y_continuous(limits = c(0, 0.006), breaks = c(0.000, 0.003, 0.006)) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  annotate("text", label = "*", x=0.5, y=0.004, size = sign_size, color=male_col) +
  annotate("text", label = "*", x=0.60, y=0.004, size = sign_size, color="wheat") +
  xlab("")

Z1_DoS <- ggplot(Z1, aes(fill=sex_bias, x=DoS)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-1,1), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = Z1_DoS_median, aes(xintercept=DoS, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  scale_y_continuous(limits = c(0, 0.006), breaks = c(0.000, 0.003, 0.006)) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  annotate("text", label = "*", x=0.5, y=0.004, size = sign_size, color=male_col) +
  annotate("text", label = "*", x=0.60, y=0.004, size = sign_size, color="wheat") +
  xlab("")

Z2_DoS <- ggplot(Z2, aes(fill=sex_bias, x=DoS)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-1,1), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = Z2_DoS_median, aes(xintercept=DoS, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  scale_y_continuous(limits = c(0, 0.006), breaks = c(0.000, 0.003, 0.006)) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlab("")

Z3_DoS <- ggplot(Z3, aes(fill=sex_bias, x=DoS)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-1,1), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = Z3_DoS_median, aes(xintercept=DoS, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  scale_y_continuous(limits = c(0, 0.006), breaks = c(0.000, 0.003, 0.006)) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title.x = element_text(size=axis_text_sizes),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlab("DoS")


#-------------------------------------------------------#
##### Tajima's D: contrast SBG on A and different Z #####
#-------------------------------------------------------#


## A ##

kruskal.test(TajD ~ sex_bias, data = Autosomes)
# p-value = 0.6894

# No difference...

Autosomes_TajD_median <- aggregate(TajD ~ sex_bias, data = Autosomes, FUN = median)


## All Z ##

kruskal.test(TajD ~ sex_bias, data = all_Z_SBG)
#p-value = 0.006544

pairwise.wilcox.test(all_Z_SBG$TajD, all_Z_SBG$sex_bias,
                     p.adjust.method = "bonferroni")

#    FBG   MBG  
#MBG 0.029 -    
#UBG 0.992 0.015

all_Z_TajD_median <- aggregate(TajD ~ sex_bias, data = all_Z_SBG, FUN = median)

## TajD all Z: FBG < MBG, MBG > UBG


## Z1 ##

kruskal.test(TajD ~ sex_bias, data = Z1)
# p-value = 0.1276

# No difference...

Z1_TajD_median <- aggregate(TajD ~ sex_bias, data = Z1, FUN = median)



## Z2 ##

kruskal.test(TajD ~ sex_bias, data = Z2)
# p-value = 0.05849

Z2_TajD_median <- aggregate(TajD ~ sex_bias, data = Z2, FUN = median)




## Z3 ##

kruskal.test(TajD ~ sex_bias, data = Z3)
# p-value = 0.0569

# No difference...

Z3_TajD_median <- aggregate(TajD ~ sex_bias, data = Z3, FUN = median)




Autosomes_TajD <- ggplot(Autosomes, aes(fill=sex_bias, x=TajD)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-2.5,2.5), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = Autosomes_TajD_median, aes(xintercept=TajD, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("")

all_Z_TajD <- ggplot(all_Z_SBG, aes(fill=sex_bias, x=TajD)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-2.5,2.5), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = all_Z_TajD_median, aes(xintercept=TajD, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", label = "*", x=1.75, y=0.004, size = sign_size, color=male_col) +
  annotate("text", label = "*", x=1.5, y=0.004, size = sign_size, color="wheat") +
  
  annotate("text", label = "*", x=1.75, y=0.0025, size = sign_size, color=male_col) +
  annotate("text", label = "*", x=1.5, y=0.0025, size = sign_size, color=fem_col) +
  xlab("")

Z1_TajD <- ggplot(Z1, aes(fill=sex_bias, x=TajD)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-2.5,2.5), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = Z1_TajD_median, aes(xintercept=TajD, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("")

Z2_TajD <- ggplot(Z2, aes(fill=sex_bias, x=TajD)) +
  geom_density(aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-2.5,2.5), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = Z2_TajD_median, aes(xintercept=TajD, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("")

Z3_TajD <- ggplot(Z3, aes(fill=sex_bias, x=TajD)) +
  geom_density( aes(y = ..count../sum(..count..)),  position="identity", alpha=.5, color = NA) +
  coord_cartesian(xlim = c(-2.5,2.5), ylim = c(0, 0.006), expand = FALSE) +
  scale_fill_manual("Sex bias", values = c(fem_col, male_col, "wheat")) +
  geom_vline(data = Z3_TajD_median, aes(xintercept=TajD, color=sex_bias, linetype=sex_bias),
             size=line_size, alpha = 0.95) +
  scale_linetype_manual("Sex bias", values = c("longdash", "dashed", "dotted")) +
  scale_color_manual("Sex bias", values = c(fem_col, male_col, "wheat4")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=axis_text_sizes),
        axis.text = element_text(size=axis_text_sizes),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=axis_text_sizes),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("Tajima's D")


Auto_label <- ggplot() +
  annotate("text", label = "Auto", x=1, y=1, size = 6) +
  theme_void()
  
all_Z_label <- ggplot() +
  annotate("text", label = "All Z", x=1, y=1, size = 6) +
  theme_void()

Z1_label <- ggplot() +
  annotate("text", label = "Z1", x=1, y=1, size = 6) +
  theme_void()

Z2_label <- ggplot() +
  annotate("text", label = "Z2", x=1, y=1, size = 6) +
  theme_void()

Z3_label <- ggplot() +
  annotate("text", label = "Z3", x=1, y=1, size = 6) +
  theme_void()




##### Analyse per Z chromosome #####



#-------------#
##### DoS #####
#-------------#

# group all autosomes...
df$chrN[df$Z_vs_A != "Z"] <- "A"

# separate by age... (SKIP!)

Lsin_anc_Z_start <- 10536855
Lsin_neo_17_end <- 10453464

anz_Z_color <- rgb(213/255,117/255,0)
neo_17_color <- rgb(102/255,141/255,60/255)

Lsin_neo_11_end <- 10209824
Lsin_neo_7_start <- 11321879

neo_11_color <- rgb(185/255,156/255,107/255)
neo_7_color <- rgb(129/255,108/255,91/255)
neo_24_color <- rgb(228/255,153/255,105/255)

neo_19_color <- rgb(131/255,146/255,159/255)
neo_15_color <- rgb(219/255,202/255,105/255)
neo_8_color <- rgb(160/255,161/255,140/255)
neo_26_color <- rgb(157/255,151/255,84/255)
neo_28_color <- rgb(183/255,166/255,173/255)

#df <- df %>% 
#  mutate(age = case_when(chrN == "Z1" & gene_position > Lsin_anc_Z_start ~ "Z1a",
#                             chrN == "Z1" & gene_position < Lsin_neo_17_end ~ "Z1b",
#                             chrN == "Z2" & gene_position > Lsin_neo_7_start ~ "Z2a",
#                             chrN == "Z2" & gene_position < Lsin_neo_11_end ~ "Z2b",
#                             chrN == "Z3" ~ "Z3",
#                             chrN != "Chr_1" & chrN != "Chr_6" & chrN != "Chr_18" ~ "A"))

kruskal.test(DoS ~ chrN, data = df)
# Kruskal-Wallis chi-squared = 44.461, df = 3, p-value = 1.204e-09
wilcox_test(df, DoS ~ chrN, p.adjust.method = "bonferroni")
#.y.   group1 group2    n1    n2 statistic            p       p.adj p.adj.signif
#* <chr> <chr>  <chr>  <int> <int>     <dbl>        <dbl>       <dbl> <chr>       
#  1 DoS   A      Z1      8784   646  2557033  0.000027     0.000162    ***         
#  2 DoS   A      Z2      8784   438  1633138. 0.0000000904 0.000000542 ****        
#  3 DoS   A      Z3      8784   265  1181898. 0.667        1           ns          
#  4 DoS   Z1     Z2       646   438   134034. 0.141        0.846       ns          
#  5 DoS   Z1     Z3       646   265    95611  0.005        0.033       *           
#  6 DoS   Z2     Z3       438   265    67796. 0.000183     0.001       **          

# A lower than Z1 & Z2 (but not Z3)
# Z1 & Z2 not different
# Z3 lower than Z1 & Z2


all_DoS_plot <- ggplot(df, aes(y=DoS, x=chrN, fill=chrN)) +
  geom_violin(color = NA) +
  geom_boxplot(width=0.1, outlier.shape = NA, color="white") +
  coord_cartesian(ylim=c(-1, 1.5)) +
  scale_fill_manual(values = c(auto_color, anz_Z_color, neo_11_color, neo_15_color)) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_text_sizes),
        axis.title.y = element_text(size=axis_text_sizes),
        axis.text.y = element_text(size=axis_text_sizes)) +
  annotate("text", label = "*", x=1.5, y=1.1, size = sign_size) +
  annotate("text", label = "*", x=2, y=1.2, size = sign_size) +
  annotate("text", label = "*", x=3, y=1.3, size = sign_size) +
  annotate("text", label = "*", x=3.5, y=1.4, size = sign_size) +
  annotate("segment", x=1, xend=2, y=1.1, yend=1.1) +
  annotate("segment", x=1, xend=3, y=1.2, yend=1.2) +
  annotate("segment", x=2, xend=4, y=1.3, yend=1.3) +
  annotate("segment", x=3, xend=4, y=1.4, yend=1.4) 

all_DoS_plot

#--------------------#
##### Tajima's D #####
#--------------------#

kruskal.test(TajD ~ chrN, data = df)
#Kruskal-Wallis chi-squared = 95.725, df = 3, p-value < 2.2e-16
wilcox_test(df, TajD ~ chrN, p.adjust.method = "bonferroni")
#.y.   group1 group2    n1    n2 statistic        p    p.adj p.adj.signif
#* <chr> <chr>  <chr>  <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
#  1 TajD  A      Z1      9709   700  4070551  2   e-18 1.20e-17 ****        
#  2 TajD  A      Z2      9709   462  2427792. 3   e- 3 1.6 e- 2 *           
#  3 TajD  A      Z3      9709   289  1598290. 5.35e- 5 3.21e- 4 ***         
#  4 TajD  Z1     Z2       700   462   143770. 1   e- 3 8   e- 3 **          
#  5 TajD  Z1     Z3       700   289    94852. 1.23e- 1 7.38e- 1 ns          
#  6 TajD  Z2     Z3       462   289    70216  2.32e- 1 1   e+ 0 ns          

# A higher than all
# Z1 lower than Z2

all_TajD_plot <- ggplot(df, aes(y=TajD, x=chrN, fill=chrN)) +
  geom_violin(color = NA) +
  geom_boxplot(width=0.1, outlier.shape = NA, color="white") +
  coord_cartesian(ylim=c(-2.5, 4.5)) +
  scale_fill_manual(values = c(auto_color, anz_Z_color, neo_11_color, neo_15_color)) +
  scale_y_continuous(name = "Tajima's D") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=axis_text_sizes),
        axis.title.y = element_text(size=axis_text_sizes),
        axis.text.y = element_text(size=axis_text_sizes)) +
  annotate("text", label = "*", x=1.5, y=3.5, size = sign_size) +
  annotate("text", label = "*", x=2, y=3.8, size = sign_size) +
  annotate("text", label = "*", x=2.5, y=4.1, size = sign_size) +
  annotate("text", label = "*", x=2.5, y=4.4, size = sign_size) +
  annotate("segment", x=1, xend=2, y=3.5, yend=3.5) +
  annotate("segment", x=1, xend=3, y=3.8, yend=3.8) +
  annotate("segment", x=1, xend=4, y=4.1, yend=4.1) +
  annotate("segment", x=2, xend=3, y=4.4, yend=4.4) 

all_TajD_plot


all_plots <- ((all_DoS_plot / all_TajD_plot) |
                (Auto_label / all_Z_label / Z1_label / Z2_label / Z3_label | 
                   Autosomes_DoS / all_Z_DoS / Z1_DoS / Z2_DoS / Z3_DoS |
                   Autosomes_TajD / all_Z_TajD / Z1_TajD / Z2_TajD / Z3_TajD) + plot_layout(guides = "collect") + 
                plot_layout(widths = c(0.5, 1.5, 1.5)) & theme(legend.position = 'bottom')) + plot_layout(widths = c(1, 2))

all_plots + plot_annotation(tag_levels = list(c("A", "B", "C", "", "", "", ""))) & 
  theme(plot.tag = element_text(size = 20))

ggsave(filename = ("figure_4.jpeg"), width = 15, height = 12, units = "in", dpi = 600, limitsize = F)
