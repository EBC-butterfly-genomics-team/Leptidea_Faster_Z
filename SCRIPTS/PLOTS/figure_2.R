
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
library(see)


setwd("~/Desktop/Plots/FASTZ")


#----------------------#
##### read in data #####
#----------------------#


df <- read.csv("lep_sin_dnds_data.txt", sep = "\t")
dunn_sign_table <- read.csv("sign_dunn_test_count.txt", sep = "\t")

# plot variables...
label_sizes=20
axis_text_sizes=12
sign_size=12
auto_color <- "grey"
anz_Z_color <- rgb(213/255,117/255,0)
neo_17_color <- rgb(102/255,141/255,60/255)


#------------------------------#
##### Test unfiltered data #####
#------------------------------#

# test but do not plot...

# dnds
dnds_plot <- ggplot(data=df, aes(y=dnds, x=Z_vs_A, fill=Z_vs_A)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  scale_y_continuous(name = "dn/ds") +
  coord_cartesian(ylim=c(0, 0.6)) +
  scale_fill_manual(values = c("skyblue3", "firebrick2")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank())

wilcox.test(dnds ~ Z_vs_A, data = df, paired=FALSE)
# p-value < 2.2e-16


# dn
dn_plot <- ggplot(data=df, aes(y=dn, x=Z_vs_A, fill=Z_vs_A)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  coord_cartesian(ylim=c(0, 0.017)) +
  scale_fill_manual(values = c("skyblue3", "firebrick2")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank())

wilcox.test(dn ~ Z_vs_A, data = df, paired=FALSE)
# p-value < 2.2e-16


# ds
ds_plot <- ggplot(data=df, aes(y=ds, x=Z_vs_A, fill=Z_vs_A)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  coord_cartesian(ylim=c(0, 0.08)) +
  scale_fill_manual(values = c("skyblue3", "firebrick2")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank())

wilcox.test(ds ~ Z_vs_A, data = df, paired=FALSE)
# p-value = 0.1241

ggarrange(dnds_plot, dn_plot, ds_plot, ncol = 3, nrow = 1)



#---------------------------#
##### Filter high dN/dS #####
#---------------------------#


df_filtered <- df[df$dnds < 999, ]


#----------------------#
##### dN/dS A vs Z #####
#----------------------#


# test but do not plot...

# dn/ds #
wilcox.test(dnds ~ Z_vs_A, data = df_filtered, paired=FALSE, conf.int=TRUE)
# p-value < 2.2e-16

dn_ds_sign_pos <- boxplot(df_filtered$dnds)$stats[5, ]

dnds_filtered_plot <- ggplot(data=df_filtered, aes(y=dnds, x=Z_vs_A, fill=Z_vs_A)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  coord_cartesian(ylim=c(0, 0.7)) +
  scale_y_continuous(name = "dN/dS") +
  scale_fill_manual(values = c("grey", "firebrick2")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "grey10"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=label_sizes),
        axis.title.y = element_text(size=label_sizes),
        axis.text.y = element_text(size=axis_text_sizes)) +
  annotate("text", label = "*", x=1.5, y=dn_ds_sign_pos, size = sign_size) 


# dn #
wilcox.test(dn ~ Z_vs_A, data = df_filtered, paired=FALSE)
# p-value < 2.2e-16

dn_sign_pos <- boxplot(df_filtered$dn)$stats[5, ]

dn_filtered_plot <- ggplot(data=df_filtered, aes(y=dn, x=Z_vs_A, fill=Z_vs_A)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  coord_cartesian(ylim=c(0, 0.015)) +
  scale_y_continuous(name = "dN") +
  scale_fill_manual(values = c("grey", "firebrick2")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "grey10"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=label_sizes),
        axis.title.y = element_text(size=label_sizes),
        axis.text.y = element_text(size=axis_text_sizes)) +
  annotate("text", label = "*", x=1.5, y=dn_sign_pos, size = sign_size)


# ds #
wilcox.test(ds ~ Z_vs_A, data = df_filtered, paired=FALSE)
# p-value = 0.001892

ds_sign_pos <- boxplot(df_filtered$ds)$stats[5, ]

ds_filtered_plot <- ggplot(data=df_filtered, aes(y=ds, x=Z_vs_A, fill=Z_vs_A)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  coord_cartesian(ylim=c(0, 0.07)) +
  scale_y_continuous(name = "dS") +
  scale_fill_manual(values = c("grey", "firebrick2")) +
  theme() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "grey10"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=label_sizes),
        axis.title.y = element_text(size=label_sizes),
        axis.text.y = element_text(size=axis_text_sizes)) +
  annotate("text", label = "*", x=1.5, y=ds_sign_pos, size = sign_size)

ggarrange(dnds_filtered_plot, dn_filtered_plot, ds_filtered_plot,
          ncol = 3, nrow = 1)



#-------------------#
##### Figure 2A #####
#-------------------#


# dN/dS anc vs neo #

# dn/ds #
kruskal.test(dnds ~ anc_vs_neo, data = df_filtered)
# p-value < 2.2e-16

pairwise.wilcox.test(df_filtered$dnds, df_filtered$anc_vs_neo,
                     p.adjust.method = "bonferroni")
#      A       anc_Z
#anc_Z 2.1e-05 -    
#neo_Z 2.9e-16 0.71 

Z_dnds <- df_filtered[df_filtered$Z_vs_A == "Z", ]

dn_ds_A_ancZ_sign_pos <- boxplot(Z_dnds$dnds)$stats[5, ]

dnds_anc_neo_filtered_plot <- ggplot(data=df_filtered, aes(y=dnds, x=anc_vs_neo, fill=anc_vs_neo)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  scale_y_continuous(name = "dN/dS") +
  coord_cartesian(ylim=c(0, 0.9)) +
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
  annotate("text", label = "*", x=1.5, y=dn_ds_A_ancZ_sign_pos*1.15, size = sign_size) +
  annotate("text", label = "*", x=2, y=dn_ds_A_ancZ_sign_pos*1.25, size = sign_size) +
  annotate("segment", x=1, xend=2, y=dn_ds_A_ancZ_sign_pos*1.15, yend=dn_ds_A_ancZ_sign_pos*1.15) +
  annotate("segment", x=1, xend=3, y=dn_ds_A_ancZ_sign_pos*1.25, yend=dn_ds_A_ancZ_sign_pos*1.25) 


# dn #
kruskal.test(dn ~ anc_vs_neo, data = df_filtered)
# p-value < 2.2e-16

pairwise.wilcox.test(df_filtered$dn, df_filtered$anc_vs_neo,
                     p.adjust.method = "bonferroni")
#      A       anc_Z  
#anc_Z 0.00484 -      
#neo_Z < 2e-16 0.00063

dn_A_ancZ_sign_pos <- boxplot(Z_dnds$dn)$stats[5, ]

dn_anc_neo_filtered_plot <- ggplot(data=df_filtered, aes(y=dn, x=anc_vs_neo, fill=anc_vs_neo)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  scale_y_continuous(name = "dN") +
  coord_cartesian(ylim=c(0, 0.018)) +
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
  annotate("text", label = "*", x=2.5, y=dn_A_ancZ_sign_pos*1.25, size = sign_size) +
  annotate("text", label = "*", x=2, y=dn_A_ancZ_sign_pos*1.15, size = sign_size) +
  annotate("text", label = "*", x=1.5, y=dn_A_ancZ_sign_pos*1.05, size = sign_size) +
  annotate("segment", x=2, xend=3, y=dn_A_ancZ_sign_pos*1.25, yend=dn_A_ancZ_sign_pos*1.25) +
  annotate("segment", x=1, xend=3, y=dn_A_ancZ_sign_pos*1.15, yend=dn_A_ancZ_sign_pos*1.15) +
  annotate("segment", x=1, xend=2, y=dn_A_ancZ_sign_pos*1.05, yend=dn_A_ancZ_sign_pos*1.05)


# ds #
kruskal.test(ds ~ anc_vs_neo, data = df_filtered)
# p-value = 2.695e-10

pairwise.wilcox.test(df_filtered$ds, df_filtered$anc_vs_neo,
                     p.adjust.method = "bonferroni")
#      A       anc_Z  
#anc_Z 0.02    -      
#neo_Z 1.8e-08 2.0e-09

ds_A_ancZ_sign_pos <- boxplot(Z_dnds$ds)$stats[5, ]

ds_anc_neo_filtered_plot <- ggplot(data=df_filtered, aes(y=ds, x=anc_vs_neo, fill=anc_vs_neo)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  scale_y_continuous(name = "dS") +
  coord_cartesian(ylim=c(0, 0.085)) +
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
  annotate("text", label = "*", x=2.5, y=ds_A_ancZ_sign_pos*1.35, size = sign_size) +
  annotate("text", label = "*", x=2, y=ds_A_ancZ_sign_pos*1.25, size = sign_size) +
  annotate("text", label = "*", x=1.5, y=ds_A_ancZ_sign_pos*1.15, size = sign_size) +
  annotate("segment", x=2, xend=3, y=ds_A_ancZ_sign_pos*1.35, yend=ds_A_ancZ_sign_pos*1.35) +
  annotate("segment", x=1, xend=2, y=ds_A_ancZ_sign_pos*1.15, yend=ds_A_ancZ_sign_pos*1.15) +
  annotate("segment", x=1, xend=3, y=ds_A_ancZ_sign_pos*1.25, yend=ds_A_ancZ_sign_pos*1.25)



#-------------------#
##### Figure 2B #####
#-------------------#


# neo vs anc on Z1 #

Z1 <- df_filtered[df_filtered$chrN == "Chr 1", ]

# dnds
wilcox.test(dnds ~ anc_vs_neo, data = Z1, paired=FALSE, conf.int=TRUE)
# p-value = 0.631

Z1_dnds_ans_vs_neo <- ggplot(data=Z1, aes(y=dnds, x=anc_vs_neo, fill=anc_vs_neo)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  coord_cartesian(ylim=c(0, 0.75)) +
  scale_y_continuous(name = "dN/dS") +
  scale_x_discrete(labels = c("anc Z1", "neo Z1")) +
  scale_fill_manual(values = c(anz_Z_color, neo_17_color)) +
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
        axis.text.y = element_text(size=axis_text_sizes)) 


# dn
wilcox.test(dn ~ anc_vs_neo, data = Z1, paired=FALSE, conf.int=TRUE)
# p-value = 0.8213

Z1_dn_ans_vs_neo <- ggplot(data=Z1, aes(y=dn, x=anc_vs_neo, fill=anc_vs_neo)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  coord_cartesian(ylim=c(0, 0.015)) +
  scale_y_continuous(name = "dN") +
  scale_fill_manual(values = c(anz_Z_color, neo_17_color)) +
  scale_x_discrete(labels = c("anc Z1", "neo Z1")) +
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
        axis.text.y = element_text(size=axis_text_sizes)) 


# ds
wilcox.test(ds ~ anc_vs_neo, data = Z1, paired=FALSE, conf.int=TRUE)
# p-value = 0.1134

Z1_ds_ans_vs_neo <- ggplot(data=Z1, aes(y=ds, x=anc_vs_neo, fill=anc_vs_neo)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  coord_cartesian(ylim=c(0, 0.06)) +
  scale_y_continuous(name = "dS") +
  scale_fill_manual(values = c(anz_Z_color, neo_17_color)) +
  scale_x_discrete(labels = c("anc Z1", "neo Z1")) +
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
        axis.text.y = element_text(size=axis_text_sizes)) 



#-------------------#
##### Figure 2C #####
#-------------------#


# plot individual chromosomes...


# sort by dN/dS...
df_filtered$chrN = with(df_filtered, reorder(chrN, dnds, median))


# get medians for plotting...
A_dnds <- df_filtered[df_filtered$anc_vs_neo == "A", ]
median_dnds_A <- median(A_dnds$dnds)
Z_anc_dnds <- df_filtered[df_filtered$anc_vs_neo == "anc Z", ]
median_dnds_anc_Z <- median(Z_anc_dnds$dnds)
Z_neo_dnds <- df_filtered[df_filtered$anc_vs_neo == "neo Z", ]
median_dnds_neo_Z <- median(Z_neo_dnds$dnds)


all_chr_plot <- ggplot(data=df_filtered, aes(y=dnds, x=chrN, fill=chrN)) +
  geom_boxplot(varwidth = TRUE, outlier.shape = NA) +
  coord_flip(ylim=c(0, 0.75)) +
  scale_y_continuous(name = "dN/dS") +
  scale_fill_manual(values = c(rep(auto_color, 25), rep("firebrick2", 2), auto_color, "firebrick2")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "grey10"),
        axis.title.x = element_text(size=axis_text_sizes),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=axis_text_sizes),
        axis.text.y = element_text(size=axis_text_sizes, hjust=0.5)) +
  geom_hline(yintercept = median_dnds_neo_Z, linetype = "dashed", color = "firebrick4", size=0.75) +
  geom_hline(yintercept = median_dnds_anc_Z, linetype = "dashed", color = anz_Z_color, size=0.75) +
  geom_hline(yintercept = median_dnds_A, linetype = "dashed", color = "grey50", size=0.75) +
  annotate("text", x = 7, y = 0.675, label = "Medians") +
  annotate("text", x = 6, y = 0.7, label = "A", hjust = 0) +
  annotate("segment", x = 6, xend = 6, y = 0.625, yend = 0.68, linetype = "dashed", color=auto_color, size=0.75) +
  annotate("text", x = 5, y = 0.7, label = "anc Z", hjust = 0) +
  annotate("segment", x = 5, xend = 5, y = 0.625, yend = 0.68, linetype = "dashed", color=anz_Z_color, size=0.75) +
  annotate("text", x = 4, y = 0.7, label = "neo Z", hjust = 0) +
  annotate("segment", x = 4, xend = 4, y = 0.625, yend = 0.68, linetype = "dashed", color="firebrick4", size=0.75) 


#---------------------------#
##### Dunn test heatmap #####
#---------------------------#


dunn_sign_table <- read.csv("sign_dunn_test_count_for_plot.txt", sep = "\t")

chr_order <- aggregate(dnds ~ chrN, median, data=df_filtered)

sorted_dunn_sign_table <- merge(dunn_sign_table, chr_order, by="chrN")

sorted_dunn_sign_table$chrN = with(sorted_dunn_sign_table, reorder(chrN, dnds))


dunn_heat_plot <- ggplot(sorted_dunn_sign_table, aes(x=chrN, y=0, fill=count)) +
  geom_tile() +
  coord_flip(ylim = c(-0.5, 0.5), clip = "off") +
  theme_minimal() +
  theme(legend.position="left",
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank()) +
  scale_fill_continuous(name= "",
                       low = "white", high = "darkred",
                       limits=c(0, max(sorted_dunn_sign_table$count)),
                       breaks=c(0, 9, 18, 27)) +
  annotate("text", x = 18.5, y = -2, label = "Number of\nsignificantly\ndifferent\nchromosomes")



#-----------------------#
##### Combine plots #####
#-----------------------#


Figure_2 <- (dnds_anc_neo_filtered_plot + dn_anc_neo_filtered_plot + ds_anc_neo_filtered_plot) /
            (Z1_dnds_ans_vs_neo + Z1_dn_ans_vs_neo + Z1_ds_ans_vs_neo) |
            (plot_spacer() + dunn_heat_plot + all_chr_plot + plot_layout(widths = c(1, 1, 8)))

Figure_2 + plot_annotation(tag_levels = list(c("A", "", "", "B", "", "", "C"))) & 
  theme(plot.tag = element_text(size = 20))

ggsave(filename = ("figure_2.png"), width = 15, height = 10, units = "in", dpi = 300, limitsize = F)

