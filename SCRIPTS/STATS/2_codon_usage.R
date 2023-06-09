
#-------------------------------#
##### Test codon usage bias ##### 
#-------------------------------#

library(coRdon)

setwd("~/Desktop/Plots/FASTZ")

Lsin_auto <- readSet(file = "autosomal_genes.fasta")
Lsin_auto_ct <- codonTable(Lsin_auto)

Lsin_ancZ <- readSet(file = "Z1_anc_genes.fasta")
Lsin_ancZ_ct <- codonTable(Lsin_ancZ)

Lsin_neoZ <- readSet(file = "Z_all_neo_genes.fasta")
Lsin_neoZ_ct <- codonTable(Lsin_neoZ)

enc_p_ancZ <- ENCprime(Lsin_ancZ_ct, filtering = "hard")
enc_p_neoZ <- ENCprime(Lsin_neoZ_ct, filtering = "hard")
enc_p_auto <- ENCprime(Lsin_auto_ct, filtering = "hard")


ENC_Z <- as.data.frame(enc_p_ancZ)
ENC_Z$chr <- "anc_Z"
colnames(ENC_Z) <- c("ENC", "chr")

ENC_neoZ <- as.data.frame(enc_p_neoZ)
ENC_neoZ$chr <- "neo_Z"
colnames(ENC_neoZ) <- c("ENC", "chr")

ENC_A <- as.data.frame(enc_p_auto)
ENC_A$chr <- "A"
colnames(ENC_A) <- c("ENC", "chr")


ENC_all <- rbind(ENC_A, ENC_Z, ENC_neoZ)

ENC_median = aggregate(ENC ~ chr, data=ENC_all, median)

#chr      ENC
#    A 58.35151
#anc_Z 59.41952
#neo_Z 58.30952


kruskal.test(ENC ~ chr, data = ENC_all)

# Kruskal-Wallis chi-squared = 27.534, df = 2, p-value = 1.05e-06

pairwise.wilcox.test(ENC_all$ENC, ENC_all$chr,
                     p.adjust.method = "bonferroni")

#A       anc_Z  
#anc_Z 4.7e-07 -      
#neo_Z 1       6.4e-05
