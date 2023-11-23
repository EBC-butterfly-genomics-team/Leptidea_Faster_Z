# Convert from FPKM to normalized counts using scaling factors computed after excluding top 25% expressed genes (4th qantile) 



##### Clear all states
rm(list=ls(all=TRUE))


#v results folder
results_dir <- "/home/larshook/LarsH/REVISIT_DC/NORMALIZED_LIBS"


##### Load gene count matrix (raw counts) and labels ####

setwd("/home/larshook/LarsH/REVISIT_DC/NORMALIZED_LIBS")
countData <- read.delim("gene_count_matrix.csv", header = TRUE, sep = ",")
head(countData)

setwd("/home/larshook/LarsH/REVISIT_DC/SCRIPTS")
colData <- read.csv("sample_info.txt", sep="\t", row.names=1)
colData





##### Load info about A/Z gene assignment 

setwd("/home/larshook/LarsH/REVISIT_DC/Assign_genes_to_chromosome")

instar_V <- read.delim("instar_V-assigned_A_or_Z-filtered.txt", header = TRUE)
head(instar_V)
pupa <- read.delim("pupa-assigned_A_or_Z-filtered.txt", header = TRUE)
head(pupa)
adult <- read.delim("adult-assigned_A_or_Z-filtered.txt", header = TRUE)
head(adult)





# cd to output folder
setwd(results_dir)




#Make separate data sets for each developmental stage ####
Instar_countdata <- subset(countData, select = c("gene_id","P5052_202_S48", "P5052_210_S51", "P5052_218_S55", 
                                                 "P5052_203_S49", "P5052_211_S52", "P5052_219_S56"))
Instar_coldata <- as.data.frame(colData[c(1, 5, 9, 2, 6, 10), ])

Pupa_countdata <- subset(countData, select = c("gene_id","P5052_204_S34", "P5052_212_S36", "P5052_220_S57",
                                               "P5052_205_S35", "P5052_213_S53", "P5052_221_S58"))
Pupa_coldata <- as.data.frame(colData[c(3, 7, 11, 4, 8, 12)  , ])

Adult_countdata <- subset(countData, select = c("gene_id","P5052_233_S60", "P5052_234_S61", "P5052_235_S62",
                                                "P5052_226_S59", "P5052_227_S38", "P7553_325_S10"))
Adult_coldata <- as.data.frame(colData[c(15, 16, 17, 13, 14, 18), ])

Instar_coldata <- Instar_coldata[ , 2:3]
Pupa_coldata <- Pupa_coldata[ , 2:3]
Adult_coldata <- Adult_coldata[ , 2:3]




#Select only autosomal genes"
Instar_countdata_AZ <-merge(x = Instar_countdata, y = instar_V, by = "gene_id", all.y = TRUE)       #select only genes assigned to A or Z 
Instar_countdata_A <- Instar_countdata_AZ[Instar_countdata_AZ$chromosome == 'A', ]                  #select only autosomal genes
Instar_countdata_A <- Instar_countdata_A[,-c(8:12)]

Pupa_countdata_AZ <-merge(x = Pupa_countdata, y = pupa, by = "gene_id", all.y = TRUE)
Pupa_countdata_A <- Pupa_countdata_AZ[Pupa_countdata_AZ$chromosome == 'A', ]
Pupa_countdata_A <- Pupa_countdata_A[,-c(8:12)]

Adult_countdata_AZ <-merge(x = Adult_countdata, y = adult, by = "gene_id", all.y = TRUE)
Adult_countdata_A <- Adult_countdata_AZ[Adult_countdata_AZ$chromosome == 'A', ]
Adult_countdata_A <- Adult_countdata_A[,-c(8:12)]









#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####
Instar_countdata_A <- Instar_countdata_A[rowSums(Instar_countdata_A[2:7] > 0) >= 2 , ]
Pupa_countdata_A <- Pupa_countdata_A[rowSums(Pupa_countdata_A[2:7] > 0) >= 2 , ]
Adult_countdata_A <- Adult_countdata_A[rowSums(Adult_countdata_A[2:7] > 0) >= 2 , ]








# Compute quantiles for each sample (based on counts); include only non-zero counts when computing quantiles
# Compute scaling factors (sum of non-zero gene expression levels below 75% percentile)


## male instar

Aux_1 <- Instar_countdata_A 
# remove zero counts
Aux_1 <- Aux_1[Aux_1$P5052_202_S48 != 0, ]    
#compute quantiles
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_202_S48, quantile(P5052_202_S48, probs=0:4/4), include.lowest=TRUE)))
# Compute scaling factors (sum of non-zero gene expression levels below 75% percentile)
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_202 <- sum(Aux_1$P5052_202_S48)
sf_202

rm(Aux_1)
Aux_1 <- Instar_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_210_S51 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_210_S51, quantile(P5052_210_S51, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_210<- sum(Aux_1$P5052_210_S51)
sf_210

rm(Aux_1)
Aux_1 <- Instar_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_218_S55 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_218_S55, quantile(P5052_218_S55, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_218<- sum(Aux_1$P5052_218_S55)
sf_218


## female instar

rm(Aux_1)
Aux_1 <- Instar_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_203_S49 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_203_S49, quantile(P5052_203_S49, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_203 <- sum(Aux_1$P5052_203_S49)
sf_203

rm(Aux_1)
Aux_1 <- Instar_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_211_S52 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_211_S52, quantile(P5052_211_S52, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_211<- sum(Aux_1$P5052_211_S52)
sf_211

rm(Aux_1)
Aux_1 <- Instar_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_219_S56 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_219_S56, quantile(P5052_219_S56, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_219<- sum(Aux_1$P5052_219_S56)
sf_219



## male pupa

rm(Aux_1)
Aux_1 <- Pupa_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_204_S34 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_204_S34, quantile(P5052_204_S34, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_204 <- sum(Aux_1$P5052_204_S34)
sf_204

rm(Aux_1)
Aux_1 <- Pupa_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_212_S36 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_212_S36, quantile(P5052_212_S36, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_212 <- sum(Aux_1$P5052_212_S36)
sf_212

rm(Aux_1)
Aux_1 <- Pupa_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_220_S57 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_220_S57, quantile(P5052_220_S57, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_220 <- sum(Aux_1$P5052_220_S57)
sf_220


## female pupa

rm(Aux_1)
Aux_1 <- Pupa_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_205_S35 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_205_S35, quantile(P5052_205_S35, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_205 <- sum(Aux_1$P5052_205_S35)
sf_205

rm(Aux_1)
Aux_1 <- Pupa_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_213_S53 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_213_S53, quantile(P5052_213_S53, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_213 <- sum(Aux_1$P5052_213_S53)
sf_213

rm(Aux_1)
Aux_1 <- Pupa_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_221_S58 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_221_S58, quantile(P5052_221_S58, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_221 <- sum(Aux_1$P5052_221_S58)
sf_221


## male adult

rm(Aux_1)
Aux_1 <- Adult_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_233_S60 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_233_S60, quantile(P5052_233_S60, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_233 <- sum(Aux_1$P5052_233_S60)
sf_233

rm(Aux_1)
Aux_1 <- Adult_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_234_S61 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_234_S61, quantile(P5052_234_S61, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_234 <- sum(Aux_1$P5052_234_S61)
sf_234

rm(Aux_1)
Aux_1 <- Adult_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_235_S62 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_235_S62, quantile(P5052_235_S62, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_235 <- sum(Aux_1$P5052_235_S62)
sf_235


## female adult

rm(Aux_1)
Aux_1 <- Adult_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_226_S59 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_226_S59, quantile(P5052_226_S59, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_226 <- sum(Aux_1$P5052_226_S59)
sf_226

rm(Aux_1)
Aux_1 <- Adult_countdata_A 
Aux_1 <- Aux_1[Aux_1$P5052_227_S38 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P5052_227_S38, quantile(P5052_227_S38, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_227 <- sum(Aux_1$P5052_227_S38)
sf_227

rm(Aux_1)
Aux_1 <- Adult_countdata_A 
Aux_1 <- Aux_1[Aux_1$P7553_325_S10 != 0, ]    
Aux_1<- within(Aux_1, quartile_sample <- as.integer(cut(P7553_325_S10, quantile(P7553_325_S10, probs=0:4/4), include.lowest=TRUE)))
Aux_1 <- Aux_1[Aux_1$quartile_sample < '4', ]   
sf_325 <- sum(Aux_1$P7553_325_S10)
sf_325




# Summary scaling factors

# instar-males, below 75% percentile (excluding zero-count genes)
c(sf_202,sf_210,sf_218)

# instar-males, all genes (raw)
c(sum(countData$P5052_202_S48),sum(countData$P5052_210_S51),sum(countData$P5052_218_S55))



# instar-females, below 75% percentile (excluding zero-count genes)
c(sf_203,sf_211,sf_219)

# instar-females, all genes (raw)
c(sum(countData$P5052_203_S49),sum(countData$P5052_211_S52),sum(countData$P5052_219_S56))



# pupa-males, below 75% percentile (excluding zero-count genes)
c(sf_204,sf_212,sf_220)

# pupa-males,all genes (raw)
c(sum(countData$P5052_204_S34),sum(countData$P5052_212_S36),sum(countData$P5052_220_S57))



# pupa-females, below 75% percentile (excluding zero-count genes)
c(sf_205,sf_213,sf_221)

# pupa-females,all genes (raw)
c(sum(countData$P5052_205_S35),sum(countData$P5052_213_S53),sum(countData$P5052_221_S58))



# adult-males, below 75% percentile (excluding zero-count genes)
c(sf_233,sf_234,sf_235)

# adult-males,all genes (raw)
c(sum(countData$P5052_233_S60),sum(countData$P5052_234_S61),sum(countData$P5052_235_S62))


# adult-females, below 75% percentile (excluding zero-count genes)
c(sf_226,sf_227,sf_325)

# adult-females,all genes (raw)
c(sum(countData$P5052_226_S59),sum(countData$P5052_227_S38),sum(countData$P7553_325_S10))










#Read STAR/StringTie output files (STAR.Aligned.out.gene_abund.tab)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_202_S48/stringtie_2")
stringtie202 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie202)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_203_S49/stringtie_2")
stringtie203 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie203)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_204_S34/stringtie_2")
stringtie204 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie204)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_205_S35/stringtie_2")
stringtie205 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie205)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_210_S51/stringtie_2")
stringtie210 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie210)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_211_S52/stringtie_2")
stringtie211 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie211)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_212_S36/stringtie_2")
stringtie212 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie212)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_213_S53/stringtie_2")
stringtie213 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie213)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_218_S55/stringtie_2")
stringtie218 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie218)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_219_S56/stringtie_2")
stringtie219 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie219)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_220_S57/stringtie_2")
stringtie220 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie220)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_221_S58/stringtie_2")
stringtie221 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie221)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_226_S59/stringtie_2")
stringtie226 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie226)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_227_S38/stringtie_2")
stringtie227 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie227)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_233_S60/stringtie_2")
stringtie233 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie233)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_234_S61/stringtie_2")
stringtie234 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie234)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P5052_235_S62/stringtie_2")
stringtie235 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie235)

setwd("/home/larshook/LarsH/RNAseq_Lsinapis_DC/RNAseq_Lsinapis/2e_STAR/P7553_325_S10/stringtie_2")
stringtie325 <- read.delim("STAR.Aligned.out.gene_abund.tab", header = TRUE, sep = "\t")
head(stringtie325)





### change column names (back to original names)

colnames(stringtie202)[1] <- "Gene_ID"
colnames(stringtie203)[1] <- "Gene_ID"
colnames(stringtie204)[1] <- "Gene_ID"
colnames(stringtie205)[1] <- "Gene_ID"
colnames(stringtie210)[1] <- "Gene_ID"
colnames(stringtie211)[1] <- "Gene_ID"
colnames(stringtie212)[1] <- "Gene_ID"
colnames(stringtie213)[1] <- "Gene_ID"
colnames(stringtie218)[1] <- "Gene_ID"
colnames(stringtie219)[1] <- "Gene_ID"
colnames(stringtie220)[1] <- "Gene_ID"
colnames(stringtie221)[1] <- "Gene_ID"
colnames(stringtie226)[1] <- "Gene_ID"
colnames(stringtie227)[1] <- "Gene_ID"
colnames(stringtie233)[1] <- "Gene_ID"
colnames(stringtie234)[1] <- "Gene_ID"
colnames(stringtie235)[1] <- "Gene_ID"
colnames(stringtie325)[1] <- "Gene_ID"

colnames(stringtie202)[2] <- "Gene Name"
colnames(stringtie203)[2] <- "Gene Name"
colnames(stringtie204)[2] <- "Gene Name"
colnames(stringtie205)[2] <- "Gene Name"
colnames(stringtie210)[2] <- "Gene Name"
colnames(stringtie211)[2] <- "Gene Name"
colnames(stringtie212)[2] <- "Gene Name"
colnames(stringtie213)[2] <- "Gene Name"
colnames(stringtie218)[2] <- "Gene Name"
colnames(stringtie219)[2] <- "Gene Name"
colnames(stringtie220)[2] <- "Gene Name"
colnames(stringtie221)[2] <- "Gene Name"
colnames(stringtie226)[2] <- "Gene Name"
colnames(stringtie227)[2] <- "Gene Name"
colnames(stringtie233)[2] <- "Gene Name"
colnames(stringtie234)[2] <- "Gene Name"
colnames(stringtie235)[2] <- "Gene Name"
colnames(stringtie325)[2] <- "Gene Name"




### Read file listing StringTie's per-million scalling factor (#reads/1e6)

setwd(results_dir)
stringtie_SF <- read.delim("StringTie_scalling-factors.txt", header = TRUE, sep = "\t")

StringTieSF_202 <- stringtie_SF[stringtie_SF$sample == 'P5052_202_S48', 2]
StringTieSF_210 <- stringtie_SF[stringtie_SF$sample == 'P5052_210_S51', 2]
StringTieSF_218 <- stringtie_SF[stringtie_SF$sample == 'P5052_218_S55', 2]

StringTieSF_203 <- stringtie_SF[stringtie_SF$sample == 'P5052_203_S49', 2]
StringTieSF_211 <- stringtie_SF[stringtie_SF$sample == 'P5052_211_S52', 2]
StringTieSF_219 <- stringtie_SF[stringtie_SF$sample == 'P5052_219_S56', 2]

StringTieSF_204 <- stringtie_SF[stringtie_SF$sample == 'P5052_204_S34', 2]
StringTieSF_212 <- stringtie_SF[stringtie_SF$sample == 'P5052_212_S36', 2]
StringTieSF_220 <- stringtie_SF[stringtie_SF$sample == 'P5052_220_S57', 2]

StringTieSF_205 <- stringtie_SF[stringtie_SF$sample == 'P5052_205_S35', 2]
StringTieSF_213 <- stringtie_SF[stringtie_SF$sample == 'P5052_213_S53', 2]
StringTieSF_221 <- stringtie_SF[stringtie_SF$sample == 'P5052_221_S58', 2]

StringTieSF_233 <- stringtie_SF[stringtie_SF$sample == 'P5052_233_S60', 2]
StringTieSF_234 <- stringtie_SF[stringtie_SF$sample == 'P5052_234_S61', 2]
StringTieSF_235 <- stringtie_SF[stringtie_SF$sample == 'P5052_235_S62', 2]

StringTieSF_226 <- stringtie_SF[stringtie_SF$sample == 'P5052_226_S59', 2]
StringTieSF_227 <- stringtie_SF[stringtie_SF$sample == 'P5052_227_S38', 2]
StringTieSF_325 <- stringtie_SF[stringtie_SF$sample == 'P7553_325_S10', 2]




# re-normalize FPKM counts (replace STRINGTIE scaling factors with new scaling factors)

stringtie202_n <- stringtie202
stringtie203_n <- stringtie203
stringtie204_n <- stringtie204
stringtie205_n <- stringtie205
stringtie210_n <- stringtie210
stringtie211_n <- stringtie211
stringtie212_n <- stringtie212
stringtie213_n <- stringtie213
stringtie218_n <- stringtie218
stringtie219_n <- stringtie219
stringtie220_n <- stringtie220
stringtie221_n <- stringtie221
stringtie226_n <- stringtie226
stringtie227_n <- stringtie227
stringtie233_n <- stringtie233
stringtie234_n <- stringtie234
stringtie235_n <- stringtie235
stringtie325_n <- stringtie325


stringtie202_n$FPKM <- (stringtie202$FPKM * StringTieSF_202) *1000000 / sf_202
stringtie203_n$FPKM <- (stringtie203$FPKM * StringTieSF_203) *1000000 / sf_203
stringtie204_n$FPKM <- (stringtie204$FPKM * StringTieSF_204) *1000000 / sf_204
stringtie205_n$FPKM <- (stringtie205$FPKM * StringTieSF_205) *1000000 / sf_205
stringtie210_n$FPKM <- (stringtie210$FPKM * StringTieSF_210) *1000000 / sf_210
stringtie211_n$FPKM <- (stringtie211$FPKM * StringTieSF_211) *1000000 / sf_211
stringtie212_n$FPKM <- (stringtie212$FPKM * StringTieSF_212) *1000000 / sf_212
stringtie213_n$FPKM <- (stringtie213$FPKM * StringTieSF_213) *1000000 / sf_213
stringtie218_n$FPKM <- (stringtie218$FPKM * StringTieSF_218) *1000000 / sf_218
stringtie219_n$FPKM <- (stringtie219$FPKM * StringTieSF_219) *1000000 / sf_219
stringtie220_n$FPKM <- (stringtie220$FPKM * StringTieSF_220) *1000000 / sf_220
stringtie221_n$FPKM <- (stringtie221$FPKM * StringTieSF_221) *1000000 / sf_221
stringtie226_n$FPKM <- (stringtie226$FPKM * StringTieSF_226) *1000000 / sf_226
stringtie227_n$FPKM <- (stringtie227$FPKM * StringTieSF_227) *1000000 / sf_227
stringtie233_n$FPKM <- (stringtie233$FPKM * StringTieSF_233) *1000000 / sf_233
stringtie234_n$FPKM <- (stringtie234$FPKM * StringTieSF_234) *1000000 / sf_234
stringtie235_n$FPKM <- (stringtie235$FPKM * StringTieSF_235) *1000000 / sf_235
stringtie325_n$FPKM <- (stringtie325$FPKM * StringTieSF_325) *1000000 / sf_325










#create output folders and save normalized FPKM values
setwd(results_dir)
dir.create("STAR_norm", showWarnings = FALSE)
STAR_dir <- paste(results_dir, "/STAR_norm/", sep="")
setwd(STAR_dir)
dir.create("P5052_202_S48/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_203_S49/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_204_S34/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_205_S35/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_210_S51/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_211_S52/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_212_S36/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_213_S53/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_218_S55/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_219_S56/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_220_S57/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_221_S58/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_226_S59/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_227_S38/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_233_S60/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_234_S61/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P5052_235_S62/stringtie_2", showWarnings = FALSE, recursive = TRUE)
dir.create("P7553_325_S10/stringtie_2", showWarnings = FALSE, recursive = TRUE)

setwd(paste(STAR_dir, "P5052_202_S48/stringtie_2", sep=""))
write.table(stringtie202_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_203_S49/stringtie_2", sep=""))
write.table(stringtie203_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_204_S34/stringtie_2", sep=""))
write.table(stringtie204_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_205_S35/stringtie_2", sep=""))
write.table(stringtie205_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_210_S51/stringtie_2", sep=""))
write.table(stringtie210_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_211_S52/stringtie_2", sep=""))
write.table(stringtie211_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_212_S36/stringtie_2", sep=""))
write.table(stringtie212_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_213_S53/stringtie_2", sep=""))
write.table(stringtie213_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_218_S55/stringtie_2", sep=""))
write.table(stringtie218_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_219_S56/stringtie_2", sep=""))
write.table(stringtie219_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_220_S57/stringtie_2", sep=""))
write.table(stringtie220_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_221_S58/stringtie_2", sep=""))
write.table(stringtie221_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_226_S59/stringtie_2", sep=""))
write.table(stringtie226_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_227_S38/stringtie_2", sep=""))
write.table(stringtie227_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_233_S60/stringtie_2", sep=""))
write.table(stringtie233_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_234_S61/stringtie_2", sep=""))
write.table(stringtie234_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P5052_235_S62/stringtie_2", sep=""))
write.table(stringtie235_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)

setwd(paste(STAR_dir, "P7553_325_S10/stringtie_2", sep=""))
write.table(stringtie325_n, 
            file = "STAR.Aligned.out.gene_abund.tab", sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE)


