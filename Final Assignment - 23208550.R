new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )

setwd(new_dir)

# We will use the following files:

# data_clinical_patient.txt, data_mrna_seq_v2_rsem.txt, data_mutations.txt and data_cna.txt

clinical = read.delim("data_clinical_patient.txt")

rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

# in this assignment we will delete the genes for which there's more than one Hugo Symbol
# These are typically genes with no Hugo Symbol ("" as an entry) or pseudogenes.

# This is more for simplicity.If you keep your analysis would still be correct so no worries.

keep = !duplicated(rnaseq[,1])

rnaseq = rnaseq[keep,]

# set rownames of rnaseq to hugo symbols

rownames(rnaseq)  = rnaseq[,1]

# Read CNA Data

cna = read.delim('data_cna.txt')




# find ERBB2 in cna

erbb2_indx = which(cna[,1] == 'ERBB2')

# Plot histogram to visualize explore the data.

hist(as.numeric(cna[erbb2_indx,-c(1,2)]))

# match patients in rnaseq to patients in cna.

rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))

# select only the rna cases which have cna data.

rna_cna_sub = rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna

no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rnaseq[,2+rna_cna_id]), colnames(cna[,-c(1,2)]))) 

# sanity check.This will print an error if the result is not the same.

sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]




# Pre-allocate memory for ERBB2

meta_erbb2 = matrix(0,length(rna_cna_id),1)

#gets only ERBB2 greater than 0 

for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
}




# This are some checks you can do to make sure your code worked.
# There's some more systematic checks you can do. See unit testing.

# simple checks to make sure. 

col_i = colnames(rna_cna_sub)[1]

col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.

pos_example = which(meta_erbb2==1)[1]


col_i = colnames(rna_cna_sub)[pos_example]

col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]

# botch checks should print true.

# We will add a title to the metadata.

colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers

rna_cna_sub = round(rna_cna_sub)





# Install DESeq2.

BiocManager::install(c(
  "cli", "fansi", "matrixStats"
), update = TRUE, ask = FALSE, force = TRUE)


library(BiocManager)
BiocManager::valid()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDb",force = TRUE)

# Install DeSeq2

BiocManager::install("DESeq2")

library(DESeq2)





# This is the differential gene analysis 

rna_cna_sub[is.na(rna_cna_sub)] = 0  # Impute with zeros the NA
rna_cna_sub[rna_cna_sub<0] = 0

dds <- DESeqDataSetFromMatrix(countData = round(rna_cna_sub),
                              colData = meta_erbb2,
                              design = ~ ERBB2Amp)




# Filter

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Normalize

dds <- DESeq(dds)

# Transform the data to visualize
rld <- vst(dds, blind=FALSE)

# Do Principal Components Analysis
assay = as.matrix(rna_cna_sub)

pc = prcomp(assay(rld))

# Plot 
plot(pc$rotation[,1], pc$rotation[,2], col = 1+(meta_erbb2), pch = 19,main="Principle Component Analysis",xlab="First principle Component",ylab="Second Principle Component")

# The pca plot might not show much but it is good to show it.

hist(assay[1000,], breaks = 50)

hist(assay(rld)[1000,], breaks = 50)

# Get Results

res <- results(dds)

# Summary

summary(res)
rownames(res) = rnaseq[keep,1]

# Significantly Differentially Expressed

signif = which(res$padj<0.05)
deg = res[signif,]



# Separate them 
dup = deg[deg[,2]>0.,]

ddown = deg[deg[,2]<0.,]


#The top 10 differentially expressed genes based on absolute value 

top_10 <- deg[sort(abs(deg$log2FoldChange),decreasing=T,index.return=T)[[2]],]
head(top_10,n=20)

top_10['BNIP3',]

top_10['BNIP3L',]



entrez_ids = rnaseq[keep,2]

entrez_all = entrez_ids[signif]
entrez_up = entrez_all[signif[deg[,2]>0.]]
entrez_down = entrez_all[signif[deg[,2]<0.]]




# Step 10 Peform a Pathway Enrichment Analysis 


BiocManager::install(c(
  "clusterProfiler"
), update = TRUE, ask = FALSE, force = TRUE)

library(clusterProfiler)

# Do a KEGG pathway over-representation analysis

all_paths =   enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.05)
head(all_paths)

# Optionally you can divide between up and down.
# Both options are Ok for the assignment.

up_paths = enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.05)
head(up_paths)


down_paths = enrichKEGG(gene = entrez_down, organism = 'hsa', pvalueCutoff = 0.05)
head(down_paths)











#Survival Analysis 

#creates a new dataframe of patient data called surv_patient 
surv_patient = data_patient
cna_surv = cna

#This checks changes the patient identifier so that the cna and patient data can be compared
surv_patient$X.Patient.Identifier <- gsub("-", ".",surv_patient$X.Patient.Identifier)
surv_patient <- surv_patient |> mutate(X.Patient.Identifier=paste0(X.Patient.Identifier,c(".01",".01")))
surv_patient$X.Patient.Identifier

rownames(surv_patient)  = surv_patient[,1]




# Getting all of the patients in the cna data and matching to the patient data 
pat_id = which(is.element(rownames(surv_patient), colnames(cna[,-c(1,2)])))

# select only the patients that are in the cna 
pat_sub = surv_patient[pat_id,]

#new column for HER2 positive
pat_sub['ERBB2Amp'] = 0 




#This adds a new column to pat sub which shows if the patient is HER2 + etc
for (i in 1:length(pat_id)){
  # access the colnames of i
  col_i = rownames(pat_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  pat_sub[i,39] <- 1*(cna[erbb2_indx,col_cna]>0)
    }

install.packages("ranger")
install.packages("ggfortify")

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

#creates new dataframe with just HER2 positive patients 
pat_sub_her <- filter(pat_sub,ERBB2Amp == 1)
#Gets rid of any blank or NA rows in the disease specific survival status column 
pat_sub_her <- pat_sub_her[!(is.na(pat_sub_her$Overall.Survival.Status) | pat_sub_her$Overall.Survival.Status==""), ]
#Changes the survival status to 0 or 1 where 0 is alive or tumor free dead or dead from tumor 



pat_sub_her$Overall.Survival.Status <- factor(pat_sub_her$Overall.Survival.Status, order=TRUE, labels=c(0,1))



pat_sub_her$Overall.Survival..Months. <- as.numeric(pat_sub_her$Overall.Survival..Months.)


#pat_sub_her$Overall.Survival.Status
#pat_sub_her$Disease.specific.Survival.status

km <- with(pat_sub_her,Surv(time = as.numeric(Overall.Survival..Months.),as.numeric(Overall.Survival.Status)))

km_fit <- survfit(Surv(as.numeric(Overall.Survival..Months.),as.numeric(Overall.Survival.Status)) ~ 1, data=pat_sub_her)
summary(km_fit, times = c(1,30,60,90*(1:10)))
autoplot(km_fit,xlab="Time in months",ylab = "Survival %", main ="Survival % vs Time in Months")

