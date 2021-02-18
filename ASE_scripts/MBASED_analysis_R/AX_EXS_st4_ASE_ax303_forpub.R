#########################################################
# ASE ANALYSIS using MBASED (developped by Oleg Mayba)
# 
# Allelic expression of color gene 
# between axillaris and exserta 
# using petal (limb only) tissue from stage 4 buds
#
#########################################################

#########################################################
# Workflow was written by Michel Moser                  #
# For questions and issues: michel.moser@ips.unibe.ch   #
#########################################################

####INSTALLING packages ########
# source("https://bioconductor.org/biocLite.R")
# biocLite("MBASED")
# biocLite("ggbio")
# biocLite("VariantAnnotation")
# # detach("package:DESeq2", unload=TRUE)
# install.packages("gtools")



#######loading packages ##############
library("MBASED")
library(gtools)
library(ggbio)
library(ggplot2)
library(VariantAnnotation)
library(rtracklayer)
library(reshape2)
library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)
library(scales)
library(tidyr)



#set a seed:
set.seed(23)

#function to extract results provided by Oleg

summarizeASEResults_1s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequency=assays(MBASEDOutput)$majorAlleleFrequency[,1],
    pValueASE=assays(MBASEDOutput)$pValueASE[,1],  
    pValueHeterogeneity=assays(MBASEDOutput)$pValueHeterogeneity[,1]
  )   
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID))) 
  return(
    list(
      geneOutput=geneOutputDF,
      locusOutput=lociOutputList
    )
  )
}


#########################################################
# GENERAL INFO                                          #
#########################################################

# working with list of the candidate genes
# MAF = major allele frequency
# make sure intronic regions get extracted
# as depth and MAF varies between exons and introns


#### Ax * Exs RNAseq ####

#########################################################
#                                                       #     
#         DATA:                                         #
#                                                       #
# 2018 F1 AX x EXS stage 4 limb     3 bio. replicates   #
# per genotype                                          #
#########################################################

# non conservative SNPs in exons only
# done with a raw vcf file, minimal filtering
# remove lines with identical allele

#setwd("~/yourdirectory")
#name the headers with appropriate species names - we use Ax for axillaris alleles and Ex for exserta alleles
header = c("gene", "first", "last", "scaffold", "pos", "baseAx", "baseExs", "rep1Ax", "rep1Exs", "rep1tot", "rep2Ax", "rep2Exs", "rep2tot", "rep3Ax", "rep3Exs", "rep3tot")

axexs18limb <- read.table("F1s.strict.counts.short.sorted", sep="\t")
dim(axexs18limb)

names(axexs18limb) <- header

dim(axexs18limb)

#remove SNPs from Exserta that are also in axillaris (this shouldn't matter anymore since I did this bioinformatically a while back in the pipeline)
axexs18limb <- axexs18limb[which(axexs18limb$baseAx != axexs18limb$baseExs),]

dim(axexs18limb)


####### general dataset inspection: ###########3
head(axexs18limb)
dim(axexs18limb)
#how many genes are present in the dataset?
length(unique(axexs18limb$gene))
#[1] 9525

#how many SNPs are there for each of the present genes?
table(axexs18limb$gene)



#check class of each column in the dataframe
lapply(axexs18limb, class)

#create mean of each allelic count as new column to dataframe
axexs18limb$axmean <- (axexs18limb$rep1Ax + axexs18limb$rep2Ax + axexs18limb$rep3Ax)/3
axexs18limb$exsmean <- (axexs18limb$rep1Exs + axexs18limb$rep2Exs + axexs18limb$rep3Exs)/3

axexs18limb$gene <- as.character(axexs18limb$gene)
#order according to gene_name
axexs18limb_ord <- axexs18limb[mixedorder(axexs18limb$gene),]
dim(axexs18limb_ord)
head(axexs18limb_ord)

#v303/304, remove ".1.path1" from the gene name and ".2.path1" for JAF13 which is misnamed in the file
axexs18limb_ord$gene
axexs18limb_ord$gene <- gsub(".1.path1", "", axexs18limb_ord$gene)
axexs18limb_ord$gene <- gsub(".2.path1", "", axexs18limb_ord$gene)

#look at SNPs from a specific gene:

#DPL
axexs18limb_ord[which(axexs18limb_ord$gene == "Peaxi162Scf01210g00002"),]



#### load genes of interest here OR can enter a list of all genes in dataset. ####

phenyl_genes <- read.table("petunia_flavonoid_pathway_gene_names_and_ids.txt", sep = "\t", header=FALSE)

names(phenyl_genes) <- c("gene", "AXILLARISv4") 

#
gene_list <- as.character(phenyl_genes$AXILLARISv4)
length(gene_list) 

#extract all SNPs from the candidate genes
head(axexs18limb_ord[axexs18limb_ord$gene %in% gene_list,])
dim(axexs18limb_ord[axexs18limb_ord$gene %in% gene_list,])

#get all SNPs within the candidate genes from the scent list
candidate_loci <- axexs18limb_ord[axexs18limb_ord$gene %in% gene_list,]

dim(candidate_loci)
candidate_loci

#for which genes of the list are SNPs found?
present_candidategenes <- as.character(candidate_loci$gene)
#gives you summary of the unique values and their count
length(unique(present_candidategenes))
present <- unique(present_candidategenes)
head(present)


#Genes present in the ASE dataset
phenyl_genes[ phenyl_genes$AXILLARISv4 %in% unique(present_candidategenes), ]
#Genes NOT present in the ASE dataset
phenyl_genes[ !phenyl_genes$AXILLARISv4 %in% unique(present_candidategenes), ]

phenyl_genes

################################
# write informative SNPs and their allele counts out for each of the candidate genes
#################################

#add gene name to list
dim(candidate_loci)
phenyl_genes

candidate_loci$name <- phenyl_genes[match(candidate_loci$gene, phenyl_genes$AXILLARISv4),1]

head(candidate_loci)
dim(candidate_loci)
#filter out SNVs which have 0 coverage
candidate_loci[round(candidate_loci$axmean) + round(candidate_loci$exsmean) != 0,]

candidate_loci <- candidate_loci[round(candidate_loci$axmean) + round(candidate_loci$exsmean) != 0,]
dim(candidate_loci)


#########################
##  CREATE INPUT     ####
## allele counts     ####
#########################


#############################
## SET WHICH READS TO DO ASE WITH
################################


head(candidate_loci)


names(candidate_loci)
lapply(candidate_loci, class)
candidate_loci$scaffold <- as.character(candidate_loci$scaffold)
lapply(candidate_loci, class)


#create candidate_loci for only one gene
#candidate_loci <- candidate_loci[which(candidate_loci$gene == "Peaxi162Scf00050g00423"),]



############################
### ASE Read input data ####
############################

#define counts from which replicate to test
ca.count <- round(candidate_loci$axmean)
ce.count <- round(candidate_loci$exsmean)

atLeast1ReadSubv <- (ca.count+ce.count)>0
table(atLeast1ReadSubv)

ca.counts <- ca.count[atLeast1ReadSubv]
ce.counts <- ce.count[atLeast1ReadSubv]


ca.count
ce.count

#########################
##  CREATE Data frame ###
##                   ####
#########################

mySNVs <- GRanges(
  seqnames = candidate_loci$scaffold,
  ranges = IRanges(start = candidate_loci$pos, width = 1),
  aseID = as.character(candidate_loci$gene),
  allele1= candidate_loci$baseAx, 
  allele2= candidate_loci$baseExs
)
mySNVs


#create dataframe with this 
#how many snps per gene 

unique(present_candidategenes)

#list used not loci!
d= NULL
for (genei in unique(present_candidategenes)){
  print(genei)
  print(length(grep(genei, candidate_loci$gene)))
  snpcount <- length(grep(genei, candidate_loci$gene))
  d = rbind(d, data.frame(genei, snpcount))
}

d
#filter out zeros: 

d_good <- d[d$snpcount != 0,] 
d_good


#create list of names for each snp with format: GENE:##nth SNP in gene## 
ASEnames = NULL
for(i in 1:nrow(d_good)){
  for (e in 1:d[i,2]){
    padded <- sprintf("%02d",e)
    file_name<-paste(d[i,1],".",padded,sep='')
    #names1[e] = file_name    
    ASEnames = rbind(ASEnames, data.frame(file_name))
  }
}

ASEnames <- lapply(ASEnames, as.character)
ASEnames

names(mySNVs) <- ASEnames$file_name

names(mySNVs)

length(ca.count)

########
# specify rho
#######

#set to default right now a bit lower than in the publication of MBASED (0.004)
rho = 0.002

mySample <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts= matrix(ca.count, 
                              ncol = 1, 
                              dimnames=list(names(mySNVs), 'mySample')
    ), 
    lociAllele2Counts=matrix(ce.count, 
                             ncol = 1, 
                             dimnames= list(names(mySNVs), 'mySample')
    ), 
    lociCountsDispersions=matrix(
      rep(rho, length(ca.count)), 
      ncol= 1, 
      dimnames= list(names(mySNVs), 'mySample')
    )
  ), 
  rowRanges = mySNVs
)
traceback()
mySNVs$aseID


#############################
#############################
##  RUN ASE TEST         ####
##  optimal results with ####
##  at least             ####
##  1 mio iterations     ####
#############################
#############################


ASEtest <- runMBASED(ASESummarizedExperiment = mySample, 
                     isPhased = TRUE, 
                     numSim = 10^6, #10^4 for shorter check
                     BPPARAM = SerialParam())



###########################
##  Inspect results     ####
############################

t <- summarizeASEResults_1s(ASEtest)

ASE_2018f1axExslimb<- t$geneOutput

ASE_2018f1axExslimb$gene <- rownames(ASE_2018f1axExslimb)

ASE_2018f1axExslimb$name <- phenyl_genes[match(ASE_2018f1axExslimb$gene, phenyl_genes$AXILLARISv4),1]

head(ASE_2018f1axExslimb)

###correct ASE p-values per multiple comparisons
#Must correct p-values for multiple comparisons with BH (Bonferonni-Holm) test, as per manual p15
#"Note that p-values are not adjusted for multiple hypothesis testing, and the users should 
#carry out such an adjustment themselves, e.g. by employing the utilities in the multtest package."
# https://www.bioconductor.org/packages/release/bioc/manuals/MBASED/man/MBASED.pdf

head(ASE_2018f1axExslimb)
#p.adjust(p, method = p.adjust.methods, n = length(p))
ASE_2018f1axExslimb$pValueASEcorrected <- p.adjust(ASE_2018f1axExslimb$pValueASE, method = "BH", n = length(ASE_2018f1axExslimb$pValueASE))
ASE_2018f1axExslimb$pValueHeterogeneitycorrected <- p.adjust(ASE_2018f1axExslimb$pValueHeterogeneity, method = "BH", n = length(ASE_2018f1axExslimb$pValueHeterogeneity))

#which genes have a MAF over 0.7?
ASE_2018f1axExslimb[which(ASE_2018f1axExslimb$majorAlleleFrequency > 0.7),]

write.table(ASE_2018f1axExslimb, "ASE_2018f1axexs_ASE.txt")



#what are the real counts: ASE towards which side?
ASEmean_2018f1axExslimb_counts <- candidate_loci[which(candidate_loci$gene %in% rownames(ASE_2018f1axExslimb)),c(1,17:19)]
ASEmean_2018f1axExslimb_counts <- candidate_loci[which(candidate_loci$gene %in% rownames(ASE_2018f1axExslimb)), ]

head(ASEmean_2018f1axExslimb_counts)
tail(ASEmean_2018f1axExslimb_counts)

#look at single genes
ASEmean_2018f1axExslimb_counts[which(ASEmean_2018f1axExslimb_counts$name == "DPL"),]
ASEmean_2018f1axExslimb_counts[which(ASEmean_2018f1axExslimb_counts$name == "PH4"),]



#per gene, sum up all SNP site counts
SNPsumup <- aggregate(cbind(axmean, exsmean) ~ gene, data=ASEmean_2018f1axExslimb_counts, FUN=sum)
SNPsumup_wreps <- aggregate(cbind(rep1Ax, rep1Exs, rep2Ax, rep2Exs, rep3Ax, rep3Exs, axmean, exsmean) ~ gene, data=ASEmean_2018f1axExslimb_counts, FUN=sum)


#####
# get mean counts per site and its position
#####

head(candidate_loci)
snp_table_ASE <- candidate_loci[,c(1,4, 5, 6, 7, 17,18, 19)]
head(snp_table_ASE)

#get counts per replicate and position: 

snp_rep_ASE <- candidate_loci[, c(1,4,5,6,7,8,9,11,12,14,15,19)]
snp_rep_ASE

SNPsumup$snp_count <- table(ASEmean_2018f1axExslimb_counts$gene)

dim(SNPsumup)
head(SNPsumup)

SNPsumup$name <- phenyl_genes[match(SNPsumup$gene, phenyl_genes$AXILLARISv4),1]

SNPsumup

write.table(SNPsumup, "ASE_2018f1axexs_snpsumup.txt")


#### package data for figure making ####
SNPsumup$snp_count <- table(ASEmean_2018f1axExslimb_counts$gene)

dim(SNPsumup)
head(SNPsumup)

AXExs_SNPsum <- DataFrame()
AXExs_SNPsum_reps <- DataFrame()

for(i in 1:nrow(SNPsumup)) {
  row <- SNPsumup[i,]
  gene_id <- SNPsumup[i,]$gene
  gene_name <- phenyl_genes[which(phenyl_genes$AXILLARISv4 == gene_id),]$Name[1]
  row$name <- gene_name
  AXExs_SNPsum = rbind(AXExs_SNPsum,row)
}
for(i in 1:nrow(SNPsumup_wreps)) {
  row <- SNPsumup_wreps[i,]
  gene_id <- SNPsumup_wreps[i,]$gene
  gene_name <- phenyl_genes[which(phenyl_genes$AXILLARISv4 == gene_id),]$Name[1]
  row$name <- gene_name
  AXExs_SNPsum_reps = rbind(AXExs_SNPsum_reps,row)
}


AXExs_SNPsum
#adding the gene names for usefulness
AXExs_SNPsum$name <- phenyl_genes[match(AXExs_SNPsum$gene, phenyl_genes$AXILLARISv4),1]

AXExs_SNPsum_reps #for graphs w reps
AXExs_SNPsum_reps$name <- phenyl_genes[match(AXExs_SNPsum_reps$gene, phenyl_genes$AXILLARISv4),1]
head(AXExs_SNPsum_reps)

#### making figures #####

#first calculate the standard error for each gene/species
ASE <- as.data.frame(AXExs_SNPsum_reps)

#calculate the standard error of the mean for the axillaris and exserta replicates
st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

ASE$ax_se <- apply(cbind(ASE$rep1Ax, ASE$rep2Ax, ASE$rep3Ax),MARGIN = 1,FUN=st.err)
ASE$exs_se <- apply(cbind(ASE$rep1Exs, ASE$rep2Exs, ASE$rep3Exs),MARGIN = 1,FUN=st.err)

ASE.means <- ASE[,c(1,8:10)]
ASE.ses <- ASE[,c(1,11:12,10)]

library(reshape)
library(plyr)

ASE.mean <- melt(ASE.means, id=c("gene","name"))
ASE.mean$variable <- revalue(ASE.mean$variable, c("axmean"="axillaris", "exsmean"="exserta"))
names(ASE.mean)[3] <- c("species")
names(ASE.mean)[4] <- c("mean")


ASE.se <- melt(ASE.ses, id=c("gene","name"))
ASE.se$variable <- revalue(ASE.se$variable, c("ax_se"="axillaris", "exs_se"="exserta"))
names(ASE.se)[3] <- c("species")
names(ASE.se)[4] <- c("se")
ASE.forfig <- merge(ASE.mean, ASE.se, by=c("gene", "name", "species"))

#add back in the # of SNPs per locus
ASE.forfig <- merge(ASE.forfig, SNPsumup[,c(1,4)], by=c("gene"))
head(ASE.forfig)
#sort by gene name
ASE.forfig <- ASE.forfig[order(ASE.forfig$name),]

write.csv(ASE.forfig, "ASE_axexs_st4_forfig_flavgenes.csv", row.names = F)

