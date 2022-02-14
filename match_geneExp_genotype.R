# This script designs input files for matrix eQTL software.
# Given a gene expression and genotype matrix of individuals along with their sample information.
# I am designing two types of matrixeQTL runs (1) run Matrix eQTL separately for each colony; 
# (2) given slopes are intact in the first step (means not changing between experimental groups); run 
# a giant model for soldiers and foragers separately.

rm(list = ls())
#setting up working directory
setwd("/Volumes/GoogleDrive/My Drive/Postdoc_BellLab_GNDP/HB_eQTL/Variants_Arian")

#loading required libraries
stopifnot(require("edgeR"))
stopifnot(require(plyr))
stopifnot(require(limma))
stopifnot(require(vcfR))

#load up gene expression count files. and create covarite files. 

targets = readTargets("targets.txt")

counts = read.delim("counts_group_postfilter.csv", header = T)
count_names = strsplit2(colnames(counts), "_")
count_names = paste(count_names[,1], count_names[,2], count_names[,3], count_names[,4], sep = "_")

#matching colnames of counts matrix with targets filenames.
temp = strsplit2(targets$file, "_")
temp = paste(temp[,1], temp[,2], temp[,3], temp[,4], sep = "_")
targets$names = temp
intersect(targets$names, count_names)

#selecting counts matrix column info from target file.
counts_sampleIDs = targets[which(targets$names %in% count_names),c("names", "seqID", "group","Colony", "MDS1_Aggression")]

colnames(counts) = counts_sampleIDs$seqID

#taking log cpm of counts .. MeQTL works well with normally distributed data.
counts = cpm(counts, log = T)

#write.table(counts, file = "counts_group_postfilter.cpm.txt", quote = F, sep = "\t")

counts_sampleIDs = counts_sampleIDs[,c(-1)]
counts_sampleIDs$group = ifelse(counts_sampleIDs$group == "soldier", 1, 0)
counts_sampleIDs$colony_group = paste(counts_sampleIDs$Colony, counts_sampleIDs$group, sep = ".")

# Reading genotype file
genotype = read.delim("genotype-master.filtered.replaced.txt", sep = "\t")

######################
library("GenomicRanges")
library("limma")

variants = readRDS("/Volumes/GoogleDrive/My Drive/Postdoc_BellLab_GNDP/HB_eQTL/Variants_Arian/Arian_variants/significant_markers_colony_aggression.RDS")

variants = paste(variants@seqnames, variants@ranges@start, sep = "_")


length(intersect(variants, rownames(genotype)))

common_variants = intersect(variants, variants_eQTL)

write.table(common_variants, file = "eQTL-variants_DEGs_overlap_soldiers_forager.txt", quote = F, sep = "\t", row.names = F)

#####################
# finding genetic relatedness between members of a colony using svd. 
# removing snps with zero variance and centering the snps before doing svd.
genotype = genotype[which(rownames(genotype) %in% variants),]
subDir = paste("MeQTL_input_variants_avalos-etal", Sys.Date(), sep = ".")
ifelse(!dir.exists(file.path(getwd(), subDir)), dir.create(file.path(getwd(), subDir)), FALSE)
setwd(file.path(getwd(), subDir))


library(matrixStats)
library(useful) # to see the corner of a matrix.
counts_sampleIDs$svd_colony = 0
for(i in 1:length(unique(counts_sampleIDs$Colony)))
{
  covariates = counts_sampleIDs[which(counts_sampleIDs$Colony == unique(counts_sampleIDs$Colony)[i]),]
  genotype_covariates = genotype[,c(covariates$seqID)]
  # 1. removing SNPs with zero variance; 2. centering SNPs (feeding transposed matrix where SNPs are on the columns)
  var = rowVars(as.matrix(genotype_covariates))
  temp = genotype_covariates[which(var > 0),]
  temp = as.matrix(t(temp))
  temp = scale(temp, center = T, scale = F)
  #temp = sweep(temp, 2, colMeans(temp), "-")
  # https://genomicsclass.github.io/book/pages/pca_svd.html [see this for centering matrix for svd]
  temp = svd(temp)
  counts_sampleIDs$svd_colony[which(counts_sampleIDs$Colony == unique(counts_sampleIDs$Colony)[i])] = temp$u[,1]
}
  
##########################  
for(i in 1:length(unique(counts_sampleIDs$colony_group)))
{
  covariates = counts_sampleIDs[which(counts_sampleIDs$colony_group == unique(counts_sampleIDs$colony_group)[i]),]
  counts_covariates = counts[, c(covariates$seqID)]
  genotype_covariates = genotype[,c(covariates$seqID)]
  var = rowVars(as.matrix(genotype_covariates))
  genotype_covariates = genotype_covariates[which(var > 0),]
  temp = t(covariates[,c("svd_colony")])
  colnames(temp) = covariates$seqID
  write.table(temp, file = paste("covariates", "colony_group", unique(counts_sampleIDs$colony_group)[i],  Sys.Date(), sep = "_"), quote = F, sep = "\t")
  write.table(counts_covariates, file = paste("counts", "colony_group", unique(counts_sampleIDs$colony_group)[i], Sys.Date(), sep = "_"), quote = F, sep = "\t")
  write.table(genotype_covariates, file = paste("genotype", "colony_group", unique(counts_sampleIDs$colony_group)[i], Sys.Date(), sep = "_"), quote = F, sep = "\t")
}



# Separate run for soldiers and forager.
# before that lets also do svd for everyone in a single matrix. 
genotype = genotype[,c(counts_sampleIDs$seqID)]
counts_sampleIDs$svd = 0
var = rowVars(as.matrix(genotype))
genotype = genotype[which(var > 0),]
temp = as.matrix(t(genotype))
temp = scale(temp, center = T, scale = F)
temp = svd(temp)
counts_sampleIDs$svd = temp$u[,1]

write.table(genotype, file = "all_genotype.txt", quote = F, sep = "\t")
write.table(counts, file= "all_counts.txt", quote = F, sep = "\t")
write.table(counts_sampleIDs, file = "all_covariates.txt", quote = F, sep = "\t")
#soldiers
covariates_soldiers = counts_sampleIDs[which(counts_sampleIDs$group == "1"),]
counts_soldiers = counts[,c(covariates_soldiers$seqID)]
genotype_soldiers = genotype[,c(covariates_soldiers$seqID)]
# removing genotypes with 0 variance
var = rowVars(as.matrix(genotype_soldiers))
genotype_soldiers = genotype_soldiers[which(var > 0),]

# setting up colony dummy variables to block for colonies.
cov_mat = t(as.data.frame(model.matrix(~(factor(covariates_soldiers$Colony)))))[2:9,]
rownames(cov_mat) = paste("Colony", unique(covariates_soldiers$Colony)[2:length(unique(covariates_soldiers$Colony))], sep = "_")
colnames(cov_mat) = covariates_soldiers$seqID
cov_mat = rbind(cov_mat, covariates_soldiers$svd)
rownames(cov_mat)[9] = "Genetic_relatedness"

write.table(cov_mat, file = paste("covariates", "soldiers", Sys.Date(), sep = "_"), quote = F, sep = "\t")
write.table(counts_soldiers, file = paste("counts", "soldiers", Sys.Date(), sep = "_"), quote = F, sep = "\t")
write.table(genotype_soldiers, file = paste("genotype", "soldiers",Sys.Date(), sep = "_"), quote = F, sep = "\t")


# foragers
covariates_foragers = counts_sampleIDs[which(counts_sampleIDs$group == "0"),]
counts_foragers = counts[,c(covariates_foragers$seqID)]
genotype_foragers = genotype[,c(covariates_foragers$seqID)]
# removing genotypes with 0 variance
var = rowVars(as.matrix(genotype_foragers))
genotype_foragers = genotype_foragers[which(var > 0),]

# setting up colony dummy variables to block for colonies.
cov_mat = t(as.data.frame(model.matrix(~(factor(covariates_foragers$Colony)))))[2:9,]
rownames(cov_mat) = paste("Colony", unique(covariates_foragers$Colony)[2:length(unique(covariates_foragers$Colony))], sep = "_")
colnames(cov_mat) = covariates_foragers$seqID
cov_mat = rbind(cov_mat, covariates_foragers$svd)
rownames(cov_mat)[9] = "Genetic_relatedness"

write.table(cov_mat, file = paste("covariates", "foragers", Sys.Date(), sep = "_"), quote = F, sep = "\t")
write.table(counts_foragers, file = paste("counts", "foragers", Sys.Date(), sep = "_"), quote = F, sep = "\t")
write.table(genotype_foragers, file = paste("genotype", "foragers",Sys.Date(), sep = "_"), quote = F, sep = "\t")

