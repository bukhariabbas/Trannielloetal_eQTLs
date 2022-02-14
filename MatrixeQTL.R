# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Modified by Syed Abbas Bukhari
rm(list = ls())
library(MatrixEQTL)
#sink(args[6])
setwd("./")
args = commandArgs(trailingOnly=TRUE)
## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');


## Settings

#arg1 = choice of linear model, can be [modelANOVA, modelLINEAR, or modelLINEAR_CROSS]
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
#useModel = modelLINEAR_CROSS; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
if(args[1] == "1")
{
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
}else {
  if(args[1] == "2")
    useModel = modelANOVA
  else
  	if(args[1] == "3") 
    		useModel = modelLINEAR_CROSS
	else
		print("Use 1, 2 and 3 for ModelLinear, modelANOVA and modelLINEAR_CROSS")
}

# args2 is the choice of counts and genotype file
SNP_file_name = args[2]
snps_location_file_name = "genotype-master.filtered.snp-loc"

# args3 is the choice of gene expression counts file
expression_file_name = args[3]
gene_location_file_name = "gene_loc.txt"

# arg4 is the corresponding covariate file name, include variables of interest
# Set to character() for no covariates
covariates_file_name = args[4]

# Distance parameter for calling cis-eQTLs
cisDist = as.numeric(args[5])
# Output file name
output_file_name_cis = paste(getwd(),"/model",args[1],args[4],args[5],"cis_eQTL.txt", sep = "");
output_file_name_tra = paste(getwd(),"/model",args[1],args[4], "tra_eQTL.txt", sep = "");

file.create(output_file_name_cis)
file.create(output_file_name_tra)

# Only associations significant at this level will be saved
#arg3 = threshold for calling significance in cis.
#arg4 = threshold for calling significance in trans.
pvOutputThreshold_cis = as.numeric(args[6])
pvOutputThreshold_tra = as.numeric(args[7])

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(me)

#args6 = label to save the workspace.
save.image(file=paste("workspace","model",args[1],"distance",args[5], args[4],Sys.Date(), sep = "_"))

#unlink(args[6])
