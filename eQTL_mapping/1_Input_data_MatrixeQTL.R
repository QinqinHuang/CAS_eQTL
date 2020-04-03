#----------------------------------------------
# 2018-05-04
# Prepare input of MatrixeQTL for each condition.
#  Expression: probe-level or gene-level

# Note that MatrixeQTL requires all input files 
# to have matched columns (sample names).

# Number of PEER factors to include in the model
num_PEER <- 10
#----------------------------------------------
# For genotype data and covariates, generate a dataset for all 135 individuals.
# Load into MatrixeQTL.
# For each condition, load expression data for overlapping samples.
# And subset individuals in genotype and covariate ME objects.

library(foreach)
library(doMC)
#registerDoMC(cores = 8)
options(stringsAsFactors = F)

setwd("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/")

# Useful file names
# Prefix of plink binary files
geno_plink_fileprefix <- "/projects/qinqinhuang/CAS/Genotype_Data/clean_data/cas_michigan_imp_filtered_135"
# File name of 10 genotype PCs:
pc_filename <- "/projects/qinqinhuang/CAS/Genotype_Data/clean_data/PCs_10_cas_michigan_imp_filtered_135indiv.txt"

# Prefix of expression file names (Inverse normal transformed datasets)
expr_fileprefix <- "/projects/qinqinhuang/CAS/Expression_Data/clean_data/Rank_normal_transformation_"

# Prefix of the folder where PEER results are under:
peerdirprefix <- "/projects/qinqinhuang/CAS/Expression_Data/process_data/limma_exc19_badquality_normexp_quantile_log2_normalised_separately_2018-05-04/PEER_3genoPCs_"


#----- Overlapping samples -----
# Read sample info for microarray dataset
sampleqc <- read.table("/projects/qinqinhuang/CAS/Expression_Data/process_data/Sample_QC_information.txt", header = T)
sampleqc_filtered <- sampleqc[which(sampleqc$exclude != "badquality"),]
# Overlapping samples: 494 (135 individuals)
sampleqc_overlapping <- sampleqc_filtered[which(sampleqc_filtered$geno == TRUE),]
# Sort by subject ID
sampleqc_overlapping <- sampleqc_overlapping[order(sampleqc_overlapping$subjid),]
#overlap_ind <- unique(sampleqc_overlapping$subjid)


#----- Generate genotype dataset -----
# 135 overlapping individuals 
# In CAS plink binary files, family ID is the same with within-family ID for each kid.

# Filter SNPs with MAF <10% in these 135 individuals
# 4,326,483 SNPs
system(paste0("plink1.9 --bfile ",geno_plink_fileprefix," --maf 0.1 --allow-no-sex --make-bed --out cas_imp_135ind_filtered_maf10"))

# We want to count the dosage of the alternative allele in HRC reference.
# Generate a file with two column,
# the 1st is SNP ID and the 2nd is the alternative allele in HRC pannel.
# The SNP ID is "chr:pos"
system("tail -n +2 ~/Downloaded_data/HRC_r1.1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab | cut -f 1 > hrc_chr.txt;
tail -n +2 ~/Downloaded_data/HRC_r1.1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab | cut -f 2 > hrc_pos.txt;
tail -n +2 ~/Downloaded_data/HRC_r1.1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab | cut -f 5 > hrc_alter.txt;
paste -d : hrc_chr.txt hrc_pos.txt > hrc_variant_id.txt;
paste hrc_variant_id.txt hrc_alter.txt > HRC_alt_allele.txt;
rm hrc_*")

# SNP dosage
system("plink1.9 --bfile cas_imp_135ind_filtered_maf10 --recode A-transpose --recode-allele HRC_alt_allele.txt --out TEMP1")
# SNP genotype
system("cut -f 2,7- TEMP1.traw > SNP_135indiv.txt")
# SNP location
system("awk '{print$2,$1,$4}' TEMP1.traw > snpsloc.txt")
# Which allele is counted in the dataset
system(paste0("cut -f 1-6 TEMP1.traw > SNP_counted_allele.txt"))
# Clean-up
system("rm TEMP*")


#----- Gene expression and covariates -----
# For four conditions
# Use gene or probe level data
foreach(celltype = rep(c("m","m","t","t"), 2), treatment = rep(c("Ctrl","LPS","Ctrl","PHA"), 2), 
        whichdata = rep(c("gene","probe"), each = 4)) %dopar% {
  # Create a folder
  dircond <- paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_",whichdata)
  if(!file.exists(dircond)) {dir.create(dircond)}
  setwd(dircond)
  
  # Samples in this condition
  curr_sam <- sampleqc_overlapping[which(sampleqc_overlapping$cell == celltype &
                                           sampleqc_overlapping$treatment == treatment),]
  # Sort the subject ID
  curr_sam <- curr_sam[order(curr_sam$subjid),]
  
  #----- Generate gene expression datasets -----
  # Read normalised gene expression data
  exprall <- read.table(paste0(expr_fileprefix,whichdata,"_expression_",celltype,"_",treatment,".txt"), header = T)
  # Keep samples that have genotype data
  expr <- exprall[,c(1,match(curr_sam$SampleID, colnames(exprall)))]
  if(!identical(colnames(expr)[-1], curr_sam$SampleID)) {
    cat("  *Error - wrong order of samples in expression file. \n")
    return(NULL)
  }
  cat("  *", ncol(expr)-1, "samples in", celltype, treatment, "\n")
  
  # Use the subject ID (FID_IID) as column names instead of sample IDs.
  colnames(expr)[-1] <- paste0(curr_sam$subjid,"_",curr_sam$subjid)
  
  #------ Generate covariate datasets -----
  # Gender, first 3 genotype PCs, and PEER factors
  # Load PEER factors for this condition
  peerfactor <- readRDS(paste0(peerdirprefix,whichdata,"_level/PEER_factors_",celltype,"_",treatment,"_",nrow(curr_sam),"_samples.rds"))
  peerfactor <- peerfactor[,c("sex",paste0("PC",1:3),paste0("PEER",1:num_PEER))]
  # Replace sample id by subject id
  peerfactor <- as.data.frame(peerfactor)
  peerfactor$subjid <- curr_sam$subjid[match(rownames(peerfactor),curr_sam$SampleID)]
  peerfactor <- peerfactor[order(peerfactor$subjid),]
  rownames(peerfactor) <- paste0(peerfactor$subjid, "_", peerfactor$subjid)
  peerfactor$subjid <- NULL
  # Samples are in columns
  peerfactor <- as.data.frame(t(peerfactor))
  peerfactor$covariates <- rownames(peerfactor)
  peerfactor <- peerfactor[,c("covariates",setdiff(colnames(peerfactor), "covariates"))]
  
  # Confirm that the columns are matched in gene expression and covariate
  if(!identical(colnames(expr)[-1], colnames(peerfactor)[-1])) {
    cat("  * Error - expr and covariates don't have matched columns! \n")
    return(NULL)
  }
  
  # Save data
  write.table(expr, paste0("gene_expression_",celltype,"_",treatment,".txt"), quote = F, sep = "\t", row.names = F)
  write.table(peerfactor, paste0("Covariates_",celltype,"_",treatment,".txt"), quote = F, sep = "\t", row.names = F)
  
  return(NULL)
}





