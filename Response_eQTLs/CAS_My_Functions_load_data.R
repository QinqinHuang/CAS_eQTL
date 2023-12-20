# Created 2018-06-26
# Last modified 2018-05-10
#
# Script on the server:
# ~/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R
#
# 1. Extract microarray expression data of a gene.
# (1) gene expression normalised separately (IRN, one condition)
# extract_expression(geneID, celltype = "m", treatment = "Ctrl")
#  *Return a data frame with one row; the number of column is 
#   the number of samples + 1; the first column is "geneid".
#   Column names are subject IDs (fid_iid).
#
# (2) gene expression normalised together in a cell type (no IRN, all conditions)
# extract_expression_norm_togethercelltype(geneID)
#  *Return a data frame with one row; the number of column is 
#   the number of samples + 1; the first column is "geneid".
#   Column names are sample IDs.
# 
# (3) gene expression normalised together (no IRN, all conditions)
# extract_expression_norm_alltogether(geneID)
#  *Return a data frame with one row; the number of column is 
#   the number of samples + 1; the first column is "geneid".
#   Column names are sample IDs.
#  
#
# 2. Extract the genotype data of a SNP.
# extract_geno(SNPID, recodedfilename = NULL)
#  *Return a data frame with one row; the number of column is 
#   the number of samples + 1; the first column is "SNP".
#  Default: genotype data for 135 kids
#
# 3. Get the info on a SNP (e.g. counted allele, frequency).
# get_CAS_SNP_info(SNPID)
#  *Return a data frame with one row.
#
#
# 4. Load signficant cis-eQTLs
# load_cis_eQTLs(type = "LD")
#  *Return a list of four sets of cis-eQTLs.
#  If type == "LD", load the data with eSNPs' LD with the corresponding top eSNP (if ≥0.5). 
#
# 5. Load the results of all eQTL tests from Matrix eQTL
# load_MatrixeQTL_out()
#  *Return a list of four sets of full summary statistics.
#
# 6. Load the results of conditional analysis
# load_conditional_cis_eQTLs()
#  *Return a list of four sets of results.
# 
# 7. Independent eQTLs from the backward stage of the conditional analysis (top eSNP per signal)
# Ind_conditional_cis_eQTLs()
#  *Return a table, with a column indicating conditions.
#
# 8. Load response eQTLs
# load_reQTLs()
#  *Resturn a table, with a column indicating cell type. 
#
#
# 9. Extract GWAS summary statistics of a SNP.
# extract_GWAS_stat(targetcolumn, target, gwasfile)
#  *targetcolumn: column name
#  *target: search by it
#  *Return a data frame with one row.
#
#
# 10. Calcualte LD r2 between two SNPs.
# LDr2_twoSNPs(SNP1, SNP2)
#
#
# 11. Check a gene
# checkgene(mygene)
# It does not return any values, print info:
#  - whether it was included in our eQTL analysis;
#  - whether it had significant cis-eQTLs in any conditions;
#  - whether it had significant reQTLs;
#  - whether it had trans-eQTLs.
#
#
# 12. Load trans-eQTLs
# load_trans_eQTLs()
#  
#
#

options(stringsAsFactors = F)


# 1. Extract microarray expression data of a gene.
# (1) gene expression normalised separately (IRN, one condition)
extract_expression <- function(geneID, celltype = "m", treatment = "Ctrl") {
  exprdd <- read.table(pipe(paste0("cat ~/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_gene/gene_expression_",celltype,"_",treatment,".txt | grep ",geneID)), header = F)
  colnames(exprdd) <- read.table(pipe(paste0("head -n 1 ~/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_gene/gene_expression_",celltype,"_",treatment,".txt")), header = F)
  exprdd <- exprdd[which(exprdd$geneid == geneID),]
  return(exprdd)
}

# (2) gene expression normalised together in a cell type (no IRN, all conditions)
extract_expression_norm_togethercelltype <- function(geneID) {
  exprdd <- read.table(pipe(paste0("cat ~/CAS/Expression_Data/process_data/limma_exc19_badquality_normexp_quantile_log2_normalised_together_celltype_2018-05-01/Filtered_normalised_gene_expression_13109_noIRN_4conditions_557samples.txt | grep ",geneID)), header = F)
  colnames(exprdd) <- read.table(pipe("head -n 1 ~/CAS/Expression_Data/process_data/limma_exc19_badquality_normexp_quantile_log2_normalised_together_celltype_2018-05-01/Filtered_normalised_gene_expression_13109_noIRN_4conditions_557samples.txt"), header = F)
  exprdd <- exprdd[which(exprdd$geneid == geneID),]
  return(exprdd)
}

# (3) gene expression normalised together (no IRN, all conditions)
extract_expression_norm_alltogether <- function(geneID) {
  exprdd <- read.table(pipe(paste0("cat ~/CAS/Expression_Data/process_data/limma_exc19_badquality_normexp_quantile_log2_normalised_alltogether_2018-05-01/Filtered_normalised_gene_expression_13109_noIRN_4conditions_557samples.txt | grep ",geneID)), header = F)
  colnames(exprdd) <- read.table(pipe("head -n 1 ~/CAS/Expression_Data/process_data/limma_exc19_badquality_normexp_quantile_log2_normalised_alltogether_2018-05-01/Filtered_normalised_gene_expression_13109_noIRN_4conditions_557samples.txt"), header = F)
  exprdd <- exprdd[which(exprdd$geneid == geneID),]
  return(exprdd)
}


# 2. Get the genotype data of a SNP.
extract_geno <- function(SNPID, recodedfilename = NULL) {
  if(is.null(recodedfilename)) {
    recodedfilename <- "~/CAS/Analysis/eQTL_mapping/SNP_135indiv.txt"
  }
  genodd <- read.table(pipe(paste0("cat ",recodedfilename," | grep ", SNPID)), header = F)
  colnames(genodd) <- read.table(pipe(paste0("head -n 1 ",recodedfilename)), header = F)
  genodd <- genodd[which(genodd$SNP == SNPID),]
  return(genodd)
}

# 3. Get the info on a SNP (e.g. location, counted allele, frequency)
get_CAS_SNP_info <- function(SNPID) {
  dd <- read.table(pipe(paste0("cat ~/CAS/Analysis/eQTL_mapping/SNP_pos_rsID_counted_freq_freq1kG.txt | grep ", SNPID)), header = F, colClasses = c("character","integer","integer","character","character", "character","character","logical","numeric","character", "character","character","numeric","numeric"))
  colnames(dd) <- read.table(pipe("head -n 1 ~/CAS/Analysis/eQTL_mapping/SNP_pos_rsID_counted_freq_freq1kG.txt"), header = F)
  dd <- dd[which(dd$snps == SNPID), ]
  return(dd)
}


# 4. Load significant eQTLs
load_cis_eQTLs <- function(type = "LD") {
  if(type == "LD") {
    return(readRDS("~/CAS/Analysis/eQTL_mapping/Significant_cis_eQTLs_four_conditions_LDtopeSNP.RDS"))
  }
  
  eQTL_m_resting <- read.table("~/CAS/Analysis/eQTL_mapping/m_Ctrl_gene/sigeAsso_eigenMTBH.txt", header = T)
  eQTL_m_LPS <- read.table("~/CAS/Analysis/eQTL_mapping/m_LPS_gene/sigeAsso_eigenMTBH.txt", header = T)
  eQTL_t_resting <- read.table("~/CAS/Analysis/eQTL_mapping/t_Ctrl_gene/sigeAsso_eigenMTBH.txt", header = T)
  eQTL_t_PHA <- read.table("~/CAS/Analysis/eQTL_mapping/t_PHA_gene/sigeAsso_eigenMTBH.txt", header = T)
  
  # Make a list
  eQTL_stat_sig <- list(m_resting = eQTL_m_resting, m_LPS = eQTL_m_LPS, t_resting = eQTL_t_resting, t_PHA = eQTL_t_PHA)
  return(eQTL_stat_sig)
}


# 5. Load the results of all eQTL tests from Matrix eQTL
load_MatrixeQTL_out <- function() {
  me_m_resting <- readRDS("~/CAS/Analysis/eQTL_mapping/m_Ctrl_gene/MatrixeQTL_out_m_Ctrl.RDS")
  cis_m_resting <- me_m_resting$cis$eqtl
  me_m_LPS <- readRDS("~/CAS/Analysis/eQTL_mapping/m_LPS_gene/MatrixeQTL_out_m_LPS.RDS")
  cis_m_LPS <- me_m_LPS$cis$eqtl
  me_t_resting <- readRDS("~/CAS/Analysis/eQTL_mapping/t_Ctrl_gene/MatrixeQTL_out_t_Ctrl.RDS")
  cis_t_resting <- me_t_resting$cis$eqtl
  me_t_PHA <- readRDS("~/CAS/Analysis/eQTL_mapping/t_PHA_gene/MatrixeQTL_out_t_PHA.RDS")
  cis_t_PHA <- me_t_PHA$cis$eqtl
  
  # Make a list
  eQTL_stat <- list(m_resting = cis_m_resting, m_LPS = cis_m_LPS, t_resting = cis_t_resting, t_PHA = cis_t_PHA)
  return(eQTL_stat)
}


# 6. Load the results of conditional analysis (backward stage)
load_conditional_cis_eQTLs <- function() {
  d1 <- read.table("~/CAS/Analysis/eQTL_mapping/m_Ctrl_gene/conditional_analysis/Backward_significant_eSNPs.txt", header = T)
  d2 <- read.table("~/CAS/Analysis/eQTL_mapping/m_LPS_gene/conditional_analysis/Backward_significant_eSNPs.txt", header = T)
  d3 <- read.table("~/CAS/Analysis/eQTL_mapping/t_Ctrl_gene/conditional_analysis/Backward_significant_eSNPs.txt", header = T)
  d4 <- read.table("~/CAS/Analysis/eQTL_mapping/t_PHA_gene/conditional_analysis/Backward_significant_eSNPs.txt", header = T)
  # Return a list
  return(list(m_resting = d1, m_LPS = d2, t_resting = d3, t_PHA = d4))
}

# 7. Independent eQTLs from the backward stage of the conditional analysis (top eSNP per signal)
Ind_conditional_cis_eQTLs <- function() {
  alldd <- load_conditional_cis_eQTLs
  indeQTLs <- data.frame()
  for(ii in 1:4) {
    currdd <- alldd[[ii]]
    currdd <- currdd[order(currdd$pvalue),]
    currdd <- currdd[!duplicated(currdd[,c("gene","Rank")]),]
    rownames(currdd) <- 1:nrow(currdd)
    currdd$condition <- names(alldd)[[ii]]
    indeQTLs <- rbind(indeQTLs, currdd)
  }
  return(indeQTLs)
}

# 8. Load response eQTLs
load_reQTLs <- function() {
  reQTLs_m <- read.table("~/CAS/Analysis/Response_eQTLs_ind_specific_randomeffect/FDR_significant_reQTL_m.txt", header = T)
  reQTLs_m$celltype <- "m"
  reQTLs_t <- read.table("~/CAS/Analysis/Response_eQTLs_ind_specific_randomeffect/FDR_significant_reQTL_t.txt", header = T)
  reQTLs_t$celltype <- "t"
  reQTLs <- rbind(reQTLs_m, reQTLs_t)
  # Keep significant reQTLs
  reQTLs <- reQTLs[which(reQTLs$lmer_perm_fdr == "Significant"), c("celltype","Condition","gene","snps","LDr2_othertop","keep","genotype_beta","Interaction_term_beta","Interaction_term_pval","lmer_genotype_beta","lmer_Interaction_term_beta","lmer_Interaction_term_dpval")]
  rownames(reQTLs) <- 1:nrow(reQTLs)
  return(reQTLs)
}



# 9. Extract GWAS summary statistics of a SNP.
extract_GWAS_stat <- function(targetcolumn, target, gwasfile) {
  GWASstat <- read.table(pipe(paste0("cat ",gwasfile," | grep ",target)), header = F)
  colnames(GWASstat) <- read.table(pipe(paste0("head -n 1 ",gwasfile)), header = F)
  GWASstat <- GWASstat[which(GWASstat[,targetcolumn] == target),]
  return(GWASstat)
}



# 10. Calcualte LD r2 between two SNPs.
LDr2_twoSNPs <- function(SNP1, SNP2) {
  if(is.na(SNP1) | is.na(SNP2)) {cat("Give two SNPs!\n");return(NA)}
  if(identical(SNP1, SNP2)) {return(1)}
  # Get the current directory, and change back to it when plink finishes
  currdir <- getwd()
  setwd("~/CAS/Analysis/eQTL_mapping/")
  write.table(data.frame(snps = c(SNP1, SNP2)), "TEMP_SNP_list", row.names = F, quote = F, col.names = F, sep = "\t")
  system(paste0("plink1.9 --bfile cas_imp_135ind_filtered_maf10 --r2 square --extract TEMP_SNP_list --out TEMP_LD_r2_two_SNPs"))
  d <- read.table("TEMP_LD_r2_two_SNPs.ld", header = F)[2,1]
  system("rm TEMP_SNP_list TEMP_LD_r2_two_SNPs*")
  setwd(currdir)
  return(d)
}


# 11. Check a gene
checkgene <- function(mygene) {
  # Whether it was included in our eQTL analysis
  geneloc <- tryCatch(read.table(pipe(paste0("cat ~/CAS/Expression_Data/clean_data/Gene_location.txt |grep ",mygene)), header = F), error = function(err) NA)
  if(is.na(geneloc)[1]) {
    print("  Gene not included in our analysis!")
    return(NULL)
  }
  geneloc <- geneloc[which(geneloc[,1] == mygene),]
  if(nrow(geneloc) == 0) {
    print("  Gene not included in our analysis!")
    return(NULL)
  }
  print("  Gene was tested.")
  
  # Whether it had significant cis-eQTLs in any conditions
  sigeQTL <- readRDS("~/CAS/Analysis/eQTL_mapping/Significant_cis_eQTLs_four_conditions_LDtopeSNP.RDS")
  
  for(ii in names(sigeQTL)) {
    dd <- sigeQTL[[ii]]
    if(!mygene %in% dd$gene) {
      next
    }
    dd <- dd[which(dd$gene == mygene),]
    print(paste("  * eQTLs in", ii))
    print(dd[which(dd$R2_topeSNP == 1), -1])
  }
  
  # Whether it had significant reQTLs
  for(cc in c("monocytes","Tcells")) {
    ee <- tryCatch(read.csv(pipe(paste0("cat ~/CAS/Analysis/Response_eQTLs_ind_specific_randomeffect/ReQTLs_eGenes_",cc,".csv |grep ",mygene)), header = F), error = function(err) NA)
    if(is.na(ee)[1]) { next }
    ee <- ee[which(ee[,4] == mygene),]
    if(nrow(ee) > 0) {
      print(paste("  * Response eQTLs in", cc))
      print(ee)
    }
  }
  
  # Whether it had trans-eQTLs
  transeQTLs <- read.table("~/CAS/Analysis/Trans_eQTLs/transeQTLs/GenomewideFDR_diffchr.txt", header = T)
  if(mygene %in% transeQTLs$gene) {
    dd <- transeQTLs[which(transeQTLs$gene == mygene),]
    dd$condition <- paste0(dd$celltype, "_", dd$treatment)
    print(paste("  * Trans eQTLs in", unique(dd$condition)))
  }
  
  
  return(NULL)
}


# 12. Load trans-eQTLs
# whichfile: "allexceptcis" or "diffchr"
load_trans_eQTLs <- function(whichfile = "diffchr") {
  return(read.table(paste0("~/CAS/Analysis/Trans_eQTLs/transeQTLs/GenomewideFDR_",whichfile,".txt"), header = T))
}





