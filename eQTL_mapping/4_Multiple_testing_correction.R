#----------------------------------------------
# 2018-05-07
# Multiple testing correction: eigenMT-BH
#
# Which dataset used - gene or probe level data.
whichdata <- "gene"   #"probe"
#----------------------------------------------
library(foreach)
library(doMC)
#registerDoMC(cores = 20)
options(stringsAsFactors = F)
source("/projects/qinqinhuang/eQTL_Simulations/Scripts/001_all_functions_for_eQTL_simulations.R")


# Read gene location
if(whichdata == "gene") {
  geneloc_filename <- "~/CAS/Expression_Data/clean_data/Gene_location.txt"
} else if(whichdata == "probe") {
  geneloc_filename <- "~/CAS/Expression_Data/clean_data/Probe_location.txt"
}
geneloc <- read.table(geneloc_filename, header = T)

# Read output of eigenMT: estimated number of independent tests per gene
eigenMTout <- read.table("~/CAS/Analysis/eQTL_mapping/eigenMT/eigenMT_output_gene_level.txt", header = T)
# Probe-level data
if(whichdata == "probe") {
  # Probes for the same gene have the same number of estimated effective tests
  probe_anno <- read.table("~/CAS/Analysis/eQTL_mapping/process_data/Probes_reliable_filtered_19230_GENCODE_anno.txt", header = T)
  eigenMTout <- merge(probe_anno[,c("Probe_Id","Gene_symbol")], eigenMTout, by.x = "Gene_symbol", by.y = "gene")
  eigenMTout <- eigenMTout[,c("Probe_Id","TESTS")]
  colnames(eigenMTout)[1] <- "gene"
}
# Not all genes are tested
cat("  * Genes that were not tested:", setdiff(geneloc[,1], eigenMTout$gene), "\n")

# 4 conditions
foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA")) %do% {
  cat("  * Working on", celltype, treatment, "...\n")
  setwd(paste0("~/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_",whichdata))
  
  if(!file.exists("MTC_process")) {dir.create("MTC_process")}
  
  # Read MatrixeQTL output
  me <- readRDS(paste0("MatrixeQTL_out_",celltype,"_",treatment,".RDS"))
  cis <- me$cis$eqtls
  # Sort by absulute statistics
  cis <- cis[order(abs(cis$statistic), decreasing = T),]
  
  # Top SNP for each gene
  topSNP <- cis[!duplicated(cis$gene), c("gene","snps","pvalue","beta")]
  # Sort by gene names
  topSNP <- topSNP[order(topSNP$gene),]
  rownames(topSNP) <- 1:nrow(topSNP)
  
  #---------- eigenMT-BH ---------
  eMTBH_returnlist <- mul_corr_eigenMT(alltests = cis, eMT_output = eigenMTout, step2_met = "BH")
  sigeAsso_eigenMTBH <- eMTBH_returnlist[[1]]
  # Nominal thresholds
  nominalthresholds <- eMTBH_returnlist[[2]]
  
  
  # Save results
  #write.table(sigeAsso_eigenMTBH, "sigeAsso_eigenMTBH.txt", quote = F, sep = "\t", row.names = F)
  write.table(nominalthresholds, "Nominal_thresholds.txt", quote = F, sep = "\t", row.names = F)
  
  cis_sub <- cis[which(cis$gene %in% sigeAsso_eigenMTBH$gene & cis$snps %in% sigeAsso_eigenMTBH$snps), c("snps","gene","statistic","pvalue","beta")]
  dd <- merge(sigeAsso_eigenMTBH[,1:2], cis_sub, by = c("gene","snps"))
  dd <- dd[order(abs(dd$statistic), decreasing = T),]
  write.table(dd, "sigeAsso_eigenMTBH.txt", quote = F, sep = "\t", row.names = F)
  
  return(NULL) 
}


#* Genes that were not tested: TEKT4P2 BAGE2



