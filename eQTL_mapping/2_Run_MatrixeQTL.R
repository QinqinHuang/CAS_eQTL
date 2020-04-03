#----------------------------------------------
# 2018-05-04
# 1. Run MatrixeQTL for each condition.
# 2. Save MatrixeQTL results, which will be used 
# as eigenMT input files.

# Note that eigenMT works on each chromosome 
# separately, so input files (QTL, GEN, GENPOS) 
# all should be one for each chromosome.

# Which dataset used - gene or probe level data.
whichdata <- "gene"   #"probe"
#----------------------------------------------

library(MatrixEQTL)
library(qqman)
library(foreach)
library(doMC)
#registerDoMC(cores = 4)
options(stringsAsFactors = F)

setwd("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/")

# Load genotype data for 135 individuals 
# When running MatrixeQTL, subset the genotype
snpsall = SlicedData$new()
snpsall$fileDelimiter = "\t"
snpsall$fileOmitCharacters = "NA"
snpsall$fileSkipRows = 1
snpsall$fileSkipColumns = 1
snpsall$fileSliceSize = 2000   # read file in pieces of 2000 rows
snpsall$LoadFile("SNP_135indiv.txt")

# SNP location 
snpspos <- read.table("snpsloc.txt", header = T)

# Read gene location data
if(whichdata == "gene") {
  geneloc_filename <- "/projects/qinqinhuang/CAS/Expression_Data/clean_data/Gene_location.txt"
} else if(whichdata == "probe") {
  geneloc_filename <- "/projects/qinqinhuang/CAS/Expression_Data/clean_data/Probe_location.txt"
}
geneloc <- read.table(geneloc_filename, header = T)

# 4 conditions
foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA")) %do% {
  cat("  * Working on", celltype, treatment, ", using", whichdata, "level data...\n")
  
  # Directory of this condition
  setwd(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_",whichdata))
  
  #----- Martrix eQTL ------
  # Load gene expression data
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 2000   # read file in pieces of 2000 rows
  gene$LoadFile(paste0("gene_expression_",celltype,"_",treatment,".txt"))
  
  # Load covariates
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = "\t"
  cvrt$fileOmitCharacters = "NA"
  cvrt$fileSkipRows = 1
  cvrt$fileSkipColumns = 1
  cvrt$LoadFile(paste0("Covariates_",celltype,"_",treatment,".txt"))
  
  # Subset genotype data
  snps <- snpsall$Clone()
  snps$ColumnSubsample(match(gene$columnNames, snps$columnNames))
  
  # Double check column names
  flag <- F
  if(!identical(snps$columnNames, gene$columnNames)) {flag <- T}
  if(!identical(snps$columnNames, cvrt$columnNames)) {flag <- T}
  if(flag) {cat("  *Error: sample names are not in the same order!\n"); return(NULL)}
  
  # Run MatrixeQTL
  me <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    pvOutputThreshold = 0,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = geneloc,
    cisDist = 1000000,
    pvalue.hist = "qqplot",
    output_file_name = NULL,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  
  # Save the results
  saveRDS(me, file = paste0("MatrixeQTL_out_",celltype,"_",treatment,".RDS"))
  
  #----- QQ plot -----
  png("qqplot.png", width=6*300, height=6*300, res=300)
  plot(me)
  dev.off()
  
  #----- Manhattan plot ----- 
  cis <- me$cis$eqtls
  
  # SNP location data
  cis$chr <- sapply(cis$snps, function(x){unlist(strsplit(x, split = ":"))[1]})
  cis$pos <- sapply(cis$snps, function(x){unlist(strsplit(x, split = ":"))[2]})
  cis$chr <- as.numeric(cis$chr)
  cis$pos <- as.numeric(cis$pos)
  
  jpeg(paste0("Manhattan_plot_cis_eQTL_",celltype,"_",treatment,".jpeg"), width = 12, height = 8, units= "in", res = 300)
  manhattan(cis, chr = "chr", bp = "pos", p = "pvalue", snp = "snps",
            col = c("gray10", "gray60"), chrlabs = NULL, logp = TRUE,
            suggestiveline = FALSE, genomewideline = FALSE)#, highlight = SNPsig)
  dev.off()
  
  #----- eigenMT -----
  # Generate input files.
  if(celltype == "m" & treatment == "Ctrl" & whichdata == "gene") {
    # Each chromosome has a separate file.
    # All cis tests
    cis <- me$cis$eqtls[, c("snps","gene","pvalue")]
    
    # Create a diretory to save eigenMT input
    setwd("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/")
    if(!file.exists("eigenMT/input_QTL")) {dir.create("eigenMT/input_QTL", recursive = T)}
    
    # 22 chromosomes
    ntests <- foreach(chr = 1:22, .combine = "+") %dopar% {
      cat("  Saving MatrixeQTL results for Chr",chr,"...\n")
      # SNPs on this chromosome
      chrsnps <- snpspos[which(snpspos$CHR == chr),]
      # Tests involving these SNPs
      eigenMTinput <- cis[which(cis$snps %in% chrsnps$SNP),]
      write.table(eigenMTinput, paste0("eigenMT/input_QTL/chr",chr,"_MatrixEQTL_all_cis_tests_pval.txt"), quote = F, sep = "\t", row.names = F)
      return(nrow(eigenMTinput))  
    }
    # Check the number of tests
    if(ntests != nrow(cis)) {cat("  * Error - check eigenMT input files! \n")}
  }
  
  rm(me)
  return(NULL) 
}


