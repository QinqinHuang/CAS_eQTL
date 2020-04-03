#----------------------------------------------
# 2018-09-04
# Run MatrixeQTL for trans eQTL identification.
#
# Save only relatively significant results.
nominalpval <- 1e-5
# 
# Which dataset used - gene or probe level data.
whichdata <- "gene"   #"probe"
#----------------------------------------------
library(MatrixEQTL)
library(foreach)
library(doMC)
#registerDoMC(cores = 8)
options(stringsAsFactors = F)

# Directory:
# /projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/

# Gene location
geneloc <- read.table("/projects/qinqinhuang/CAS/Expression_Data/clean_data/Gene_location.txt", header = T)

# Load genotype data for 135 individuals 
# When running MatrixeQTL, subset the genotype
snpsall = SlicedData$new()
snpsall$fileDelimiter = "\t"
snpsall$fileOmitCharacters = "NA"
snpsall$fileSkipRows = 1
snpsall$fileSkipColumns = 1
snpsall$fileSliceSize = 2000   # read file in pieces of 2000 rows
snpsall$LoadFile("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/SNP_135indiv.txt")


# Four conditions
foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA")) %do% {
  cat("  * Working on", celltype, treatment, ", using", whichdata, "level data...\n")
  
  # Directory of this condition
  dirname <- paste0("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/",celltype,"_",treatment,"_",whichdata)
  if(!file.exists(dirname)) {dir.create(dirname)} 
  setwd(dirname)
  
  #----- Martrix eQTL ------
  # Load gene expression data
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 2000   # read file in pieces of 2000 rows
  gene$LoadFile(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_",whichdata,"/gene_expression_",celltype,"_",treatment,".txt"))
  
  # Load covariates
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = "\t"
  cvrt$fileOmitCharacters = "NA"
  cvrt$fileSkipRows = 1
  cvrt$fileSkipColumns = 1
  cvrt$LoadFile(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_",whichdata,"/Covariates_",celltype,"_",treatment,".txt"))

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
    pvOutputThreshold = nominalpval,
    output_file_name = NULL,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  
  # Save the results of MatrixeQTL
  saveRDS(me, file = paste0("MatrixeQTL_out_",celltype,"_",treatment,".RDS"))
  
  
  #--- Add location for SNPs and genes ---- 
  asso <- me$all$eqtls
  
  # Location of each gene
  asso <- merge(asso, geneloc[,1:3], by.x = "gene", by.y = "geneid")
  colnames(asso)[which(colnames(asso) == "Chr")] <- "gene_chr"
  colnames(asso)[which(colnames(asso) == "left")] <- "gene_TSS"
  
  # Location of each SNP
  asso$snp_chr <- as.numeric(sapply(asso$snps, function(x) {unlist(strsplit(x, split = ":"))[1]}))
  asso$snp_pos <- as.numeric(sapply(asso$snps, function(x) {unlist(strsplit(x, split = ":"))[2]}))
  
  # Label tests where SNP and gene were not on the same chromosome
  asso$label_chr <- ifelse(asso$gene_chr == asso$snp_chr, "same", "different")
  asso$label_posdist <- (asso$snp_pos - asso$gene_TSS)/1e6
  asso$label_posdist[which(asso$label_chr == "different")] <- NA
  
  # Labe tests that were cis
  asso$label_cis <- ifelse((asso$label_chr == "same" & asso$label_posdist <= 1), "cis", "not")
  
  # Sort by significance
  asso <- asso[order(asso$pvalue),]
  rownames(asso) <- 1:nrow(asso[order(asso$pvalue),])
  write.table(asso, paste0("All_tests_sig_",as.character(nominalpval),".txt"), quote = F, sep = "\t", row.names = F)

  rm(me)
  return(NULL) 
}








