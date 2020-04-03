#----------------------------------------------
# 2018-05-14
# Run conditional analysis to identify 
# independent eQTL signals.
# Nominal thresholds determined by eigenMT-BH.

# Which dataset used - gene or probe level data.
whichdata <- "gene"   #"probe"
#----------------------------------------------
# Output: 1. Top SNPs in forward stage
#         2. All significant SNPs in backward stage using leave-out-one model

library(foreach)
library(doMC)
#registerDoMC(cores = 20)
library(MatrixEQTL)
options(stringsAsFactors = F)

# Working directory
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
nothing <- foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA")) %do% {
  # Directory of this condition
  setwd(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_",whichdata))
  cat("  * Conditional analysis on", celltype, treatment, ", using", whichdata, "level data...\n")
  if(!file.exists("conditional_analysis")) { dir.create("conditional_analysis") }
  
  # Focus on significant eGenes
  eAsso <- read.table("sigeAsso_eigenMTBH.txt", header = T)
  eAsso <- eAsso[order(abs(eAsso$statistic), decreasing = T),]
  # Top SNP for eGenes
  topSNP <- eAsso[!duplicated(eAsso$gene),]
  
  # Read nominal thresholds
  nomthre <- read.table("Nominal_thresholds.txt", header = T)
  nomthre <- nomthre[,c("gene","eigenMT.BH")]
  
  # Load gene expression data
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 2000   # read file in pieces of 2000 rows
  gene$LoadFile(paste0("gene_expression_",celltype,"_",treatment,".txt"))
  
  # Subset genotype data
  snps <- snpsall$Clone()
  snps$ColumnSubsample(match(gene$columnNames, snps$columnNames))
  
  # Double check column names
  flag <- F
  if(!identical(snps$columnNames, gene$columnNames)) {flag <- T}
  if(flag) {write("  *Error: sample names are not in the same order!", "error.txt", append = T); return(NULL)}
  
  
  #----- Forward stage of conditional analysis for each gene -----
  ind_forward <- foreach(mygene = topSNP$gene, .combine = rbind) %dopar% {
    # Gene expression data
    gene_mygene <- gene$Clone()
    gene_mygene$RowReorder(which(gene$GetAllRowNames() == mygene))
    
    # Nominal threshold for this gene, which is used to decide whether we identify an additional independent eSNP.
    nomthre_mygene <- nomthre[which(nomthre$gene == mygene),2]
    
    # A file to keep covariates and SNPs will be added
    # File name
    cov_filename <- paste0("./conditional_analysis/Cov_corrected_SNPs_",gsub(";","_",mygene),".txt")
    system(paste0("cp Covariates_*.txt ",cov_filename))
    
    #--- Iterations of forward stage ---
    # All significant independent eSNPs for this gene in forward stage
    ind_forward_mygene <- topSNP[which(topSNP$gene == mygene), ]
    ind_forward_mygene$Rank <- 0
    
    # Repeat conditional analysis until no additional significant signal can be identified
    while(ind_forward_mygene$pvalue[which.max(ind_forward_mygene$Rank)] <= nomthre_mygene) {
      # SNP to be adjusted
      SNPcov <- ind_forward_mygene$snps[which.max(ind_forward_mygene$Rank)]
      # Add its genotype data to covariate file
      write(paste(c(SNPcov, snps$FindRow(SNPcov)$row), collapse = "\t"), file = cov_filename, append = T)
      
      # Load the new covariate file
      cvrt = SlicedData$new()
      cvrt$fileDelimiter = "\t"
      cvrt$fileOmitCharacters = "NA"
      cvrt$fileSkipRows = 1
      cvrt$fileSkipColumns = 1
      cvrt$LoadFile(cov_filename)
      
      # Double check rows and columns
      if(!all(cvrt$nCols() == gene_mygene$nCols(), cvrt$nCols() == snps$nCols(), 
              cvrt$nRows() == nrow(ind_forward_mygene)+14, 
              all(cvrt$GetAllRowNames()[-1:-14] %in% ind_forward_mygene$snps) )) {
        writelog("Error - Coverate file!", IdxRep, SS, MAF)
      }
      
      # Run linear regression
      me_curr <- Matrix_eQTL_main(
        snps = snps,
        gene = gene_mygene,
        cvrt = cvrt,
        pvOutputThreshold = 0,
        output_file_name.cis = NULL,
        pvOutputThreshold.cis = 1,
        snpspos = snpspos,
        genepos = geneloc,
        cisDist = 1000000,
        output_file_name = NULL,
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        verbose = TRUE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = FALSE)
      
      # New top SNP
      cis_curr <- me_curr$cis$eqtls
      cis_curr <- cis_curr[order(abs(cis_curr$statistic), decreasing = T),]
      top_curr <- cis_curr[1, -5]
      
      # Add the number of iterations
      top_curr$Rank <- max(ind_forward_mygene$Rank) + 1
      
      # Add this new top SNP into forward list,
      # even if it didn't pass significance threshold.
      ind_forward_mygene <- rbind(ind_forward_mygene, top_curr)
      
    } # End of forward while loop
    
    # The while loop stopped because the last SNP was not significant
    # Remove it
    if(which(ind_forward_mygene$pvalue > nomthre_mygene) == which.max(ind_forward_mygene$Rank)) {
      ind_forward_mygene <- ind_forward_mygene[-which.max(ind_forward_mygene$Rank), ]
    } else {write(paste0("Error - last SNP in forward stage, gene ",mygene), "error.txt", append = T)}
    
    # Remove covariate file
    system(paste0("rm ", cov_filename))
    
    return(ind_forward_mygene)
  } # End of forward stage for all eGenes
  
  #-- Save forward results --
  write.table(ind_forward, file = "./conditional_analysis/Forward_top_SNPs.txt", quote = F, sep = "\t", row.names = F)
  
  
  #----- Backward stage of conditional analysis for each gene -----
  ind_backward <- foreach(mygene = topSNP$gene, .combine = rbind) %dopar% {
    # Gene expression data
    gene_mygene <- gene$Clone()
    gene_mygene$RowReorder(which(gene$GetAllRowNames() == mygene))
    
    # Nominal threshold
    nomthre_mygene <- nomthre[which(nomthre$gene == mygene), 2]
    
    # Top SNPs identified in forward stage
    ind_forward_mygene <- ind_forward[which(ind_forward$gene == mygene),]
    
    #-- If only one independent SNP, no additional analysis should be done in backward stage --
    if(nrow(ind_forward_mygene) == 1) {
      ind_backward_mygene <- eAsso[which(eAsso$gene == mygene),]
      # Rank
      ind_backward_mygene$Rank <- 0
      return(ind_backward_mygene)
    }
    
    # If >1 independent SNPs, run leave-out-one model
    ind_backward_mygene <- foreach(jj = 1:nrow(ind_forward_mygene), .combine = rbind) %do% {
      # Covairate file name
      cov_filename <- paste0("./conditional_analysis/Cov_corrected_SNPs_",gsub(";","_",mygene),".txt")
      system(paste0("cp Covariates_*.txt ",cov_filename))
      
      # SNP to be adjusted - all except one
      for(SNPcov in ind_forward_mygene$snps[-jj]) {
        # Add its genotype data to covariate file
        write(paste(c(SNPcov, snps$FindRow(SNPcov)$row), collapse = "\t"), file = cov_filename, append = T)
      }
      
      # Load the new covariate file
      cvrt = SlicedData$new()
      cvrt$fileDelimiter = "\t"
      cvrt$fileOmitCharacters = "NA"
      cvrt$fileSkipRows = 1
      cvrt$fileSkipColumns = 1
      cvrt$LoadFile(cov_filename)
      
      # Double check rows and columns
      if(!all(cvrt$nCols() == gene_mygene$nCols(), cvrt$nRows() == nrow(ind_forward_mygene)-1+14, 
              all(cvrt$GetAllRowNames()[-1:-14] %in% ind_forward_mygene$snps[-jj]) )) {
        write(paste0("Error in coverate file for gene ",mygene), "error.txt", append = T)
      }
      
      # Run linear regression
      me_curr <- Matrix_eQTL_main(
        snps = snps,
        gene = gene_mygene,
        cvrt = cvrt,
        pvOutputThreshold = 0,
        output_file_name.cis = NULL,
        pvOutputThreshold.cis = 1,
        snpspos = snpspos,
        genepos = geneloc,
        cisDist = 1000000,
        output_file_name = NULL,
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        verbose = TRUE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = FALSE)
      
      # Only return significant SNPs
      cis_curr <- me_curr$cis$eqtls
      cis_curr <- cis_curr[order(abs(cis_curr$statistic), decreasing = T), -5]
      cis_curr <- cis_curr[which(cis_curr$pvalue <= nomthre_mygene),]
      
      # If no significant SNPs
      if(nrow(cis_curr) == 0) { return(NULL) }
      
      # Rank
      cis_curr$Rank <- ind_forward_mygene$Rank[jj]
      
      # Remove covariate file
      system(paste0("rm ", cov_filename))
      
      return(cis_curr)
    } # End of the loop for all top SNPs identified in forward stage
    
    return(ind_backward_mygene)
  } # End of backward stage for all eGenes
  
  #-- Save backward results --
  write.table(ind_backward, file = "./conditional_analysis/Backward_significant_eSNPs.txt", quote = F, sep = "\t", row.names = F)
  
  
  return(NULL) 
}

