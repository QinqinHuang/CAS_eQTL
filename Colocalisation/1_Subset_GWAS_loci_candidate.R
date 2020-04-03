#----------------------------------------------
# 2018-07-25
# Last GWAS added: 2018-09-18
#
# Considering only significant eQTLs from the 
# initial unconditional eQTL scan. 
#
# Aims:
# 1. Find GWAS loci that are likely to share a 
# causal variant with eQTL signals:
# Select all signficant GWAS SNPs that passed
# 1e-5, or 5e-8, and look for overlap with 
# significant eQTL SNPs in our dataset.
gwas_thre <- 1e-5
#
# 2. If any, prepare data for coloc test
# Window is defined as +/- 200kb around the 
# top SNP in each eQTL locus. 
windowsize_oneside <- 200
#
# Also generate files for locuszoom plot.
# Note that SNP ID is: chr[chr]:[pos]
#
#----------------------------------------------
# GWAS summary statistics downloaded from GWAS catalog are under: 
downloaded_gwas_folder <- "/projects/qinqinhuang/CAS/Analysis/Downloaded_GWAS_Summary_data_clean/"

library(foreach)
library(doMC)
#registerDoMC(cores = 20)
library(coloc)
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

# My function to keep errors
writelog <- function(newline) {
  write(newline, "/projects/qinqinhuang/CAS/Analysis/Colocalisation/Log_subset_coloc_input.txt", append = T)
}

# Working directory
setwd("/projects/qinqinhuang/CAS/Analysis/Colocalisation/coloc_input/")

# Read phenotypes that have GWAS summary statistics
pheno_list <- read.table("../GWAS_phenotypes_files.txt", header = T, sep = "\t", quote = "")

# Make sure all GWAS files are available
all(pheno_list$file_name %in% list.files(downloaded_gwas_folder))

# Each study will have a subdirectory, the name of the folder:
pheno_list$subdir <- gsub(".txt", "", pheno_list$file_name)

# Load info of SNPs tested in CAS
CAS_SNPs <- read.table("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/SNP_pos_rsID_counted_freq_freq1kG.txt", header = T)
# Add new SNP ID: "chr[chr]:[pos]"
CAS_SNPs$SNPID <- paste0("chr", CAS_SNPs$snps)

colnames(CAS_SNPs)[c(which(colnames(CAS_SNPs) == "COUNTED"), which(colnames(CAS_SNPs) == "ALT"), which(colnames(CAS_SNPs) == "freq_counted"))] <- c("effect_allele","other_allele","eaf")


#----- Load eQTL data -----
# Load significant eQTLs from the initial unconditional eQTL scan
eQTL_stat_sig <- load_cis_eQTLs()

# Load all eQTL summary statistics
eQTL_stat <- load_MatrixeQTL_out()

# Add a new column for SNP ID: "chr[chr]:[pos]"
for(tt in 1:4) {
  eQTL_stat_sig[[tt]]$SNPID <- paste0("chr",eQTL_stat_sig[[tt]]$snps)
  eQTL_stat[[tt]]$FDR <- NULL
  eQTL_stat[[tt]]$SNPID <- paste0("chr",eQTL_stat[[tt]]$snps)
}


#------ Candidate loci for coloc -----
# Go through all diseases 
nothing <- foreach(ii = 1:nrow(pheno_list), .combine = rbind) %do% {
  setwd("/projects/qinqinhuang/CAS/Analysis/Colocalisation/coloc_input/")
  
  # If folder exists
  if(file.exists(pheno_list$subdir[ii])) {
    writelog("Folder exists!"); writelog("---")
    return(NULL)
  }
  
  # Create a new subdirectory
  dir.create(pheno_list$subdir[ii])
  setwd(pheno_list$subdir[ii])
  
  # Load gwas summary statistics
  writelog(paste0("  *GWAS phenotype: ", pheno_list$Disease[ii]))
  gwas_stat <- read.table(paste0(downloaded_gwas_folder,pheno_list$file_name[ii]), header = T, sep = "\t")
  writelog(paste0("  *Number of SNPs tested on autosomes: ", nrow(gwas_stat)/1e6, " million."))
  
  # New SNP ID: "chr[chr]:[pos]"
  gwas_stat$SNPID <- paste0("chr",gwas_stat$chr,":",gwas_stat$position)
  writelog(paste0("  *among which ", length(intersect(gwas_stat$SNPID, CAS_SNPs$SNPID))/1e6, " million SNPs were tested in our eQTL analysis."))
  
  # Significant GWAS SNPs
  gwas_stat_sig <- gwas_stat[which(gwas_stat$pval <= gwas_thre),]
  writelog(paste0("  *GWAS on ", pheno_list$Disease[ii], " had ", nrow(gwas_stat_sig), " SNPs on autosomes significant at ",gwas_thre," significance level."))
  writelog(paste0("  *among which ", length(intersect(gwas_stat_sig$SNPID, CAS_SNPs$SNPID)), " significant SNPs were tested in our eQTL analysis."))
  
  # Look for overlap with our significant eQTLs
  overlap_eQTL <- foreach(tt = 1:4, .combine = rbind) %dopar% {
    dd <- eQTL_stat_sig[[tt]]
    dd_overlap <- dd[which(dd$SNPID %in% gwas_stat_sig$SNPID),]
    if(nrow(dd_overlap) == 0) {return(NULL)}
    dd_overlap$condition <- names(eQTL_stat_sig)[tt]
    return(dd_overlap)
  }
  if(is.null(overlap_eQTL)) {
    writelog(paste0("#### No overlapping eQTL signals. ###")); writelog("---")
    return(NULL)
  }
  
  # eGenes that we are going to test in coloc
  test_eGene_list <- unique(overlap_eQTL[,c("condition","gene")])
  rownames(test_eGene_list) <- 1:nrow(test_eGene_list)
  
  # Go through all possible eGenes and prepare the data for testing
  test_eGene_list_topeSNP <- foreach(gg = 1:nrow(test_eGene_list), .combine = rbind) %dopar% {
    # eGene
    mygene <- test_eGene_list$gene[gg]
    
    #--- eQTL ---
    # Summary statistics for this gene
    curr_eQTL <- eQTL_stat[[test_eGene_list$condition[gg]]]
    curr_eQTL <- curr_eQTL[which(curr_eQTL$gene == mygene),]
    colnames(curr_eQTL)[match("pvalue",colnames(curr_eQTL))] <- "pval"
    
    # Top eSNP
    alltopSNP <- curr_eQTL[which(curr_eQTL$pval == min(curr_eQTL$pval)),"SNPID"]
    # If multiple eSNPs in perfect LD had the min pval, choose the one with the lowest pval in GWAS.
    if(any(alltopSNP %in% gwas_stat_sig$SNPID)) {
      tempdd <- gwas_stat_sig[which(gwas_stat_sig$SNPID %in% alltopSNP),]
      tempdd <- tempdd[order(tempdd$pval),]
      topeSNP <- tempdd$SNPID[1]
    } else {
      topeSNP <- alltopSNP[1]
    }
    # SNP ID for file name
    topeSNP_fn <- gsub(":", "_", topeSNP)
    
    # Chromosome
    chr <- CAS_SNPs[which(CAS_SNPs$SNPID == topeSNP),]$chr
    # Genomic location
    topeSNPpos <- CAS_SNPs[which(CAS_SNPs$SNPID == topeSNP),]$position
    
    
    #--- eQTL summary statistics within the test window ---
    ## Save all summary data within this 2MB window ##
    write.table(curr_eQTL[,c("SNPID","pval")], paste0("LocusZoom_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_eQTL_2Mb.txt"), row.names = F, sep = "\t", quote = F)
    
    # SNP info
    curr_eQTL <- merge(curr_eQTL, CAS_SNPs[,c("SNPID","chr","position","rsid","effect_allele","other_allele","eaf")], by = "SNPID")
    
    # Subset SNPs within +/- 200kb window
    curr_eQTL <- curr_eQTL[which(curr_eQTL$position >= topeSNPpos-windowsize_oneside*1000 & curr_eQTL$position <= topeSNPpos+windowsize_oneside*1000),]
    
    # se of beta estimates
    curr_eQTL$se <- curr_eQTL$beta/curr_eQTL$statistic
    
    ## Save the data for locuszoom plotting - all SNPs tested in eQTL analysis in +/-200kb window ##
    write.table(curr_eQTL[,c("SNPID","pval")], paste0("LocusZoom_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_eQTL_all.txt"), row.names = F, sep = "\t", quote = F)
    
    
    #--- GWAS ---
    # Summary statistics in this window
    curr_gwas <- gwas_stat[which(gwas_stat$chr == chr),]
    curr_gwas <- curr_gwas[which(curr_gwas$position >= topeSNPpos-windowsize_oneside*1000 & curr_gwas$position <= topeSNPpos+windowsize_oneside*1000),]
    
    ## Save the data for locuszoom plotting - all SNPs tested in GWAS ##
    write.table(curr_gwas[,c("SNPID","pval")], paste0("LocusZoom_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_GWAS_all.txt"), row.names = F, sep = "\t", quote = F)
    
    
    #--- Overlapping SNPs are tested in coloc ---
    # eQTL summary stat
    curr_eQTL_overlap <- curr_eQTL[which(curr_eQTL$SNPID %in% curr_gwas$SNPID),]
    # Order by genomic location
    curr_eQTL_overlap <- curr_eQTL_overlap[order(curr_eQTL_overlap$position),]
    rownames(curr_eQTL_overlap) <- 1:nrow(curr_eQTL_overlap)
    
    curr_eQTL_overlap <- curr_eQTL_overlap[,c("SNPID","rsid","pval","beta","se","effect_allele","other_allele","eaf")]

    ## Save the data for locuszoom plotting and coloc test ##
    write.table(curr_eQTL_overlap, paste0("coloc_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_eQTL_overlapping.txt"), row.names = F, sep = "\t", quote = F)
    
    
    # GWAS summary stat
    curr_gwas_overlap <- curr_gwas[match(curr_eQTL_overlap$SNPID, curr_gwas$SNPID),]
    rownames(curr_gwas_overlap) <- 1:nrow(curr_gwas_overlap)
    curr_gwas_overlap$snps <- NULL
    
    ## Save the data for locuszoom plotting and coloc test ##
    write.table(curr_gwas_overlap, paste0("coloc_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_GWAS_overlapping.txt"), row.names = F, sep = "\t", quote = F)
    
    
    # Top eSNP
    curr_gene_topeSNP <- test_eGene_list[gg,]
    curr_gene_topeSNP$topeSNP <- topeSNP
    return(curr_gene_topeSNP)
  }
  
  write.table(test_eGene_list_topeSNP, "Tested_eGene_list.txt", row.names = F, sep = "\t", quote = F)
  
  writelog("---")
  
  return(NULL)
}




