#----------------------------------------------
# 2018-10-04
# Aim: using Mendelian randomisation analysis 
# to test causal associations between cis-eGenes 
# and immune related diseases.
#   Use cis-eQTLs as instrumental variables.
#
# IV: cis-eQTLs (r2 0.1)
# GWAS: I downloaded summary stat
#
# 2018-10-12
# In GWAS data, allele values of some variants 
# are not A/T/C/G, exclude those.
#----------------------------------------------
library(foreach)
library(doMC)#; registerDoMC(cores = 10)
library(TwoSampleMR)
library(ggplot2)
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

# Working directory
setwd("/projects/qinqinhuang/CAS/Analysis/Mendelian_Randomisation_ciseGene_diseases/")


# ----- Load data -----
# Load significant eQTLs
eQTLs_all <- load_cis_eQTLs()

# Sample size
# Exposure
eQTL_samplesize <- data.frame(condition = c("m_resting","m_LPS","t_resting","t_PHA"),
                              samplesize = c(116, 125, 126, 127))
# Outcome
gwas_phenotype <- read.table("/projects/qinqinhuang/CAS/Analysis/Colocalisation/GWAS_phenotypes_files.txt", header = T, sep = "\t", quote = "")
gwas_phenotype$samplesize <- gwas_phenotype$n_cases + gwas_phenotype$n_controls
# Keep those with reported beta and se
gwas_phenotype <- gwas_phenotype[which(gwas_phenotype$coloc_input == "se"),]

# Read SNP info
CAS_SNPs <- read.table("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/SNP_pos_rsID_counted_freq_freq1kG.txt", header = T)
colnames(CAS_SNPs)[c(6,7,9)] <- c("effect_allele","other_allele","eaf")

# Gene location
geneloc <- read.table("/projects/qinqinhuang/CAS/Expression_Data/clean_data/Gene_location.txt", header = T)
# ---------- 


# ----- MR -----
# Go through each disease
nothing <- foreach(gwasfilename = gwas_phenotype$file_name) %do% {
  # Whether all output files are generated
  flag <- F
  for(ii in 1:4) {
    if(!file.exists(paste0("./output/",gsub(".txt","",gwasfilename),"_Celltype_",eQTL_samplesize$condition[ii],".RDS"))) {
      flag <- T
      cat("  Missing: ", gwasfilename, eQTL_samplesize$condition[ii], "\n")
    }
  }
  if(!flag) {return(NULL)}
  
  cat("  * Working on GWAS:", gwasfilename, "\n")
  
  # Read gwas summary stat
  gwas_full <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/Downloaded_GWAS_Summary_data_clean/",gwasfilename), header = T, sep = "\t")
  if(!"snps" %in% colnames(gwas_full)) {
    gwas_full$snps <- paste0(gwas_full$chr,":",gwas_full$position)
  }
  
  # Go through four conditions
  nothing2 <- foreach(ii = 1:4) %do% {
    cat("  * Cell type condition:", eQTL_samplesize$condition[ii], "\n")
    
    # Cis-eGenes in this condition
    eQTLs <- eQTLs_all[[eQTL_samplesize$condition[ii]]]
    eGenelist <- eQTLs[!duplicated(eQTLs$gene),]
    
    # Go through each cis-eGene
    MR_results <- foreach(mygene = eGenelist$gene, .combine = c) %dopar% {
      cat(mygene, "\t")
      
      #-- Exposure data --
      # Significant cis-eQTLs
      sigeQTL <- eQTLs[which(eQTLs$gene == mygene),]
      sigeQTL <- merge(sigeQTL, CAS_SNPs[,c("snps","rsid","effect_allele","other_allele","eaf")], by = "snps")
      
      # Calculate se
      sigeQTL$se <- sigeQTL$beta / sigeQTL$statistic
      
      # Sample size
      sigeQTL$samplesize <- eQTL_samplesize$samplesize[ii]
      
      # Update some column names
      colnames(sigeQTL)[which(colnames(sigeQTL)=="rsid")] <- "SNP"
      colnames(sigeQTL)[which(colnames(sigeQTL)=="gene")] <- "Phenotype"
      colnames(sigeQTL)[which(colnames(sigeQTL)=="pvalue")] <- "pval"
      
      #-- Outcome data --
      gwas_stat <- gwas_full[which(gwas_full$snps %in% sigeQTL$snps),]
      # Remove variants whose allele valus are not A/T/C/G
      gwas_stat <- gwas_stat[which(gwas_stat$effect_allele %in% c("A","T","C","G") & gwas_stat$other_allele %in% c("A","T","C","G")),]
      if(nrow(gwas_stat) == 0) {
        #write(paste0(" No overlapping SNPs for ",gwasfilename," and ",mygene," in ",eQTL_samplesize$condition[ii]), file = "Log_MR_inc_trans.txt", append = T)
        cat(paste0(" No overlapping SNPs for ",mygene," in ",eQTL_samplesize$condition[ii]," and ",gwasfilename,"\n"))
        return(NULL)
      }
      
      # Overlapping
      sigeQTL_overlapping <- sigeQTL[match(gwas_stat$snps, sigeQTL$snps),]
      
      # RSID
      gwas_stat$SNP <- sigeQTL_overlapping$SNP
      
      # Phenotype
      gwas_stat$Phenotype <- gsub(".txt","",gwasfilename)
      
      # Frequency
      # Use the allele frequency based on the GWAS if provided;
      # some studies also provide frequency calculated from 1000 Genomes project.
      if("AF_effect" %in% colnames(gwas_stat)) {
        gwas_stat$eaf <- gwas_stat$AF_effect
      } else if("AF_effect_1000G" %in% colnames(gwas_stat)) {
        gwas_stat$eaf <- gwas_stat$AF_effect_1000G
      }
      
      # Sample size
      if(!"samplesize" %in% colnames(gwas_stat)) {
        gwas_stat$samplesize <- gwas_phenotype$samplesize[which(gwas_phenotype$file_name == gwasfilename)]
      }
      
      # Reformat
      exposure_dat <- format_data(sigeQTL_overlapping, type = "exposure")
      outcome_dat <- format_data(gwas_stat, type = "outcome")
      
      # Harmonisation
      dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
      
      # Clumping data
      dat_clump <- clump_data(dat[which(dat$mr_keep),], clump_r2 = 0.1)
      
      # MR
      res <- mr(dat_clump)
      
      # Direction
      directiontest <- directionality_test(dat_clump)
      
      ## If there are multiple instrumental variables ##
      if(nrow(dat_clump) > 2) {
        p1 <- mr_scatter_plot(res, dat_clump)
        ggsave(paste0("./Scatter_plot/",mygene,"_",eQTL_samplesize$condition[ii],"_",gwasfilename,".pdf"), p1[[1]], height = 5, width = 5)
      }
      
      #-- Save input data --
      mylist <- list(phenotype = c(gwasfilename, mygene, eQTL_samplesize$condition[ii]),
                     input = dat_clump, results = res, direction = directiontest)
      return(mylist)
      
    } # End of all cis-eGenes 
    
    saveRDS(MR_results, paste0("./output/",gsub(".txt","",gwasfilename),"_Celltype_",eQTL_samplesize$condition[ii],".RDS"))
    cat("  \n")
    
    return(NULL)
  } # End of four conditions
  
  return(NULL)
} # End of all GWASs




