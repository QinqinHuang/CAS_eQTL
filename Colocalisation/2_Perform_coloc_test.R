#----------------------------------------------
# 2018-07-05
# Run coloc to test whether the same cuasal 
# variant is underlying eQTL and GWAS signals.
#
# 2018-09-18
# Add GWAS on allergic rhinitis and allergic 
# sensitisation.
# And use an updated reQTL list.
#
# 2018-10-02
# Use an updated reQTL list.
# Also add back eGenes whose test with the top 
# eSNP was removed because of high LD with the
# other top eSNP.
#
# 2018-10-09
# Add GWAS on childhood/adult onset asthma.
# 
windowsize_oneside <- 200
#----------------------------------------------
library(foreach)
library(coloc)
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

# My function to keep errors
writelog <- function(newline) {
  write(newline, "/projects/qinqinhuang/CAS/Analysis/Colocalisation/Log_run_coloc.txt", append = T)
}

# List of GWAS diseases
setwd("/projects/qinqinhuang/CAS/Analysis/Colocalisation/")
pheno_list <- read.table("GWAS_phenotypes_files.txt", header = T, sep = "\t")

# Each study has a subdirectory
pheno_list$subdir <- gsub(".txt", "", pheno_list$file_name)
all(pheno_list$subdir %in% list.files("./coloc_input/"))

# SNPs tested in eQTL data and rs ID
CAS_SNPs <- read.table("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/SNP_tested_info_maf.txt", header = T)

# Add new SNP ID: "chr[chr]:[pos]"
CAS_SNPs$SNPID <- paste0("chr", CAS_SNPs$snps)


#------------------------------------------------------------------------
# Run coloc
#------------------------------------------------------------------------
# Go through each disease phenotype
allresults <- foreach(ii = c(1:33,35:nrow(pheno_list)), .combine = rbind) %do% {
  
  # Working directory
  setwd("/projects/qinqinhuang/CAS/Analysis/Colocalisation/coloc_input/")
  setwd(pheno_list$subdir[ii])
  
  # If job has finished
  if(file.exists("coloc_results.txt")) {
    results_curr_pheno <- read.table("coloc_results.txt", header = T)
    results_curr_pheno$Study <- pheno_list$subdir[ii]
    return(results_curr_pheno)
  }
  
  writelog(paste0("  Working on: ",pheno_list$subdir[ii]))
  
  # If no overlapping with eQTL signals
  if(!file.exists("Tested_eGene_list.txt")) {
    writelog(paste0("  * No coloc tests for ",pheno_list$Disease[ii]))
    return(NULL)
  }
  
  # Load eGenes that we are going to test in coloc
  test_eGene_list <- read.table("Tested_eGene_list.txt", header = T)
  
  # GWAS sample size & proportion of cases
  gwas_ss <- pheno_list$n_cases[ii] + pheno_list$n_controls[ii]
  prop_cases <- pheno_list$n_cases[ii]/gwas_ss
  
  
  # Go through each possible eGene
  results_curr_pheno <- foreach(gg = 1:nrow(test_eGene_list), .combine = rbind) %dopar% {
    # eGene
    mygene <- test_eGene_list$gene[gg]
    # Top eSNP
    topeSNP <- test_eGene_list$topeSNP[gg]
    chr <- gsub("chr", "", unlist(strsplit(topeSNP, split = ":"))[1])
    topeSNPpos <- unlist(strsplit(topeSNP, split = ":"))[2]
    topeSNP_fn <- gsub(":","_",topeSNP)
    
    
    # Read GWAS summary statistics
    # Overlapping SNPs
    curr_gwas_overlap <- read.table(paste0("coloc_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_GWAS_overlapping.txt"), header = T, sep = "\t")
    # All SNPs in the window
    curr_gwas <- read.table(paste0("LocusZoom_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_GWAS_all.txt"), header = T, sep = "\t")
    
    # Read eQTL summary statistics
    # Overlapping SNPs
    curr_eQTL_overlap <- read.table(paste0("coloc_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_eQTL_overlapping.txt"), header = T)
    # All SNPs in the window
    curr_eQTL <- read.table(paste0("LocusZoom_",test_eGene_list$condition[gg],"_",mygene,"_",topeSNP_fn,"_eQTL_all.txt"), header = T)
    
    # SNPs should be in the same order
    if(!identical(curr_eQTL_overlap$SNPID, curr_gwas_overlap$SNPID)) {
      writelog(paste0(" Error - SNPs are not in the same order - gene ", mygene))
      return(NULL)
    }
    
    #--- Run coloc ---
    # Are allele frequency available in downloaded GWAS summary statistics
    # (It is the same to use frequency or 1-frequency)
    if("AF_effect" %in% colnames(curr_gwas_overlap)) {
      if(gg == 1) {
        writelog(paste0("  * Use the allele frequency provided by the study: ",pheno_list$Disease[ii]))
      }
      currAF <- curr_gwas_overlap$AF_effect
    } else if("AF_effect_1000G" %in% colnames(curr_gwas_overlap)) {  
      if(gg == 1) {
        writelog(paste0("  * Use the allele frequency (1000 Genomes) provided by the study: ",pheno_list$Disease[ii]))
      }
      currAF <- curr_gwas_overlap$AF_effect_1000G
    } else {
      # Otherwise, get allele frequency from 1000 Genomes phase3 European population
      currAF <- CAS_SNPs$EUR_alt_freq[match(curr_gwas_overlap$SNPID, CAS_SNPs$SNPID)]
    }
    # Minor allele frequency
    currAF[which(currAF>0.5)] <- 1-currAF[which(currAF>0.5)]
    
    # Are beta & variance of beta available in GWAS summary 
    flag_var <- F
    if(pheno_list$coloc_input[ii] == "se" & "beta" %in% colnames(curr_gwas_overlap) & "se" %in% colnames(curr_gwas_overlap)) {
      # If any variance (se) is NA, coloc doesnt work
      if(sum(is.na(curr_gwas_overlap$se)) == 0) {
        # If any variance (se) is 0, coloc doesn't work
        if(all(curr_gwas_overlap$se >0)) {flag_var <- T} 
      }
    }
    if(flag_var) {
      gwasdataset <- list(beta = curr_gwas_overlap$beta, varbeta = (curr_gwas_overlap$se)^2,
                          MAF = currAF,
                          N = gwas_ss, type = "cc", s = prop_cases)
    } else {
      gwasdataset <- list(pvalues = curr_gwas_overlap$pval, 
                          MAF = currAF,
                          N = gwas_ss, type = "cc", s = prop_cases)
    }
    # eQTL summary stat
    eQTLdataset <- list(beta = curr_eQTL_overlap$beta, varbeta = (curr_eQTL_overlap$se)^2,
                        sdY = 1, type = "quant")
    
    # Run coloc
    colocresult_0 <- coloc.abf(dataset1 = gwasdataset, dataset2 = eQTLdataset,
                               p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
    
    colocresult_1 <- coloc.abf(dataset1 = gwasdataset, dataset2 = eQTLdataset,
                               p1 = 1e-04, p2 = 1e-04, p12 = 1e-06)
    
    returndd <- data.frame(gene = mygene, condition = test_eGene_list$condition[gg],
                           topeSNP = topeSNP, 
                           topeSNP_tested_in_coloc = ifelse(topeSNP %in% curr_eQTL_overlap$SNPID, T, F),
                           min_gwas_pval = min(curr_gwas$pval),
                           min_gwas_pval_tested = min(curr_gwas_overlap$pval),
                           Chr = chr, topeSNPpos = topeSNPpos,
                           n_SNP_eQTL = nrow(curr_eQTL),
                           n_SNP_GWAS = nrow(curr_gwas),
                           n_SNP_overlap = nrow(curr_eQTL_overlap),
                           p12 = c(1e-05, 1e-06),
                           n_SNP_tested = c(colocresult_0$summary[1], colocresult_1$summary[1]), 
                           PP0 = c(colocresult_0$summary[2], colocresult_1$summary[2]),
                           PP1 = c(colocresult_0$summary[3], colocresult_1$summary[3]),
                           PP2 = c(colocresult_0$summary[4], colocresult_1$summary[4]),
                           PP3 = c(colocresult_0$summary[5], colocresult_1$summary[5]),
                           PP4 = c(colocresult_0$summary[6], colocresult_1$summary[6]) )
    
    return(returndd)
  }
  
  # Save results
  write.table(results_curr_pheno, "coloc_results.txt", quote = F, sep = "\t", row.names = F)
  
  
  results_curr_pheno$Study <- pheno_list$subdir[ii]
  
  writelog("---")
  return(results_curr_pheno)
}


setwd("/projects/qinqinhuang/CAS/Analysis/Colocalisation/coloc_results/")

# PP3+PP4
allresults$sum_pp3_pp4 <- allresults$PP3 + allresults$PP4
# PP4/PP3
allresults$ratio_pp4_pp3 <- allresults$PP4 / allresults$PP3
# Order by PP4
allresults <- allresults[order(allresults$p12, allresults$PP4, decreasing = T),]

# Save
write.table(allresults, "Coloc_output_added_childhood_onset_asthma.txt", quote = F, sep = "\t", row.names = F)



#------------------------------------------------------------------------
# Plot to compare default prior and a lower prior
#------------------------------------------------------------------------
library(ggplot2)
library(reshape)
options(stringsAsFactors = F)
setwd("/projects/qinqinhuang/CAS/Analysis/Colocalisation/coloc_results/")

# Read coloc results
allresults <- read.table("Coloc_output_added_childhood_onset_asthma.txt", header = T)

# Keep tests where sum of PP3 and PP4 ≥0.8
ddplot <- allresults[which(allresults$sum_pp3_pp4 >= 0.8), c("n_SNP_tested","PP3","PP4","sum_pp3_pp4","ratio_pp4_pp3","p12")]

# y axis: pp4/(pp3+pp4)
ddplot$posterior <- ddplot$PP4/ddplot$sum_pp3_pp4
# The expected value
ddplot$prior <- ddplot$p12/(ddplot$p12 + (ddplot$n_SNP_tested - 1)*1e-8)

ddplot_2 <- melt(ddplot, measure.vars = c("posterior","prior"))

# Plot all data
pp <- ggplot(ddplot_2, aes(n_SNP_tested, value, color = variable)) + 
  geom_point(size = 0.2) + geom_smooth() +
  facet_wrap(~p12) +
  theme(legend.title = element_blank()) +
  xlab("Number of SNPs tested in a window") + ylab("PP4/(PP3+PP4)")
ggsave("Coloc_compare_p12.pdf", pp, width = 6, height = 4)

# >50, <1500 SNPs
pp <- ggplot(ddplot_2[which(ddplot_2$n_SNP_tested >= 50 & ddplot_2$n_SNP_tested <= 1500),], aes(n_SNP_tested, value, color = variable)) + 
  geom_point(size = 0.2) + geom_smooth() +
  facet_wrap(~p12) +
  theme(legend.title = element_blank()) +
  xlab("Number of SNPs tested in a window") + ylab("PP4/(PP3+PP4)")
ggsave("Coloc_compare_p12_more_than_50_SNPs_smallerthan1500.pdf", pp, width = 6, height = 4)



