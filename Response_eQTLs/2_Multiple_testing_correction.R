#----------------------------------------------
# 2018-08-25
# 2018-09-28: added lmer estimated beta.
# 2018-10-02: use updated LD r2 between top eSNPs.
# Aim: to correct for multiple testing and get a 
# final list of eGenes with significant eQTLs.
#   BH FDR correction
# 
# For eGenes whose top eSNPs in two conditions
# are not independent (LD r2 ≥0.8), test on the
# more significant eQTL would be kept.
#----------------------------------------------
library(foreach)
library(ggplot2)
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

# Working directory
setwd("/projects/qinqinhuang/CAS/Analysis/Response_eQTLs_ind_specific_randomeffect/")
if(!file.exists("Plot_sanitycheck")) {dir.create("Plot_sanitycheck")}


#----- Functions -----
# 1. P-value histogram
plot_pval_histogram <- function(pval, filename) {
  dd <- data.frame(p = pval)
  p <- ggplot(dd, aes(p)) + geom_histogram() + xlab("P-value") + ylab("Frequency") +
    expand_limits(x = c(0,1))
  ggsave(filename, p, width = 3, height = 2.3)
}

# 2. Volcano plot
# x-axis: beta estimates, y-axis: -log10(pval)
plot_volcano <- function(pval, est, significant, filename) {
  dd <- data.frame(y = -log10(pval), x = est, 
                   sig = factor(significant, levels = c("Significant", "Not")))
  dd <- dd[order(dd$y),]
  p <- ggplot(dd, aes(x, y, colour = sig)) + geom_point(shape = 1) +
    scale_color_manual(values = c("red","grey")) +
    theme(legend.title = element_blank()) +
    xlab("Beta estimate for SNP:condition") + ylab("-log10 p-value")
  ggsave(filename, p, width = 3.6, height = 2.3)
}


# Two cell types
nothing <- foreach(celltype = c("m","t")) %do% {
  # Cell type
  cat("  * Working on", celltype, "cells...\n")

  # Interaction test
  dd <- read.table(paste0("Interaction_test_topSNP_perm_pval_",celltype,".txt"), header = T)
  
  # For top eSNPs in two conditions are in high LD, keep the more significant eQTL
  dd_labelled <- foreach(mygene = unique(dd$gene), .combine = rbind) %do% {
    mygene_tests <- dd[which(dd$gene == mygene),]
    # Significant in only one condition
    if(nrow(mygene_tests) == 1) {
      mygene_tests$keep <- "Single"
    } else if(nrow(mygene_tests) == 2) {
      # Independent top SNPs and keep both tests
      if(mygene_tests$LDr2_othertop[1] < 0.8) {
        mygene_tests$keep <- "Independent"
      } else if(mygene_tests$LDr2_othertop[1] >= 0.8) {
        # Keep the more significant condition, and will remove the less significant condition
        mygene_tests$keep <- "Remove"
        mygene_tests$keep[which.min(mygene_tests$pvalue)] <- "Keep"
      }
    }
    return(mygene_tests)
  }
  print(table(dd_labelled$keep))
  
  # Apply BH FDR on interaction tests
  dd_tests <- dd_labelled[which(dd_labelled$keep != "Remove"),]
  FDRdd <- dd_tests[,c("gene","Condition","snps","pvalue","LDr2_othertop","genotype_beta","Interaction_term_beta","Interaction_term_pval","lmer_genotype_beta","lmer_Interaction_term_beta","lmer_Interaction_term_dpval","keep")]
  
  # (1) Linear model
  # Interation term
  FDRdd$Intterm_fdr <- ifelse(p.adjust(dd_tests$Interaction_term_pval, "BH") <= 0.05, "Significant", "Not")
  plot_pval_histogram(pval = dd_tests$Interaction_term_pval, filename = paste0("Plot_sanitycheck/",celltype,"_Histgram_pval_InteractionTerm.pdf"))
  plot_volcano(pval = dd_tests$Interaction_term_pval, est = dd_tests$Interaction_term_beta, significant = FDRdd$Intterm_fdr, filename = paste0("Plot_sanitycheck/",celltype,"_Volcano_InteractionTerm.pdf"))
  # ANOVA
  FDRdd$ANOVA_fdr <- ifelse(p.adjust(dd_tests$ANOVA_pval, "BH") <= 0.05, "Significant", "Not")
  plot_pval_histogram(pval = dd_tests$ANOVA_pval, filename = paste0("Plot_sanitycheck/",celltype,"_Histgram_pval_ANOVA.pdf"))
  plot_volcano(pval = dd_tests$ANOVA_pval, est = dd_tests$Interaction_term_beta, significant = FDRdd$ANOVA_fdr, filename = paste0("Plot_sanitycheck/",celltype,"_Volcano_ANOVA.pdf"))
  # Score
  FDRdd$score_fdr <- ifelse(p.adjust(dd_tests$score_pval, "BH") <= 0.05, "Significant", "Not")
  plot_pval_histogram(pval = dd_tests$score_pval, filename = paste0("Plot_sanitycheck/",celltype,"_Histgram_pval_score.pdf"))
  plot_volcano(pval = dd_tests$score_pval, est = dd_tests$Interaction_term_beta, significant = FDRdd$score_fdr, filename = paste0("Plot_sanitycheck/",celltype,"_Volcano_score.pdf"))
  # RobWald
  FDRdd$robWald_fdr <- ifelse(p.adjust(dd_tests$robWald_intterm_pval, "BH") <= 0.05, "Significant", "Not")
  plot_pval_histogram(pval = dd_tests$robWald_intterm_pval, filename = paste0("Plot_sanitycheck/",celltype,"_Histgram_pval_robWald.pdf"))
  plot_volcano(pval = dd_tests$robWald_intterm_pval, est = dd_tests$Interaction_term_beta, significant = FDRdd$robWald_fdr, filename = paste0("Plot_sanitycheck/",celltype,"_Volcano_robWald.pdf"))
  # LR
  FDRdd$LR_fdr <- ifelse(p.adjust(dd_tests$LR_pval, "BH") <= 0.05, "Significant", "Not")
  plot_pval_histogram(pval = dd_tests$LR_pval, filename = paste0("Plot_sanitycheck/",celltype,"_Histgram_pval_LR.pdf"))
  plot_volcano(pval = dd_tests$LR_pval, est = dd_tests$Interaction_term_beta, significant = FDRdd$LR_fdr, filename = paste0("Plot_sanitycheck/",celltype,"_Volcano_LR.pdf"))
  # Permutation
  FDRdd$lm_perm_fdr <- ifelse(p.adjust(dd_tests$Interaction_term_dpval, "BH") <= 0.05, "Significant", "Not")
  plot_pval_histogram(pval = dd_tests$Interaction_term_dpval, filename = paste0("Plot_sanitycheck/",celltype,"_Histgram_pval_lm_permutation.pdf"))
  plot_volcano(pval = dd_tests$Interaction_term_dpval, est = dd_tests$Interaction_term_beta, significant = FDRdd$lm_perm_fdr, filename = paste0("Plot_sanitycheck/",celltype,"_Volcano_lm_permutation.pdf"))
  
  # (2) lmer
  # Permutation
  FDRdd$lmer_perm_fdr <- ifelse(p.adjust(dd_tests$lmer_Interaction_term_dpval, "BH") <= 0.05, "Significant", "Not")
  plot_pval_histogram(pval = dd_tests$lmer_Interaction_term_dpval, filename = paste0("Plot_sanitycheck/",celltype,"_Histgram_pval_lmer_permutation.pdf"))
  plot_volcano(pval = dd_tests$lmer_Interaction_term_dpval, est = dd_tests$Interaction_term_beta, significant = FDRdd$lmer_perm_fdr, filename = paste0("Plot_sanitycheck/",celltype,"_Volcano_lmer_permutation.pdf"))
  
  # (3) gee
  # Permutation
  FDRdd$gee_perm_fdr <- ifelse(p.adjust(dd_tests$gee_Interaction_term_dpval, "BH") <= 0.05, "Significant", "Not")
  plot_pval_histogram(pval = dd_tests$gee_Interaction_term_dpval, filename = paste0("Plot_sanitycheck/",celltype,"_Histgram_pval_gee_permutation.pdf"))
  plot_volcano(pval = dd_tests$gee_Interaction_term_dpval, est = dd_tests$Interaction_term_beta, significant = FDRdd$gee_perm_fdr, filename = paste0("Plot_sanitycheck/",celltype,"_Volcano_gee_permutation.pdf"))
  
  write.table(FDRdd, paste0("FDR_significant_reQTL_",celltype,".txt"), quote = F, sep = "\t", row.names = F)
  return(NULL)
}




