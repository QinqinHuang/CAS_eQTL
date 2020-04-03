#----------------------------------------------
# 2018-08-24
# Aim of this analysis is to identify reQTLs 
# using interaction test combining data from 
# resting and stimulated samples.
#   - Run interaction tests on top eSNPs of 
#   eGenes that are significant (eigenMT-BH) in
#   at least one condition. 
#   - Run permutations to get empirical pvalues.
# 
# Interaction test:
# - Full model, including interaction terms 
# between covariates and condition.
# - Given that the majority of the individuals
# have both resting and stimualted samples, we
# also adjust for that.
#   (1) Allow individual specific random effect.
#   (2) Or use GEE model, which was designed for 
# longitudinal analysis and it accounts for 
# repeated measures for the same individual.
#
# Use the following methods to compare two models
# when possible.
#   ANOVA, score test, robWald test, LR test
#
# Also use permutations to estimate empirical
# pvalues. 
#   Permutate condition within each individual.
#
# We will restrict our tests to independent top
# eSNPs:
# - If an eGene has different top eSNPs in two
# conditions, keep both tests if in high LD
# (r2 ≥0.8).
# - If it has the same top eSNP or top eSNPs 
# are in high LD, keep one test.
#----------------------------------------------
library(foreach)
library(doMC)#; registerDoMC(cores = 20)
library(lme4)
library(gee)
source("/projects/qinqinhuang/CAS/Analysis/Scripts/score_rwald_LR_test_Agus.R")
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

# Number of permutations
npermutations <- 1000

# Working directory
setwd("/projects/qinqinhuang/CAS/Analysis/Response_eQTLs_ind_specific_randomeffect/")


#----- Functions -----
# 1. Get top eSNP and SNPs in perfect LD
get_topeSNP <- function(eAsso) {
  eAsso <- eAsso[order(eAsso$pvalue),]
  topSNP <- eAsso[!duplicated(eAsso$gene),]
  # All SNPs in perfect LD
  topSNP$snps_perfLD <- sapply(1:nrow(topSNP), function(ii) {
    paste(sort(merge(topSNP[ii,c("gene","statistic","pvalue","beta")], eAsso, by = c("gene","statistic","pvalue","beta"))$snps), collapse = ";")
  })
  topSNP <- topSNP[order(topSNP$gene),]
  rownames(topSNP) <- 1:nrow(topSNP)
  return(topSNP)
}

# 2. Permute condition within individual
perm_withinsubj <- function(dd, cond_coln = "condition", ind_coln = "subjid") {
  permdd <- foreach(subj = unique(dd[,ind_coln]), .combine = rbind) %do% {
    dd_subj <- dd[which(dd[,ind_coln] == subj),]
    if(nrow(dd_subj) > 1) {
      dd_subj[,cond_coln] <- sample(dd_subj[,cond_coln], 2, replace = F)
    }
    return(dd_subj)
  }
  return(permdd)
}


#----- Load data -----
# Load cis-eQTLs
# A list of cis eQTLs in m_resting, m_LPS, t_resting, t_PHA
ciseQTLs <- load_cis_eQTLs(type = "NoLD")


#----- Interaction test within cell type -----
nothing <- foreach(celltype = c("m","t")) %do% {
  # Cell type and treatment
  cat("  * Working on", celltype, "cells...\n")
  if(celltype == "m") {treatment <- "LPS"} else {treatment <- "PHA"}
  
  # eQTLs in two conditions
  eAsso_resting <- ciseQTLs[[1+ifelse(celltype == "m",0,2)]]
  eAsso_sti <- ciseQTLs[[2+ifelse(celltype == "m",0,2)]]

  #-- Top eSNP --
  topSNP0 <- get_topeSNP(eAsso_resting)
  topSNP0$Condition <- "Resting"
  topSNP1 <- get_topeSNP(eAsso_sti)
  topSNP1$Condition <- "Stimulated"
  #cat("  * Number of eGenes:", nrow(topSNP0), "in resting cells, and", nrow(topSNP1), "after stimulation. In total", length(unique(c(topSNP0$gene, topSNP1$gene))), "eGenes. \n")
  
  # All eGene-topSNP pairs we will test:
  top_pairs <- rbind(topSNP0, topSNP1)
  top_pairs <- top_pairs[order(top_pairs$gene, top_pairs$Condition, top_pairs$snps),]
  rownames(top_pairs) <- 1:nrow(top_pairs)
  top_pairs <- top_pairs[,c("gene","Condition",setdiff(colnames(top_pairs), c("gene","Condition")))]
  
  #-- LD r2 between top eSNPs from two conditions --
  # Read LD r2 between eGenes' top eSNPs in two conditions
  # NA if the eGene was significant in only one condition
  topeSNP_LD <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/LD_top_SNP_conditions/Top_eSNP_eGene_LD_r2_",celltype,".txt"), header = T)
  top_pairs$LDr2_othertop <- topeSNP_LD$LDr2[match(top_pairs$gene, topeSNP_LD$gene)]
  
  
  #-- Load gene expression and covariate data --
  # gene expression
  expr_resting <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_Ctrl_gene/gene_expression_",celltype,"_Ctrl.txt"), header = T)
  expr_sti <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_gene/gene_expression_",celltype,"_",treatment,".txt"), header = T)
  
  # covariates
  cov_resting <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_Ctrl_gene/Covariates_",celltype,"_Ctrl.txt"), header = T)
  cov_sti <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_gene/Covariates_",celltype,"_",treatment,".txt"), header = T)
  
  # Make sure samples are in the same order
  if(!(identical(colnames(expr_resting)[-1], colnames(cov_resting)[-1]) & identical(colnames(expr_sti)[-1], colnames(cov_sti)[-1]))) {
    cat("  * Samples are not in the same order in express matrix and in cov matrix in",celltype,"cells.\n")
    return(NULL)
  }
  
  
  #-- Run interaction test --
  # Covariate name
  cov_name <- cov_resting$covariates
  
  # (1) No individual info
  # H0: No interaction term for genotype
  formula_qtl <- as.formula(paste("expression ~ genotype + condition ", 
                                  paste(cov_name, collapse = " + "), 
                                  paste(paste0(cov_name, ":condition"), collapse = " + "), 
                                  sep = "+ "))
  # H1: add interaction term to model the interaction between genotype and condition
  formula_interaction <- as.formula(paste("expression ~ genotype + condition:genotype + condition ", 
                                          paste(cov_name, collapse = " + "), 
                                          paste(paste0(cov_name, ":condition"), collapse = " + "), 
                                          sep = "+ "))
  
  # (2) Use lem4 package to model individual specific random effect
  formula_qtl_ind <- as.formula(paste("expression ~ genotype + condition + (1|subjid) ", 
                                      paste(cov_name, collapse = " + "), 
                                      paste(paste0(cov_name, ":condition"), collapse = " + "), 
                                      sep = "+ "))
  formula_interaction_ind <- as.formula(paste("expression ~ genotype + condition + (1|subjid) + condition:genotype ", 
                                              paste(cov_name, collapse = " + "), 
                                              paste(paste0(cov_name, ":condition"), collapse = " + "), 
                                              sep = "+ "))

  
  # Go through each eGene-top SNP pair
  interaction_test_results <- foreach(ii = 1:nrow(top_pairs), .combine = rbind) %dopar% {
    
    # The eGene-topeSNP pair we focus on:
    mySNP <- top_pairs$snps[ii]
    mygene <- top_pairs$gene[ii]
    
    # Prepare data - rows are samples and columns are variables
    # Genotype data
    genodata <- extract_geno(SNPID = mySNP)
    tgeno <- as.data.frame(t(genodata[,-1]))
    colnames(tgeno) <- "genotype"
    tgeno$subjid <- sapply(rownames(tgeno), function(x) {unlist(strsplit(x, split = "_"))[2]})
    
    # Gene expression, covariate, and genotype data in resting cells
    texpr_resting <- as.data.frame(t(expr_resting[which(expr_resting$geneid == mygene),-1]))
    colnames(texpr_resting) <- "expression"
    texpr_resting$subjid <- sapply(rownames(texpr_resting), function(x) {unlist(strsplit(x, split = "_"))[2]})
    # Covariate resting
    tcov_resting <- as.data.frame(t(cov_resting[,-1]))
    colnames(tcov_resting) <- cov_resting$covariates
    tcov_resting$subjid <- sapply(rownames(tcov_resting), function(x) {unlist(strsplit(x, split = "_"))[2]})
    # Merge
    t_resting <- merge(texpr_resting, tcov_resting, by = "subjid")
    t_resting <- merge(t_resting, tgeno, by = "subjid")
    # Condition
    t_resting$condition <- 0
    
    # Gene expression, covariate, and genotype data in stimulated cells
    texpr_sti <- as.data.frame(t(expr_sti[which(expr_sti$geneid == mygene),-1]))
    colnames(texpr_sti) <- "expression"
    texpr_sti$subjid <- sapply(rownames(texpr_sti), function(x) {unlist(strsplit(x, split = "_"))[2]})
    # Covariate stimulated
    tcov_sti <- as.data.frame(t(cov_sti[,-1]))
    colnames(tcov_sti) <- cov_sti$covariates
    tcov_sti$subjid <- sapply(rownames(tcov_sti), function(x) {unlist(strsplit(x, split = "_"))[2]})
    # Merge
    t_sti <- merge(texpr_sti, tcov_sti, by = "subjid")
    t_sti <- merge(t_sti, tgeno, by = "subjid")
    # Condition
    t_sti$condition <- 1
    
    # Combined two datasets
    ddall <- rbind(t_resting, t_sti)
    
    
    #--------------------------------------------------------------------------------
    # (1) Test two linear regression models without considering individual info
    lmH0 <- lm(formula_qtl, ddall)
    lmH1 <- lm(formula_interaction, ddall)
    
    # Compare two models
    # ANOVA
    anovatest <- anova(lmH0, lmH1)
    # Score test
    scoretest <- score.test(lmH0, lmH1)
    # robWald test
    robWaldtest <- robWald.test(lmH1)
    # LR test
    LRtest <- LR.test(lmH0, lmH1)
    
    # Return results from all tests
    returndd <- top_pairs[ii,]
    
    # Statistics from interaction test
    # Genotype term
    returndd$genotype_beta <- summary(lmH1)$coefficients["genotype",1]
    returndd$genotype_se <- summary(lmH1)$coefficients["genotype",2]
    returndd$genotype_tvalue <- summary(lmH1)$coefficients["genotype",3]
    returndd$genotype_pval <- summary(lmH1)$coefficients["genotype",4]
    # Interaction term
    returndd$Interaction_term_beta <- summary(lmH1)$coefficients["genotype:condition",1]
    returndd$Interaction_term_se <- summary(lmH1)$coefficients["genotype:condition",2]
    returndd$Interaction_term_tvalue <- summary(lmH1)$coefficients["genotype:condition",3]
    returndd$Interaction_term_pval <- summary(lmH1)$coefficients["genotype:condition",4]
    
    # ANOVA
    returndd$ANOVA_F <- anovatest[[5]][2]
    returndd$ANOVA_pval <- anovatest[[6]][2]
    
    # Score test
    returndd$score_stat <- scoretest$stat[1,1]
    returndd$score_pval <- scoretest$pval[1,1]
    
    # robWald test
    returndd$robWald_genotype_stat <- robWaldtest$stat["genotype"]
    returndd$robWald_genotype_pval <- robWaldtest$pval["genotype"]
    returndd$robWald_intterm_stat <- robWaldtest$stat["genotype:condition"]
    returndd$robWald_intterm_pval <- robWaldtest$pval["genotype:condition"]
    
    # LR test
    returndd$LR_stat <- LRtest$stat
    returndd$LR_pval <- LRtest$pval
    
    
    #--------------------------------------------------------------------------------
    # (2) Test two linear mixed effect models adjusting for repeated individual data (lme4)
    lmerH0 <- lmer(formula_qtl_ind, ddall)
    lmerH1 <- lmer(formula_interaction_ind, ddall)
    
    # Compare two models
    # ANOVA
    anovatest_lmer <- anova(lmerH0, lmerH1)
    # Score test
    #scoretest_lmer <- score.test(lmerH0, lmerH1)
    # robWald test
    #robWaldtest <- robWald.test(lmerH1)
    # LR test
    LRtest_lmer <- LR.test(lmerH0, lmerH1)
    
    # Statistics from interaction test adjusting for repeated measures for the same individual
    # Genotype term
    returndd$lmer_genotype_beta <- summary(lmerH1)$coefficients["genotype",1]
    returndd$lmer_genotype_se <- summary(lmerH1)$coefficients["genotype",2]
    returndd$lmer_genotype_tvalue <- summary(lmerH1)$coefficients["genotype",3]
    # Interaction term
    returndd$lmer_Interaction_term_beta <- summary(lmerH1)$coefficients["genotype:condition",1]
    returndd$lmer_Interaction_term_se <- summary(lmerH1)$coefficients["genotype:condition",2]
    returndd$lmer_Interaction_term_tvalue <- summary(lmerH1)$coefficients["genotype:condition",3]
    
    # ANOVA
    returndd$lmer_ANOVA_Chisq <- anovatest_lmer$Chisq[2]
    returndd$lmer_ANOVA_pval <- anovatest_lmer$Pr[2]
    
    # LR test
    returndd$lmer_LR_stat <- LRtest_lmer$stat
    returndd$lmer_LR_pval <- LRtest_lmer$pval
    
    
    #--------------------------------------------------------------------------------
    # (3) Test two models adjusting for repeated individual data (gee)
    # Sort the data by subject
    ddgee <- ddall[order(ddall$subjid),]
    ddgee$subjid <- factor(ddgee$subjid)
    #geeH0 <- gee(formula_qtl, id = subjid, data = ddgee)
    geeH1 <- gee(formula_interaction, id = subjid, data = ddgee)
    g <- summary(geeH1)
    
    # Statistics from interaction test adjusting for repeated measures for the same individual using GEE model
    # Genotype term
    returndd$gee_genotype_beta <- g$coefficients["genotype",1]
    returndd$gee_genotype_se <- g$coefficients["genotype",4]
    returndd$gee_genotype_zvalue <- g$coefficients["genotype",5]
    # Interaction term
    returndd$gee_Interaction_term_beta <- g$coefficients["genotype:condition",1]
    returndd$gee_Interaction_term_se <- g$coefficients["genotype:condition",4]
    returndd$gee_Interaction_term_zvalue <- g$coefficients["genotype:condition",5]
    
    
    #--------------------------------------------------------------------------------
    # Permutation
    cat("  *Permutation analysis...")
    permuted_pval <- foreach(pp = 1:npermutations, .combine = rbind) %dopar% {
      cat(pp,"")
      
      # Permute condition data within each subject
      perm_dd <- ddall[order(ddall$subjid),]
      perm_dd <- perm_withinsubj(perm_dd)
      #perm_dd$condition <- sample(perm_dd$condition, nrow(perm_dd), replace = F)
      
      # Interaction test
      # Linear model
      permlmH1 = lm(formula_interaction, perm_dd)
      # Linear mixed effect model
      permlmerH1 <- lmer(formula_interaction_ind, perm_dd)
      # GEE model
      perm_dd <- perm_dd[order(perm_dd$subjid),]; perm_dd$subjid <- factor(perm_dd$subjid)
      permgeeH1 <- gee(formula_interaction, id = subjid, data = perm_dd)
      
      # Permuted pvalues
      curr_pval <- data.frame(lmH1_pval = summary(permlmH1)$coefficients["genotype:condition",4],
                              lmerH1_tvalue = summary(permlmerH1)$coefficients["genotype:condition",3],
                              geeH1_zvalue = summary(permgeeH1)$coefficients["genotype:condition",5] )
      
      # Remove some R objects
      rm(perm_dd); rm(permlmH1); rm(permlmerH1); rm(permgeeH1)
      
      return(curr_pval)
    } # End of permutations
    
    # Calculate empirical pvalues from permutations
    returndd$Interaction_term_dpval <- (sum(permuted_pval$lmH1_pval < returndd$Interaction_term_pval) + 1) / (nrow(permuted_pval) + 1)
    returndd$lmer_Interaction_term_dpval <- (sum(abs(permuted_pval$lmerH1_tvalue) > abs(returndd$lmer_Interaction_term_tvalue)) + 1) / (nrow(permuted_pval) + 1)
    returndd$gee_Interaction_term_dpval <- (sum(abs(permuted_pval$geeH1_zvalue) > abs(returndd$gee_Interaction_term_zvalue)) + 1) / (nrow(permuted_pval) + 1)
    
    rm(permuted_pval)  # clean up
    
    
    # Save the result for the current eGene-top SNP pair
    savecoln <- ifelse(ii == 1, T, F)
    write.table(returndd, paste0("Temp_Interaction_test_topSNP_perm_pval_",celltype,".txt"), quote = F, sep = "\t", row.names = F, col.names = savecoln, append = T)
    
    return(NULL)
  } # End of the loop for all tests
  
  
  # Run in parallel so column name may not be in the first line
  dd <- read.table(paste0("Temp_Interaction_test_topSNP_perm_pval_",celltype,".txt"), header = F)
  colnames(dd) <- dd[which(dd$V1 == "gene"),]
  dd <- dd[-which(dd$gene == "gene"),]
  dd <- dd[order(dd$gene, dd$Condition),]
  write.table(dd, paste0("Interaction_test_topSNP_perm_pval_",celltype,".txt"), row.names = F, sep = "\t", quote = F)
  
  return(NULL) 
}


# 20 parallel, 9.4 hr


# Monocytes
# * Number of eGenes: 971 in resting cells, and 1347 after stimulation. In total 1749 eGenes.

# T cells
# * Number of eGenes: 136 in resting cells, and 376 after stimulation. In total 398 eGenes.


