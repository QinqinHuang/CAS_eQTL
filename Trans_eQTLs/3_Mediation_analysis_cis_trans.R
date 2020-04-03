#----------------------------------------------
# 2018-09-06
# Aim: 
# 1. to get trans-eQTLs that have cis-eGenes
# 2. to perform mediation analysis to identify
# cis-eGenes that act as cis-mediators of 
# trans-associations.
#
# Methods: "mediation" R package
# (Bootstrap method based on simulations)
#
# Which SNP is tested in mediation analysis?
#  - From all trans-eQTLs that have cis eGenes,
# the SNP with the minimum p-value.
#  
# Number of simulations for bootstrap and Monte 
# Carlo:
ndraws <- 10000
#----------------------------------------------
library(foreach)
library(mediation)
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

setwd("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/Mediation_analysis/")
if(!file.exists(paste0("Plots_mediation_n",ndraws))) {
  dir.create(paste0("Plots_mediation_n",ndraws))
}


#----------------------------------------------
# 1. List all eQTL-ciseGene-transeGene trios that will be tested for mediation 

# Load significant cis eQTLs
allciseQTLs <- load_cis_eQTLs()
names(allciseQTLs) <- c("m_Ctrl","m_LPS","t_Ctrl","t_PHA")

trios <- foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA"), .combine = rbind) %do% {
  cat("  * Working on", celltype, treatment, "...\n")

  # Read trans-eQTLs (all tests except for cis tests)
  transeQTLs <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/", celltype, "_", treatment, "_gene/Trans_allexceptcis_genomewideFDR_sig.txt"), header = T)
  transeQTLs$label_cis <- NULL; transeQTLs$FDR <- NULL
  cat("    Number of trans-eGenes:", length(unique(transeQTLs$gene)), "\n")
  
  # Cis eQTLs for this condition
  ciseQTLs <- allciseQTLs[[paste0(celltype, "_", treatment)]]
  colnames(ciseQTLs)[-2] <- paste0(colnames(ciseQTLs)[-2], "_cis")
  
  # Trans-eQTLs that also have cis-eGenes
  trans_cis_overlap <- merge(transeQTLs, ciseQTLs, by = "snps")
  cat("    Number of trans-eGenes whose trans-eQTLs also have cis-eGenes:", length(unique(trans_cis_overlap$gene)), "\n")
  
  # Canditate mediation trios
  candi_trios <- unique(trans_cis_overlap[,c("gene","gene_cis")])

  # We test the best SNP from SNPs that had significant cis eGenes
  # (Keep all SNPs in perfect LD, they may have different stats in cis mapping)
  candi_trios_topeSNP <- foreach(ii = 1:nrow(candi_trios), .combine = rbind) %do% {
    
    # All trans associations that have cis eGenes for this pair
    dd <- merge(candi_trios[ii,], trans_cis_overlap, by = c("gene","gene_cis"))
    # SNPs with the minimum pvalue
    ddtop <- dd[which(dd$pvalue == min(dd$pvalue)),]
    
    # Label whether the SNPs that will be tested in mediation analysis are top eSNPs
    # Trans associations
    curr_transeQTLs <- transeQTLs[which(transeQTLs$gene == candi_trios$gene[ii]),]
    ddtop$top_trans <- ifelse(ddtop$snps %in% curr_transeQTLs$snps[which(curr_transeQTLs$pvalue == min(curr_transeQTLs$pvalue))], TRUE, FALSE)
    # Cis associations
    curr_ciseQTLs <- ciseQTLs[which(ciseQTLs$gene == candi_trios$gene_cis[ii]),]
    ddtop$top_cis <- ifelse(ddtop$snps %in% curr_ciseQTLs$snps[which(curr_ciseQTLs$pvalue == min(curr_ciseQTLs$pvalue))], TRUE, FALSE)
    
    return(ddtop)
  }
  
  rownames(candi_trios_topeSNP) <- 1:nrow(candi_trios_topeSNP)
  candi_trios_topeSNP$celltype <- celltype
  candi_trios_topeSNP$treatment <- treatment
  candi_trios_topeSNP <- candi_trios_topeSNP[,c("celltype","treatment",setdiff(colnames(candi_trios_topeSNP), c("celltype","treatment")))]
  
  return(candi_trios_topeSNP) 
}

write.table(trios, "eQTL_ciseGene_transeGene_trios_mediation_analysis.txt", quote = F, sep = "\t", row.names = F)


#----------------------------------------------
# 2. Mediation analysis
#trios <- read.table("eQTL_ciseGene_transeGene_trios_mediation_analysis.txt", header = T)

# Four conditions
mediation_results <- foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA"), .combine = rbind) %do% {
  cat("  * Working on", celltype, treatment, "...\n")
  
  # Read covariates
  cov <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_gene/Covariates_",celltype,"_",treatment,".txt"), header = T)
  rownames(cov) <- cov$covariates
  cov <- as.data.frame(t(cov[,-1]))
  cov$subjid <- sapply(rownames(cov), function(x) {unlist(strsplit(x, split = "_"))[2]})
  covnames <- setdiff(colnames(cov), "subjid")
  
  # Go through each trio in this condition
  currtrios <- trios[which(trios$celltype == celltype & trios$treatment == treatment),]
  currresults <- foreach(ii = 1:nrow(currtrios), .combine = rbind) %do% {
    # Mediation trio
    mySNP <- currtrios$snps[ii]
    mycisgene <- currtrios$gene_cis[ii]
    mytransgene <- currtrios$gene[ii]
    cat("  Trans-eGene:", mytransgene, "cis-eGene:", mycisgene, "SNP:", mySNP, "\n")
    
    # Read genotype data
    geno_X <- extract_geno(mySNP)
    geno_X <- as.data.frame(t(geno_X[,-1]))
    colnames(geno_X) <- "SNP"
    geno_X$subjid <- sapply(rownames(geno_X), function(x) {unlist(strsplit(x, split = "_"))[2]})
    
    # Cis eGene expression
    expr_M <- extract_expression(mycisgene, celltype = celltype, treatment = treatment)
    expr_M <- as.data.frame(t(expr_M[,-1]))
    colnames(expr_M) <- "ciseGene"
    expr_M$subjid <- sapply(rownames(expr_M), function(x) {unlist(strsplit(x, split = "_"))[2]})
    
    # Trans eGene expression
    expr_Y <- extract_expression(mytransgene, celltype = celltype, treatment = treatment)
    expr_Y <- as.data.frame(t(expr_Y[,-1]))
    colnames(expr_Y) <- "transeGene"
    expr_Y$subjid <- sapply(rownames(expr_Y), function(x) {unlist(strsplit(x, split = "_"))[2]})
    
    # Merge data
    currdd <- merge(expr_Y, expr_M, by = "subjid")
    currdd <- merge(currdd, geno_X, by = "subjid")
    currdd <- merge(currdd, cov, by = "subjid")
    if(nrow(currdd) != nrow(cov)) {
      cat("  ** Error - number of subjects is not correct!\n"); return(NULL)
    }
    
    # Step 1: Trans-eGene ~ eQTL
    lm1 <- lm(formula(paste0("transeGene ~ SNP + ", paste(covnames, collapse = " + "))), currdd)
    
    # Step 2: Cis-eGene ~ eQTL
    lm2 <- lm(formula(paste0("ciseGene ~ SNP + ", paste(covnames, collapse = " + "))), currdd)
    
    # Step 3: Trans-eGene ~ eQTL + cis-eGene
    lm3 <- lm(formula(paste0("transeGene ~ SNP + ciseGene + ", paste(covnames, collapse = " + "))), currdd)
    
    # Method by default: the quasi-Bayesian Monte Carlo method based on normal approximation is used for CI.
    # sims: number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.
    # robustSE TRUE: heteroskedasticity-consistent standard errors will be used in quasi-Bayesian simulations.
    medout_quasiBayesian <- mediate(model.m = lm2, model.y = lm3, treat = "SNP", mediator = "ciseGene", robustSE = TRUE, sims = ndraws)
    summ_quasiBayesian <- summary(medout_quasiBayesian)
    pdf(paste0("Plots_mediation_n",ndraws,"/",celltype,"_",treatment,"_transeGene_",mytransgene,"_ciseGene_",mycisgene,"_SNP_",mySNP,"_quasiBaye.pdf"), 5, 5)
    plot(medout_quasiBayesian)
    dev.off()
    
    # Nonparametric bootstrap is used for CI; default "perc" basic percentile intervals
    medout_bootstrap <- mediate(model.m = lm2, model.y = lm3, boot = TRUE, treat = "SNP", mediator = "ciseGene", sims = ndraws)
    summ_bootstrap <- summary(medout_bootstrap)
    pdf(paste0("Plots_mediation_n",ndraws,"/",celltype,"_",treatment,"_transeGene_",mytransgene,"_ciseGene_",mycisgene,"_SNP_",mySNP,"_bootstrap.pdf"), 5, 5)
    plot(medout_bootstrap)
    dev.off()
    
    # Bootstrap method: "bca" bias-corrected and accelerated CIs
    medout_bootstrap_bca <- mediate(model.m = lm2, model.y = lm3, boot = TRUE, boot.ci.type = "bca", treat = "SNP", mediator = "ciseGene", sims = ndraws)
    summ_bootstrap_bca <- summary(medout_bootstrap_bca)
    pdf(paste0("Plots_mediation_n",ndraws,"/",celltype,"_",treatment,"_transeGene_",mytransgene,"_ciseGene_",mycisgene,"_SNP_",mySNP,"_bootstrap_bca.pdf"), 5, 5)
    plot(medout_bootstrap_bca)
    dev.off()
    
    # Summary statistics from the third linear regression
    returnstat <- currtrios[ii,]
    returnstat$SNP_beta <- summary(lm3)$coefficients["SNP",1]
    returnstat$SNP_se <- summary(lm3)$coefficients["SNP",2]
    returnstat$SNP_pvalue <- summary(lm3)$coefficients["SNP",4]
    returnstat$ciseGene_beta <- summary(lm3)$coefficients["ciseGene",1]
    returnstat$ciseGene_se <- summary(lm3)$coefficients["ciseGene",2]
    returnstat$ciseGene_pvalue <- summary(lm3)$coefficients["ciseGene",4]
    
    # Statistics from summary of the "mediation" output
    # ACME (the average causal mediation effects): $d.avg
    # ADE (the average direct effects): $z.avg
    # Total Effect: $tau
    # Prop. Mediated: $n.avg
    # Quasi-Bayesian Monte Carlo
    returnstat$ACME_quasiBayeMonteCarlo <- summ_quasiBayesian$d.avg
    returnstat$ACME_pval_quasiBayeMonteCarlo <- summ_quasiBayesian$d.avg.p
    returnstat$ADE_quasiBayeMonteCarlo <- summ_quasiBayesian$z.avg
    returnstat$ADE_pval_quasiBayeMonteCarlo <- summ_quasiBayesian$z.avg.p
    returnstat$TE_quasiBayeMonteCarlo <- summ_quasiBayesian$tau.coef
    returnstat$TE_pval_quasiBayeMonteCarlo <- summ_quasiBayesian$tau.p
    returnstat$PropMed_quasiBayeMonteCarlo <- summ_quasiBayesian$n.avg
    returnstat$PropMed_pval_quasiBayeMonteCarlo <- summ_quasiBayesian$n.avg.p
    # Bootstrap
    returnstat$ACME_boot <- summ_bootstrap$d.avg
    returnstat$ACME_pval_boot <- summ_bootstrap$d.avg.p
    returnstat$ADE_boot <- summ_bootstrap$z.avg
    returnstat$ADE_pval_boot <- summ_bootstrap$z.avg.p
    returnstat$TE_boot <- summ_bootstrap$tau.coef
    returnstat$TE_pval_boot <- summ_bootstrap$tau.p
    returnstat$PropMed_boot <- summ_bootstrap$n.avg
    returnstat$PropMed_pval_boot <- summ_bootstrap$n.avg.p
    # Bias corrected and accelerated bootstrap
    returnstat$ACME_bootbca <- summ_bootstrap_bca$d.avg
    returnstat$ACME_pval_bootbca <- summ_bootstrap_bca$d.avg.p
    returnstat$ADE_bootbca <- summ_bootstrap_bca$z.avg
    returnstat$ADE_pval_bootbca <- summ_bootstrap_bca$z.avg.p
    returnstat$TE_bootbca <- summ_bootstrap_bca$tau.coef
    returnstat$TE_pval_bootbca <- summ_bootstrap_bca$tau.p
    returnstat$PropMed_bootbca <- summ_bootstrap_bca$n.avg
    returnstat$PropMed_pval_bootbca <- summ_bootstrap_bca$n.avg.p
    
    return(returnstat)
  }
  
  return(currresults) 
}

write.table(mediation_results, paste0("Mediation_analysis_results_n",ndraws,".txt"), quote = F, sep = "\t", row.names = F)

  
 
 

