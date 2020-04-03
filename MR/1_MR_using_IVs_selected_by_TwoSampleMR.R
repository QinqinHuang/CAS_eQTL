#----------------------------------------------
# 2019-05-06 last modified
# Aim: to perform Mendelian Randomisation 
# analysis and make effect size plots using
# Scott's codes.
#
# The genetic IVs were harmanised and selected
# using the TwoSampleMR package:
#   "MR_IVs.txt"
# Focus on genes with at least three IVs.
#----------------------------------------------
library(MendelianRandomization)
library(data.table)
library(ggplot2)
library(foreach)
options(stringsAsFactors = F)

# Functions defined by Scott
source("~/CAS/Analysis/Scripts/MR/MR_functions_ScottRitchie_2019_04_17.R")

# Working directory
setwd("~/CAS/Analysis/Mendelian_Randomisation_ciseGene_diseases/MendelianRandomization_package/")


#----- Load the genetic IVs -----
# Load the results calcualted by TwoSampleMR
MRresults <- read.table("../All_MR_results.txt", header = T, sep = "\t")
# Most of them have only one IV
table(MRresults$nsnp)

# Load the instrumental variables harmanised and seleted by TwoSampleMR
allMRIVs <- read.table("../MR_IVs.txt", header = T)

# The list of external GWAS
gwas_phenotype <- read.table("~/CAS/Analysis/Colocalisation/GWAS_phenotypes_files.txt", header = T, sep = "\t", quote = "")
gwas_phenotype$outcome <- gsub(pattern = ".txt", replacement = "", x = gwas_phenotype$file_name)


#----- Select exposure-disease pairs that ≥3 IVs are available -----
MR_cand <- MRresults[which(MRresults$nsnp >= 3),]
MR_cand <- unique(MR_cand[,c("outcome","exposure","condition")])

# For 19 GWAS we focused on,
# 52 genes had ≥3 IVs in at least one condition;
length(unique(MRresults$exposure[which(MRresults$nsnp >= 3 & MRresults$outcome %in% gwas_phenotype$outcome[which(gwas_phenotype$keep == 1)])]))
# 15 genes had ≥4 IVs in at least one condition;
length(unique(MRresults$exposure[which(MRresults$nsnp >= 4 & MRresults$outcome %in% gwas_phenotype$outcome[which(gwas_phenotype$keep == 1)])]))
# 9 genes had ≥5 IVs in at least one condition;
length(unique(MRresults$exposure[which(MRresults$nsnp >= 5 & MRresults$outcome %in% gwas_phenotype$outcome[which(gwas_phenotype$keep == 1)])]))
# 3 genes had ≥6 IVs: "BTN3A2" "HLA-C"  "MICB"
unique(MRresults$exposure[which(MRresults$nsnp >= 6 & MRresults$outcome %in% gwas_phenotype$outcome[which(gwas_phenotype$keep == 1)])])
# 1 gene had ≥7 IVs: BTN3A2
unique(MRresults$exposure[which(MRresults$nsnp >= 7 & MRresults$outcome %in% gwas_phenotype$outcome[which(gwas_phenotype$keep == 1)])])

# 963 pairs across four conditions
dim(MR_cand)
MR_cand_IVs <- merge(allMRIVs, MR_cand, by = c("outcome","exposure","condition"))
dim(MR_cand_IVs); dim(unique(MR_cand_IVs[,c("outcome","exposure","condition")]))


#----- Perform MR tests using Scott's codes -----
# Test all genes with multiple IVs
mr_results_all <- foreach(ii = 1:nrow(MR_cand), .combine = rbind) %do% {
  # Parameters
  mygene <- MR_cand$exposure[ii]
  mycondition <- MR_cand$condition[ii]
  myoutcome <- MR_cand$outcome[ii]
  mydiseasename <- gwas_phenotype[which(gwas_phenotype$outcome == myoutcome), "Disease_name"]
  
  testdd <- MR_cand_IVs[which(MR_cand_IVs$exposure == mygene & MR_cand_IVs$outcome == myoutcome & MR_cand_IVs$condition == mycondition),]
  nIVs <- nrow(testdd)
  if(nIVs <3) {cat(ii,"\n")}
  
  # Construct a data.table containing information for each instrument
  iv_table <- testdd[,c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome", "SNP", 
                        "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")]
  iv_table$type <- "cis"
  iv_table$colocalises <- F
  colnames(iv_table)[-9:-10] <- c("effect.pqtl", "se.pqtl", "effect.gwas", "se.gwas", "ALT_ID",
                                  "EA", "OA", "EAF.pqtl")
  iv_table <- as.data.table(iv_table)
  
  # Use the minor allele
  iv_table[EAF.pqtl > 0.5, c("EA", "OA", "EAF.pqtl", "effect.pqtl", "effect.gwas") :=
             .(OA, EA, 1 - EAF.pqtl, -effect.pqtl, -effect.gwas)]
  
  # Construct MR input object:
  mri <- mr_input(
    bx = iv_table$effect.pqtl,
    bxse = iv_table$se.pqtl,
    by = iv_table$effect.gwas,
    byse = iv_table$se.gwas)
  
  # Run all available MR methods constructing a data.table of results
  mr_results <- mr_tryall(mri)
  
  # Generate plot of dose response curve
  if(mycondition == "m_resting") {
     mycondition_label <- "Resting M"
   } else if(mycondition == "m_LPS") {
     mycondition_label <- "LPS M"
   } else if(mycondition == "t_resting") {
     mycondition_label <- "Resting T"
   } else if(mycondition == "t_PHA") {
     mycondition_label <- "PHA T"
   }
  g <- gg_dose_response(iv_table, mr_results, label=FALSE) + # label=TRUE needs fixing to work on new iv_table
    ylab(paste0("Log odds of ", mydiseasename)) +
    xlab(paste0("SD effect on ", mygene, " (", mycondition_label, ")"))
  ggsave(g, file = paste0("./Dose_response_curve/",myoutcome,"_",mygene,"_",mycondition,".pdf"), width = 6.8, height = 4.8)
  
  mr_results$gene <- mygene
  mr_results$condition <- mycondition_label
  mr_results$outcome <- myoutcome
  mr_results$nIVs <- nIVs
  
  # Write the results for each loop
  write.table(mr_results, "mr_results_append.txt", quote = F, sep = "\t", row.names = F, col.names = ifelse(ii == 1, T, F), append = T)
  cat("  ",ii,"finished! \n")
  
  return(as.data.frame(mr_results))
}

mr_results <- read.table("mr_results_append.txt", header = T, sep = "\t")

mr_results <- mr_results[,c("outcome","gene","condition","nIVs","method","mr_estimate","mr_se","mr_L95","mr_U95","mr_pval")]
write.table(mr_results, "MR_16methods_statistics.txt", quote = F, sep = "\t", row.names = F)



