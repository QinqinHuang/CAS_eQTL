#----------------------------------------------
# 2019-05-08 Last modified
# Focus on these four methods: 
#  "IVW","Weighted median",
#  "Weighted mode (simple SE)", and "MR-Egger"
#
# Plot dose response curve showing these 4
# methods for significant results.
#   (1) Effect sizes of the minor allele
#whichallele <- "MinorAllele"
#   (2) IV-exposure associations are positive,
#       so that the MR-Egger is easier to show.
whichallele <- "IV_expo_Pos"
#
# Get the supplemental table for all significant
# MR results.
#----------------------------------------------
library(foreach)
library(ggplot2)
options(stringsAsFactors = F)
source("~/CAS/Analysis/Scripts/MR/MR_functions_ScottRitchie_2019_04_17.R")

setwd("~/CAS/Analysis/Mendelian_Randomisation_ciseGene_diseases/MendelianRandomization_package/")


#----- Load data -----
# The list of external GWAS
gwas_phenotype <- read.table("~/CAS/Analysis/Colocalisation/GWAS_phenotypes_files.txt", header = T, sep = "\t", quote = "")
gwas_phenotype$outcome <- gsub(pattern = ".txt", replacement = "", x = gwas_phenotype$file_name)

# Load the IVs harmanised and seleted by TwoSampleMR
allMRIVs <- read.table("../MR_IVs.txt", header = T)
allMRIVs$condition[which(allMRIVs$condition == "m_resting")] <- "Resting M"
allMRIVs$condition[which(allMRIVs$condition == "m_LPS")] <- "LPS M"
allMRIVs$condition[which(allMRIVs$condition == "t_resting")] <- "Resting T"
allMRIVs$condition[which(allMRIVs$condition == "t_PHA")] <- "PHA T"

# Read MR results
mr_results <- read.table("MR_16methods_statistics.txt", header = T, sep = "\t")
# Disease name
mr_results$Disease_label <- gwas_phenotype$Disease_label[match(mr_results$outcome, gwas_phenotype$outcome)]
# Label significant methods
mr_results$Significant05 <- ifelse(mr_results$mr_pval <= 0.05, "Significant", "Not")

# Focus on 19 GWAS datasets
mr_results_sub <- mr_results[which(mr_results$outcome %in% gwas_phenotype$outcome[which(gwas_phenotype$keep == 1)]),]

# Tested pairs
mr_tested <- unique(mr_results_sub[,c("outcome","gene","condition")])
rownames(mr_tested) <- 1:nrow(mr_tested)
mr_tested$n_met <- NA
mr_tested$n_sig <- NA

# Causal associations with ≥3 significant out of 4 methods
mr_sig <- foreach(ii = 1:nrow(mr_tested), .combine = rbind) %do% {
  dd <- mr_results_sub[which(mr_results_sub$gene == mr_tested$gene[ii] 
                             & mr_results_sub$outcome == mr_tested$outcome[ii] 
                             & mr_results_sub$condition == mr_tested$condition[ii]),]
  # Focus on 4 method
  dd_sub <- dd[c(which(dd$method %in% c("IVW","Weighted median","Weighted mode (simple SE)")),
                 which(dd$method == "MR-Egger"), which(dd$method == "MR-Egger")+1), ]
  # Number of methods succeed
  mr_tested$n_met[ii] <- length(which(!is.na(dd_sub$mr_estimate)))
  # Number of significant methods out of all 4
  mr_tested$n_sig[ii] <- length(which(dd_sub$Significant05[which(dd_sub$method != "(intercept)")] == "Significant"))
  
  # Return significant associations if 3 out of 4 significiant
  if(mr_tested$n_sig[ii] >= 3) {
    # IVs
    testdd <- allMRIVs[which(allMRIVs$exposure == mr_tested$gene[ii]
                             & allMRIVs$outcome == mr_tested$outcome[ii]
                             & allMRIVs$condition == mr_tested$condition[ii]),]
    
    # Construct a data.table containing information for each instrument
    iv_table <- testdd[,c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome", "SNP", 
                          "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")]
    colnames(iv_table) <- c("effect.pqtl", "se.pqtl", "effect.gwas", "se.gwas", "ALT_ID",
                            "EA", "OA", "EAF.pqtl")
    iv_table$type <- "cis"
    iv_table$colocalises <- F
    iv_table <- as.data.table(iv_table)
    
    # Use the minor allele
    if(whichallele == "MinorAllele") {
      iv_table[EAF.pqtl > 0.5, c("EA", "OA", "EAF.pqtl", "effect.pqtl", "effect.gwas") :=
                 .(OA, EA, 1 - EAF.pqtl, -effect.pqtl, -effect.gwas)]
    } else if(whichallele == "IV_expo_Pos") {       # Make the IV-exposure associations positive
      iv_table[effect.pqtl < 0, c("EA", "OA", "EAF.pqtl", "effect.pqtl", "effect.gwas") :=
                 .(OA, EA, 1 - EAF.pqtl, -effect.pqtl, -effect.gwas)]
    }
    
    # MR results
    mr_results <- as.data.table(dd_sub)
    
    # Disease name (on y lab)
    mydiseasename <- unique(dd_sub$Disease_label)
    
    # Plotting effect size (minor allele)
    g <- gg_dose_response(iv_table, mr_results, label=FALSE) + 
      ylab(paste0("Log odds of ", mydiseasename)) +
      xlab(paste0("SD effect on ", mr_tested$gene[ii], " (", mr_tested$condition[ii], ")"))
    # Effect sizes of the minor allele
    if(whichallele == "MinorAllele") {
      ggsave(g, file = paste0("./Sig_MR_3_out_of_4/Dose_response_curve_Sig_MR_MinorAllele/",mr_tested$outcome[ii],"_",mr_tested$gene[ii],"_",mr_tested$condition[ii],".pdf"), width = 4.5, height = 2.7)
    } else if(whichallele == "IV_expo_Pos") {  # make IV-exposure associations positive
      ggsave(g, file = paste0("./Sig_MR_3_out_of_4/Dose_response_curve_Sig_MR_IVexpo_positive/",mr_tested$outcome[ii],"_",mr_tested$gene[ii],"_",mr_tested$condition[ii],".pdf"), width = 4.5, height = 2.7)
    }
    
    # Return significant results
    return(dd_sub)
  } else {
    return(NULL)
  }
  
}

# Save the gene-disease pair that were tested in MR and the number of significant methods
write.table(mr_tested, "./Sig_MR_3_out_of_4/Number_significant_MR_methods_out_of_4.txt", quote = F, sep = "\t", row.names = F)

# Save the significant MR tests (≥3 out of 4)
mr_sig$Disease_label <- NULL
write.table(mr_sig, "./Sig_MR_3_out_of_4/MR_significant_3_out_of_4_results.txt", quote = F, sep = "\t", row.names = F)



#----- Supplemental table for all significant results -----
# The list of external GWAS
gwas_phenotype <- read.table("~/CAS/Analysis/Colocalisation/GWAS_phenotypes_files.txt", header = T, sep = "\t", quote = "")
gwas_phenotype$outcome <- gsub(pattern = ".txt", replacement = "", x = gwas_phenotype$file_name)

# Read all significant results
mr_sig <- read.table("./Sig_MR_3_out_of_4/MR_significant_3_out_of_4_results.txt", header = T, sep = "\t", quote = "")

# 28 eGenes had significant causal associations
length(unique(mr_sig$gene))

# Disease name
mr_sig$Disease <- gwas_phenotype$Disease_name[match(mr_sig$outcome, gwas_phenotype$outcome)]
# Sort by disease and gene names
mr_sig <- mr_sig[order(mr_sig$Disease, mr_sig$gene),]

mr_sig <- mr_sig[,c("Disease","gene","condition","nIVs","method","mr_estimate","mr_se","mr_L95","mr_U95","mr_pval","Significant05","outcome")]
colnames(mr_sig) <- c("Disease","Exposure","Condition","N_IVs","Method","Causal_estimate","SE","L95","U95","P-value","Significant","GWAS_stat")
write.table(mr_sig, "./Sig_MR_3_out_of_4/MR_stat_significant_3_out_of_4_results_suppletable.txt", quote = F, row.names = F, sep = "\t")




