#----------------------------------------------
# 2018-10-02
# Aim: to generate supplemental tables for 
# reQTLs.
# Use SNP rsids.
#
# 2018-11-13
# For eGenes whose two top eSNPs were in high
# LD and only one was tested, condition change
# to "Stimulated; Resting".
#----------------------------------------------
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

# Working directory
setwd("/projects/qinqinhuang/CAS/Analysis/Response_eQTLs_ind_specific_randomeffect/")

# Load annotations for SNPs tested in CAS
SNPanno <- read.table("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/SNP_pos_rsID_counted_freq_freq1kG.txt", header = T)

# Load reQTLs
reQTLs <- load_reQTLs()

# Two cell types
for(celltype in c("m", "t")) {
  myreQTLs <- reQTLs[which(reQTLs$celltype == celltype),]
  
  # Add other info on rsid, effective allele, MAF, etc.
  myreQTLs <- merge(myreQTLs, SNPanno[,c("snps","rsid","chr","position","COUNTED","ALT","freq_counted","Genotyped")], by = "snps")
  
  # For eGene - "Keep"
  myc <- myreQTLs[which(myreQTLs$keep == "Keep"),]$Condition
  myreQTLs[which(myreQTLs$keep == "Keep"),]$Condition <- ifelse(myc == "Stimulated", yes = "Stimulated;Resting", no = "Resting;Stimulated")
  
  # Sort by significance level
  myreQTLs <- myreQTLs[order(myreQTLs$lmer_Interaction_term_dpval, -abs(myreQTLs$lmer_Interaction_term_beta)),]
  
  myreQTLs <- myreQTLs[,c("Condition","gene","snps","rsid","chr","position","LDr2_othertop","keep","genotype_beta","Interaction_term_beta","Interaction_term_pval","lmer_genotype_beta","lmer_Interaction_term_beta","lmer_Interaction_term_dpval","COUNTED","ALT","freq_counted")]
  colnames(myreQTLs) <- c("Condition","Gene","SNP","RSID","Chr","Position","LD_othertop","Keep","Beta_resting","Beta_interaction","Interaction_term_pval","Beta_resting_lmer","Beta_interaction_lmer","PermPval","Counted_allele","Other_allele","Freq_counted")
  write.csv(myreQTLs, paste0("ReQTLs_eGenes_",ifelse(celltype == "m","monocytes","Tcells"),"_Table.csv"), row.names = F)
}







  
  