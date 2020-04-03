#----------------------------------------------
# 2018-10-02
# Aim: to calculate LD r2 between an eGene's 
# top eSNPs from two conditions. 
# (For interaction test for reQTLs.)
#----------------------------------------------
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

setwd("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/LD_top_SNP_conditions/")


# ---- Use the top SNPs identified from the initial (unconditional) eQTL scan -----
for(celltype in c("m","t")) {
  treatment <- ifelse(celltype == "m", "LPS", "PHA")
  
  # Read significant eGenes
  Ctrl <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_Ctrl_gene/sigeAsso_eigenMTBH.txt"), header = T)
  Sti <- read.table(paste0("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_gene/sigeAsso_eigenMTBH.txt"), header = T)
  
  # One top eSNP for each gene
  Ctrl <- Ctrl[order(Ctrl$pvalue),]
  Ctrl_top <- Ctrl[!duplicated(Ctrl$gene),]
  Sti <- Sti[order(Sti$pvalue),]
  Sti_top <- Sti[!duplicated(Sti$gene),]
  
  # eGenes that are significant in either condition
  # NA: not significant in this condition
  eGene_list <- merge(Ctrl_top[,c("gene","snps","pvalue","beta")], Sti_top[,c("gene","snps","pvalue","beta")], by = "gene", all = T)
  colnames(eGene_list) <- c("gene","snps_Ctrl","pvalue_Ctrl","beta_Ctrl","snps_Sti","pvalue_Sti","beta_Sti")
  
  # For eGenes that have significant eQTLs in both conditions,
  # calculate LD between the two top eSNPs.
  eGene_list$LDr2 <- NA
  for(ii in 1:nrow(eGene_list)) {
    if(!is.na(eGene_list$snps_Ctrl[ii]) & !is.na(eGene_list$snps_Sti[ii])) {
      # Run Plink to calculate pairwise LD 
      top_eSNP_list <- c(eGene_list$snps_Ctrl[ii], eGene_list$snps_Sti[ii])
      # If the same top eSNP
      if(length(unique(top_eSNP_list)) == 1) {
        eGene_list$LDr2[ii] <- 1
      } else {
        write.table(data.frame(snps = top_eSNP_list), "Temp_SNP_list.txt", row.names = F, quote = F, col.names = F)
        system(paste0("plink1.9 --bfile /projects/qinqinhuang/CAS/Analysis/eQTL_mapping/cas_imp_135ind_filtered_maf10 --r2 square --extract Temp_SNP_list.txt --out Temp_LD_r2_top_eSNPs"))
        LD_mat <- read.table("Temp_LD_r2_top_eSNPs.ld", header = F)
        if(nrow(LD_mat) != 2) {cat("Wrong", ii, "\n")}
        eGene_list$LDr2[ii] <- LD_mat[1,2]
      }
    }
  }
  eGene_list$LDr2 <- as.numeric(eGene_list$LDr2)
  
  write.table(eGene_list, paste0("Top_eSNP_eGene_LD_r2_",celltype,".txt"), quote = F, sep = "\t", row.names = F)
  
  system("rm Temp_*")
}


