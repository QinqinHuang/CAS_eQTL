#----------------------------------------------
# 2018-11-19
# Aim: to generate a supplemental table 
#  Calculate LD between each eSNP and the 
# corresponding top trans-eSNP.
#  Also add SNP info.
#  Add cis-eGene if the trans-eSNP is also a 
# cis-eQTL.
#
# Sort by condition, and top eSNP P-value.
#----------------------------------------------
library(foreach)
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")
setwd("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/transeQTLs")

# Load trans-eQTLs
transeQTLs <- read.table("GenomewideFDR_diffchr.txt", header = T)

# Top eSNPs in four conditions
topeSNPs <- transeQTLs[order(transeQTLs$pvalue),]
topeSNPs <- topeSNPs[!duplicated(topeSNPs[,c("gene","celltype","treatment")]),]

# A list of all top eSNPs across four conditions
write.table(data.frame(snps = unique(topeSNPs$snps)), "Unique_topeSNP_list.txt", quote = F, row.names = F, col.names = F)

# Keep all SNPs that are in LD r2 ≥0.8 with the top eSNP
system("plink1.9 --bfile /projects/qinqinhuang/CAS/Analysis/eQTL_mapping/cas_imp_135ind_filtered_maf10 --ld-snp-list Unique_topeSNP_list.txt --r2 --ld-window 2000 --ld-window-r2 0.01 --out Top_eSNPs_LD001")

# LD r2 0.01
ld001 <- read.table("Top_eSNPs_LD001.ld", header = T)

# Add LD R2 with top eSNP
colnames(topeSNPs)[2] <- "topeSNP"
# Top eSNP for each gene
transeQTLs <- merge(transeQTLs, topeSNPs[,c("gene","celltype","treatment","topeSNP")], by = c("gene","celltype","treatment"))
transeQTLs <- merge(transeQTLs, ld001[,c("SNP_A","SNP_B","R2")], by.x = c("topeSNP","snps"), by.y = c("SNP_A","SNP_B"), all.x = T)
colnames(transeQTLs)[ncol(transeQTLs)] <- "R2_topeSNP"


# Load annotations for SNPs tested in CAS
SNPanno <- read.table("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/SNP_pos_rsID_counted_freq_freq1kG.txt", header = T)
transeQTLs <- merge(transeQTLs, SNPanno[,c("snps","rsid","chr","position","COUNTED","ALT","freq_counted")], by = "snps")

# Load cis-eQTLs
ciseQTLs_all <- load_cis_eQTLs()
names(ciseQTLs_all) <- c("m_Ctrl","m_LPS","t_Ctrl","t_PHA" )
transeQTLs_cis <- foreach(ii = 1:nrow(transeQTLs), .combine = rbind) %do% {
  dd <- transeQTLs[ii,]
  ciscurr <- ciseQTLs_all[[paste0(dd$celltype,"_",dd$treatment)]]
  dd$cis_eGene <- paste(unique(ciscurr[which(ciscurr$snps == dd$snps), "gene"]), collapse = ";")
  return(dd)
}

# Sort by condition and P-value 
transeQTLs_cis <- transeQTLs_cis[order(transeQTLs_cis$celltype, transeQTLs_cis$treatment, transeQTLs_cis$pvalue),]
# Trans-eSNPs for the same gene should stay together
unique_gene <- unique(transeQTLs_cis[,c("gene","celltype","treatment")])
dd <- foreach(ii = 1:nrow(unique_gene), .combine = rbind) %do% {
  unique_gene[ii,]
  returnme <- transeQTLs_cis[which(transeQTLs_cis$gene == unique_gene$gene[ii] & transeQTLs_cis$celltype == unique_gene$celltype[ii] & transeQTLs_cis$treatment == unique_gene$treatment[ii]),]
  returnme <- returnme[order(returnme$pvalue),]
  return(returnme)
}

# Column names
dd <- dd[,c("celltype","treatment","rsid","snp_chr","snp_pos","gene","gene_chr","gene_TSS","pvalue","beta","R2_topeSNP","COUNTED","ALT","freq_counted","cis_eGene")]
colnames(dd) <- c("Cell","Treatment","RSID","SNP_chr","SNP_pos","Trans_eGene","Chr","TSS","Pvalue","Beta","R2_topeSNP","Counted_allele","Other_allele","Freq_counted","Cis_eGene")
dd$Cell <- ifelse(dd$Cell == "m", yes = "Monocyte", no = "T cell")
dd$Treatment[which(dd$Treatment == "Ctrl")] <- "Resting"
write.csv(dd, "GenomewideFDR_diffchr_SNPinfo_ciseGene_table.csv", row.names = F)








