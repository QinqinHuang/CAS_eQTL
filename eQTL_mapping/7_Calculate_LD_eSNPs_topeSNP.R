#----------------------------------------------
# 2018-11-09
# Aim: calculate the LD between each eSNP and 
# the top eSNP of the gene.
#----------------------------------------------
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")
library(foreach)
setwd("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/")

# Load cis-eQTLs
ciseQTLs <- load_cis_eQTLs(type = "other")

# Top eSNPs in four conditions
topeSNPs <- lapply(ciseQTLs, function(x) {
  x[which(!duplicated(x$gene)),]
})

# A list of all top eSNPs across four conditions
topeSNPlist <- foreach(ii = 1:4, .combine = c) %do% {
  return(topeSNPs[[ii]]$snps)
}
topeSNPlist <- unique(topeSNPlist)
write.table(data.frame(snps = topeSNPlist), "Unique_topeSNP_list.txt", quote = F, row.names = F, col.names = F)

# Keep all SNPs that are in LD r2 ≥0.8 with the top eSNP
system("plink1.9 --bfile cas_imp_135ind_filtered_maf10 --ld-snp-list Unique_topeSNP_list.txt --r2 --ld-window 2000 --ld-window-r2 0.01 --out Top_eSNPs_LD001")

# LD r2 0.01
ld001 <- read.table("Top_eSNPs_LD001.ld", header = T)


# Add LD R2 with top eSNP
ciseQTLs_LD <- lapply(ciseQTLs, function(x) {
  # Top eSNP for each eGene
  top <- x[which(!duplicated(x$gene)), c("gene","snps")]
  colnames(top)[2] <- "topeSNP"
  x <- merge(x, top, by = "gene")
  
  # LD r2 if ≥0.01
  x <- merge(x, ld001[,c("SNP_A","SNP_B","R2")], by.x = c("topeSNP","snps"), by.y = c("SNP_A","SNP_B"), all.x = T)
  x <- x[order(x$pvalue), ]
  colnames(x)[ncol(x)] <- "R2_topeSNP"
  rownames(x) <- 1:nrow(x)
  return(x)
})
saveRDS(ciseQTLs_LD, "Significant_cis_eQTLs_four_conditions_LDtopeSNP.RDS")





