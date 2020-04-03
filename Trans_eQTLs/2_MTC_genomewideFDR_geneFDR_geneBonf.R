#----------------------------------------------
# 2018-09-04
# Aim: to apply multiple testing correction.
#
# Trans-eQTLs: on different chromosomes.
# 
# Multiple testing:
# (1) Genome-wide FDR correction
# (2) Gene-level FDR correction
# (3) Gene-level Bonferroni correction
#----------------------------------------------
library(foreach)
library(doMC)
#registerDoMC(cores = 4)
options(stringsAsFactors = F)

setwd("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/")

# Gene location
geneloc <- read.table("/projects/qinqinhuang/CAS/Expression_Data/clean_data/Gene_location.txt", header = T)
nrow(geneloc)   # 13,109 genes
# SNP location
snploc <- read.table("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/snpsloc.txt", header = T)
nrow(snploc)   # 4,326,483 SNPs


# Number of tests:
# number of SNPs on each chr
nSNP_chr <- data.frame(table(snploc$CHR))
# number of genes on each chr
ngene_chr <- data.frame(table(geneloc$Chr))
# Number of tests
ngene_chr$ntests <- nrow(snploc) - nSNP_chr$Freq
ntests_diffchr <- sum(as.numeric(ngene_chr$Freq) * as.numeric(ngene_chr$ntests))    # 53,835,949,427


# Four conditions
foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA")) %do% {
  
  # Directory for this condition
  cat("  * Working on", celltype, treatment, "...\n")
  setwd(paste0("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/",celltype,"_",treatment,"_gene/"))
  
  # Read associations with nominal pval <= 1e-5
  asso <- read.table("All_tests_sig_1e-05.txt", header = T)
  
  
  #---------- Trans eQTL definition: (B) on a different chromosome ----------
  tests_diffchr <- asso[which(asso$label_chr == "different"),]

  # (1) Genome-wide FDR
  tests_diffchr$FDR_diffchr <- p.adjust(tests_diffchr$pvalue, "BH", n = ntests_diffchr)
  tests_diffchr_gwFDRsig <- tests_diffchr[which(tests_diffchr$FDR_diffchr <= 0.05),]
  cat("    (1) Gemone-wide FDR:", length(unique(tests_diffchr_gwFDRsig$gene)), "eGenes with", nrow(tests_diffchr_gwFDRsig), "associations:", paste(unique(tests_diffchr_gwFDRsig$gene), collapse = ", "), " \n")
  write.table(tests_diffchr_gwFDRsig, "Trans_diffchr_genomewideFDR_sig.txt", quote = F, sep = "\t", row.names = F)
  
  # (2) Gene-level FDR correction
  # Multiply p-values by 1e6 to adjust multiple SNPs tested.
  tests_diffchr$adjpval <- tests_diffchr$pvalue * (1e6)
  # Best association for each gene
  tests_diffchr_topSNP <- tests_diffchr[!duplicated(tests_diffchr$gene),]
  tests_diffchr_topSNP$gene_FDR <- p.adjust(tests_diffchr_topSNP$adjpval, "BH", n = nrow(geneloc))
  
  # (3) Gene-level Bonferroni correction
  tests_diffchr_topSNP$gene_Bonf <- tests_diffchr_topSNP$adjpval * nrow(geneloc)
  
  tests_diffchr_topSNP_sig <- tests_diffchr_topSNP[which(tests_diffchr_topSNP$gene_FDR <= 0.05),]
  cat("    (2) Gene-level FDR:", nrow(tests_diffchr_topSNP_sig), "eGenes:", paste(tests_diffchr_topSNP_sig$gene, collapse = ", "), " \n")
  cat("    (3) Gene-level Bonferroni:", sum(tests_diffchr_topSNP_sig$gene_Bonf <= 0.05), "eGenes:", paste(tests_diffchr_topSNP_sig$gene[which(tests_diffchr_topSNP_sig$gene_Bonf <= 0.05)], collapse = ", "), " \n")
  write.table(tests_diffchr_topSNP_sig, "Trans_diffchr_genelevelFDR_Bonf_sig.txt", quote = F, sep = "\t", row.names = F)
  
  cat("----------- \n"); return(NULL)
}



# Merge all result from 4 conditions
# Genome-wide FDR
# Different chromosomes:
gwfdr <- foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA"), .combine = rbind) %do% {
  setwd(paste0("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/",celltype,"_",treatment,"_gene"))
  asso_diffchr <- read.table("Trans_diffchr_genomewideFDR_sig.txt", header = T)
  asso_diffchr$FDR <- NULL
  asso_diffchr$celltype <- celltype
  asso_diffchr$treatment <- treatment
  return(asso_diffchr) 
}
colnames(gwfdr)[which(colnames(gwfdr) == "FDR_diffchr")] <- "FDR_diffchr_per_condition"
write.table(gwfdr, "/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/transeQTLs/GenomewideFDR_diffchr.txt", quote = F, sep = "\t", row.names = F)






