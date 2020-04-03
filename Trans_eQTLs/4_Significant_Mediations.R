#----------------------------------------------
# 2018-09-12
# Aim: to get the list of significant mediations
#----------------------------------------------
library(foreach)
library(ggplot2)
options(stringsAsFactors = F)
setwd("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/Mediation_analysis/")

# Read mediation analysis results
med_out <- read.table("Mediation_analysis_results_n10000.txt", header = T)
head(med_out)

# Move important columns in the front
cn <- c("celltype","treatment","gene","gene_cis","ACME_pval_boot","ACME_pval_bootbca","ACME_pval_quasiBayeMonteCarlo","ACME_boot","ACME_bootbca","ACME_quasiBayeMonteCarlo")
med_out <- med_out[,c(cn, setdiff(colnames(med_out), cn))]
med_out <- med_out[order(med_out$celltype, med_out$treatment, med_out$gene_cis, med_out$gene),]

# Keep unique mediation tests
med_out <- med_out[!duplicated(med_out[,1:4]),]


#----------------------------------------------
# Define trans-eQTLs as those located on different chromosomes
transeQTLs <- read.table("/projects/qinqinhuang/CAS/Analysis/Trans_eQTLs/transeQTLs/GenomewideFDR_diffchr.txt", header = T)

# Keep significant eGenes identified when focusing on tests where SNPs and genes were located on different chromosomes.
transeGene <- unique(transeQTLs[,c("gene","celltype","treatment")])
transeGene$trans <- T
dd <- merge(med_out, transeGene, by = c("celltype","treatment","gene"))
dd <- dd[order(dd$celltype, dd$treatment, dd$gene_cis),]


# Apply BH-FDR on all tests
FDR <- dd
FDR$FDR <- p.adjust(FDR$ACME_pval_boot, "BH")
cat("  Significant at FDR 5% - ACME_pval_boot:", sum(FDR$FDR <= 0.05), "\n")

# Significant results
sigresults <- FDR[, c("celltype","treatment","snps","gene_cis","gene","beta","beta_cis","ciseGene_beta","ACME_boot","PropMed_boot","ACME_pval_boot","FDR")]
write.table(sigresults, "Sig_Medeiation_bootstrap_FDR_transeQTLs_diffchr.txt", quote = F, sep = "\t", row.names = F)



