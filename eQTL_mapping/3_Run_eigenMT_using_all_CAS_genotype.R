#----------------------------------------------
# 2018-05-05
# 1. Prepare Genotype/SNP pos file for each 
# chromosome (eigenMT input)
# 2. Run eigenMT.
#
# Davis et al. recommended using genotype data
# for all individuals and running eigenMT once
# to get the estimated number of effective tests;
# don't have to run it for each of the condition.
# We use the same set of SNPs that were tested 
# in eQTL mapping (MAF ≥10% in 135 individuals).
#
# Note that eigenMT works on each chromosome 
# separately, so input files (QTL, GEN, GENPOS) 
# should be splitted.
#----------------------------------------------

library(foreach)
library(doMC)
#registerDoMC(cores = 10)
options(stringsAsFactors = F)

# Working directory
setwd("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/eigenMT")

# Get the genotyp data for SNPs that were tested in eQTL mapping in all 215 CAS individuals.
# The list of SNPs
system("cut -f 2 /projects/qinqinhuang/CAS/Analysis/eQTL_mapping/cas_imp_135ind_filtered_maf10.bim > SNPs_tested_eQTL.txt")
# Extract the list of SNPs from the filtered genotype dataset of all 215 individuals
system("plink1.9 --bfile /projects/qinqinhuang/CAS/Genotype_Data/process_data/Michigan_imputation_server_results/cas_auto_sub_filt-updated_auto_chr_2018_01_11_01_15/filtered_geno/cas_michigan_imp_filtered --extract SNPs_tested_eQTL.txt --allow-no-sex --make-bed --out genotype_in_all_cas_imp_215ind")
system("rm SNPs_tested_eQTL.txt")

# Genotype and SNP position files
# Split into 22 chromosomes
if(!file.exists("input_Geno")) {dir.create("input_Geno")}; setwd("input_Geno")
nothing <- foreach(chr = 1:22) %dopar% {
  system(paste0("plink1.9 --bfile ../genotype_in_all_cas_imp_215ind --chr ",chr," --allow-no-sex --make-bed --out TEMP_chr",chr))
  # SNP dosage
  system(paste0("plink1.9 --bfile TEMP_chr",chr," --recode A-transpose --recode-allele ../../HRC_alt_allele.txt --out TEMP2_chr",chr))
  system(paste0("cut -f 2,7- TEMP2_chr",chr,".traw > SNP_chr",chr,".txt"))
  # SNP location
  system(paste0("awk '{print$2,$1,$4}' TEMP2_chr",chr,".traw > snpsloc_chr",chr,".txt"))
  return(NULL)
}
system("rm TEMP*")


#----- run eigenMT -----
# Run eigenMT to estimate the empirical number of independent tests.
setwd("/projects/qinqinhuang/CAS/Analysis/eQTL_mapping/eigenMT")
if(!file.exists("output_eigenMT")) {dir.create("output_eigenMT")}

# Run eigenMT for 22 chromosomes
eigenMTout <- foreach(chr = 1:22, .combine = rbind) %dopar% {
  cat("  Running eigenMT for Chr",chr,"...\n")
  system(paste0("time python ~/software/eigenMT/eigenMT_QH.py --CHROM ", chr, " --QTL ./input_QTL/chr", chr, "_MatrixEQTL_all_cis_tests_pval.txt --GEN ./input_Geno/SNP_chr", chr, ".txt --GENPOS ./input_Geno/snpsloc_chr",chr,".txt --PHEPOS /projects/qinqinhuang/CAS/Expression_Data/clean_data/Gene_location.txt --OUT output_eigenMT/chr", chr, "_eigenMT_output.txt"))
  
  # Read the output
  dd <- read.table(paste0("output_eigenMT/chr",chr,"_eigenMT_output.txt"), header = T)
  dd <- dd[which(complete.cases(dd$TESTS)), c("gene","TESTS")]
  return(dd)
}

write.table(eigenMTout, file = paste0("eigenMT_output_gene_level.txt"), quote = F, sep = "\t", row.names = F)


# 57min

