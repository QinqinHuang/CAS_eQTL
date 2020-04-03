#----------------------------------------------
# 2018-12-17
# Aim: to extract the results from a list and 
# convert to data frames.
#----------------------------------------------
library(foreach)
source("/projects/qinqinhuang/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

# Working directory
setwd("/projects/qinqinhuang/CAS/Analysis/Mendelian_Randomisation_ciseGene_diseases/")

# A list of outcome diseases tested in MR
gwas_phenotype <- read.table("/projects/qinqinhuang/CAS/Analysis/Colocalisation/GWAS_phenotypes_files.txt", header = T, sep = "\t", quote = "")
gwas_phenotype <- gwas_phenotype[which(gwas_phenotype$coloc_input == "se"),]


# Each MR test had four elements in the list:
# [1] exposure, outcome, condition
# [2] input data frame
# [3] MR output
# [4] Direction test

#------ MR test results -----
allMR <- foreach(gwasfilename = gwas_phenotype$file_name, .combine = rbind) %do% {
  
  # Four cell conditions
  asso_gwas <- foreach(cc = c("m_resting","m_LPS","t_resting","t_PHA"), .combine = rbind) %do% {
    
    filename <- paste0("./output/",gsub(".txt","",gwasfilename),"_Celltype_",cc,".RDS")
    if(!file.exists(filename)) {cat(" *Output not available for", gwasfilename, "\n"); return(NULL)}
    MR_results <- readRDS(paste0("./output/",gsub(".txt","",gwasfilename),"_Celltype_",cc,".RDS"))
    if(is.null(MR_results)) {return(NULL)}
    
    curr <- foreach(ii = 1:(length(MR_results)/4), .combine = rbind) %do% {
      dd <- MR_results[[ii*4-1]][,-1:-2]
      if(nrow(dd) == 0) {return(NULL)}
      dd$condition <- MR_results[[ii*4-3]][3]
      return(dd)
    } # End of one RDS file
    
    if(!all(curr$condition == cc)) {cat(" * wrong cell condition.\n     ", cc, "\n")}
    return(curr)
    
  } # End of four cell conditions
  
  return(asso_gwas)
}

# Number of genetic IVs
table(allMR$nsnp)
#    1     2     3     4     5     6     7     8 
#56939  5449  3440   860   300   120    30    65 
write.table(allMR, "All_MR_results.txt", quote = F, sep = "\t", row.names = F)


#------ Results of the direction test -----
directiontest  <- foreach(gwasfilename = gwas_phenotype$file_name, .combine = rbind) %do% {
  
  # Four cell conditions
  dir_gwas <- foreach(cc = c("m_resting","m_LPS","t_resting","t_PHA"), .combine = rbind) %do% {
    
    filename <- paste0("./output/",gsub(".txt","",gwasfilename),"_Celltype_",cc,".RDS")
    if(!file.exists(filename)) {cat(" *Output not available for", gwasfilename, "\n"); return(NULL)}
    MR_results <- readRDS(paste0("./output/",gsub(".txt","",gwasfilename),"_Celltype_",cc,".RDS"))
    if(is.null(MR_results)) {return(NULL)}
    
    curr <- foreach(ii = 1:(length(MR_results)/4), .combine = rbind) %do% {
      dd <- MR_results[[ii*4]][,-1:-2]
      if(nrow(dd) == 0) {return(NULL)}
      dd$condition <- MR_results[[ii*4-3]][3]
      if(nrow(dd) != 1) {cat(" * N of rows of direction test is not 1.\n")}
      return(dd)
    } # End of one RDS file
    
    if(!all(curr$condition == cc)) {cat(" * wrong cell condition.\n     ", cc, "\n")}
    return(curr)
    
  } # End of four cell conditions
  
  return(dir_gwas)
}

table(directiontest$correct_causal_direction)
# All direction tests show that the direction is correct


