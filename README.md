# CAS_eQTL
This repository contains the scripts I used for the CAS (r)eQTL project. 

Preprint: https://www.biorxiv.org/content/10.1101/683086v1

Huang, Q.Q., Tang, H.H.F., Teo, S.M., Ritchie, S.C., Nath, A.P. et al. Neonatal genetics of gene expression reveal the origins of autoimmune and allergic disease risk. Under review. 

## eQTL mapping
###	1_Input_data_MatrixeQTL.R
Generate input files for Matrix eQTL.
###	2_Run_MatrixeQTL.R
Perform eQTL mapping using Matrix eQTL.
###	3_Run_eigenMT_using_all_CAS_genotype.R
Perform eigenMT to estimate the number of independent SNPs/tests in the 2-Mb window for each gene.
###	4_Multiple_testing_correction.R
Correct for multiple testing (local correction: eigenMT; global correction: BH FDR) and get the list of significant cis-eQTLs.
###	5_LD_top_eSNPs_two_conditions.R
For eGenes that were significant in both resting and stimulated cells (monocytes or T cells), calculate the LD correlation between the two top eSNPs from these two conditions. (This is actually for response eQTL mapping)
###	6_Conditional_analysis.R
Two stage conditional analysis.
###	7_Calculate_LD_eSNPs_topeSNP.R
For each significant eQTL SNP, calculate the LD correlation with the top eSNP of the corresponding gene.

## Response eQTL mapping
###	1_Interaction_test_topeSNPs_permutation.R
Perform interaction tests on top eSNPs and run permutations.
###	2_Multiple_testing_correction.R
Apply BH-FDR on permutation adjusted P-values to correct for multiple testing and get significant reQTLs.
###	3_ReQTL_SuppleTable.R
Gather some more information on the reQTLs for supplemental table 7 and 8.
## Trans-eQTL mapping and mediation analysis
###	1_Run_MatrixeQTL_genomewide.R
Perform genome-wide eQTL mapping using Matrix eQTL (output and write tests with P-value ≤1e-5).
###	2_MTC_genomewideFDR_geneFDR_geneBonf.R
Multiple testing using three ways: genome-wide FDR, gene-level FDR, gene-level Bonferroni.
###	3_Mediation_analysi_cis_trans.R
Get all mediation trios (eQTL–cis-eGene–trans-eGene) and perform mediation analysis.
###	4_Significant_Mediations.R
Apply BH FDR controlling procedure to correct for multiple testing.
###	5_LD_eSNPs_top_transeSNP_SNPinfo_SuppleTable.R
Calculate the LD correlation between each trans-eQTL SNP and the corresponding top SNP, and gather information for supplemental table 10.

## Colocalisation analysis
###	1_Subset_GWAS_loci_candidate.R
Find loci where GWAS and eQTL signals have overlap, and prepare input data for coloc.
###	2_Perform_coloc_test.R
Perform colocalisation analysis using the R package coloc.

## Mendelian randomisation analysis
###	0_1_MRbase_ciseQTLs_eGene_immunediseases.R
Get genetic instrumental variables and perform MR analysis using the R package TwoSampleMR.
###	0_2_MR_results_convert_list_to_dataframe.R
Extract MR test statistics.
###	1_MR_using_IVs_selected_by_TwoSampleMR.R
Perform the MR analysis again using Scott’s codes. The genetic IVs were harmonised and selected by the TwoSampleMR package.
###	2_Significant_MR_results_3_out_of_4_SuppleTable.R
Select significant causal associations and generate a supplementary table.
###	MR_functions_ScottRitchie_2019_04_17.R
I used this script from Scott Ritchie (https://github.com/sritchie73) to perform MR analysis and generate dose responsive curves.

