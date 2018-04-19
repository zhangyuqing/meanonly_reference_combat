########  Copied everything from /restricted/projectnb/pathsig/work/dfj/20160215_combat_meanonlypaper/results
########  Missing a script to calculate correlation, so this script calculates correlations in table 2 of paper

#### Correlations
rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/refcombat_review_201803/compare_sva_ruv")
source("/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/scripts/Key_ASSIGN_functions_balancedsig_v2.R")

test_name_lst <- c("TCGA", "ICBP")
#batch_methods <- c("original", "meanonly", "reference") 
batch_methods <- c("original", "meanonly", "reference", "SVA", "RUV")
corr_res <- matrix(NA, nrow=length(batch_methods), ncol=4,
                   dimnames=list(batch_methods, c("prote_icbp", "prote_tcga", "Erlotinib", "GSK1120212")))

for(test_name in test_name_lst){
  for(method_name in batch_methods){
    ####  Read in ASSIGN prediction
    pred_res_dir <- sprintf("%s/%s/egfr_50_gene_list/adap_adap_single/pathway_activity_testset.csv", test_name, method_name)
    pred_res <- read.csv(pred_res_dir, row.names="X")
    
    if(test_name=="TCGA"){
      ####  TCGA protein expression
      ##  Read in proteomics data
      tcga_prote_dir <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/Datasets/TCGA-BRCA-RBN.csv"
      tcga_prote <- read.csv(tcga_prote_dir)[, c("TCGA_patient_barcode", "EGFR")]
      
      ##  Match samples
      # sample list to be considered
      sample_dics <- read.csv("../tcga_sample_dics.csv", as.is=TRUE)
      # take subset corresponding to sample list
      tcga_prote_subset <- tcga_prote[match(sample_dics$short_name, tcga_prote[,1]), ]
      pred_res_subset <- pred_res[sample_dics$long_name, ]
      
      ##  Compute correlation
      corr_res[method_name, "prote_tcga"] <- cor(pred_res_subset, tcga_prote_subset[,"EGFR"], method="spearman")
    }else if(test_name=="ICBP"){
      ####  ICBP protein expression
      ##  Read in protein expression 
      icbp_prote_dir <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/correlation_data/proteomics.txt"
      icbp_prote <- read.table(icbp_prote_dir)
      
      ##  Match samples
      for(i in grep("^[0-9]",rownames(pred_res))){rownames(pred_res)[i] <- paste("X",rownames(pred_res)[i], sep="")}
      # take subset
      overlap_cellnames <- intersect(rownames(icbp_prote), rownames(pred_res))
      pred_res_subset <- pred_res[overlap_cellnames, ]
      icbp_prote_subset <- icbp_prote[overlap_cellnames, "EGFR"]
      
      ##  Compute correlation 
      corr_res[method_name, "prote_icbp"] <- cor(pred_res_subset, icbp_prote_subset, method="spearman")
      
      
      ####  ICBP drug response
      ##  Read in drug response
      icbp_drug_dir <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/correlation_data/ICBP_drugs.txt"
      icbp_drug <- read.delim(icbp_drug_dir, header=1, sep='\t',row.names=1)
      pred_res <- read.csv(pred_res_dir, row.names="X")
      ##  Match samples: take subset
      merged_dat <- merge_drop(pred_res, icbp_drug)
      
      # overlap_cellnames <- intersect(rownames(icbp_drug), rownames(pred_res))
      # icbp_drug_subset <- icbp_drug[overlap_cellnames, c("Erlotinib", "GSK1120212")]
      # pred_res_subset <- pred_res[overlap_cellnames, 1]
      
      ##  Compute correlation 
      corr_res[method_name, "Erlotinib"] <- cor(merged_dat[,"egfr"], merged_dat[,"Erlotinib"], use="complete.obs", method="spearman")
      corr_res[method_name, "GSK1120212"] <- cor(merged_dat[,"egfr"], merged_dat[,"GSK1120212"], use="complete.obs", method="spearman")
    }
  }
}
print(round(corr_res,3))



#### Intersect of signature genes
#batch_methods <- c("original", "meanonly", "reference")
batch_methods <- c("original", "meanonly", "reference", "SVA", "RUV")
n_overlap <- rep(NA, length(batch_methods)); names(n_overlap) <- batch_methods

for(method_name in batch_methods){
  tcga_sig_dir <- sprintf("TCGA/%s/egfr_50_gene_list/adap_adap_single/posterior_delta.csv", method_name)
  icbp_sig_dir <- sprintf("ICBP/%s/egfr_50_gene_list/adap_adap_single/posterior_delta.csv", method_name)

  tcga_sigs <- read.csv(tcga_sig_dir, as.is=TRUE)[,1]
  icbp_sigs <- read.csv(icbp_sig_dir, as.is=TRUE)[,1]
  n_overlap[method_name] <- length(intersect(tcga_sigs, icbp_sigs))
}
print(n_overlap)
print(round(n_overlap/50, 2))

table2 <- data.frame(n_overlap, round(corr_res,3))
write.csv(table2, file="corr_summary.csv", 
          quote=FALSE, row.names=TRUE)
