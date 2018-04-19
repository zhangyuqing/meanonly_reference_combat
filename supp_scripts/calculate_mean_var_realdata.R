######## Calculate mean and variance from real data: signature, TCGA, bladderbatch
rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/refcombat_review_201803/")
lapply(c("BatchQC", "data.table", "ggplot2", "reshape"), require, character.only = TRUE)


###  EGFR signature data
data_dir <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/signatures/18_GFP_EGFR_TPMlog2.txt"
sig_dir <- "/restricted/projectnb/pathsig/work/dfj/20160215_combat_meanonlypaper/results/TCGA/original/egfr_50_gene_list/adap_adap_single/posterior_delta.csv"

EGFR_dat <- as.matrix(read.table(data_dir, sep='\t', header=1, row.names=1))
## normalize
# EGFR_dat_norm <- gnormalize(EGFR_dat)
# rownames(EGFR_dat_norm) <- rownames(EGFR_dat); colnames(EGFR_dat_norm) <- colnames(EGFR_dat)
# EGFR_dat_norm <- EGFR_dat_norm[apply(EGFR_dat_norm,1,function(x){!any(is.na(x))}), ]
## don't normalize
EGFR_dat_norm <- EGFR_dat[apply(EGFR_dat, 1, var)!=0, ]  # remove genes with 0 variance

EGFR_assignres <- read.csv(sig_dir, as.is=TRUE)
EGFR_siglst <- EGFR_assignres[,1]  ## DE genes
EGFR_siglst_up <- EGFR_siglst[EGFR_assignres[,4]>=0]
EGFR_siglst_down <- EGFR_siglst[EGFR_assignres[,4]<0]

sigUP_off <- EGFR_dat_norm[EGFR_siglst_up, grep("Control",colnames(EGFR_dat_norm))]
sigUP_on <- EGFR_dat_norm[EGFR_siglst_up, grep("EGFR",colnames(EGFR_dat_norm))]
sigDOWN_off <- EGFR_dat_norm[EGFR_siglst_down, grep("Control",colnames(EGFR_dat_norm))]
sigDOWN_on <- EGFR_dat_norm[EGFR_siglst_down, grep("EGFR",colnames(EGFR_dat_norm))]
nonsig_off <- EGFR_dat_norm[setdiff(rownames(EGFR_dat_norm),EGFR_siglst), grep("Control",colnames(EGFR_dat_norm))]
nonsig_on <- EGFR_dat_norm[setdiff(rownames(EGFR_dat_norm),EGFR_siglst), grep("EGFR",colnames(EGFR_dat_norm))]

mean1_up_off <- mean(sigUP_off); mean1_up_on <- mean(sigUP_on)
mean1_down_off <- mean(sigDOWN_off); mean1_down_on <- mean(sigDOWN_on) 
mean1_nonsig_off <- mean(nonsig_off); mean1_nonsig_on <- mean(nonsig_on)

print(round(mean1_up_on - mean1_up_off, 2))
print(round(mean1_down_on - mean1_down_off, 2))
print(round(mean1_nonsig_on - mean1_nonsig_off, 2))
## for normalized:
## 0.25 change between sig genes
## 0.12 change between non-sig genes

var1_up_off <- mean(apply(sigUP_off,1,var)); var1_up_on <- mean(apply(sigUP_on,1,var))
var1_down_off <- mean(apply(sigDOWN_off,1,var)); var1_down_on <- mean(apply(sigDOWN_on,1,var))
var1_nonsig_off <- mean(apply(nonsig_off,1,var)); var1_nonsig_on <- mean(apply(nonsig_on,1,var))
print(round(var1_up_on - var1_up_off, 2))
print(round(var1_down_on - var1_down_off, 2))
print(round(var1_nonsig_on - var1_nonsig_off, 2))
## for normalized:
## 0.08 change between sig genes
## 0.18 change between non-sig genes

save(mean1_up_off, mean1_up_on, mean1_down_off, mean1_down_on, mean1_nonsig_off, mean1_nonsig_on,
     var1_up_off, var1_up_on, var1_down_off, var1_down_on, var1_nonsig_off, var1_nonsig_on,
     file="meanvar_egfr.rda")



###  TCGA 
rm(list=ls())
tcga_dir <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/test_data/PANCAN24_BRCA_1119_TPMlog2.txt"
tcga_dat <- data.frame(fread(tcga_dir), check.names=F, row.names=1)
tcga_prote_dir <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/Datasets/TCGA-BRCA-RBN.csv"
tcga_prote_egfr <- read.csv(tcga_prote_dir, as.is=TRUE)[,c("TCGA_patient_barcode", "EGFR")]

## normalize
# tcga_dat_norm <- gnormalize(as.matrix(tcga_dat))
# rownames(tcga_dat_norm) <- rownames(tcga_dat); colnames(tcga_dat_norm) <- colnames(tcga_dat)
## don't normalize
tcga_dat_norm <- as.matrix(tcga_dat[apply(tcga_dat,1,var)!=0, ])

##selected genes
# selected from heatmap 
#sel_genes <- c("E2F1", "MELK", "MYBL2", "MCM4", "MCM2", "E2F8", "GBP5")  
#heatmap of selected genes
# plt_tcga <- as.matrix(tcga_dat[sel_genes, ])
# plt_tcga_mlt <- melt(plt_tcga)
# plt_tcga_mlt$X2 <- factor(plt_tcga_mlt$X2, levels=tcga_prote_egfr[order(tcga_prote_egfr$EGFR), 1])
# ggplt_tcga <- ggplot(data = plt_tcga_mlt, aes(x=X2, y=X1, fill=value)) + 
#   geom_tile() +
#   scale_fill_gradientn(colors=colorRampPalette(c("white", "steelblue"))(3), breaks=c(0,5,10,15)) +
#   ggtitle("Selected genes in EGFR signature from TCGA")
# ggsave(filename=paste0("./heatmap_selgenes_TCGA.pdf"), 
#        plot=ggplt_tcga, width=340, height=80, units="mm", dpi=300)
# #dev.off()

##match proteomics barcode with sample names
duplicated_id <- names(table(tcga_prote_egfr[, 1]))[table(tcga_prote_egfr[, 1])>1]
tcga_prote_egfr <- tcga_prote_egfr[tcga_prote_egfr[,1] %in% setdiff(tcga_prote_egfr[, 1], duplicated_id), ]   # remove duplicate barcodes from proteomics
rownames(tcga_prote_egfr) <- tcga_prote_egfr[, 1]   # set the barcodes as rownames
short_name <- sapply(colnames(tcga_dat_norm), substr, start=1, stop=12)  # shorten sample names to barcodes
tcga_id_dics <- data.frame(short_name=as.character(short_name), long_name=names(short_name))
duplicated_shortnames <- names(table(tcga_id_dics$short_name))[table(tcga_id_dics$short_name)>1]
tcga_id_dics <- tcga_id_dics[tcga_id_dics$short_name %in% setdiff(tcga_id_dics$short_name, duplicated_shortnames), ]  # remove duplicate barcodes in sample names
# unique overlapping barcodes between sample name in gene expression and barcodes from proteomics
uni_overlap_barcodes <- intersect(tcga_id_dics$short_name, rownames(tcga_prote_egfr))
tcga_prote_egfr <- tcga_prote_egfr[uni_overlap_barcodes, ]
tcga_id_dics <- tcga_id_dics[match(uni_overlap_barcodes, as.character(tcga_id_dics$short_name)), ]
identical(as.character(tcga_id_dics$short_name), rownames(tcga_prote_egfr))
tcga_dat_norm <- tcga_dat_norm[, as.character(tcga_id_dics$long_name)]
# final check of sample names
dim(tcga_dat_norm); dim(tcga_prote_egfr)
identical(substr(colnames(tcga_dat_norm), 1, 12), rownames(tcga_prote_egfr))
identical(substr(colnames(tcga_dat_norm), 1, 12), tcga_prote_egfr[, "TCGA_patient_barcode"])
write.csv(tcga_id_dics, file="tcga_sample_dics.csv", row.names=FALSE, quote=FALSE)

##EGFR signature 50-gene list
sig_dir <- "/restricted/projectnb/pathsig/work/dfj/20160215_combat_meanonlypaper/results/TCGA/original/egfr_50_gene_list/adap_adap_single/posterior_delta.csv"
EGFR_assignres <- read.csv(sig_dir, as.is=TRUE)
EGFR_siglst <- EGFR_assignres[,1]  ## DE genes
EGFR_siglst_up <- EGFR_siglst[EGFR_assignres[,4]>=0]
EGFR_siglst_down <- EGFR_siglst[EGFR_assignres[,4]<0]

##proteomics data, bin EGFR protein expression to find patient groups
# bin EGFR protein expression in 5 groups
egfr_breaks <- c(-Inf, quantile(tcga_prote_egfr$EGFR, probs=1/6*(1:5)), Inf)   #c(-Inf, -1, -0.5, 0, 0.5, Inf)
group_ind <- findInterval(tcga_prote_egfr$EGFR, egfr_breaks)

# take subset of TCGA and calculate mean, variance
tcga_subset <- tcga_dat_norm[EGFR_siglst_up, ]  #tcga_dat_norm[sel_genes, ]
mean2_up_vec <- sapply(1:6, function(i){mean(as.matrix(tcga_subset[, group_ind==i]))})
var2_up_vec <- sapply(1:6, function(i){mean(apply(as.matrix(tcga_subset[, group_ind==i]), 1, var))})
rm(tcga_subset)

tcga_subset <- tcga_dat_norm[EGFR_siglst_down, ]  #tcga_dat_norm[sel_genes, ]
mean2_down_vec <- sapply(1:6, function(i){mean(as.matrix(tcga_subset[, group_ind==i]))})
var2_down_vec <- sapply(1:6, function(i){mean(apply(as.matrix(tcga_subset[, group_ind==i]), 1, var))})
rm(tcga_subset)

tcga_subset <- tcga_dat_norm[setdiff(rownames(tcga_dat_norm), EGFR_siglst), ]  #tcga_dat_norm[sel_genes, ]
mean2_nonsig_vec <- sapply(1:6, function(i){mean(as.matrix(tcga_subset[, group_ind==i]))})
var2_nonsig_vec <- sapply(1:6, function(i){mean(apply(as.matrix(tcga_subset[, group_ind==i]), 1, var))})
rm(tcga_subset)

png("qqplot.png", width=5, height=5, units="in", res=300)
qqnorm(tcga_subset[1000, group_ind==1])
qqline(tcga_subset[1000, group_ind==1], col="red")
dev.off()

save(mean2_up_vec, var2_up_vec, mean2_down_vec, var2_down_vec, mean2_nonsig_vec, var2_nonsig_vec,
     file="meanvar_tcga.rda")



###  Write out mean and variance estimates from signature and TCGA
rm(list=ls())
load("meanvar_egfr.rda")
load("meanvar_tcga.rda")
egfr_df <- data.frame(mean_ctrl=c(mean1_up_off, mean1_down_off, mean1_nonsig_off),
                      mean_egfr=c(mean1_up_on, mean1_down_on, mean1_nonsig_on),
                      var_ctrl=c(var1_up_off, var1_down_off, var1_nonsig_off),
                      var_egfr=c(var1_up_on, var1_down_on, var1_nonsig_on))
rownames(egfr_df) <- c("up", "down", "nonsig")
egfr_df <- t(egfr_df)
write.csv(round(egfr_df,3), file="meanvar_egfr.csv", quote=F)

tcga_df <- data.frame(up=c(mean2_up_vec, var2_up_vec),
                      down=c(mean2_down_vec, var2_down_vec),
                      nonsig=c(mean2_nonsig_vec, var2_nonsig_vec))
rownames(tcga_df) <- c(paste0("mean",1:6), paste0("var",1:6))
write.csv(round(tcga_df,3), file="meanvar_tcga.csv", quote=F)



###  Bladdarbatch
rm(list=ls())
setwd("~/Dropbox/Work/EvanJohnson/refComBat_forSubmittion_BMCbioinfo/review20180320/modsel_sensitivity")
lapply(c("bladderbatch", "sva", "limma"), require, character.only = TRUE)

## load data
data(bladderdata)
dat_mat <- exprs(bladderEset)
condition <- bladderEset$cancer
batch <- bladderEset$batch

## ComBat batch adjustment (mean and variance)
dat_adj <- ComBat(dat_mat, batch=batch, mod=model.matrix(~condition))

## DE analysis on combat adjusted data with limma
design <- model.matrix(~ 0 + condition)
colnames(design) <- c("Biopsy", "Cancer", "Normal")
contrast.matrix <- makeContrasts(Cancer-Normal, levels=design)
fit <- lmFit(dat_adj, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, coef=1, adjust.method="BH", number=nrow(dat_mat))
DE_up_genes <- rownames(res[res$logFC > 1, ])
DE_down_genes <- rownames(res[res$logFC < (-1), ])
nonDE_genes <- rownames(res[abs(res$logFC) <= 1, ])  
# log2 fold change at least 1 -> at least 2 fold changes
# length(DE_up_genes)+length(DE_down_genes)+length(nonDE_genes) == nrow(dat_mat)

## Condition difference: use batch corrected data to estimate
# mean
mean_up_case <- mean(dat_adj[DE_up_genes, condition=="Cancer"])
mean_down_case <- mean(dat_adj[DE_down_genes, condition=="Cancer"])
mean_nonDE_case <- mean(dat_adj[nonDE_genes, condition=="Cancer"])
mean_up_ctrl <- mean(dat_adj[DE_up_genes, condition=="Normal"])
mean_down_ctrl <- mean(dat_adj[DE_down_genes, condition=="Normal"])
mean_nonDE_ctrl <- mean(dat_adj[nonDE_genes, condition=="Normal"])
# var
var_up_case <- mean(apply(dat_adj[DE_up_genes, condition=="Cancer"], 1, var))
var_down_case <- mean(apply(dat_adj[DE_down_genes, condition=="Cancer"], 1, var))
var_nonDE_case <- mean(apply(dat_adj[nonDE_genes, condition=="Cancer"], 1, var))
var_up_ctrl <- mean(apply(dat_adj[DE_up_genes, condition=="Normal"], 1, var))
var_down_ctrl <- mean(apply(dat_adj[DE_down_genes, condition=="Normal"], 1, var))
var_nonDE_ctrl <- mean(apply(dat_adj[nonDE_genes, condition=="Normal"], 1, var))

## Batch difference:
# Ctrl, Batch 2 vs Batch 3
#mean
batch2_mean_up_ctrl <- mean(dat_mat[DE_up_genes, batch==2&condition=="Normal"])
batch3_mean_up_ctrl <- mean(dat_mat[DE_up_genes, batch==3&condition=="Normal"])
batch2_mean_down_ctrl <- mean(dat_mat[DE_down_genes, batch==2&condition=="Normal"])
batch3_mean_down_ctrl <- mean(dat_mat[DE_down_genes, batch==3&condition=="Normal"])
batch2_mean_nonDE_ctrl <- mean(dat_mat[nonDE_genes, batch==2&condition=="Normal"])
batch3_mean_nonDE_ctrl <- mean(dat_mat[nonDE_genes, batch==3&condition=="Normal"])
#var
mean(apply(dat_mat[, batch==2&condition=="Normal"], 1, var))
mean(apply(dat_mat[, batch==3&condition=="Normal"], 1, var))
# Case, Batch 1 vs 2 vs 5
#mean
batch1_mean_up_case <- mean(dat_mat[DE_up_genes, batch==1&condition=="Cancer"])
batch5_mean_up_case <- mean(dat_mat[DE_up_genes, batch==5&condition=="Cancer"])
batch1_mean_down_case <- mean(dat_mat[DE_down_genes, batch==1&condition=="Cancer"])
batch5_mean_down_case <- mean(dat_mat[DE_down_genes, batch==5&condition=="Cancer"])
batch1_mean_nonDE_case <- mean(dat_mat[nonDE_genes, batch==1&condition=="Cancer"])
batch5_mean_nonDE_case <- mean(dat_mat[nonDE_genes, batch==5&condition=="Cancer"])
#var
mean(apply(dat_mat[, batch==1&condition=="Cancer"], 1, var))
mean(apply(dat_mat[, batch==2&condition=="Cancer"], 1, var))
mean(apply(dat_mat[, batch==5&condition=="Cancer"], 1, var))



###  Oncogenic signature
rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/refcombat_review_201803/")

## load data
data_name <- "Sig"
filePaths <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/signatures/"
batch1 <- read.table(file=paste(filePaths, "GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog", sep=""))
batch2 <- read.table(file=paste(filePaths, "GFP30_KRAS-GV_KRAS-QH_KRAS-WT_tpmlog.txt", sep=""))
batch3 <- read.table(file=paste(filePaths, "18_GFP_EGFR_TPMlog2.txt", sep=""))
## intersect genes
dim(batch1); dim(batch2); dim(batch3)
setdiff(rownames(batch2), rownames(batch1)) # Extra gene: DDX11L1
batch2 <- batch2[rownames(batch1), ]
identical(rownames(batch1), rownames(batch2)); identical(rownames(batch1), rownames(batch3))
## Remove ERK.1-6 in batch 1 according to paper (P11)
batch1 <- batch1[, -grep("ERK", colnames(batch1))]
## Combine data
dat_mat <- cbind(batch1, batch2, batch3)
batch <- c(rep(1, ncol(batch1)), rep(2, ncol(batch2)), rep(3, ncol(batch3)))
cond1 <- rep(0, ncol(batch1)); cond2 <- rep(0, ncol(batch2)); cond3 <- rep(0, ncol(batch3))
# condition in batch 1
cond1[grep("GFP", colnames(batch1))] <- "GFP"
cond1[grep("AKT", colnames(batch1))] <- "AKT"
cond1[grep("BAD", colnames(batch1))] <- "BAD"
cond1[grep("IGF1R", colnames(batch1))] <- "IGF1R"
cond1[grep("RAF", colnames(batch1))] <- "RAF"
cond1[grep("HER2", colnames(batch1))] <- "HER2"
# condition in batch 2
cond2[grep("GFP", colnames(batch2))] <- "GFP"
cond2[grep("KRAS_WT", colnames(batch2))] <- "KRAS_WT"
cond2[grep("KRAS_GV", colnames(batch2))] <- "KRAS_GV"
cond2[grep("KRAS_QH", colnames(batch2))] <- "KRAS_QH"
# condition in batch 3
cond3[grep("Control", colnames(batch3))] <- "GFP"
cond3[grep("EGFR", colnames(batch3))] <- "EGFR"
condition <- as.factor(c(cond1, cond2, cond3))
# remove genes with 0 variance in any of batches
keep1 <- apply(dat_mat[, batch==1], 1, sd)!=0
keep2 <- apply(dat_mat[, batch==2], 1, sd)!=0
keep3 <- apply(dat_mat[, batch==3], 1, sd)!=0
dat_mat <- dat_mat[keep1 & keep2 & keep3, ]
dat_mat <- as.matrix(dat_mat)

## load signatures
sig_dir <- "/restricted/projectnb/pathsig/work/dfj/20160215_combat_meanonlypaper/results/TCGA/original/egfr_50_gene_list/adap_adap_single/posterior_delta.csv"
EGFR_assignres <- read.csv(sig_dir, as.is=TRUE)
EGFR_siglst <- EGFR_assignres[,1]  ## DE genes
EGFR_siglst_up <- EGFR_siglst[EGFR_assignres[,4]>=0]
EGFR_siglst_down <- EGFR_siglst[EGFR_assignres[,4]<0]


## Batch differences: compare between GFP samples
# mean
mean1_up <- mean(dat_mat[EGFR_siglst_up, batch==1&condition=="GFP"])
mean2_up <- mean(dat_mat[EGFR_siglst_up, batch==2&condition=="GFP"])
mean3_up <- mean(dat_mat[EGFR_siglst_up, batch==3&condition=="GFP"])
mean1_down <- mean(dat_mat[EGFR_siglst_down, batch==1&condition=="GFP"])
mean2_down <- mean(dat_mat[EGFR_siglst_down, batch==2&condition=="GFP"])
mean3_down <- mean(dat_mat[EGFR_siglst_down, batch==3&condition=="GFP"])
mean1_nonsig <- mean(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==1&condition=="GFP"])
mean2_nonsig <- mean(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==2&condition=="GFP"])
mean3_nonsig <- mean(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==3&condition=="GFP"])
# var
var1_up <- mean(apply(dat_mat[EGFR_siglst_up, batch==1&condition=="GFP"], 1, var))
var2_up <- mean(apply(dat_mat[EGFR_siglst_up, batch==2&condition=="GFP"], 1, var))
var3_up <- mean(apply(dat_mat[EGFR_siglst_up, batch==3&condition=="GFP"], 1, var))
var1_down <- mean(apply(dat_mat[EGFR_siglst_down, batch==1&condition=="GFP"], 1, var))
var2_down <- mean(apply(dat_mat[EGFR_siglst_down, batch==2&condition=="GFP"], 1, var))
var3_down <- mean(apply(dat_mat[EGFR_siglst_down, batch==3&condition=="GFP"], 1, var))
var1_nonsig <- mean(apply(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==1&condition=="GFP"], 1, var))
var2_nonsig <- mean(apply(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==2&condition=="GFP"], 1, var))
var3_nonsig <- mean(apply(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==3&condition=="GFP"], 1, var))

sigBatch_df <- data.frame(up=c(mean1=mean1_up, mean2=mean2_up, mean3=mean3_up,
                               var1=var1_up, var2=var2_up, var3=var3_up),
                          down=c(mean1=mean1_down, mean2=mean2_down, mean3=mean3_down,
                                 var1=var1_down, var2=var2_down, var3=var3_down),
                          nonsig=c(mean1=mean1_nonsig, mean2=mean2_nonsig, mean3=mean3_nonsig,
                                   var1=var1_nonsig, var2=var2_nonsig, var3=var3_nonsig)) 
write.csv(round(sigBatch_df, 2), file="meanvar_sigBatch.csv", quote=F)
save(mean1_up, mean2_up, mean3_up, mean1_down, mean2_down, mean3_down, mean1_nonsig, mean2_nonsig, mean3_nonsig,
     var1_up, var2_up, var3_up, var1_down, var2_down, var3_down, var1_nonsig, var2_nonsig, var3_nonsig,
     file="meanvar_sigBatch.rda")


## Condition differences: compare within EGFR (batch 3)
#mean
mean_ctrl_up <- mean(dat_mat[EGFR_siglst_up, batch==3&condition=="GFP"])
mean_case_up <- mean(dat_mat[EGFR_siglst_up, batch==3&condition=="EGFR"])
mean_ctrl_down <- mean(dat_mat[EGFR_siglst_down, batch==3&condition=="GFP"])
mean_case_down <- mean(dat_mat[EGFR_siglst_down, batch==3&condition=="EGFR"])
mean_ctrl_nonsig <- mean(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==3&condition=="GFP"])
mean_case_nonsig <- mean(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==3&condition=="EGFR"])
#var
var_ctrl_up <- mean(apply(dat_mat[EGFR_siglst_up, batch==3&condition=="GFP"],1,var))
var_case_up <- mean(apply(dat_mat[EGFR_siglst_up, batch==3&condition=="EGFR"],1,var))
var_ctrl_down <- mean(apply(dat_mat[EGFR_siglst_down, batch==3&condition=="GFP"],1,var))
var_case_down <- mean(apply(dat_mat[EGFR_siglst_down, batch==3&condition=="EGFR"],1,var))
var_ctrl_nonsig <- mean(apply(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==3&condition=="GFP"],1,var))
var_case_nonsig <- mean(apply(dat_mat[setdiff(rownames(dat_mat), EGFR_siglst), batch==3&condition=="EGFR"],1,var))

sigBio_df <- data.frame(up=c(mean_gfp=mean_ctrl_up, mean_egfr=mean_case_up, 
                             var_gfp=var_ctrl_up, var_egfr=var_case_up),
                        down=c(mean_gfp=mean_ctrl_down, mean_egfr=mean_case_down, 
                             var_gfp=var_ctrl_down, var_egfr=var_case_down),
                        nonsig=c(mean_gfp=mean_ctrl_nonsig, mean_egfr=mean_case_nonsig, 
                               var_gfp=var_ctrl_nonsig, var_egfr=var_case_nonsig))
write.csv(round(sigBio_df, 2), file="meanvar_sigBio.csv", quote=F)
save(mean_ctrl_up, mean_case_up, mean_ctrl_down, mean_case_down, mean_ctrl_nonsig, mean_case_nonsig,
     var_ctrl_up, var_case_up, var_ctrl_down, var_case_down, var_ctrl_nonsig, var_case_nonsig,  
     file="meanvar_sigBio.rda")



###  Nitric Oxide
rm(list=ls())
setwd("~/Dropbox/Work/EvanJohnson/refComBat_forSubmittion_BMCbioinfo/review20180320/modsel_sensitivity")
lapply(c("limma", "sva", "BatchQC"), require, character.only = TRUE)

## load data
fdir <- "~/Dropbox/Work/EvanJohnson/refComBat_paper_092116/BatchQC/"
dat_mat <- as.matrix(read.delim(paste0(fdir, "combat_paper_nitric_oxide_dataset/arielGeneric.txt")))
metdat <- read.delim(paste0(fdir, "combat_paper_nitric_oxide_dataset/arielGenericSampleInfo.txt"))
batch <- metdat[, "Batch"]
condition <- metdat[, "Treatment"]

## Preprocessing
#remove 0 variance genes in either batch
# keep1 <- apply(dat_mat[, batch==1], 1, var) != 0 
# keep2 <- apply(dat_mat[, batch==2], 1, var) != 0 
# dat_mat <- dat_mat[keep1 & keep2, ]
#log transform
# log_dat <- log2(dat_mat+1)
#normalize data
dat_mat <- gnormalize(dat_mat)

## ComBat adjustment
dat_adj <- ComBat(dat_mat, batch=batch, mod=model.matrix(~condition))

## DE analysis on ComBat adjusted data using limma
design <- model.matrix(~ 0 + condition)
colnames(design) <- c("Control", "NO")
contrast.matrix <- makeContrasts(NO-Control, levels=design)
fit <- lmFit(dat_adj, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, coef=1, adjust.method="BH", number=nrow(dat_mat))
# at least 2 fold change
DE_up_genes <- as.numeric(rownames(res[res$logFC > 1, ]))
DE_down_genes <- as.numeric(rownames(res[res$logFC < (-1), ]))
nonDE_genes <- as.numeric(rownames(res[abs(res$logFC) <= 1, ]))  


## Condition difference: estimate from ComBat adjusted data
mean_ctrl_up <- mean(dat_adj[DE_up_genes, condition=="Control"])
mean_case_up <- mean(dat_adj[DE_up_genes, condition=="NO"])
mean_ctrl_down <- mean(dat_adj[DE_down_genes, condition=="Control"])
mean_case_down <- mean(dat_adj[DE_down_genes, condition=="NO"])
mean_ctrl_nonsig <- mean(dat_adj[nonDE_genes, condition=="Control"])
mean_case_nonsig <- mean(dat_adj[nonDE_genes, condition=="NO"])

var_ctrl_up <- mean(apply(dat_adj[DE_up_genes, condition=="Control"],1,var))
var_case_up <- mean(apply(dat_adj[DE_up_genes, condition=="NO"],1,var))
var_ctrl_down <- mean(apply(dat_adj[DE_down_genes, condition=="Control"],1,var))
var_case_down <- mean(apply(dat_adj[DE_down_genes, condition=="NO"],1,var))
var_ctrl_nonsig <- mean(apply(dat_adj[nonDE_genes, condition=="Control"],1,var))
var_case_nonsig <- mean(apply(dat_adj[nonDE_genes, condition=="NO"],1,var))

NO_bio_df <- data.frame(up=c(mean_ctrl=mean_ctrl_up, mean_case=mean_case_up, 
                             var_ctrl=var_ctrl_up, var_case=var_case_up),
                        down=c(mean_ctrl=mean_ctrl_down, mean_case=mean_case_down, 
                               var_ctrl=var_ctrl_down, var_case=var_case_down),
                        nonsig=c(mean_ctrl=mean_ctrl_nonsig, mean_case=mean_case_nonsig, 
                                 var_ctrl=var_ctrl_nonsig, var_case=var_case_nonsig))
write.csv(round(NO_bio_df, 2), file="meanvar_NO_bio.csv", quote=F)
save(mean_ctrl_up, mean_case_up, mean_ctrl_down, mean_case_down, mean_ctrl_nonsig, mean_case_nonsig,
     var_ctrl_up, var_case_up, var_ctrl_down, var_case_down, var_ctrl_nonsig, var_case_nonsig,  
     file="meanvar_NO_bio.rda")


## Batch difference: estimated from data subtracting biological 
## component (function from BatchQC)
batch_mat <- batchQC_condition_adjusted(dat_mat, batch=batch, condition=condition)

# mean
mean1_up <- mean(batch_mat[DE_up_genes, batch==1])
mean2_up <- mean(batch_mat[DE_up_genes, batch==2])
mean3_up <- mean(batch_mat[DE_up_genes, batch==3])
mean1_down <- mean(batch_mat[DE_down_genes, batch==1])
mean2_down <- mean(batch_mat[DE_down_genes, batch==2])
mean3_down <- mean(batch_mat[DE_down_genes, batch==3])
mean1_nonsig <- mean(batch_mat[nonDE_genes, batch==1])
mean2_nonsig <- mean(batch_mat[nonDE_genes, batch==2])
mean3_nonsig <- mean(batch_mat[nonDE_genes, batch==3])
# var
var1_up <- mean(apply(batch_mat[DE_up_genes, batch==1], 1, var))
var2_up <- mean(apply(batch_mat[DE_up_genes, batch==2], 1, var))
var3_up <- mean(apply(batch_mat[DE_up_genes, batch==3], 1, var))
var1_down <- mean(apply(batch_mat[DE_down_genes, batch==1], 1, var))
var2_down <- mean(apply(batch_mat[DE_down_genes, batch==2], 1, var))
var3_down <- mean(apply(batch_mat[DE_down_genes, batch==3], 1, var))
var1_nonsig <- mean(apply(batch_mat[nonDE_genes, batch==1], 1, var))
var2_nonsig <- mean(apply(batch_mat[nonDE_genes, batch==2], 1, var))
var3_nonsig <- mean(apply(batch_mat[nonDE_genes, batch==3], 1, var))

NO_batch_df <- data.frame(up=c(mean1=mean1_up, mean2=mean2_up, mean3=mean3_up,
                               var1=var1_up, var2=var2_up, var3=var3_up),
                          down=c(mean1=mean1_down, mean2=mean2_down, mean3=mean3_down,
                                 var1=var1_down, var2=var2_down, var3=var3_down),
                          nonsig=c(mean1=mean1_nonsig, mean2=mean2_nonsig, mean3=mean3_nonsig,
                                   var1=var1_nonsig, var2=var2_nonsig, var3=var3_nonsig)) 
write.csv(round(NO_batch_df, 2), file="meanvar_NO_batch.csv", quote=F)
save(mean1_up, mean2_up, mean3_up, mean1_down, mean2_down, mean3_down, mean1_nonsig, mean2_nonsig, mean3_nonsig,
     var1_up, var2_up, var3_up, var1_down, var2_down, var3_down, var1_nonsig, var2_nonsig, var3_nonsig,
     file="meanvar_NO_batch.rda")
