library(devtools)
library(ruv)
library(RUVnormalize)
library(limma)
library(data.table)

# Name:    ASSIGN_merge_and_combat.R
#
# Purpose: Merge together the signature data, and a test dataset and perform
#          the mean only version of ComBat and save a session that can be used
#          with ASSIGN_run_predictions.R to run ASSIGN
#
# Usage:   Rscript ASSIGN_merge_and_combat.R
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-09-29
#
################################################################################

#----------------------------------------------------#
#Input Files (modify these locations for your system)#
#----------------------------------------------------#
signatures_dir      <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/signatures"
expr_file           <- paste(signatures_dir,"GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog",sep="/")
control_egfr_l_file <- paste(signatures_dir,"18_GFP_EGFR_TPMlog2.txt",sep="/")
gfp_kras_file       <- paste(signatures_dir,"GFP30_KRAS-GV_KRAS-QH_KRAS-WT_tpmlog.txt",sep="/")
key_assign_file     <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/scripts/Key_ASSIGN_functions_balancedsig_v2.R"
testFile            <- "/restricted/projectnb/pathsig/work/dfj/20150929_bild_paper_new_ASSIGN/data/test_data/icbp_Rsubread_tpmlog.txt"

#--------------------------------------#
#Output Files (modify these every time)#
#--------------------------------------#
working_dir         <- "/restricted/projectnb/combat/work/yuqingz/refcombat_review_201803/compare_sva_ruv/ICBP/RUV"
output_rda          <- "icbp_RUV.rda"

#---------#
#Load Data#
#---------#
source(key_assign_file)
setwd(working_dir)
expr<-as.matrix(read.table(expr_file,sep='\t',row.names=1,header=1))
control<-subset(expr, select=GFP.1:GFP.12)
her2<-subset(expr, select=HER2.1:HER2.6)
akt<-subset(expr,select=AKT.1:AKT.6)
bad<-subset(expr,select=BAD.1:BAD.6)
igf1r<-subset(expr,select=IGF1R.1:IGF1R.6)
raf<-subset(expr,select=RAF.1:RAF.6)
expr_all<-cbind(control,akt,bad,her2,igf1r,raf)
expr_all_f <-expr_all[apply(expr_all[,1:41]==0,1,mean) < 0.85,]
control_egfr_l<-read.table(control_egfr_l_file, sep='\t', header=1, row.names=1)
gfp_egfr_multi_f <- merge_drop(control_egfr_l,expr_all_f)
gfp_kras<-read.table(gfp_kras_file, sep='\t', header=1, row.names=1)
gfp_egfr_kras_multi_f<-merge_drop(gfp_egfr_multi_f,gfp_kras)
#load in test data frame
test<-data.frame(fread(testFile), check.names=F,row.names=1)

#------#
#RUV#
#------#
bat<-as.data.frame(cbind(c(rep(1,12),rep(2,41),rep(3,36)),
                         c(rep("GFP",6),rep("EGFR",6),
                           rep("GFP",12),rep("AKT",6),rep("BAD",6),rep("HER2",5),rep("IGF1R",6),rep("RAF",6),
                           rep("GFP",9),rep("KRAS_GV",9),rep("KRAS_QH",9),rep("KRAS_WT",9))))
colnames(bat)<-c("Batch","Model")
rownames(bat)<-colnames(gfp_egfr_kras_multi_f)
path_names <- c("GFP", "EGFR", "AKT", "BAD", "HER2", "IGF1R", "RAF", "KRAS_GV", "KRAS_QH", "KRAS_WT")

## Search for control genes: least DE between any pathway on and controls in signature data
mod <- model.matrix(~ 0 + factor(bat$Model, levels=path_names), data=bat)
colnames(mod) <- path_names
leastDE_genes <- list()
for(ii in 2:10){
  contrast_matrix <- makeContrasts(paste0(path_names[ii],"-GFP"), levels=mod)
  fit <- lmFit(as.matrix(gfp_egfr_kras_multi_f), mod)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, coef=1, adjust.method="BH", number=nrow(gfp_egfr_kras_multi_f))
  leastDE_genes[[ii-1]] <- rownames(tail(res, n=3000))
}
ctrl_geneID <- Reduce(intersect, leastDE_genes)
ctl <- match(ctrl_geneID, rownames(gfp_egfr_kras_multi_f))

## ruv::RUVIII - Globally adjust data matrix using both negative controls and replicates.
Y <- t(gfp_egfr_kras_multi_f)
M <- replicate.matrix(factor(bat$Model, levels=path_names))
ruv_expr <- RUVIII(Y, M, ctl=ctl, k=length(ctl))

## RUVnormalize
test <- test[colnames(Y), ]
ruv_expr1 <- t(naiveRandRUV(rbind(Y, t(test)), ctl))

c_gfp<-subset(ruv_expr1, select=GFP.1:GFP.12)
c_akt<-subset(ruv_expr1, select=AKT.1:AKT.6)
c_bad<-subset(ruv_expr1, select=BAD.1:BAD.6)
c_her2<-subset(ruv_expr1, select=HER2.1:HER2.6)
c_igf1r<-subset(ruv_expr1, select=IGF1R.1:IGF1R.6)
c_raf<-subset(ruv_expr1, select=RAF.1:RAF.6)
train_egfr<-ruv_expr1[,1:12]
c_egfr_gfp <- train_egfr[,1:6]
c_egfr <- train_egfr[,7:12]
c_kras_gfp<-subset(ruv_expr1,select=GFP30.1:GFP30.9)
c_kraswt<-subset(ruv_expr1,select=KRAS_WT.1:KRAS_WT.9)
c_krasqh<-subset(ruv_expr1,select=KRAS_QH.1:KRAS_QH.9)
c_krasgv<-subset(ruv_expr1,select=KRAS_GV.1:KRAS_GV.9)
c_test<-ruv_expr1[,(ncol(gfp_egfr_kras_multi_f)+1):ncol(ruv_expr1)]

save.image(file=output_rda)
