library(devtools)
library(sva)
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
working_dir         <- "/restricted/projectnb/combat/work/yuqingz/refcombat_review_201803/compare_sva_ruv/ICBP/SVA"
output_rda          <- "icbp_SVA.rda"

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
#SVA#
#------#
bat<-as.data.frame(cbind(c(rep(1,12),rep(2,41),rep(3,36)),c(rep(1,6),rep(2,6),rep(1,12),rep(3,6),rep(4,6),rep(5,5),rep(6,6),rep(7,6),rep(1,9),rep(8,9),rep(9,9),rep(10,9))))
colnames(bat)<-c("Batch","Model")
rownames(bat)<-colnames(gfp_egfr_kras_multi_f)
mod <- model.matrix(~as.factor(bat$Model))
mod0 <- model.matrix(~1, data=bat)
trainSv <- sva(as.matrix(gfp_egfr_kras_multi_f), mod, mod0)
test <- test[rownames(gfp_egfr_kras_multi_f), ]
fsvaobj <- fsva(as.matrix(gfp_egfr_kras_multi_f), mod, trainSv, as.matrix(test))
sva_expr1 <- merge_drop(fsvaobj$db, fsvaobj$new)

c_gfp<-subset(sva_expr1, select=GFP.1:GFP.12)
c_akt<-subset(sva_expr1, select=AKT.1:AKT.6)
c_bad<-subset(sva_expr1, select=BAD.1:BAD.6)
c_her2<-subset(sva_expr1, select=HER2.1:HER2.6)
c_igf1r<-subset(sva_expr1, select=IGF1R.1:IGF1R.6)
c_raf<-subset(sva_expr1, select=RAF.1:RAF.6)
train_egfr<-sva_expr1[,1:12]
c_egfr_gfp <- train_egfr[,1:6]
c_egfr <- train_egfr[,7:12]
c_kras_gfp<-subset(sva_expr1,select=GFP30.1:GFP30.9)
c_kraswt<-subset(sva_expr1,select=KRAS_WT.1:KRAS_WT.9)
c_krasqh<-subset(sva_expr1,select=KRAS_QH.1:KRAS_QH.9)
c_krasgv<-subset(sva_expr1,select=KRAS_GV.1:KRAS_GV.9)
c_test<-sva_expr1[,(ncol(gfp_egfr_kras_multi_f)+1):ncol(sva_expr1)]

save.image(file=output_rda)
