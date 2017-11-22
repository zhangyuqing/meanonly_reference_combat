rm(list=ls())
setwd("C:/Users/zhang/Dropbox/Work/EvanJohnson/refComBat_paper_092116/BatchQC/")
library(BatchQC)


######## nitric oxide ########
dat <- read.delim("combat_paper_nitric_oxide_dataset/arielGeneric.txt")
metdat <- read.delim("combat_paper_nitric_oxide_dataset/arielGenericSampleInfo.txt")
batch <- metdat[, "Batch"]
condition <- metdat[, "Treatment"]
batchQC(dat, batch=batch, condition=condition,
        report_file="nitric_oxide_report.html", report_dir=".", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)



######## oncogenic signature (Real signature dataset / human mammary epithelial cells) ########
#### batchQC included data ####
# data(example_batchqc_data)
# batch <- batch_indicator$V1
# condition <- batch_indicator$V2
# batchQC(signature_data, batch=batch, condition=condition, 
#         report_file="oncogenic_signature_report.html", report_dir=".", 
#         report_option_binary="111111111",
#         view_report=FALSE, interactive=TRUE)

#### oncogenic signature data from David ####
filePaths <- "C:/Users/zhang/Documents/Work/Evan/oncogenic_signature_ABild/"
batch1 <- read.table(file=paste(filePaths, "GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog", sep=""))
batch2 <- read.table(file=paste(filePaths, "GFP30_KRAS-GV_KRAS-QH_KRAS-WT_tpmlog.txt", sep=""))
batch3 <- read.table(file=paste(filePaths, "18_GFP_EGFR_TPMlog2.txt", sep=""))
# intersect genes
dim(batch1)
dim(batch2)
setdiff(rownames(batch2), rownames(batch1)) # Extra gene: DDX11L1
dim(batch3)
identical(rownames(batch1), rownames(batch3))
matchID <- c()
for(i in 1:nrow(batch1)){
  matchID[i] <- which(rownames(batch2)==rownames(batch1)[i])
}
batch2 <- batch2[matchID, ]
identical(rownames(batch1), rownames(batch2))
# Remove ERK.1-6 in batch 1 according to paper (P11)
batch1 <- batch1[, -grep("ERK", colnames(batch1))]
# Combine data
edata <- cbind(batch1, batch2, batch3)
batch <- c(rep(1, ncol(batch1)), rep(2, ncol(batch2)), rep(3, ncol(batch3)))
cond1 <- rep(0, ncol(batch1)); cond2 <- rep(0, ncol(batch2)); cond3 <- rep(0, ncol(batch3))
# condition in batch 1
cond1[grep("GFP", colnames(batch1))] <- 0
cond1[grep("AKT", colnames(batch1))] <- 1
cond1[grep("BAD", colnames(batch1))] <- 2
cond1[grep("IGF1R", colnames(batch1))] <- 3
cond1[grep("RAF", colnames(batch1))] <- 4
cond1[grep("HER2", colnames(batch1))] <- 5
# condition in batch 2
cond2[grep("GFP", colnames(batch2))] <- 0
cond2[grep("KRAS_WT", colnames(batch2))] <- 6
cond2[grep("KRAS_GV", colnames(batch2))] <- 7
cond2[grep("KRAS_QH", colnames(batch2))] <- 8
# condition in batch 3
cond3[grep("Control", colnames(batch3))] <- 0
cond3[grep("EGFR", colnames(batch3))] <- 9
condition <- c(cond1, cond2, cond3);condition <- as.factor(condition)
rm(batch1, batch2, batch3, cond1, cond2, cond3, matchID, filePaths)

batchQC(edata, batch=batch, condition=condition,
        report_file="oncogenic_signature_report.html", report_dir=".",
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)



######## lung cancer ########
lnames <- load("stripiness/stripiness_data.rda")
mod = model.matrix(~ 1, data=ann)
# Set up batch variable (=batch)
batch.rma <- as.character(ann$batch[ match(colnames(edata.rma), rownames(ann))])
batch.scan <- as.character(ann$batch[ match(colnames(edata.scan), rownames(ann))])

batchQC(edata.rma, batch=batch.rma, condition=mod, 
        report_file="lung_cancer_rma_report.html", report_dir=".", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)    ## works on desktop

batchQC(edata.scan, batch=batch.scan, condition=mod, 
        report_file="lung_cancer_scan_report.html", report_dir=".", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)



######## bladder cancer ########
library(bladderbatch)
data(bladderdata)
pheno <- pData(bladderEset)
edata <- exprs(bladderEset)
batch <- pheno$batch  
condition <- pheno$cancer
tmp <- batchQC(edata, batch=batch, condition=condition, 
        report_file="baldder_cancer_report.html", report_dir=".", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)
