############## Simulations on reference batch ComBat ##############
rm(list=ls())
setwd("C:/Users/zhang/Dropbox/Work/EvanJohnson/refComBat_forSubmission_plosone/")
set.seed(1)
x <- c("ggplot2", "reshape2", "plyr", "moments", "sva", "gridExtra", "Biobase", "RColorBrewer")
lapply(x, require, character.only = TRUE)


########  Simulate datasets  ########
## batch 1
batch1 <- matrix(rep(0, 6*200), nrow=200, ncol=6)
for(i in 1:nrow(batch1)){
  if(i<=100){
    # signature genes
    batch1[i,1:3] <- rnorm(3, mean=0, sd=sqrt(0.1)) # control
    batch1[i,4:6] <- rnorm(3, mean=1, sd=sqrt(0.1)) # case
  }else{
    # NULL genes
    batch1[i,1:6] <- rnorm(6, mean=0, sd=sqrt(0.1))
  }
}
## batch 2
batch2 <- matrix(rep(0,600*200), nrow=200, ncol=600)
for(i in 1:nrow(batch2)){
  if(i<=100){
    batch2[i,1:100] <- rnorm(100, mean=0.5, sd=sqrt(10))
    batch2[i,101:200] <- rnorm(100, mean=0.7, sd=sqrt(10))
    batch2[i,201:300] <- rnorm(100, mean=0.9, sd=sqrt(10))
    batch2[i,301:400] <- rnorm(100, mean=1.1, sd=sqrt(10))
    batch2[i,401:500] <- rnorm(100, mean=1.3, sd=sqrt(10))
    batch2[i,501:600] <- rnorm(100, mean=1.5, sd=sqrt(10))
  }else{
    batch2[i,1:600] <- rnorm(600, mean=0.5, sd=sqrt(10))
  }
}
## Construct ExpressionSet with combined data
dat <- cbind(batch1, batch2)
colnames(dat) <- 1:ncol(dat); rownames(dat) <- 1:nrow(dat)
batch <- c(rep(1, ncol(batch1)), rep(2, ncol(batch2)))
treatment <- c(rep(c(0, 1), each=3),
               rep(seq(from=0, to=1, by=0.2), each=100))
treatment <- as.data.frame(treatment)
colnames(treatment) <- "Treatment"; rownames(treatment) <- 1:nrow(treatment)
pheno <- new("AnnotatedDataFrame", data=treatment)
eset <- ExpressionSet(assayData=dat, phenoData=pheno)  



########  Batch Adjustment  ########
mod = model.matrix(~as.factor(Treatment), data=pData(eset))
### Combat - mod
res1 <- ComBat(exprs(eset), batch=batch, mod=mod, prior.plot=FALSE)
combat_mod_batch1 <- res1[, 1:ncol(batch1)]
combat_mod_batch2 <- res1[, (ncol(batch1)+1):ncol(dat)]
### Combat - NULL
res2 <- ComBat(exprs(eset), batch=batch, mod=NULL, prior.plot=FALSE)
combat_null_batch1 <- res2[, 1:ncol(batch1)]
combat_null_batch2 <- res2[, (ncol(batch1)+1):ncol(dat)]
### refCombat - mod
res_ref <- ComBat(exprs(eset), batch=batch, mod=mod, prior.plot=FALSE, ref.batch=1)
refcombat_mod_batch1 <- res_ref[, 1:ncol(batch1)]
refcombat_mod_batch2 <- res_ref[, (ncol(batch1)+1):ncol(dat)]
### refCombat - NULL
res_ref_2 <- ComBat(exprs(eset), batch=batch, mod=NULL, prior.plot=FALSE, ref.batch=1)
refcombat_null_batch1 <- res_ref_2[, 1:ncol(batch1)]
refcombat_null_batch2 <- res_ref_2[, (ncol(batch1)+1):ncol(dat)]



######## Fig4: heatmap for batch 1 vs batch 2, original vs combat vs refcombat (full mod) ########
breaks_seq <- c(-14,-5,-1,-0.7,-0.5,-0.3,0,0.2,0.5,0.8,1,1.2,1.5,1.7,2,5,15)
pal <- colorRampPalette(c("darkblue", "white", "orange"))(n =length(breaks_seq)-1)

datlst <- list(batch1, batch2,
               combat_mod_batch1, combat_mod_batch2,
               refcombat_mod_batch1, refcombat_mod_batch2)
names(datlst) <- c("Batch 1 - Original", "Batch2 - Original",
                   "Batch 1 - ComBat", "Batch2 - ComBat",
                   "Batch 1 - Ref ComBat", "Batch2 - Ref ComBat")

hmaps_lst <- lapply(1:length(datlst), function(ii){
  dat <- datlst[[ii]]
  rownames(dat) <- paste0("gene", 1:nrow(dat))
  colnames(dat) <- paste0("sample", 1:ncol(dat))
  
  dat_m <- melt(dat)
  colnames(dat_m) <- c("Genes", "Samples", "Expression")
  dat_m$Genes <- factor(dat_m$Genes, levels=rev(levels(dat_m$Genes)))
  
  gene_annot <- c(0, nrow(dat_m))
  gene_annot[dat_m$Genes %in% paste0("gene", 1:100)] <- "Signature genes"
  gene_annot[dat_m$Genes %in% paste0("gene", 101:200)] <- "Control genes"
  dat_m$gene_annot <- factor(gene_annot, levels=c("Signature genes", "Control genes"))
  
  sample_annot <- c(0, nrow(dat_m))
  if(ii %in% c(1,3,5)){
    sample_annot[dat_m$Sample %in% paste0("sample", 1:3)] <- "Pathway off"
    sample_annot[dat_m$Sample %in% paste0("sample", 4:6)] <- "Pathway on"
    dat_m$sample_annot <- factor(sample_annot, levels=c("Pathway off", "Pathway on"))
  }else{
    sample_annot <- "Pathway activity increasing (from left to right)"
    dat_m$sample_annot <- factor(sample_annot, levels=c("Pathway activity increasing (from left to right)"))
  }
  
  plt <- ggplot(dat_m, aes(x=Samples, y=Genes, fill=Expression)) +
    facet_grid(gene_annot ~ sample_annot, scales="free", switch="y") + 
    geom_tile() +
    scale_fill_gradientn(colors=pal, breaks=breaks_seq, limits=c(-5,5),
                         labels = c('',-5,'','','','',0,'','','','','','','','',5,'')) +
    labs(x=NULL, y=NULL, title=names(datlst[ii])) +
    theme(strip.text=element_text(size=4),
          panel.background = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.title=element_text(hjust=0, size=5),
          panel.spacing.x=unit(0.03, "line"),
          panel.spacing.y=unit(0.03, "line"),
          legend.title=element_blank(),
          legend.text=element_text(size=3),
          legend.position="bottom",
          legend.key.size=unit(0.2, "line"),
          legend.key.width=unit(0.8, "line"))
  return(plt)
})
fig4 <- grid.arrange(hmaps_lst[[1]], hmaps_lst[[2]], hmaps_lst[[3]], 
                     hmaps_lst[[4]], hmaps_lst[[5]], hmaps_lst[[6]], ncol=2)
ggsave("./figures/Fig4/Fig4.png", plot=fig4, width=4.5, height=6, units="in", dpi=550)
dev.off()



######## S4 Fig: heatmap for batch 1 vs batch 2, combat vs refcombat (null mod) ########
breaks_seq <- c(-14,-5,-1,-0.7,-0.5,-0.3,0,0.2,0.5,0.8,1,1.2,1.5,1.7,2,5,15)
pal <- colorRampPalette(c("darkblue", "white", "orange"))(n =length(breaks_seq)-1)

datlst_nullmod <- list(combat_null_batch1, combat_null_batch2,
                       refcombat_null_batch1, refcombat_null_batch2)
names(datlst_nullmod) <- c("Batch 1 - ComBat", "Batch2 - ComBat",
                           "Batch 1 - Ref ComBat", "Batch2 - Ref ComBat")

hmaps_lst_nullmod <- lapply(1:length(datlst_nullmod), function(ii){
  dat <- datlst_nullmod[[ii]]
  rownames(dat) <- paste0("gene", 1:nrow(dat))
  colnames(dat) <- paste0("sample", 1:ncol(dat))
  
  dat_m <- melt(dat)
  colnames(dat_m) <- c("Genes", "Samples", "Expression")
  dat_m$Genes <- factor(dat_m$Genes, levels=rev(levels(dat_m$Genes)))
  
  gene_annot <- c(0, nrow(dat_m))
  gene_annot[dat_m$Genes %in% paste0("gene", 1:100)] <- "Signature genes"
  gene_annot[dat_m$Genes %in% paste0("gene", 101:200)] <- "Control genes"
  dat_m$gene_annot <- factor(gene_annot, levels=c("Signature genes", "Control genes"))
  
  sample_annot <- c(0, nrow(dat_m))
  if(ii %in% c(1,3)){
    sample_annot[dat_m$Sample %in% paste0("sample", 1:3)] <- "Pathway off"
    sample_annot[dat_m$Sample %in% paste0("sample", 4:6)] <- "Pathway on"
    dat_m$sample_annot <- factor(sample_annot, levels=c("Pathway off", "Pathway on"))
  }else{
    sample_annot <- "Pathway activity increasing (from left to right)"
    dat_m$sample_annot <- factor(sample_annot, levels=c("Pathway activity increasing (from left to right)"))
  }
  
  plt <- ggplot(dat_m, aes(x=Samples, y=Genes, fill=Expression)) +
    facet_grid(gene_annot ~ sample_annot, scales="free", switch="y") + 
    geom_tile() +
    scale_fill_gradientn(colors=pal, breaks=breaks_seq, limits=c(-5,5),
                         labels = c('',-5,'','','','',0,'','','','','','','','',5,'')) +
    labs(x=NULL, y=NULL, title=names(datlst_nullmod[ii])) +
    theme(strip.text=element_text(size=4),
          panel.background = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.title=element_text(hjust=0, size=5),
          panel.spacing.x=unit(0.03, "line"),
          panel.spacing.y=unit(0.03, "line"),
          legend.title=element_blank(),
          legend.text=element_text(size=3),
          legend.position="bottom",
          legend.key.size=unit(0.2, "line"),
          legend.key.width=unit(0.8, "line"))
  return(plt)
})
S4_fig <- grid.arrange(hmaps_lst_nullmod[[1]], hmaps_lst_nullmod[[2]], 
                       hmaps_lst_nullmod[[3]], hmaps_lst_nullmod[[4]], ncol=2)
ggsave("./support_info/S4_fig/S4_fig.png", plot=S4_fig, width=4.5, height=4, units="in", dpi=550)
dev.off()



############## K-means biomarker detection ###############
origin_dat <- cbind(batch1, batch2)
combat_mod_dat <- cbind(combat_mod_batch1, combat_mod_batch2)
combat_null_dat <- cbind(combat_null_batch1, combat_null_batch2)
refcombat_mod_dat <- cbind(refcombat_mod_batch1, refcombat_mod_batch2)
refcombat_null_dat <- cbind(refcombat_null_batch1, refcombat_null_batch2)

### cluster on genes ###
set.seed(1)
origin_kmeans <- kmeans(origin_dat, 2, algorithm="Lloyd", iter.max=50)
combat_mod_kmeans <- kmeans(combat_mod_dat, 2, algorithm="Lloyd", iter.max=50)
combat_null_kmeans <- kmeans(combat_null_dat, 2, algorithm="Lloyd", iter.max=50)
refcombat_mod_kmeans <- kmeans(refcombat_mod_dat, 2, algorithm="Lloyd", iter.max=50)
refcombat_null_kmeans <- kmeans(refcombat_null_dat, 2, algorithm="Lloyd", iter.max=50)
 
#### only cluster on batch1 ####
batch1_origin_kmeans <- kmeans(batch1, 2, algorithm="Lloyd",iter.max=50)
batch1_combat_mod_kmeans <- kmeans(combat_mod_batch1, 2, algorithm="Lloyd",iter.max=50)
batch1_combat_null_kmeans <- kmeans(combat_null_batch1, 2, algorithm="Lloyd",iter.max=50)
batch1_refcombat_mod_kmeans <- kmeans(refcombat_mod_batch1, 2, algorithm="Lloyd",iter.max=50)
batch1_refcombat_null_kmeans <- kmeans(refcombat_null_batch1, 2, algorithm="Lloyd",iter.max=50)
#### only cluster on batch2 ####
batch2_origin_kmeans <- kmeans(batch2, 2, algorithm="Lloyd",iter.max=50)
batch2_combat_mod_kmeans <- kmeans(combat_mod_batch2, 2, algorithm="Lloyd",iter.max=50)
batch2_combat_null_kmeans <- kmeans(combat_null_batch2, 2, algorithm="Lloyd",iter.max=50)
batch2_refcombat_mod_kmeans <- kmeans(refcombat_mod_batch2, 2, algorithm="Lloyd",iter.max=50)
batch2_refcombat_null_kmeans <- kmeans(refcombat_null_batch2, 2, algorithm="Lloyd",iter.max=50)

#### Cluster accuracy ####
source("scripts/clustAccu.R")
# batch 1+2
accu_origin_comb <- clustAccu(origin_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_combat_mod_comb <- clustAccu(combat_mod_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_combat_null_comb <- clustAccu(combat_null_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_refcombat_mod_comb <- clustAccu(refcombat_mod_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_refcombat_null_comb <- clustAccu(refcombat_null_kmeans$cluster, n_gene=200, n_biomarker=100)
# batch 1
accu_origin_1 <- clustAccu(batch1_origin_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_combat_mod_1 <- clustAccu(batch1_combat_mod_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_combat_null_1 <- clustAccu(batch1_combat_null_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_refcombat_mod_1 <- clustAccu(batch1_refcombat_mod_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_refcombat_null_1 <- clustAccu(batch1_refcombat_null_kmeans$cluster, n_gene=200, n_biomarker=100)
# batch 2
accu_origin_2 <- clustAccu(batch2_origin_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_combat_mod_2 <- clustAccu(batch2_combat_mod_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_combat_null_2 <- clustAccu(batch2_combat_null_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_refcombat_mod_2 <- clustAccu(batch2_refcombat_mod_kmeans$cluster, n_gene=200, n_biomarker=100)
accu_refcombat_null_2 <- clustAccu(batch2_refcombat_null_kmeans$cluster, n_gene=200, n_biomarker=100)


########  Fig 5: kmeans biomarker separation (original + full mod)  ########
#png("./figures/Fig5/Fig5_nocaps.png", width=4, height=3, units="in", res=550)
setEPS()
postscript("./figures/Fig5/Fig5_nocaps.eps", width=4, height=3)
layout(matrix(c(0, 0, 1:4, 0, 0, 0, 5:11), nrow=2, byrow=TRUE))
par(xpd=T, mar=c(1,1,2,1))
## part 1: original
# True separation
image(t(as.matrix(c(rep(1,100), rep(2,100)))),
      col=c("blue","red"), axes=FALSE)
text(x=0, y=1.15, labels=c("True separation"), cex=0.6)
# Batch 1 original
image(t(as.matrix(batch1_origin_kmeans$cluster)),
      col=c("blue","red"), axes=FALSE)
text(x=0, y=1.15, labels=c("Batch 1"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_origin_1*100, "%", sep="")), cex=0.5)
# Batch 2 original
image(t(as.matrix(batch2_origin_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.15, labels=c("Batch 2"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_origin_2*100, "%", sep="")), cex=0.5)
# combined original
image(t(as.matrix(origin_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.15, labels=c("Whole dataset"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_origin_comb*100, "%", sep="")), cex=0.5)

## part 2: full models
image(t(as.matrix(c(rep(1,100), rep(2,100)))),
      col=c("blue","red"), axes=FALSE)
text(x=0, y=1.15, labels=c("True separation"), cex=0.6)
image(t(as.matrix(batch1_combat_mod_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.15, labels=c("Batch 1\nComBat"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_combat_mod_1*100, "%", sep="")), cex=0.5)
image(t(as.matrix(batch2_combat_mod_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.15, labels=c("Batch 2\nComBat"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_combat_mod_2*100, "%", sep="")), cex=0.5)
image(t(as.matrix(combat_mod_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.15, labels=c("Whole dataset\nComBat"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_combat_mod_comb*100, "%", sep="")), cex=0.5)
image(t(as.matrix(batch1_refcombat_mod_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.15, labels=c("Batch 1\nrefComBat"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_refcombat_mod_1*100, "%", sep="")), cex=0.5)
image(t(as.matrix(batch2_refcombat_mod_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.15, labels=c("Batch 2\nrefComBat"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_refcombat_mod_2*100, "%", sep="")), cex=0.5)
image(t(as.matrix(refcombat_mod_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.15, labels=c("Whole dataset\nrefComBat"), cex=0.6)
text(x=0, y=1.05, labels=c(paste(accu_refcombat_mod_comb*100, "%", sep="")), cex=0.5)
dev.off()



########  S5 Fig: kmeans biomarker separation (null mod)  ########
par(mfrow=c(1,7))
par(xpd=T, mar=c(2,2,12,2))

image(t(as.matrix(c(rep(1,100), rep(2,100)))),
      col=c("blue","red"), axes=FALSE)
text(x=0, y=1.2, labels=c("True separation"), cex=2.5)

image(t(as.matrix(batch1_combat_null_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.2, labels=c("Batch 1 \n ComBat"), cex=2.5)
text(x=0, y=1.1, labels=c(paste(accu_combat_null_1*100, "%", sep="")), cex=2)

image(t(as.matrix(batch2_combat_null_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.2, labels=c("Batch 2 \n ComBat"), cex=2.5)
text(x=0, y=1.1, labels=c(paste(accu_combat_null_2*100, "%", sep="")), cex=2)

image(t(as.matrix(combat_null_kmeans$cluster)),
      col=c("blue","red"), axes=FALSE)
text(x=0, y=1.2, labels=c("Whole dataset \n ComBat"), cex=2.5)
text(x=0, y=1.1, labels=c(paste(accu_combat_null_comb*100, "%", sep="")), cex=2)

image(t(as.matrix(batch1_refcombat_null_kmeans$cluster)),
      col=c("blue","red"), axes=FALSE)
text(x=0, y=1.2, labels=c("Batch 1 \n refComBat"), cex=2.5)
text(x=0, y=1.1, labels=c(paste(accu_refcombat_null_1*100, "%", sep="")), cex=2)

image(t(as.matrix(batch2_refcombat_null_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.2, labels=c("Batch 2 \n refComBat"), cex=2.5)
text(x=0, y=1.1, labels=c(paste(accu_refcombat_null_2*100, "%", sep="")), cex=2)

image(t(as.matrix(refcombat_null_kmeans$cluster)),
      col=c("red","blue"), axes=FALSE)
text(x=0, y=1.2, labels=c("Whole dataset \n refComBat"), cex=2.5)
text(x=0, y=1.1, labels=c(paste(accu_refcombat_null_comb*100, "%", sep="")), cex=2)

par(mar=c(5, 4, 4, 2) + 0.1)
