rm(list=ls()); dev.off()
setwd("C:/Users/zhang/Dropbox/Work/EvanJohnson/refComBat_forSubmission_plosone/")
x <- c("ggplot2", "reshape2", "plyr", "moments", "sva", "BatchQC", "gridExtra")
lapply(x, require, character.only = TRUE)
source("scripts/simulation-batchQC-helper.R")


############   Load Data   ############
## Bladder Cancer
library(bladderbatch)
data(bladderdata)
pheno <- pData(bladderEset)
edata <- exprs(bladderEset)
batch <- pheno$batch  
condition <- pheno$cancer
nbatch <- length(unique(batch))

## Nitric Oxide
filePaths <- "C:/Users/zhang/Dropbox/Work/EvanJohnson/refComBat_paper_092116/BatchQC/"
edata <- as.matrix(read.delim(paste0(filePaths, "combat_paper_nitric_oxide_dataset/arielGeneric.txt")))
metdat <- read.delim(paste0(filePaths, "combat_paper_nitric_oxide_dataset/arielGenericSampleInfo.txt"))
batch <- metdat[, "Batch"]
condition <- metdat[, "Treatment"]
nbatch <- length(unique(batch))

## Oncogenic signature data from David 
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


############   batchQC preprocess: filter genes   ############ 
tmp <- batchQC_filter_genes(as.matrix(edata), batch, condition)
identical(tmp, edata)
if(!identical(tmp, edata)){edata <- tmp}
rm(tmp)


############   Apply Mean-only and Regular ComBat   ############
edata.meanonly <- ComBat(edata, batch=batch, mod=model.matrix(~condition), mean.only=TRUE)
edata.combat <- ComBat(edata, batch=batch, mod=model.matrix(~condition))


############   Condition adjust: remove biological component   ############
edata.adj <- batchQC_condition_adjusted(edata, batch, condition)
edata.meanonly.adj <- batchQC_condition_adjusted(edata.meanonly, batch, condition)
edata.combat.adj <- batchQC_condition_adjusted(edata.combat, batch, condition)
rm(edata.meanonly, edata.combat)


############   shapeAnalysis preprocess: sort sample and gnormalize; get moments   ############
groupsorder <- order(batch)
sortbatch <- batch[groupsorder]
sortcond <- condition[groupsorder]

sortdata <- gnormalize(edata.adj[, groupsorder])
sortdata.meanonly <- gnormalize(edata.meanonly.adj[, groupsorder])
sortdata.combat <- gnormalize(edata.combat.adj[, groupsorder])
rm(edata.adj, edata.meanonly.adj, edata.combat.adj)

## Get Sample-wise moments
Y <- fitMoments(sortdata)
Y.meanonly <- fitMoments(sortdata.meanonly)
Y.combat <- fitMoments(sortdata.combat)
# check anova
sapply(colnames(Y), function(i){
  ff <- lm(Y.combat[,i] ~ as.factor(sortbatch))  #modify Y/Y.meanonly/Y.combat
  return(round(anova(ff)[1,5],4))
})

## Get Gene-wise moments
X <- genewise_stats(sortdata, sortbatch)
X.meanonly <- genewise_stats(sortdata.meanonly, sortbatch)
X.combat <- genewise_stats(sortdata.combat, sortbatch)
# check anova
res <- do.call(cbind, lapply(X.combat, c))    #modify X/X.meanonly/X.combat
prdctrs <- factor(rep(unique(sortbatch), each=nrow(sortdata)))
gg <- delta_f.pvalue(t(res), model.matrix(~prdctrs), matrix(rep(1, nrow(res)), ncol=1))
round(gg$p,4)

## sortdata (.meanonly, .combat are the final adjusted data to be plotted)



############  Fig1: bladder cancer, sample-wise mean and sample-wise variance boxplots  ############
allY <- data.frame(rbind(Y, Y.meanonly, Y.combat), 
                   combat_type=factor(c(rep("No\nAdjustment", nrow(Y)),
                                        rep("Mean-only\nComBat", nrow(Y.meanonly)),
                                        rep("Mean/variance\nComBat", nrow(Y.combat))),
                                      levels=c("No\nAdjustment", "Mean-only\nComBat", "Mean/variance\nComBat")))
pltdat <- data.frame(Sample=rep(colnames(edata)[groupsorder], 3), 
                     allY[, c("Mean", "Variance", "combat_type")], 
                     Batch=as.factor(rep(paste("Batch", sortbatch, sep=""), 3)))
# sample-wise mean moments
mean_plt <- ggplot(data=pltdat, aes(x=Batch, y=Mean, fill=Batch)) + 
  geom_boxplot(outlier.size=0.3) + 
  facet_grid(~combat_type) + 
  stat_summary(fun.y = "mean", colour="darkred", geom="point", 
               shape=18, size=0.8, show.legend = FALSE) +
  labs(y="Sample-wise Mean") +
  scale_y_continuous(breaks = seq(from=-0.25, to=0.3, 0.05)) +
  coord_cartesian(ylim=c(-0.25, 0.3)) +
  theme_bw() +
  theme(text = element_text(size=5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5),
        legend.key.size=unit(0.65,"line")) 
# sample-wise var moments
var_plt <- ggplot(data=pltdat, aes(x=Batch, y=Variance, fill=Batch)) + 
  geom_boxplot(outlier.size=0.3) + 
  facet_grid(~combat_type) + 
  stat_summary(fun.y = "mean", colour="darkred", geom="point", 
               shape=18, size=0.8, show.legend = FALSE) +
  labs(y="Sample-wise Variance") +
  scale_y_continuous(breaks = seq(from=0.25, to=2, 0.25)) +
  coord_cartesian(ylim=c(0.25, 2)) +
  theme_bw() +
  theme(text = element_text(size=5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5),
        legend.key.size=unit(0.65,"line")) 
# plot in one figure
fig1 <- grid.arrange(mean_plt, var_plt, nrow=2)
# save plot
ggsave("./figures/Fig1/Fig1_nocaps.eps", plot=fig1, width=3, height=5, units="in", dpi=550)
dev.off()



############  Fig2: bladder cancer, gene-wise variance boxplots  ############
genetest_varmoments <- rbind(melt(X[["Variance"]])[,-1], 
                             melt(X.meanonly[["Variance"]])[,-1], 
                             melt(X.combat[["Variance"]])[,-1])
genetest_combat_type <- factor(c(rep("No\nAdjustment", length(X[["Variance"]])),
                                 rep("Mean-only\nComBat", length(X.meanonly[["Variance"]])),
                                 rep("Mean/variance\nComBat", length(X.combat[["Variance"]]))),
                               levels=c("No\nAdjustment", "Mean-only\nComBat", "Mean/variance\nComBat"))
pltdat <- data.frame(genetest_varmoments, genetest_combat_type)
colnames(pltdat) <- c("Batch", "Stats", "combat_type")
# create plot
ggplot(data=pltdat, aes(x=Batch, y=Stats, fill=Batch)) + 
  geom_boxplot(outlier.size=0.2) + 
  facet_grid(~combat_type) + 
  stat_summary(fun.y="mean", colour="darkred", geom="point", 
               shape=18, size=0.8, show.legend = FALSE) +
  labs(y="Gene-wise Variance") +
  scale_y_continuous(breaks = seq(from=0, to=10, 1)) +
  coord_cartesian(ylim=c(0, 10)) +
  theme_bw() +
  theme(text = element_text(size=5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5),
        legend.key.size=unit(0.65,"line"))
ggsave("./figures/Fig2/Fig2.eps", width=3, height=3, units="in", dpi=550)
dev.off()



############  Fig3: bladder cancer, skewness and kurtosis density plots  ############
# sample-wise kurtosis
pltdat_1 <- data.frame(Sample=colnames(edata)[groupsorder], 
                       Stats=Y.combat[, "Kurtosis"], 
                       Batch=as.factor(paste("Batch", sortbatch, sep="")))
cdat1 <- ddply(pltdat_1, "Batch", summarise, Stats.mean=mean(Stats))
s_kurt <- ggplot(pltdat_1, aes(x=Stats, fill=Batch)) + 
  geom_density(alpha=.3) +
  facet_grid(Batch ~ .) +
  geom_vline(data=cdat1, aes(xintercept=Stats.mean, colour=Batch),
             linetype="dashed", size=0.5) +
  ggtitle("(A)") +
  labs(x="Sample-wise Kurtosis", y="Density") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(text = element_text(size=5),
        legend.position="none")
# gene-wise skewness
pltdat_2 <- melt(X.combat[["Skew"]])[,-1]
colnames(pltdat_2) <- c("Batch", "Stats")
cdat2 <- ddply(pltdat_2, "Batch", summarise, Stats.mean=mean(Stats))
g_skew <- ggplot(pltdat_2, aes(x=Stats, fill=Batch)) + 
  geom_density(alpha=.3) +
  facet_grid(Batch ~ .) +
  geom_vline(data=cdat2, aes(xintercept=Stats.mean, colour=Batch),
             linetype="dashed", size=0.5) +
  ggtitle("(B)") +
  labs(x="Gene-wise Skewness", y="") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(text = element_text(size=5),
        legend.position="none")
# gene-wise kurtosis
pltdat_3 <- melt(X.combat[["Kurtosis"]])[,-1]
colnames(pltdat_3) <- c("Batch", "Stats")
cdat3 <- ddply(pltdat_3, "Batch", summarise, Stats.mean=mean(Stats))
g_kurt <- ggplot(pltdat_3, aes(x=Stats, fill=Batch)) + 
  geom_density(alpha=.3) +
  facet_grid(Batch ~ .) +
  geom_vline(data=cdat3, aes(xintercept=Stats.mean, colour=Batch),
             linetype="dashed", size=0.5) +
  ggtitle("(C)") +
  labs(x="Gene-wise Kurtosis", y="") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(text = element_text(size=5),
        legend.position="none")
# plot in one figure
fig3 <- grid.arrange(s_kurt, g_skew, g_kurt, ncol=3)
# save plot
ggsave("./figures/Fig3/Fig3.png", plot=fig3, width=4, height=3, units="in", dpi=550)
dev.off()



############  S1 Fig: bladder cancer, plot Batch-condition-adjusted expression per sample  ############
for(method.name in c("Original", "Meanonly", "ComBat")){
  if(method.name == "Original"){pltmat <- sortdata
  }else if(method.name == "Meanonly"){pltmat <- sortdata.meanonly
  }else if(method.name == "ComBat"){pltmat <- sortdata.combat}
  
  sdata <- data.frame(Sample=factor(colnames(edata)[groupsorder], levels=colnames(edata)[groupsorder]),
                      Batch=as.factor(sortbatch), Condition=sortcond,
                      t(pltmat))
  mlt.sdata <- melt(sdata, id.vars=c("Sample","Batch","Condition"), 
                    variable_name="Gene")
  #for(i in 10000:10010){print(sum(pltmat[i,]!=mlt.sdata[(1+(i-1)*57):(57*i), "value"]))}
  colnames(mlt.sdata)[5] <- "Expression"
  png(paste("./support_info/S1_fig/bladder_Raw_", method.name, ".png", sep=""),
      width=400, height=700)
  p_gg <- ggplot(mlt.sdata, aes(x=Sample, y=Expression, fill=Batch)) + 
    geom_boxplot(notch=F, alpha=0.7) + 
    scale_y_continuous(breaks = seq(from=-5,#floor(min(pltmat)),
                                    to=7, 0.5), #ceiling(max(pltmat)), 0.5),
                       limits=c(-5, 7)) +
    theme_bw() +
    theme(text = element_text(size=20), 
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=15),
          legend.position = "bottom") 
  print(p_gg)
  dev.off()
}



############  S2 Fig: oncogenic signature, sample-wise mean and variance boxplots  ############
y_start <- c(Mean=-0.35, Variance=0.5)
y_end <- c(Mean=0.3, Variance=3.1)
y_step <- c(Mean=0.1, Variance=0.5)

for(stat.name in colnames(Y)[1:2]){
  for(method.name in c("Original", "Meanonly", "ComBat")){
    
    if(method.name == "Original"){chosenY <- Y
    }else if(method.name == "Meanonly"){chosenY <- Y.meanonly
    }else if(method.name == "ComBat"){chosenY <- Y.combat}
    
    pltdat <- data.frame(Sample=colnames(edata)[groupsorder], 
                         Statistics=chosenY[, stat.name], 
                         Batch=as.factor(paste("Batch", sortbatch, sep="")))

    png(paste("./support_info/S2_fig/scaledSig_sample", stat.name, "_box_", method.name, ".png",sep=""), 
        width=400, height=700)
    p_gg <- ggplot(data=pltdat, aes(x=Batch, y=Statistics, fill=Batch)) + 
      geom_boxplot() + 
      stat_summary(fun.y = "mean", colour="darkred", geom="point", 
                   shape=18, size=3, show.legend = FALSE) +
      labs(y=paste("Sample-wise", stat.name, sep=" ")) +
      scale_y_continuous(breaks = seq(from=y_start[stat.name], to=y_end[stat.name], y_step[stat.name]),
                         limits=c(y_start[stat.name], y_end[stat.name])) +
      theme_bw() +
      theme(text = element_text(size=20),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=15)) 
    print(p_gg)
    dev.off()
  }
}



############  S3 Fig: nitric oxide, gene-wise variance boxplot  ############
for(method.name in c("Original", "Meanonly", "ComBat")){
  if(method.name == "Original"){chosenX <- X
  }else if(method.name == "Meanonly"){chosenX <- X.meanonly
  }else if(method.name == "ComBat"){chosenX <- X.combat}
  
  prdctrs <- factor(paste("Batch", rep(unique(sortbatch), each=nrow(sortdata)), sep=""))
  pltdat <- data.frame(Batch=prdctrs,
                       Statistics=c(chosenX[["Variance"]])) 
  
  png(paste("./support_info/S3_fig/NO_gene",stat.name,"_box_", method.name, ".png", sep=""), 
      width=400, height=700)
  p_gg <- ggplot(data=pltdat, aes(x=Batch, y=Statistics, fill=Batch)) + 
    geom_boxplot() + 
    stat_summary(fun.y = "mean", colour="darkred", geom="point", 
                 shape=18, size=3, show.legend = FALSE) +
    labs(y=paste("Gene-wise", stat.name, sep=" ")) +
    scale_y_continuous(breaks = seq(from=0, to=3.5, 0.5),
                       limits=c(0, 3.5)) +
    theme_bw() +
    theme(text = element_text(size=20),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=15)) 
  print(p_gg)
  dev.off()
}    
