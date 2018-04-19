rm(list=ls())
#setwd("~/Dropbox/Work/EvanJohnson/refComBat_forSubmittion_BMCbioinfo/review20180320/modsel_sensitivity")
setwd("/restricted/projectnb/combat/work/yuqingz/refcombat_review_201803/modsel_sensitivity")
lapply(c("limma", "sva", "BatchQC", "ggplot2", "gridExtra", "reshape"), require, character.only = TRUE)

exp_name <- "sim2" #exp_name <- "nitricOxide"

## read in result files - tradeoff
tradeoff_none <- read.csv(sprintf('%s/tradeoff_none.csv', exp_name))
tradeoff_meanonly <- read.csv(sprintf('%s/tradeoff_meanonly.csv', exp_name))
tradeoff_meanvar <- read.csv(sprintf('%s/tradeoff_meanvar.csv', exp_name))

tradeoff_merge <- cbind(rbind(tradeoff_meanonly, tradeoff_meanvar),
                        data=c(rep("meanonly", nrow(tradeoff_meanonly)), rep("meanvar", nrow(tradeoff_meanvar))))
tradeoff_merge_mlt <- melt(tradeoff_merge, id.vars="data", variable_name="batchAdj")

png(sprintf('%s/tradeoff_%s.png', exp_name, exp_name),
    width=5, height=4, units="in", res=300)
ggplot(tradeoff_merge_mlt, aes(x=batchAdj, y=value)) +
  geom_boxplot() +
  facet_grid(~data) +
  scale_x_discrete(name="Batch adjustment model") +
  scale_y_continuous(name="Power VS type-I error rate trade-off") +
  #scale_y_continuous(name="Power VS type-I error rate trade-off", limits=c(0,200)) +   
  # for real data (nitric Oxide)
  ggtitle("Type of batch effect") +
  theme(plot.title = element_text(hjust = 0.5, size=10),
        axis.title = element_text(size=10))
dev.off()
  