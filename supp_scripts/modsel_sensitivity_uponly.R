rm(list=ls())
#setwd("C:/Users/zhang/Dropbox/Work/EvanJohnson/refComBat_forSubmittion_BMCbioinfo/review20180320/modsel_sensitivity")
#setwd("~/Dropbox/Work/EvanJohnson/refComBat_forSubmittion_BMCbioinfo/review20180320/modsel_sensitivity")
setwd("/restricted/projectnb/combat/work/yuqingz/refcombat_review_201803/modsel_sensitivity")
lapply(c("limma", "sva", "BatchQC", "ggplot2", "gridExtra", "reshape"), require, character.only = TRUE)


####  Parameters manually set
command_args <- commandArgs(trailingOnly=TRUE)
exp_name <- command_args[1]
type <- command_args[2]
par_type <- command_args[3]
if(!dir.exists(exp_name)){dir.create(exp_name)}
#exp_name <- "test"; type="none"; par_type="sim"

if(par_type=="real"){
  ## Using batch and bio estimates from real data
  # use parameters estimated from oncogenic signature data 
  # load("meanvar_sigBio.rda")
  # load("meanvar_sigBatch.rda")
  # use parameters estimated from nitric oxide data
  load("meanvar_NO_bio.rda")
  load("meanvar_NO_batch.rda")
  batch_mean_diff <- c(mean3_up-mean1_up, mean3_down-mean1_down, mean3_nonsig-mean1_nonsig)
  batch_var_diff <-c(var3_up/var1_up, var3_down/var1_down, var3_nonsig/var1_nonsig)
}else if(par_type=="sim"){
  ## Use manually set batch and bio estimates
  # bio
  mean_ctrl_up <- 0; mean_case_up <- 1
  mean_ctrl_down <- 1; mean_case_down <- 0
  mean_ctrl_nonsig <- mean_case_nonsig <- 0
  var_ctrl_up <- var_case_up <- 0.3
  var_ctrl_down <- var_case_down <- 0.3
  var_ctrl_nonsig <- var_case_nonsig <- 0.3
  # batch
  batch_mean_diff <- rep(-0.5,3) # c(-0.25, -0.25, -0.25)
  batch_var_diff <- rep(3,3) #c(1.25, 1.25, 1.25)
}else{
  stop("Wrong parameter type!")
}


## Other parameters
N_sims <- 100
N_genes <- 10000
N_up <- 100
N_nonsig <- N_genes - N_up 
N_samples_vec <- c(5, 5, 5, 5)

# use bio parameters for batch 1 as baseline, then add batch differences to batch 2
# batch 1 parameter
batch1_pars <- data.frame(up=c(mean_ctrl=mean_ctrl_up, mean_case=mean_case_up, 
                               var_ctrl=var_ctrl_up, var_case=var_case_up),
                          down=c(mean_ctrl=mean_ctrl_down, mean_case=mean_case_down, 
                                 var_ctrl=var_ctrl_down, var_case=var_case_down),
                          nonsig=c(mean_ctrl=mean_ctrl_nonsig, mean_case=mean_case_nonsig, 
                                   var_ctrl=var_ctrl_nonsig, var_case=var_case_nonsig))
batch1_pars <- t(batch1_pars)

# batch 2 parameter
if(type=="none"){
  batch_mean_diff <- 0
  batch_var_diff <- 1
}else if(type=="meanonly"){
  batch_var_diff <- 1
}
batch2_pars <- cbind(batch1_pars[,1] + batch_mean_diff, batch1_pars[,2] + batch_mean_diff,
                     batch1_pars[,3] * batch_var_diff, batch1_pars[,4] * batch_var_diff)
colnames(batch2_pars) <- colnames(batch1_pars)

# sanity check
if(type=="none" & (!identical(batch2_pars, batch1_pars))){stop("Wrong batch parameter calculation!")}
if(type=="meanonly" & (!identical(batch2_pars[, 3:4], batch1_pars[, 3:4]))){stop("Wrong batch parameter calculation!")}

# split into mean par matrix and variance par matrix
mean_pars <- cbind(batch1_pars[, 1:2], batch2_pars[, 1:2])
var_pars <- cbind(batch1_pars[, 3:4], batch2_pars[, 3:4])
sd_pars <- sqrt(var_pars)

mean_pars <- mean_pars[c("up", "nonsig"), ]
var_pars <- var_pars[c("up", "nonsig"), ]
sd_pars <- sd_pars[c("up", "nonsig"), ]



####  Objects computed
genes_id <- c(rep("up", N_up), rep("nonsig", N_nonsig))
batch <- c(rep(1, sum(N_samples_vec[1:2])), rep(2, sum(N_samples_vec[3:4])))
condition <- c(rep(0, N_samples_vec[1]), rep(1, N_samples_vec[2]),
               rep(0, N_samples_vec[3]), rep(1, N_samples_vec[4]))


####  Run simulations
for(ID in 1:N_sims){
  print(paste("Simulation", ID))
  
  ####  Simulate dataset 
  dat_mat <- matrix(NA, nrow=N_genes, ncol=sum(N_samples_vec))
  for(ii in 1:N_genes){
    dat_mat[ii,batch==1&condition==0]<-rnorm(N_samples_vec[1], mean_pars[genes_id[ii],1], sd_pars[genes_id[ii],1])
    dat_mat[ii,batch==1&condition==1]<-rnorm(N_samples_vec[2], mean_pars[genes_id[ii],2], sd_pars[genes_id[ii],2])
    dat_mat[ii,batch==2&condition==0]<-rnorm(N_samples_vec[3], mean_pars[genes_id[ii],3], sd_pars[genes_id[ii],3])
    dat_mat[ii,batch==2&condition==1]<-rnorm(N_samples_vec[4], mean_pars[genes_id[ii],4], sd_pars[genes_id[ii],4])
  }
  
  
  ####  Normalization
  dat_mat <- gnormalize(dat_mat)
  
  
  ####  Adjust batch effect
  # dat_meanonly <- ComBat(dat_mat, batch=batch, mod=model.matrix(~factor(condition)), mean.only=TRUE)
  # dat_meanvar <- ComBat(dat_mat, batch=batch, mod=model.matrix(~factor(condition)))
  dat_meanonly <- ComBat(dat_mat, batch=batch, mod=NULL, mean.only=TRUE)
  dat_meanvar <- ComBat(dat_mat, batch=batch, mod=NULL)
  
  data_lst <- list(original=dat_mat, meanonly=dat_meanonly, meanvar=dat_meanvar)
  
  
  ####  DE analysis with limma
  DE_res <- lapply(data_lst, function(dat, group, N_nonsig, N_up){
    design <- model.matrix(~factor(group))
    fit <- lmFit(dat, design)
    fit2 <- eBayes(fit)
    res <- topTable(fit2, coef=1, adjust.method="BH", number=nrow(dat))
    called_ind <- as.numeric(rownames(res[res$P.Value < 0.05, ]))
    # false positive rate
    fp <-  sum(called_ind %in% which(genes_id=="nonsig")) / N_nonsig    
    # power 
    pwr <- sum(called_ind %in% which(genes_id!="nonsig")) / N_up    
    return(c(fp=fp, pwr=pwr))
  }, group=condition, N_nonsig=N_nonsig, N_up=N_up)
  DE_res <- as.data.frame(do.call(cbind, DE_res))
  #print(round(DE_res,3))
  
  
  ####  Type 1 error VS power trade-off
  pwr_fp_trade <- data.frame(meanonly=DE_res$meanonly - DE_res$original,
                             meanvar=DE_res$meanvar - DE_res$original)
  rownames(pwr_fp_trade) <- rownames(DE_res)
  
  
  ####  Output results
  first.file <- !file.exists(sprintf('%s/fp_%s.csv', exp_name, type))
  # type 1 error
  write.table(DE_res["fp", ], sprintf('%s/fp_%s.csv', exp_name, type), 
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",")
  # power
  write.table(DE_res["pwr", ], sprintf('%s/pwr_%s.csv', exp_name, type), 
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",")
  # trade-off
  write.table(as.data.frame(pwr_fp_trade["pwr",]/pwr_fp_trade["fp",]), 
              sprintf('%s/tradeoff_%s.csv', exp_name, type), 
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",")
}
