clustAccu <- function(cluster_vec, n_gene, n_biomarker){
  if(length(cluster_vec)!=n_gene){print("Wrong clustering assingment vector!")}
  tab1 <- table(cluster_vec[1:n_biomarker])
  tab2 <- table(cluster_vec[(n_biomarker+1):n_gene])
  if(length(tab1)==1){
    tmp_table <- matrix(0, nrow=1, ncol=2)
    colnames(tmp_table) <- c("1","2")
    tmp_table[, names(tab1)] <- tab1[1]
    tab1 <- as.table(tmp_table)
    names(tab1) <- c("1","2")
  }
  if(length(tab2)==1){
    tmp_table <- matrix(0, nrow=1, ncol=2)
    colnames(tmp_table) <- c("1","2")
    tmp_table[, names(tab2)] <- tab2[1]
    tab2 <- as.table(tmp_table)
    names(tab2) <- c("1","2")
  }
  accu1 <- (tab1["1"] + tab2["2"]) / n_gene
  accu2 <- (tab1["2"] + tab2["1"]) / n_gene
  accu <- max(accu1, accu2)  
  return(accu)
}