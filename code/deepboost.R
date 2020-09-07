library(data.table)
library(xgboost)


deep_dir="data/DeepSEA2018/SAD_dump"
# deep_dir="data/Basenji/SAD_dump"
pos_file = "data/positive_snps.txt"
control_file = "data/control_snps.txt"
bimpath = "data/BIMS"
output_path = "data/deepboost_testpred"

xgboost_train = function(train_chr){
  train_data = c()
  train_labels = c()
  control_snps = read.delim(control_file, header=F)[,1]
  finemap_snps = read.delim(pos_file, header=F)[,1]

  for(numchr in train_chr){

    bimfile = read.delim(paste0(bimpath, "/", "1000G.EUR.QC.", numchr, ".bim"), header=F)
    colnames(bimfile) = c("CHR", "SNP", "CM", "BP", "REF", "ALT")
    fine_indices = which(!is.na(match(bimfile$SNP, finemap_snps)))
    control_indices = which(!is.na(match(bimfile$SNP, control_snps)))
    control_labels = rep(0, nrow(bimfile))
    control_labels[control_indices] = 1
    finemap_labels = rep(0, nrow(bimfile))
    finemap_labels[fine_indices] = 1

    num_snps = 0
    files = list.files(paste0(deep_dir, "/chr", numchr))
    ll = vector(mode="list", length = length(files))

    for(cc in 1:length(files)){
      temp = data.frame(fread(paste0(deep_dir, "/chr", numchr, "/", files[cc])))
      control_idx = which(!is.na(match(temp[,2], control_snps)))
      finemap_idx = which(!is.na(match(temp[,2], finemap_snps)))
      train_data = rbind(train_data, rbind(temp[control_idx, -(1:6)], temp[finemap_idx, -(1:6)]))
      train_labels = c(train_labels, c(rep(0, length(control_idx)),
                                       rep(1, length(finemap_idx))))
      num_snps = num_snps + nrow(temp)
      cat("Read file:", files[cc], "\n")
    }
    cat("Number of SNPs investigated in chr:", numchr, " is ", num_snps, "\n")
  }
  train_set = list("data" = as.matrix(train_data), "labels" = train_labels)
  return(train_set)
}


xgboost_test = function(bst, test_chr){
  for(numchr in test_chr){
    bimfile = read.delim(paste0(bimpath, "/", "1000G.EUR.QC.", numchr, ".bim"), header=F)
    colnames(bimfile) = c("CHR", "SNP", "CM", "BP", "REF", "ALT")
    files = list.files(paste0(deep_dir, "/chr", numchr))
    pred = c()
    for(cc in 1:length(files)){
      temp = data.frame(fread(paste0(deep_dir, "/chr", numchr, "/", files[cc])))
      test_data = as.matrix(temp[,-(1:6)])
      pred <- c(pred, predict(bst, as.matrix(test_data)))
      cat("Read file:", files[cc], "\n")
    }
    newdf = cbind.data.frame(bimfile[,1:4], pred)
    colnames(newdf) = c("CHR", "BP", "SNP", "CM", "AN")

    if(!dir.exists(paste0(output_path))){
      dir.create(paste0(output_path))
    }
    write.table(newdf,
                file = gzfile(paste0(output_path, "/", "test", ".", numchr, ".annot.gz")),
                quote=FALSE, row.names=FALSE)
    cat("Prediction finished for  chrom:", numchr, "\n")
  }
}

##############################  Train on Odd chr and Predict on Even chr  ####################################

res = xgboost_train(seq(1,22,2))
bstSparse <-  xgboost(data = res$data, label = res$labels,
                      n_estimators = c(200, 250, 300),
                      max_depth = c(25, 30, 35),
                      learning_rate = 0.05,
                      gamma = 10,
                      min_child_weight = c(6, 8, 10),
                      nthread = 2,
                      scale_pos_weight = 1,
                      subsample = c(0.6, 0.8, 1),
                      nrounds = 200,
                      objective = "binary:logistic",
                      eval_metric = "auc")
xgboost_test(bstSparse, seq(2,22,2))


##############################  Train on Even chr and Predict on Odd chr  ####################################

res = xgboost_train(seq(2,22,2))
bstSparse <- xgboost(data = res$data, label = res$labels,
                     n_estimators = c(200, 250, 300),
                     max_depth = c(25, 30, 35),
                     learning_rate = 0.05,
                     gamma = 10,
                     min_child_weight = c(6, 8, 10),
                     nthread = 2,
                     scale_pos_weight = 1,
                     subsample = c(0.6, 0.8, 1),
                     nrounds = 500,
                     objective = "binary:logistic",
                     eval_metric = "auc")
xgboost_test(bstSparse, seq(1,22,2))

