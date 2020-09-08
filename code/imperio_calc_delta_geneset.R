options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
numchr <- as.numeric(toString(args[1])) ## 1, 2, ..., 22
method <- toString(args[2]) ## DeepSEA or Basenji
outname <- toString(args[3]) #outname="deep_imperios2g_geneset_annot"
geneset <- toString(args[4])

library(rhdf5)
library(testit)
library(data.table)
library(xgboost)

xgboost_path = "../data/Xgboost_saved_models"
deep_dir = "../data/Basenji/SAD_dump"
snp_s2g_dir = paste0("../data/Imperio_SNP_S2G_MAP_", geneset)

model_name = paste0(xgboost_path, "/", "imperio_model", "_", geneset, ".save")

if(method == "DeepSEA"){
  files = list.files(paste0(deep_dir, "/chr", numchr))
  ll = vector(mode="list", length = length(files))
  idx_list = vector(mode="list", length = length(files))

  num_snps = 0
  for(cc in 1:length(files)){
    temp = data.frame(fread(paste0(deep_dir, "/chr", numchr, "/", files[cc])))
    ll[[cc]] = t(temp[,-c(1:6)])
    idx_list[[cc]] = num_snps + 1:ncol(ll[[cc]])
    num_snps = num_snps + ncol(ll[[cc]])
  }
  print(num_snps)
  df = data.frame(fread(paste0(snp_s2g_dir, "/",
                               "ALL_S2G_BLOOD", "/",
                               "ALL_S2G_BLOOD", ".", numchr, ".annot.gz")))

  nrow(df)

  chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
  split_indices = chunk(1:nrow(df), ceiling(nrow(df)/50000))
  sapply(split_indices, function(x) length(x))
  sapply(ll, function(x) ncol(x))

  assert(length(ll) == length(split_indices))

}else if (method == "Basenji"){
  files = list.files(paste0(deep_dir, "/chr", numchr))
  ll = vector(mode="list", length = length(files))

  num_snps = 0
  for(cc in 1:length(files)){
    ll[[cc]] = h5read(paste0(deep_dir, "/chr", numchr, "/1000G.SAD.", numchr, ".chunk.", cc, "/sad.h5"), "/SAD")
    num_snps = num_snps + ncol(ll[[cc]])
  }
  print(num_snps)
  df = data.frame(fread(paste0(snp_s2g_dir, "/",
                               "ALL_S2G_BLOOD", "/",
                               "ALL_S2G_BLOOD", ".", numchr, ".annot.gz")))

  nrow(df)

  chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
  split_indices = chunk(1:nrow(df), ceiling(nrow(df)/50000))
  sapply(split_indices, function(x) length(x))
  sapply(ll, function(x) ncol(x))

  assert(length(ll) == length(split_indices))
}else{
  stop("Method must be DeepSEA or Basenji")
}

boost = xgb.load(model_name)

total_preds = c()
predList = vector(mode="list", length = length(ll))
for(cc in 1:length(ll)){
  temp = cbind(rep(0, ncol(ll[[cc]])), t(ll[[cc]]))
  try = cbind(matrix(0, nrow(temp), ncol(temp)), matrix(0, nrow(temp), ncol(temp)), matrix(0, nrow(temp), ncol(temp)),
              matrix(0, nrow(temp), ncol(temp)), matrix(0, nrow(temp), ncol(temp)))
  baseline = min(predict(boost, try))  ### try is a vector of same value

  try = cbind(temp, matrix(0, nrow(temp), ncol(temp)), matrix(0, nrow(temp), ncol(temp)),
              matrix(0, nrow(temp), ncol(temp)), matrix(0, nrow(temp), ncol(temp)))
  preds1 = predict(boost, try) - baseline

  try = cbind(matrix(0, nrow(temp), ncol(temp)), temp, matrix(0, nrow(temp), ncol(temp)),
              matrix(0, nrow(temp), ncol(temp)), matrix(0, nrow(temp), ncol(temp)))
  preds2 = predict(boost, try) - baseline

  try = cbind(matrix(0, nrow(temp), ncol(temp)),  matrix(0, nrow(temp), ncol(temp)), temp,
              matrix(0, nrow(temp), ncol(temp)), matrix(0, nrow(temp), ncol(temp)))
  preds3 = predict(boost, try) - baseline

  try = cbind(matrix(0, nrow(temp), ncol(temp)),  matrix(0, nrow(temp), ncol(temp)),
              matrix(0, nrow(temp), ncol(temp)), temp, matrix(0, nrow(temp), ncol(temp)))
  preds4 = predict(boost, try) - baseline

  try = cbind(matrix(0, nrow(temp), ncol(temp)),  matrix(0, nrow(temp), ncol(temp)),
              matrix(0, nrow(temp), ncol(temp)),  matrix(0, nrow(temp), ncol(temp)), temp)
  preds5 = predict(boost, try) - baseline


  predList[[cc]] = cbind(preds1, preds2, preds3, preds4, preds5)
  cat("Predicted for chunk:", cc, "of", length(ll), " chunks \n")
}

total_preds = do.call("rbind", predList)

newdf = cbind.data.frame(df[,1:4], total_preds)
colnames(newdf) = c(colnames(df)[1:4], "5kb-specific", "100kb-specific",
                    "TSS-specific", "ABC9-specific", "Roadmap-specific")


if(!dir.exists(paste0("../data/", outname))){
  dir.create(paste0("../data/", outname))
}
write.table(newdf,
            file = gzfile(paste0("../data/",
                                 outname, "/",
                                 outname, ".",
                                 numchr, ".annot.gz")),
            quote=FALSE, row.names=FALSE)




