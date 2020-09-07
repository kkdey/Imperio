
ll=list.files("../data/allgenes_s2g/")
library(data.table)
for(numchr in 1:22){
  df = data.frame(fread(paste0("../data/deepboost_annot/","deepboost_annot", ".", numchr, ".annot.gz")))
  binannot = rep(0, nrow(df))
  binannot[which(df[,5] > quantile(df[,5], 0.85))] = 1
  newdf = cbind.data.frame(df[,1:4], binannot)
  colnames(newdf) = c(colnames(df)[1:4], "AN")

  if(!dir.exists(paste0("../data/deepboost_annot_output/", "deepboost_annot_binary"))){
    dir.create(paste0("../data/deepboost_annot_output/", "deepboost_annot_binary"))
  }

  write.table(newdf,
              file = gzfile(paste0("../data/deepboost_annot_output/",
                                   "deepboost_annot_binary", "/",
                                   "deepboost_annot_binary", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)


  for(numl in 1:length(ll)){
    temp = data.frame(fread(paste0("../data/allgenes_s2g/", ll[numl], "/", ll[numl], ".", numchr, ".annot.gz")))
    newdf2 = cbind.data.frame(temp[,1:4], temp[,5]*newdf[,5])
    colnames(newdf2) = c(colnames(temp)[1:4], "AN")
    if(!dir.exists(paste0("../data/deepboost_annot_output/",
                          "deepboost_annot_binary_", ll[numl]))){
      dir.create(paste0("../data/deepboost_annot_output/",
                        "deepboost_annot_binary_", ll[numl]))}
    write.table(newdf2,
                file = gzfile(paste0("../data/deepboost_annot_output/",
                                     "deepboost_annot_binary_", ll[numl], "/",
                                     "deepboost_annot_binary_", ll[numl], ".", numchr, ".annot.gz")),
                quote=FALSE, row.names=FALSE)
  }


  cat("We are at chr:", numchr, "\n")

}
