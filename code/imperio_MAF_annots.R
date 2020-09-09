outname = "deep_imperios2g_annot"
snp_s2g_dir = "../data/Imperio_SNP_S2G_MAP"

if(!dir.exists(paste0("../data/", outname, "_MAF"))){
  dir.create(paste0("../data/", outname, "_MAF"))
}

prioritize_vec = function(xx){
  y = 1-ecdf(xx)(xx)
  z = -2*log(y+1e-08)
  PR = (z - min(z))/(max(z) - min(z))
  return(PR)
}

prioritize_vec2 = function(xx){
  PR = (xx - min(xx))/(max(xx) - min(xx))
  return(PR)
}


for(numchr in 1:22){
  maf_df = data.frame(fread(paste0("../data/1000G_Phase3_frq/",
                                   "1000G.EUR.QC.", numchr, ".frq")))
  df = data.frame(fread(paste0("../data/",
                               outname, "/",
                               outname, ".",
                               numchr, ".annot.gz")))
  df2 = data.frame(fread(paste0(snp_s2g_dir, "/",
                                "ALL_S2G_BLOOD", "/",
                                "ALL_S2G_BLOOD", ".", numchr, ".annot.gz")))
  annot = apply(temp_mat2, 1, sum)
  vec = prioritize_vec2((maf_df$MAF)*(1 - maf_df$MAF)*annot)
  newdf = cbind.data.frame(df[,1:4], vec)
  colnames(newdf) = c(colnames(df)[1:4], "AN")
  if(!dir.exists(paste0("../data/", outname, "_MAF"))){
    dir.create(paste0("../data/", outname, "_MAF"))
  }
  write.table(newdf,
              file = gzfile(paste0("../data/", outname, "_MAF", "/",
                                   outname, "_MAF", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chr:", numchr, "\n")
}



