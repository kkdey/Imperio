
geneset="pLI" # PPI_enhancer
snp_s2g_dir = paste0("../data/Imperio_SNP_S2G_MAP_", geneset)

for(numchr in 1:22){
  df1 = data.frame(fread(paste0(snp_s2g_dir, "/",
                                "100kb", "/",
                                "100kb", ".",
                                numchr, ".annot.gz")))
  df2 = data.frame(fread(paste0(snp_s2g_dir, "/",
                                "5kb", "/",
                                "5kb", ".",
                                numchr, ".annot.gz")))
  df3 = data.frame(fread(paste0(snp_s2g_dir, "/",
                                "TSS", "/",
                                "TSS", ".",
                                numchr, ".annot.gz")))
  df4 = data.frame(fread(paste0(snp_s2g_dir, "/",
                                "Roadmap_Blood", "/",
                                "Roadmap_Blood", ".",
                                numchr, ".annot.gz")))
  df5 = data.frame(fread(paste0(snp_s2g_dir, "/",
                                "ABC9_Blood", "/",
                                "ABC9_Blood", ".",
                                numchr, ".annot.gz")))
  newdf = cbind.data.frame(df1[,1:4], df1[,5], df2[,5], df3[,5], df4[,5], df5[,5])
  colnames(newdf) = c(colnames(df1)[1:4], "100kb", "5kb", "TSS", "Roadmap", "ABC9")

  if(!dir.exists(paste0(snp_s2g_dir, "/", "ALL_S2G_BLOOD"))){
    dir.create(paste0(snp_s2g_dir, "/", "ALL_S2G_BLOOD"))
  }

  write.table(newdf,
              file = gzfile(paste0(snp_s2g_dir, "/",
                                   "ALL_S2G_BLOOD", "/",
                                   "ALL_S2G_BLOOD", ".",
                                   numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chrom:", numchr, "\n")

}

