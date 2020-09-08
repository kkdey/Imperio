library(rhdf5)
library(testit)
library(data.table)
outpath="../data/Imperio_models"


###########################################  5kb  ##########################################

total_df = c()
for(numchr in 1:22){
  temp = get(load(paste0(outpath, "/",
                       "5kb", "/", "chr_", numchr, ".rda")))
  total_df = rbind(total_df, temp)
  cat("We are at chr:", numchr, "\n")
}

total_df2 = total_df[which(!is.na(rowMeans(total_df))), ]
total_df3 = total_df2
rownames(total_df3) = rownames(total_df2)

save(total_df3, file = paste0(outpath, "/",
      "5kb_full.rda"))


###########################################  100kb  ##########################################

total_df = c()
for(numchr in 1:22){
  temp = get(load(paste0(outpath, "/",
                         "100kb", "/", "chr_", numchr, ".rda")))
  total_df = rbind(total_df, temp)
  cat("We are at chr:", numchr, "\n")
}

total_df2 = total_df[which(!is.na(rowMeans(total_df))), ]

total_df3 = total_df2
rownames(total_df3) = rownames(total_df2)

save(total_df3, file = paste0(outpath, "/",
                              "100kb_full.rda"))

###########################################  TSS  ##########################################

total_df = c()
for(numchr in 1:22){
  temp = get(load(paste0(outpath, "/",
                         "TSS", "/", "chr_", numchr, ".rda")))
  total_df = rbind(total_df, temp)
  cat("We are at chr:", numchr, "\n")
}

total_df2 = total_df[which(!is.na(rowMeans(total_df))), ]
total_df3 = total_df2
rownames(total_df3) = rownames(total_df2)

save(total_df3, file = paste0(outpath, "/",
                              "TSS_full.rda"))



###########################################  ABC  ##########################################

total_df = c()
for(numchr in 1:22){
  temp = get(load(paste0(outpath, "/",
                         "ABC", "/", "chr_", numchr, ".rda")))
  total_df = rbind(total_df, temp)
  cat("We are at chr:", numchr, "\n")
}

total_df2 = total_df[which(!is.na(rowMeans(total_df))), ]
total_df3 = total_df2
rownames(total_df3) = rownames(total_df2)

save(total_df3, file = paste0(outpath, "/",
                              "ABC9_full.rda"))


###########################################  Roadmap  ##########################################

total_df = c()
for(numchr in 1:22){
  temp = get(load(paste0(outpath, "/",
                         "Roadmap", "/", "chr_", numchr, ".rda")))
  total_df = rbind(total_df, temp)
  cat("We are at chr:", numchr, "\n")
}

total_df2 = total_df[which(!is.na(rowMeans(total_df))), ]
total_df3 = total_df2
rownames(total_df3) = rownames(total_df2)

save(total_df3, file = paste0(outpath, "/",
                              "Roadmap_full.rda"))

