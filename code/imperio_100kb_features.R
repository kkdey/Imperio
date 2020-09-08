options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
numchr <- as.numeric(toString(args[1])) ## 1, 2, ...., 22
method <- toString(args[2])  ## DeepSEA, Basenji
outpath <- toString(args[3])

library(rhdf5)
library(testit)
library(data.table)


get_unlinked_from_linked_snps = function(vec){
  v1 = min(vec)
  v2 = max(vec)
  ff = c()
  all_seqs = seq(v1, v2, by = 1000)
  all_seqs = c(all_seqs, max(all_seqs, v2))
  for(kk in 1:(length(all_seqs)-1)){
    idx2 = which(vec >= all_seqs[kk] & vec <= all_seqs[kk+1])
    ff = c(ff, vec[idx2[ceiling(length(idx2)/2)]])
  }
  return(ff)
}

if(!dir.exists(paste0(output_dir, "/100kb"))){
  dir.create(paste0(output_dir, "/100kb"))
}

geneanno = data.table::fread("../data/geneanno.csv")
gene_beds = read.table("../data/Gene_100kb.txt", header=T)

common_genes = intersect(gene_beds$gene, geneanno$symbol)
geneanno2 = geneanno[match(common_genes, geneanno$symbol), ]
gene_beds2 = gene_beds[match(common_genes, gene_beds$gene), ]



if(method == "DeepSEA"){
  deep_dir = "../data/DeepSEA2018/REF_dump"
  output_dir = outpath
  files = list.files(paste0(deep_dir, "/chr", numchr), pattern = ".csv.gz")
  ll = vector(mode="list", length = length(files))

  num_snps = 0
  tt = c()
  for(num in 1:length(ll)){
    temp = data.frame(fread(paste0(deep_dir, "/chr", numchr,
                                   "/", files[num])))

    ll[[num]] <- t(temp[,-(1:6)])
    cat("Read chunk:", num, "\n")
    num_snps = num_snps + ncol(ll[[num]])
    tt = rbind(tt, cbind.data.frame(num, 1:ncol(ll[[num]])))
  }
}else if (method == "Basenji"){
  deep_dir = "../data/Basenji/REF_dump"
  output_dir = outpath
  files = list.files(paste0(deep_dir, "/chr", numchr), pattern = ".h5")
  ll = vector(mode="list", length = length(files))

  num_snps = 0
  tt = c()
  for(num in 1:length(ll)){
    ll[[num]] <- h5read(paste0("/n/groups/price/tier2/kushal/ALL_DEEP/Basenji/REF_dump/chr", numchr,
                               "/", "chr", numchr, "_chunk_", num, ".h5"), "/preds")
    cat("Read chunk:", paste0("/n/groups/price/tier2/kushal/ALL_DEEP/Basenji/REF_dump/chr", numchr,
                              "/", "chr", numchr, "_chunk_", num, ".h5"), "\n")
    num_snps = num_snps + ncol(ll[[num]])
    tt = rbind(tt, cbind.data.frame(num, 1:ncol(ll[[num]])))
  }
  h5closeAll()
}else{
  stop("The method option must be DeepSEA or Basenji")
}

base = data.frame(fread(paste0("zcat ../data/Baselines/",
                               "baselineLD_v2.1", "/",
                               "baselineLD", ".", numchr, ".annot.gz")))


assert(nrow(base) == num_snps)  ## test if total snps in basenji ref files matches baselineLD
tt = cbind.data.frame(tt, base$BP)
colnames(tt) = c("chr_chunk", "id_chunk", "base_bp")

genes_in_chr = common_genes[which(geneanno2$seqnames == paste0("chr", numchr))]
temp = c()

merge_df = c()
for(gg in genes_in_chr){

  ###  Get starts and ends for 100kb regions for each gene  gg  #######################
  start = gene_beds2$start[which(gene_beds2$gene == gg)]
  end =  gene_beds2$end[which(gene_beds2$gene == gg)]
  snp_locs = which(base$BP > start & base$BP < end)
  if(length(snp_locs) >0){
    unlinked_loci = get_unlinked_from_linked_snps(base$BP[snp_locs])
    tt2 = tt[match(unlinked_loci, tt$base_bp), ]
    uu = unique(tt2[,1])
    kk = c()
    for(num in 1:length(uu)){
      kk = rbind(kk, t(ll[[uu[num]]][, tt2$id_chunk[which(tt2$chr_chunk == uu[num])], drop=F]))
    }
    merge_df = rbind(merge_df, c(nrow(kk), apply(kk, 2, sum)))
  }else{
    merge_df = rbind(merge_df, c(0, rep(NA, nrow(ll[[1]]))))
  }
  cat("Recorded data for gene:", gg, '\n')
}
rownames(merge_df) = genes_in_chr

save(merge_df, file = paste0(output_dir, "/100kb/chr_", numchr, ".rda"))









