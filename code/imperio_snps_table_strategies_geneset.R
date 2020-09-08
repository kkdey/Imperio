options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
numchr <- as.numeric(toString(args[1]))
geneset <- toString(args[2])

library(rhdf5)
library(testit)
library(data.table)

if(!dir.exists(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset))){
  dir.create(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset))
}


# geneanno = data.table::fread("/n/groups/price/kushal/ExPecto/resources/geneanno.csv")
# gene_set = read.table("/n/groups/price/kushal/Enhancer_MasterReg/data/Gene_Scores2/RWR_Master_Reg_TF.txt")
# imp_genes = gene_set[,1]
# geneanno2 = geneanno[match(intersect(imp_genes, geneanno$symbol), geneanno$symbol), ]
# write.csv(geneanno2, "/n/groups/price/kushal/ExPecto/resources/geneanno_PPI_master.csv")
#


library(rhdf5)
library(testit)
library(data.table)

geneanno = data.table::fread(paste0("/n/groups/price/kushal/ExPecto/resources/geneanno_", geneset, ".csv"))
base = data.frame(fread(paste0("zcat /n/groups/price/kushal/DiseaseNet/data/ANNOTATIONS/Baselines/",
                               "baselineLD_v2.1", "/",
                               "baselineLD", ".", numchr, ".annot.gz")))



########################################   ABC    ##################################################

abc_per_gene = get(load("/n/groups/price/kushal/Imperio/data/ABC9_per_gene_beds_Blood.rda"))

common_genes = intersect(names(abc_per_gene), geneanno$symbol)
abc_per_gene2 = abc_per_gene[match(common_genes, names(abc_per_gene))]
geneanno2 = geneanno[match(common_genes, geneanno$symbol), ]



genes_in_chr = common_genes[which(geneanno2$seqnames == paste0("chr", numchr))]
temp = c()

summ = rep(0, nrow(base))
for(gg in genes_in_chr){

  ###  Get starts and ends for ABC regions for each gene  gg  #######################
  starts = abc_per_gene[[gg]][,2]
  ends =  abc_per_gene[[gg]][,3]

  idx2 = c()
  for(cc in 1:length(starts)){
    idx2 = c(idx2, which(base$BP > starts[cc] & base$BP < ends[cc]))
  }
  snp_locs = as.numeric(names(table(idx2)))
  summ[snp_locs] = summ[snp_locs]+1
}

newdf = cbind.data.frame(base[,1:4], summ)
colnames(newdf) = c(colnames(base)[1:4], "AN")
rm(summ)

if(!dir.exists(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/",  "ABC9_Blood"))){
  dir.create(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "ABC9_Blood"))
}

write.table(newdf,
            file = gzfile(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/",
                                 "ABC9_Blood", "/",
                                 "ABC9_Blood", ".",
                                 numchr, ".annot.gz")),
            quote=FALSE, row.names=FALSE)
cat("We are at chrom:", numchr, "\n")


##########################   Roadmap    #########################################################

geneanno = data.table::fread(paste0("/n/groups/price/kushal/ExPecto/resources/geneanno_", geneset, ".csv"))
road_per_gene = get(load("/n/groups/price/kushal/Imperio/data/Roadmap_per_gene_beds_Blood.rda"))

common_genes = intersect(names(road_per_gene), geneanno$id)
road_per_gene2 = road_per_gene[match(common_genes, names(road_per_gene))]
geneanno2 = geneanno[match(common_genes, geneanno$id), ]

genes_in_chr = common_genes[which(geneanno2$seqnames == paste0("chr", numchr))]

summ = rep(0, nrow(base))

for(gg in genes_in_chr){

  ###  Get starts and ends for Roadmap regions for each gene  gg  #######################
  starts = road_per_gene[[gg]][,2]
  ends =  road_per_gene[[gg]][,3]

  idx2 = c()
  for(cc in 1:length(starts)){
    idx2 = c(idx2, which(base$BP > starts[cc] & base$BP < ends[cc]))
  }
  snp_locs = as.numeric(names(table(idx2)))
  summ[snp_locs] = summ[snp_locs]+1
}

newdf = cbind.data.frame(base[,1:4], summ)
colnames(newdf) = c(colnames(base)[1:4], "AN")
rm(summ)

if(!dir.exists(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "Roadmap_Blood"))){
  dir.create(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "Roadmap_Blood"))
}

write.table(newdf,
            file = gzfile(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/",
                                 "Roadmap_Blood", "/",
                                 "Roadmap_Blood", ".",
                                 numchr, ".annot.gz")),
            quote=FALSE, row.names=FALSE)
cat("We are at chrom:", numchr, "\n")


##########################   Promoter   #########################################################

geneanno = data.table::fread(paste0("/n/groups/price/kushal/ExPecto/resources/geneanno_", geneset, ".csv"))

gene_beds = cbind.data.frame(geneanno$seqnames, geneanno$TSS-5000, geneanno$TSS+5000, geneanno$symbol)
colnames(gene_beds) = c("chr", 'start', 'end', 'gene')
#gene_beds = read.table("/n/groups/price/kushal/Mouse_Humans/data/Gene_5kb.txt", header=T)

common_genes = intersect(gene_beds$gene, geneanno$symbol)
geneanno2 = geneanno[match(common_genes, geneanno$symbol), ]
gene_beds2 = gene_beds[match(common_genes, gene_beds$gene), ]

genes_in_chr = common_genes[which(geneanno2$seqnames == paste0("chr", numchr))]

summ = rep(0, nrow(base))

for(gg in genes_in_chr){

  ###  Get starts and ends for ABC regions for each gene  gg  #######################
  start = gene_beds2$start[which(gene_beds2$gene == gg)]
  end =  gene_beds2$end[which(gene_beds2$gene == gg)]
  snp_locs = which(base$BP > start & base$BP < end)
  summ[snp_locs] = summ[snp_locs]+1
}

newdf = cbind.data.frame(base[,1:4], summ)
colnames(newdf) = c(colnames(base)[1:4], "AN")
rm(summ)

if(!dir.exists(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "Promoter"))){
  dir.create(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "Promoter"))
}

write.table(newdf,
            file = gzfile(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/",
                                 "Promoter", "/",
                                 "Promoter", ".",
                                 numchr, ".annot.gz")),
            quote=FALSE, row.names=FALSE)
cat("We are at chrom:", numchr, "\n")


##########################   5kb   #########################################################

geneanno = data.table::fread(paste0("/n/groups/price/kushal/ExPecto/resources/geneanno_", geneset, ".csv"))
gene_beds = read.table("/n/groups/price/kushal/Enhancer_MasterReg/data/Gene_5kb.txt", header=T)

common_genes = intersect(gene_beds$gene, geneanno$symbol)
geneanno2 = geneanno[match(common_genes, geneanno$symbol), ]
gene_beds2 = gene_beds[match(common_genes, gene_beds$gene), ]

genes_in_chr = common_genes[which(geneanno2$seqnames == paste0("chr", numchr))]

summ = rep(0, nrow(base))

for(gg in genes_in_chr){

  ###  Get starts and ends for ABC regions for each gene  gg  #######################
  start = gene_beds2$start[which(gene_beds2$gene == gg)]
  end =  gene_beds2$end[which(gene_beds2$gene == gg)]
  snp_locs = which(base$BP > start & base$BP < end)
  summ[snp_locs] = summ[snp_locs]+1
}

newdf = cbind.data.frame(base[,1:4], summ)
colnames(newdf) = c(colnames(base)[1:4], "AN")
rm(summ)

if(!dir.exists(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "5kb"))){
  dir.create(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "5kb"))
}

write.table(newdf,
            file = gzfile(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/",
                                 "5kb", "/",
                                 "5kb", ".",
                                 numchr, ".annot.gz")),
            quote=FALSE, row.names=FALSE)
cat("We are at chrom:", numchr, "\n")


##########################   100kb   #########################################################

geneanno = data.table::fread(paste0("/n/groups/price/kushal/ExPecto/resources/geneanno_", geneset, ".csv"))
gene_beds = read.table("/n/groups/price/kushal/Enhancer_MasterReg/data/Gene_100kb.txt", header=T)

common_genes = intersect(gene_beds$gene, geneanno$symbol)
geneanno2 = geneanno[match(common_genes, geneanno$symbol), ]
gene_beds2 = gene_beds[match(common_genes, gene_beds$gene), ]

genes_in_chr = common_genes[which(geneanno2$seqnames == paste0("chr", numchr))]

summ = rep(0, nrow(base))

for(gg in genes_in_chr){

  ###  Get starts and ends for ABC regions for each gene  gg  #######################
  start = gene_beds2$start[which(gene_beds2$gene == gg)]
  end =  gene_beds2$end[which(gene_beds2$gene == gg)]
  snp_locs = which(base$BP > start & base$BP < end)
  summ[snp_locs] = summ[snp_locs]+1
}

newdf = cbind.data.frame(base[,1:4], summ)
colnames(newdf) = c(colnames(base)[1:4], "AN")
rm(summ)

if(!dir.exists(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "100kb"))){
  dir.create(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/", "100kb"))
}

write.table(newdf,
            file = gzfile(paste0("/n/groups/price/kushal/Imperio/data/Imperio_SNP_S2G_MAP_", geneset, "/",
                                 "100kb", "/",
                                 "100kb", ".",
                                 numchr, ".annot.gz")),
            quote=FALSE, row.names=FALSE)
cat("We are at chrom:", numchr, "\n")



