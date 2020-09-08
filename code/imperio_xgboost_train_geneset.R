##################  Gene sets (top 10%)  - pLI, PPI-enhancer, PPI-master  ##############

library(xgboost)
library(data.table)
outpath = "../data/Imperio_models"
xgboost_path = "../data/Xgboost_saved_models"

#geneset="pLI"

if(geneset == "pLI"){
  gene_set = read.table("../data/pLI_genes2.txt")
  imp_genes = as.character(gene_set[order(gene_set[,2], decreasing = T)[1:2200], 1])
}else if (geneset == "PPI_enhancer"){
  gene_set = read.table("../data/PPI_enhancer.txt")
  imp_genes = gene_set[,1]
}

score1 = get(load(paste0(outpath, "/5kb_full.rda")))
score2 = get(load(paste0(outpath, "/100kb_full.rda")))
score3 = get(load(paste0(outpath, "/TSS_full.rda")))
score4 = get(load(paste0(outpath, "/ABC9_full.rda")))
score5 = get(load(paste0(outpath, "/Roadmap_full.rda")))

geneanno = data.table::fread("../data/geneanno.csv")

common_names = Reduce(intersect, list(rownames(score1),rownames(score2), rownames(score3),
                                      rownames(score4), rownames(score5)))

XX = cbind(score1[match(common_names, rownames(score1)), ],
           score2[match(common_names, rownames(score2)), ],
           score3[match(common_names, rownames(score3)), ],
           score4[match(common_names, rownames(score4)), ],
           score5[match(common_names, rownames(score5)), ])

imp_genes_common = intersect(imp_genes, rownames(XX))

XX2 = XX[match(imp_genes_common, rownames(XX)), ]

gene_exp = data.frame(fread('../data/geneanno.exp.csv'))
pseudocount = 0.0001
YY = log(gene_exp['Whole_Blood'][match(imp_genes_common, geneanno$symbol), , drop=F]+pseudocount)

chrnames = geneanno$seqnames[match(rownames(XX2), geneanno$symbol)]


library(xgboost)
bstSparse <- xgboost(data = XX2,
                     label = YY[,1],
                     learning_rate = 0.05,
                     max_depth = c(5, 10, 20),
                     eta = 0.01,
                     nthread = 16,
                     gamma = 10,
                     min_child_weight = c(6, 8, 10),
                     nrounds = 200,
                     objective = "reg:linear",
                     base_score = 2,
                     subsample = c(0.6, 0.8, 1),
                     booster = 'gblinear',
                     early_stopping_rounds=10)

xgb.save(bstSparse, paste0(xgboost_path, "/", "imperio_model", "_", geneset, ".save"))
