outpath = "../data/Imperio_models"
xgboost_path = "../data/Xgboost_saved_models"


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

gene_exp = data.frame(fread('../data/geneanno.exp.csv'))
pseudocount = 0.0001
YY = log(gene_exp['Whole_Blood'][match(common_names, geneanno$symbol), , drop=F]+pseudocount)

chrnames = geneanno$seqnames[match(rownames(XX), geneanno$symbol)]

test_ids = which(chrnames == paste0("chr", 8))
train_ids = setdiff(1:length(chrnames), test_ids)

YY_train = YY[train_ids, ]
YY_test = YY[test_ids, ]

XX_train = XX[train_ids, ]
XX_test = XX[test_ids, ]

library(xgboost)
bstSparse <- xgboost(data = XX_train,
                     label = YY_train,
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

eval = predict(bstSparse, XX_test)
temp = cbind(eval, YY_test)
cor(temp, method="spearman")

xgb.save(bstSparse, paste0(xgboost_path, "/", "imperio_model.save"))
