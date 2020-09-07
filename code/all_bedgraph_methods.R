
PCHiC_bedgraph_calc <- function(scores,
                                output_cell,
                                output_bed = "temp.bed"){
  library(data.table)
  pchic_data = data.frame(fread("../data/PCHiC_peak_matrix_cutoff5.tsv", header=T))

  names1 = lapply(pchic_data$baitName, function(x) return(strsplit(as.character(x), ";")[[1]]))
  num_genes_per_connect = unlist(lapply(names1, function(x) return(length(x))))

  genes_all_connect = unlist(names1)
  indices_per_connect = unlist(lapply(1:nrow(pchic_data), function(x) return(rep(x, num_genes_per_connect[x]))))

  common_genes = intersect(unique(genes_all_connect), names(scores))
  idx = which(!is.na(match(genes_all_connect, common_genes)))

  new_grades = scores[match(genes_all_connect[idx], names(scores))]
  df3 = cbind(pchic_data[indices_per_connect[idx], c("oeChr", "oeStart", "oeEnd")], new_grades)
  df3[,1] = paste0("chr", df3[,1])

  write.table(df3, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

}


ABC_intergenic_bedgraph_calc <- function(scores,
                                         output_cell,
                                         output_bed = "temp.bed"){

  df_pre = data.frame(fread(paste0("../data/",
                                   "ABCpaper_NasserFulcoEngreitz2020_Blood_AvgHiC.txt.gz")))
  df = df_pre[unique(c(grep("intergenic", df_pre$name))),]
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
  colnames(df2) = c("chr", "start", "end", "TargetGene")
  common_genes = intersect(names(scores), df2[,4])
  idx = which(!is.na(match(df2[,4], common_genes)))
  df3 = df2[idx,]
  grades = scores[match(df3[,4], names(scores))]
  temp = cbind(df3, grades)
  final_bed = temp[, c(1:3, 5)]

  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

Yoshida_bedgraph_calc <- function(scores,
                                  output_cell,
                                  output_bed = "temp.bed"){

  df = data.frame(fread(paste0("../data",
                               "/", "Yoshida_correlated_summits_per_gene_inside100KB_sort.bed")))
  unique_genes = as.character(unique(df[,4]))
  gene_orthologs = read.table("../data/Orthologs_Yoshida_Mouse_Human.txt", header=T)

  temp = gene_orthologs[match(df[,4], gene_orthologs[,1]), 2]
  df2 = cbind.data.frame(df[which(!is.na(temp)), 1:3], temp[which(!is.na(temp))])
  colnames(df2) = c("chr", "start", "end", "gene")

  common_genes = intersect(names(scores), df2[,4])
  idx = which(!is.na(match(df2[,4], common_genes)))
  df3 = df2[idx,]
  grades = scores[match(df3$gene, names(scores))]
  bed1 = cbind(df3, grades)
  bed2 = bed1[,c(1:3, 5)]
  write.table(bed2, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

eQTLCPP_bedgraph_calc = function(scores,
                                 output_cell,
                                 output_bed = "temp.bed"){
  options(scipen = 999)

  Enhancer = read.table("../data/eQTL_Blood_CPP.txt",
                        header=F)
  matched_ids = match(Enhancer[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  df3 = cbind(Enhancer[,1:3], temp)

  write.table(df3, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

Roadmap_Enhancer_bedgraph_calc = function(scores,
                                          output_cell,
                                          output_bed = "temp.bed"){
  options(scipen = 999)

  Enhancer = read.table("../data/Roadmap_Enhancers_Blood.txt",
                        header=F)
  matched_ids = match(Enhancer[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  df3 = cbind(Enhancer[,1:3], temp)

  write.table(df3, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


