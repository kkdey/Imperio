# Imperio
Repository for performing DeepBoost, S2G linking and Imperio models linking deep learning models with complex traits

## Citation

If you use the data or code from this repository, please cite our paper 

Dey, K.K. et al bioRxiv. 2020. Integrative approaches to improve the informativeness of deep learning models for human complex diseases.

If you use the ABC S2G strategies please cite Nasser, J., Engreitz, J. et al. unpublished data. 2020. and [Fulco et al, 2019](https://www.nature.com/articles/s41588-019-0538-0). If you use PC-HiC S2G strategies, please cite [Javierre et al 2016 Cell](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123897/). If you use ATAC (Yoshida), cite [Yoshida et al 2019 Cell](https://www.cell.com/cell/pdf/S0092-8674(18)31650-7.pdf). If you use eQTL S2G strategies, cite [Hormozdiari et al 2018 NG](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6030458/). If you use Roadmap S2G links, cite the references [here](https://ernstlab.biolchem.ucla.edu/roadmaplinking/). 


## Description

The repository provides a set of tools to integrate functional information on top of trained genomic deep learning models. The codes can be classified into following categories 

- **DeepBoost: Boosted deep learning annotations using functionally finemapped SNPs + S2G strategies + gene sets**
      - *code/deepboost.R* : Run classification model to separate fine-mapped SNPs from control using deep learning annotations on
         even (odd) chromosomes and predict labels on SNPs from odd (even) chromosomes.
      - *code/deepboost_s2g.R*: Combine boosted deep learning annotations with S2G strategies 
      - *code/geneset_to_bed.R* : Create S2G bedgraph files from a gene set file (check `code/README_GSSSG.txt`)
      - *code/clean_bed.sh*: Postprocess S2G bedgraph files from a gene set from previous step.
      - *code/bedgraph_to_annot.py*: Convert S2G bedgraph files for a geneset to SNP annotations.

      check `code/README_GSSSG.txt` for an illustration of how to generate SNP level annotations starting from a gene set file.

- **Imperio: Predicted blood expression using deep learning models combined with S2G strategies**
      - *code/imperio_100kb_features.R*: Build Imperio per chr features for 100kb S2G strategy
      - *code/imperio_5kb_features.R*: Build Imperio per chr features for 5kb S2G strategy
      - *code/imperio_TSS_features.R*: Build Imperio per chr features for TSS S2G strategy
      - *code/imperio_ABC_features.R*: Build Imperio per chr features for ABC(blood) S2G strategy
      - *code/imperio_Roadmap_features.R*: Build Imperio per chr features for Roadmap(blood) S2G strategy
      - *code/imperio_aggregate_s2g.R*: Aggregate features across chromosomes for each S2G strategy
      - *code/imperio_xgboost_train.R*: Run Imperio regression model 
      - *code/imperio_s2g_mapping.R*: Perform mapping of SNPs to genes using different S2G strategies
      - *code/imperio_snp_table_merge.R*: Merge S2G strategies as a step to build SNP level annotation.
      - *code/imperio_calc_delta.R*: Generate Imperio SNP level annotation per S2G strategy 
      - *code/imperio_MAF_annots.R*: Generate MAF corrected Imperio SNP level annotation
      
## Annotations

All DeepBoost, restricted DeepBoost (by S2G), restricted DeepBoost (by S2G and gene set), Imperio, Imperio-geneset annotations 
are available [here]()

## How to use these annotations?

1) Download the LDSC from git (https://github.com/bulik/ldsc/wiki/Partitioned-Heritability)
2) Download the baselineLD_v2.1 annotations from Broad webpage (https://data.broadinstitute.org/alkesgroup/LDSCORE/)
3) Download deep learning predictive annotations (see above) from https://data.broadinstitute.org/alkesgroup/LDSCORE/Dey_DeepLearning
4) Use your GWAS summary statistics formatted in LDSC details is available (https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format)
5) Download the baseline frq file and weights for 1000G available (https://data.broadinstitute.org/alkesgroup/LDSCORE/)
6) Run S-LDSC with these annotations conditional on baselineLD_v2.1 (see https://github.com/bulik/ldsc/)

```
ANNOT FILE header (*.annot):

CHR -- chromosome
BP -- physical position (base pairs)
SNP -- SNP identifier (rs number)
CM -- genetic position (centimorgans)
all additional columns -- Annotations
```

NOTE: Although one would expect the genetic position to be non-negative for all 1000G SNPs, we have checked that
in fact the genetic position is negative for a handful of 1000G SNPs that have a physical position that is smaller
than the smallest physical position in the genetic map. The genetic positions were obtained by running PLINK on
the Oxford genetic map (http://www.shapeit.fr/files/genetic_map_b37.tar.gz).

MORE DETAIL CAN BE OBTAINED FROM https://github.com/bulik/ldsc/wiki/LD-File-Formats


```
LD SCORE FILE header (*.l2.ldscore):

CHR -- chromosome
BP -- physical position (base pairs)
SNP -- SNP identifier (rs number)
all additional columns -- LD Scores

```

## Contact 

In case of any questions, please open an issue or send an email to me at `kdey@hsph.harvard.edu`.





      
      
      
      







