=====================================================
Steps for converting gene set (GS) to S2G annotations
=====================================================

- First step:
  Start with a gene score file. This is a standard .txt file (no column name, no row name), with gene name in the first column, and
  the score in the second column (it does not automatically work with Ensembl and Entrez ids, so please convert to gene names for the 
  Ensembl ids)
  
 - See an example gene score file 
- INPUT_FOLDER=~/Imperio/code; geneset_dir=$INPUT_FOLDER/Test_GeneSets;  head -n 5 $INPUT_FOLDER/Test_GeneSets/geneset_example.txt`

RITA1						   0
FAM225B						   0.00320672144802546
TBC1D23						   0.0716823315695179
SIRPG						   0
L3MBTL4-AS1					   0

- Next create a folder where you want to save the bed (graph) format files created from this gene score file
- bed_dir=$INPUT_FOLDER/Test_Beds`

- Run the first script to save the bed(graph) files for ABC, Roadmap, and ABC U Roadmap S2G strategies for the given gene set. If you have multiple
  gene set files in the same geneset_dir, then also the following code should work
- bash geneset_to_bed.sh $geneset_dir $bed_dir 'BLD'

This script will create a number of folders as number of gene sets in $bed_dir with bed files corresponding to different S2G strategies (not part of baselineLD model)
- In this case, it will look like the following
- Test_Beds/
  - geneset_example/
    -100kb.bed
    -5kb.bed
    -ABC.bed
    -FinemapBloodeQTL.bed
    -PCHiC.bed
    -Roadmap_Enhancer.bed
    -Yoshida.bed
 
- Now take for the gene set you are interested (in case of multiple gene sets), run the following cleaning operation that uses bedops to clean the 
  bedgraph files from previous step.
- geneset=geneset_example
- bash clean_bed.sh $bed_dir/$geneset

As a result of running this script, the bed files in the previous step will get modified, no extra files will be created.

- Now, you are all set to generate annotations from the cleaned bed (graph) files. First, create  folder where you want to save your annotations 
  data.
- annot_path=$INPUT_FOLDER/Test_Annots

- Then run the following code 
- bimfile_path=~/Imperio/data/BIMS
- bash bed_to_annot.sh $bed_dir $bimfile_path $annot_path $geneset

- This script will create annotation files of the following form
- Test_Annots/
  - geneset_example
    -100kb
    -5kb
    -ABC
    -FinemapBloodeQTL
    -PCHiC
    -Roadmap_Enhancer
    -Yoshida
    
Each of the sub-directories contains the S2G annotations for the gene set linked by some S2G strategy.

The full pipeline:


INPUT_FOLDER=~/Documents/Imperio/code
geneset_dir=$INPUT_FOLDER/Test_GeneSets
bed_dir=$INPUT_FOLDER/Test_Beds

bash geneset_to_bed.sh $geneset_dir $bed_dir

geneset=geneset_example
bash clean_bed.sh $bed_dir/$geneset

bimfile_path=~/Documents/Imperio/data/BIMS
annot_path=$INPUT_FOLDER/Test_Annots
bash bed_to_annot.sh $bed_dir $bimfile_path $annot_path $geneset


 
