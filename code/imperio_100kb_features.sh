module load gcc/6.2.0
module load R

IFS="
"

method=DeepSEA #Basenji
outpath=../data/Imperio_models                                                                                                                       

for chrom in {1..22}
do
cmd="Rscript imperio_100kb_features.R $chrom $method $outpath"
sbatch --time=80:00 --mem=40000 --output=100kb.out --error=100kb.err -p short -c 1 --wrap="$cmd"
done
