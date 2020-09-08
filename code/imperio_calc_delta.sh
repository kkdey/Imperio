module load gcc/6.2.0
module load R

IFS="
"

method=DeepSEA #Basenji
outname=deep_imperios2g_annot

for chrom in {1..22}
do
cmd="Rscript imperio_calc_delta.R $chrom $method $outname"
sbatch --time=400:00 --mem=160000 --output=imperio.out --error=imperio.err -p short -c 1 --wrap="$cmd"
done

