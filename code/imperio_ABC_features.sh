module load gcc/6.2.0
module load R

IFS="
"
method=DeepSEA #Basenji
outpath=../data/Imperio_models

for chrom in {1..22}
do
cmd="Rscript imperio_abc_deepsea_blood_prep.R $chrom $method $outpath"
sbatch --time=120:00 --mem=80000 --output=abc.out --error=abc.err -p short -c 1 --wrap="$cmd"
done

