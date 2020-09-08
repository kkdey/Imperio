module load gcc/6.2.0
module load R

IFS="
"
method=DeepSEA #Basenji
outpath=../data/Imperio_models

for chrom in {1..22}
do
cmd="Rscript imperio_roadmap_basenji_blood_prep.R $chrom $method $outpath"
sbatch --time=120:00 --mem=80000 --output=roadmap.out --error=roadmap.err -p short -c 1 --wrap="$cmd"
done
