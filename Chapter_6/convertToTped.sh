module load vcftools/0.1.14
module load plink/1.90beta

for i in {5..20}
do
	plink --vcf simGene$i.vcf --maf 0.01 --hwe 0.001 --recode 12 transpose --out simGene$i
done
