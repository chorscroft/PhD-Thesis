module load plink/1.90beta
for i in {1..38}
do
	mkdir chr$i
	cd chr$i
	plink --freq --dog --tfile /temp/hgig/EXOME_DATA/Clare/Dogs/zalpha/DataForZalpha/chr$i --recode
	cd ..

done
