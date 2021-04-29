for i in {1..38}
do

	cd /temp/hgig/EXOME_DATA/Clare/Dogs/cluster2

	## make a directory
	mkdir /temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/zalpha_all/chr$i

	## copy zalpha .R and .sh code to folder
	cp runRcode.sh /temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/zalpha_all/chr$i/
	cp runZalpha.R /temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/zalpha_all/chr$i/

	##Run code
	cd /temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/zalpha_all/chr$i
	qsub /temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/zalpha_all/chr$i/runRcode.sh -v j=$i

done
