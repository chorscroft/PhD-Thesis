for i in {5..20}
do
	k=0

	while read f; do
		k=$(($k+1))
		newName="snp"$k
		echo ${f/./$newName} >> simNewGene$i.tped
	done < simGene$i.tped
done

