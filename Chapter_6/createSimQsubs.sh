for i in {5..20}
do
	cp baseSimQsub.sh simQsubGene$i.sh
	sed -i "s/SIMTORUN/simGene$i/g" simQsubGene$i.sh
	qsub simQsubGene$i.sh
done
