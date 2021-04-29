for i in {5..20}
do
	cp baseSimulation.txt simGene$i.txt
	sed -i "s/GENESPERMB/$i/g" simGene$i.txt
	sed -i "s/OUTPUTFILENAME/simGene$i/g" simGene$i.txt
done
