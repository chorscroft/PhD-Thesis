for i in {5..20}
do
	cp baseRunLDMAP.sh runLDMAP$i.sh
	sed -i "s/cpX/cp$i/g" runLDMAP$i.sh
	qsub runLDMAP$i.sh
done

