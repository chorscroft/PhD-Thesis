
for i in {1..100}
do

	## Create the folder
	mkdir sim_$i

	## copy in the SLiM code
	cp selectedsim.txt ./sim_$i/

	## copy in the qsub code
	cp runSim.sh ./sim_$i/

	## run the qsub code
	cd ./sim_$i
	qsub runSim.sh
	cd ..

done
