for i in {5..20}
do
	cp Xcp "cp$i"
	sed -i "s/X/$i/g" "cp$i"
done

