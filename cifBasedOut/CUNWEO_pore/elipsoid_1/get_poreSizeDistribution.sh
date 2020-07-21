
	i=0
	while [ $i -lt 10 ]; do
		network -ha -psd 1.2 1.2 50000 $i.cif
		i=$((i+1))
	done



