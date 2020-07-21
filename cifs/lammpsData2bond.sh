for i in *
do
echo $i
l_start="$(awk '/Bonds/{ print NR; exit }' Data_lammps/data.$i)"
l_end="$(awk '/Angles/{ print NR; exit }' Data_lammps/data.$i)"

lstart=$(($l_start+2))
lend=$(($l_end-2))

sed -n ''$lstart','$lend'p' Data_lammps/data.$i > $i.bond

done

rm lammpsData2bond.sh.bond
