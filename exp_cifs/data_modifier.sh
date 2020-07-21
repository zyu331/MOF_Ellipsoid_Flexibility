for i in *.cif
do

# sys file
line_start="$(grep -n "_symmetry_equiv_pos_as_xyz" $i | head -n 1 | cut -d: -f1)"
line_end="$(grep -n "_cell_length_a" $i | head -n 1 | cut -d: -f1)" 

lstart=$(($line_start+1))
lend=$(($line_end-1))

sed -n ''$lstart','$lend'p' $i   > temp
sed -e "s/'//" temp > temp2
sed -e "s/'//" temp2 > ${i}sym
rm temp temp2


# pos file
line_s="$(grep -n "_atom_site_disorder_group" $i | head -n 1 | cut -d: -f1)"
line_e="$(grep -n "_atom_site_aniso_label" $i | head -n 1 | cut -d: -f1)"

l_s=$(($line_s+1))
l_e=$(($line_e-3))

sed -n ''$l_s','$l_e'p' $i | awk '{print substr($3,0,6),substr($4,0,6),substr($5,0,6)}' | cut -d: -f3  > ${i}pos
sed -n ''$l_s','$l_e'p' $i | awk '{print substr($8,0,6)' | cut -d: -f3  > ${i}occupancy
sed -n ''$l_s','$l_e'p' $i | awk '{print $2}' | cut -d: -f3  > ${i}atomType

# atomname
sed -n ''$l_s','$l_e'p' $i | awk '{print $1}' > ${i}name

# head file

grep "_cell_length_a" $i | awk '{print $2}'| head -c 7 > ${i}cellPara
echo -e "" >> ${i}cellPara
grep "_cell_length_b" $i | awk '{print $2}'| head -c 7 >> ${i}cellPara
echo -e "" >> ${i}cellPara
grep "_cell_length_c" $i | awk '{print $2}'| head -c 7 >> ${i}cellPara
echo -e "" >> ${i}cellPara
grep "_cell_angle_a" $i | awk '{print $2}'| head -c 7 >> ${i}cellPara
grep "_cell_angle_b" $i | awk '{print $2}'| head -c 7 >> ${i}cellPara
grep "_cell_angle_g" $i | awk '{print $2}'| head -c 7 >> ${i}cellPara



done


