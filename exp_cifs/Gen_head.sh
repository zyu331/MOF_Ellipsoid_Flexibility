for i in *.cif
do
sed "s/(//g" $i > temp
sed "s/)//g" temp > $i
rm temp

j=`printf "%s" $i | cut -c -6`
mkdir $j

# cell file

grep "_cell_length_a" $i > $j/${i}cellPara
#echo -e "" >> $j/${i}cellPara
grep "_cell_length_b" $i >> $j/${i}cellPara
#echo -e "" >> $j/${i}cellPara
grep "_cell_length_c" $i >> $j/${i}cellPara
#echo -e "" >> $j/${i}cellPara
grep "_cell_angle_a" $i >> $j/${i}cellPara
grep "_cell_angle_b" $i >> $j/${i}cellPara
grep "_cell_angle_g" $i >> $j/${i}cellPara

# sys file
line_start="$(grep -n "_symmetry_equiv_pos_as_xyz" $i | head -n 1 | cut -d: -f1)"
line_end="$(grep -n "_cell_length_a" $i | head -n 1 | cut -d: -f1)" 

lstart=$(($line_start-1))
lend=$(($line_end-1))

sed -n ''$lstart','$lend'p' $i   > ${i}sysm


# pos file 
data_s="$(grep -n "_atom_site_fract_z" $i | head -n 1 | cut -d: -f1)"
data_e="$(grep -n "_atom_site_aniso_label" $i | head -n 1 | cut -d: -f1)"

# pos head
s_head_T=$(($data_s-8))
e_head_T=$(($data_s+15))
sed -n ''$s_head_T','$e_head_T'p' $i > $j/temp
s_head="$(grep -n "loop_" $j/temp | head -n 1 | cut -d: -f1)"
e_head="$(grep -n "_" $j/temp | tail -n 1 | cut -d: -f1)"
sed -n ''$s_head','$e_head'p' $j/temp > $j/${i}poshead

# head
cat $j/${i}cellPara $j/${i}sysm $j/${i}poshead > $j/${i}head

# pos num
s_num=$(($e_head+$s_head_T))
echo $e_head
echo $data_s
echo $s_num
e_num=$(($data_e-2))
echo $e_num
sed -n ''$s_num','$e_num'p' $i | awk '{print $3,$4,$5}' | cut -d: -f3  > ${i}pos
sed -n ''$s_num','$e_num'p' $i | awk '{print $1,$2}' | cut -d: -f3  > ${i}prePos
sed -n ''$s_num','$e_num'p' $i | awk '{print $6,$7,$8,$9,$10,$11,$12,$13}' | cut -d: -f3  > ${i}afterPos
sed -n ''$s_num','$e_num'p' $i | awk '{print $1,$6}' | cut -d: -f3  > ${i}pattern_T
sed "s,',,g" ${i}pattern_T > ${i}pattern

# atomname
sed -n ''$s_num','$e_num'p' $i | awk '{print $1}' > ${i}label


done


