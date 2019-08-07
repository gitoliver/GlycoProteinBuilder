

for file in `ls Rotamer*.pdb`
do
    echo "$file"
    grep "NLN" $file > $file.glycan.pdb
    grep -v "fat" $file | tail -978 >> $file.glycan.pdb
done
