python patho_ensembl.py
python gnomad_download.py

# Continuar aqui
uniq gnomad_mutations.txt > unique_mutations
sed 's/\.//' unique_mutations >o
sed 's/unique_mutations//g' o>l
sed '/Ter/d' l > pp
sed '/_/d' pp >l
sed '/Met1/d' l >pp
sed '/del/d' pp >l
sed '/0\.0/d' l>pp
mv pp l
awk '{print $1}' l >a
awk '{print $2}' l >b
awk '{print $3}' l >c
grep -Eo '[0-9]{1,9}' b > o
cut -c1-3 b >r 
cut -c4-18 b>i
grep -Eo '[abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPRSTUVXWYZ]{1,18}'  i>j 
paste a o r j c > gnomad2.txt
sed -i.bak $'s/\t/,/g' gnomad2.txt
awk -F, '!/^$/ {if($3==$4) {} else {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}}' gnomad2.txt >gnomad.txt
rm  o  l a b r i j unique_mutations  gnomad2.txt
uniq gnomad.txt>gnomad2.txt
mv gnomad2.txt gnomad.txt