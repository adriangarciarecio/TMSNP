#From uniprot download uniprot_list.txt with uniprot_list accessopn code gene name and subcellular location transmem
tail -n +2 uniprot_list.txt  > 1
tail -n +2 uniprot_list.txt  > 7
sed $'s/\t/GENE/' 1 > 2
grep -o 'TRANSMEM[^Helical]\+\|GENE' 2 > 3
tail -n +2 3 > 4
tr '\n' ' ' < 4 >5
sed $'s/GENE/\\\n/g' 5 > 6
awk '{print $1 "\t" $2}' 7 > 8
paste 8 6 > uniprot_list_tm.txt
rm 1 2 3 4 5 6 7 8
#clinvar_result.txt conte llista amb missense likely patho and patho human de Clin Var download text format
cp ClinVar.txt result_clinvar.txt
fgrep "(p"  result_clinvar.txt>clin
fgrep ">" clin>clin2
sed '/delins/d' clin2>clin
sed '/?/d' clin>clin2
sed '/=/d' clin2>clin
sed $'/\[/d' clin>clin2
sed 's/|/ /g' clin2>clin
fgrep "Pathogenic" clin> clinvarpatho
fgrep "Likely" clin > clinvarlp
awk '{print $2}' clinvarpatho>clinvarix
awk '{print $3}' clinvarpatho>clinvargene
sed 's/\.//g' clinvarix>clinvar2
sed 's/)//g' clinvar2>clinvar3
sed 's/(p//g' clinvar3>clinvar4
sed 's/(//g' clinvar4>pp
mv pp clinvar4
grep -Eo '[0-9]{1,9}' clinvar4 > clinvar5
cut -c1-3 clinvar4 >clinvar4x
cut -c4-11 clinvar4>clinvar10
grep -Eo '[abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPRSTUVXWYZ]{1,18}' clinvar10 > clinvar11
grep -Eo '[abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPRSTUVXWYZ]{1,3}' clinvar4x > clinvar6
paste clinvar5 clinvar6 clinvar11 clinvargene  >patho_clinvar.txt
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "1" }' patho_clinvar.txt >patho_clinvar2.txt
rm patho_clinvar.txt
awk '{print $2}' clinvarlp>clinvarix
awk '{print $3}' clinvarlp>clinvargene
sed 's/\.//g' clinvarix>clinvar2
sed 's/)//g' clinvar2>clinvar3
sed 's/(p//g' clinvar3>clinvar4
sed 's/(//g' clinvar4>pp
mv pp clinvar4
grep -Eo '[0-9]{1,9}' clinvar4 > clinvar5
cut -c1-3 clinvar4 >clinvar4x
cut -c4-11 clinvar4>clinvar10
grep -Eo '[abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPRSTUVXWYZ]{1,18}' clinvar10 > clinvar11
grep -Eo '[abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPRSTUVXWYZ]{1,3}' clinvar4x > clinvar6
paste clinvar5 clinvar6 clinvar11 clinvargene >lpatho_clinvar.txt
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "2" }' lpatho_clinvar.txt >lpatho_clinvar2.txt
rm lpatho_clinvar.txt
cat patho_clinvar2.txt lpatho_clinvar2.txt > patho_clinvar.txt
rm clin*txt
rm patho_clinvar2.txt lpatho_clinvar2.txt
mv result_clinvar.txt clinvar_result.txt
python ./clinvar_clean.py > clinvar_clean.txt
python ./clinvar_reclean.py > clinvar_output.txt
rm clinvar_clean.txt clinvar_result.txt
