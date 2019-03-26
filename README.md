# TMSNP

This is a collection of scripts to set up or update the TMSNP database. They should be executed in the order listed here.

- Get SNPs from the Uniprot and also the TM helix annotations (mainly Iker’s stuff). This starts with the list of all membrane proteins in Uniprot at each time the script runs. Creates the MySQL tables receptor_pfam (contains Uniprot ID - multiple PFAM ID relations), tm_segments (start end of each TM segment and linked to Uniprot ID), and snp (all info).
  - ./main/uniprot_tables.py (runs in seconds)
  - ./main/pfam_update.py (runs in hours)

- Get SNPs from Clinvar (all programs run in seconds)
  - Prepare uniprot_list.txt and ClinVar.txt (both files from web pages)
  - ./clinvar/ClinVar.sh
  - ./clinvar/clinvar2tmsnp.py
  
- Remove proteins swith non-pathogenic SNPs
  -  ./main/nonpathogenic_proteins.py (runs in seconds)

- Get SNPs from Clinvar and GNOMAD 
  - Generate ensembl_list.txt (*PENDING*)
  - ./gnomad/Gnomad.sh (runs in a few hours)
  - ./gnomad/gnomad2tmsnp.py (runs in minutes)

- Check and remove duplicates after combining Uniprot ClinVar and GNOMAD (we keep Uniprot > ClinVar > GNOMAD)
  - ./main/duplicates_3db.py (runs in seconds)

- Complete snps table with id and gene entries (needs file ensembl_list.txt)
  - ./main/get_id.py (runs in minutes)
  
- Download the PFAM alignments
  - ./main/download_pfam_aln.py (runs in hours) 
  - ./main/pfam_download.py (use this better!) --> returns a list of no downloaded pfams. Use the list in the next script as miss_pfam
      - ./main/from_full_aln.py (to complete the alignments that failed in the previous step; runs in hours)

- Evaluate entropy, substitution…
  - ./main/evaluate_snps.py (runs in hours)

      THE NEXT STEP IS NOT NECESSARY!!!!
- A second round of remove duplicates (whith duplicate PFAM records ...)
  - ./main/duplicates_pfam.py 

- Add Uniprot annotations relevant to positions
  - ./uniprot_annotations/Annotations.sh (runs in seconds)
  - ./uniprot_annotations/add_anotations.py (runs in minutes)

- Get SIFT and POLYPHEN2 scores
  - ./sift_polyphen/send_sift_polyphen.py 
  - Use the created file as input for the two web servers
  - The output should be saved and processed with the following scipts.
  - ./sift_polyphen/addtotable_sift.py
  - ./sift_polyphen/addtotable_polyphen.py
  - 
- Export table to .csv file
  - $ mysql -u lmcdb -h alf03.uab.cat -p$LMCDB_PASS tmsnp -B -e "select * from snp_eval;" > tmsnps_12-2018.csv