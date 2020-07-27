# TMSNP

TMSNP is a web application tool that was constructed using a Python backend (v.3.7) with the Flask framework  (v.1.0.2). Both the application and the associated datasets used for training and testing the predictor were built automatically using Python/Bash scripts that collected the required data and stored it in a MySQL database (v.8.0.18), facilitating regular updates.

It is divided into two parts: 
- TMSNP Database 
- TMSNP Predictor (WebApp)

## TMSNP DATABASE

It is necessary to install the next: 

- Python v3.7
- SQLAlchemy (or another mysql connector package, this package is useful because is a more directly connection in one line) 
- Pandas 
- Numpy  

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
  - ./gnomad/Gnomad.sh (runs in a few hours) *Note: If you want to get all annotations (Nonpathogenic and pathogenic) use ensembl_list_all.txt on gnomad_download.py insted of ensembl_list_model.txt 
  - ./gnomad/gnomad2tmsnp.py 0 or 1(runs in minutes/ hours) *Note: select 0 if you want evaluate_snps for training the model or 0 if you want evaluate_snps for all annotations.

- Check and remove duplicates after combining Uniprot ClinVar and GNOMAD (we keep Uniprot > ClinVar > GNOMAD)
  - ./main/duplicates_3db.py 0 or 1 (runs in seconds) *Note: select 0 if you want evaluate_snps for training the model or 0 if you want evaluate_snps for all 

- Complete snps table with id and gene entries (needs file ensembl_list.txt)
  - ./main/get_id.py 0 or 1 (runs in few minutes/ minutes) *Note: select 0 if you want evaluate_snps for training the model or 0 if you want evaluate_snps for all 
  
- Download the PFAM alignments
  - ./main/download_pfam_aln.py (runs in hours) 
  - ./main/pfam_download.py (use this better!) --> returns a list of no downloaded pfams. Use the list in the next script as miss_pfam
      - ./main/from_full_aln.py (to complete the alignments that failed in the previous step; runs in hours)

- Evaluate entropy, substitution…
  - ./main/evaluate_snps.py 0 or 1 (runs in hours) *Note: select 0 if you want evaluate_snps for training the model or 0 if you want evaluate_snps for all annotations. 
  - ./main/remove_acc.py (ONLY IN CASE IF YOU SELECT THE MODE 0)
      THE NEXT STEP IS NOT NECESSARY!!!!
- A second round of remove duplicates (whith duplicate PFAM records ...)
  - ./main/duplicates_pfam.py 

- Add Uniprot annotations relevant to positions
  - ./uniprot_annotations/Annotations.sh (runs in seconds)
  - ./uniprot_annotations/add_anotations.py (runs in minutes) *Note: select 0 if you want evaluate_snps for training the model or 0 if you want evaluate_snps for all annotations. 

- Get SIFT and POLYPHEN2 scores
  - ./sift_polyphen/send_sift_polyphen.py 
  - Use the created file as input for the two web servers
  - The output should be saved and processed with the following scipts.
  - ./sift_polyphen/addtotable_sift.py
  - ./sift_polyphen/addtotable_polyphen.py
  - 
- Export table to .csv file
  - $ mysql -u lmcdb -h alf03.uab.cat -p$LMCDB_PASS tmsnp -B -e "select * from snp_eval;" > tmsnps_12-2018.csv
  

## TMSNP PREDICTOR

Machine-learning models were built using Flame (https://github.com/phi-grib/flame; a Python modeling framework which wraps scikit-learn (http://scikit-learn.sourceforge.net)). Unfortunately, the size of the files that contain the models are higher than the maximum upload size that GitHub give.

We have also created a Virtual Environment (webapp/venv_tmsnp) used by the webapp to perform the predictions using the script webapp/predictor.py. 
