# TMSNP

This is a collection of scripts to set up or update the TMSNP database. They should be executed in the order listed here.

- Get SNPs from the Uniprot and also the TM helix annotations (mainly Iker’s stuff)
  - ./main/uniprot_tables.py (runs in seconds)

- Get SNPs from Clinvar and GNOMAD 
  - Run Mireia's scripts which will create the inputs for **clinvar2tmsnp.py** and **gnomad2tmsnp.py**
  - clinvar/clinvar2tmsnp.py (runs in seconds)
  - gnomad/gnomad2tmsnp.py (runs in minutes)

- Check and remove duplicates after combining Uniprot ClinVar and GNOMAD (we keep Uniprot > ClinVar > GNOMAD)
  - ./main/duplicates_3db.py (runs in seconds)
  
- Download the PFAM alignments
  - ./main/download_pfam_aln.py (hours) 
      - ./main/from_full_aln.py (to complete the alignments that failed in the previous step; runs in hours)

- Evaluate entropy, substitution…
  - ./main/evaluate_snps.py (hours)

- A second round of remove duplicates (whith duplicate PFAM records ...)
  - ./main/duplicates_pfam.py
  
- Add Uniprot annotations relevant to positions
  - ./uniprot_annotations/add_anotations.py

- Get SIFT and POLYPHEN2 scores
  - ./sift_polyphen/send_sift_polyphen.py
  - Use the created file as input for the two web servers
  - The output should be saved and processed with the following scipts.
  - ./sift_polyphen/addtotable_sift.py
  - ./sift_polyphen/addtotable_polyphen.py