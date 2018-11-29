awk '{for(i=1;i<=NF;i++) if ($i=="ACT_SITE") print $1 "\t" $(i+1)}' uniprot_anno.txt > act_site.txt

awk '{for(i=1;i<=NF;i++) if ($i=="BINDING") print $1 "\t" $(i+1)}' uniprot_anno.txt > binding_site.txt
awk '{for(i=1;i<=NF;i++) if ($i=="MOD_RES") print $1 "\t" $(i+1)}' uniprot_anno.txt > mod_res.txt
awk '{for(i=1;i<=NF;i++) if ($i=="CARBOHYD") print $1 "\t" $(i+1)}' uniprot_anno.txt > carbohyd.txt
awk '{for(i=1;i<=NF;i++) if ($i=="MUTAGEN") print $1 "\t" $(i+1)}' uniprot_anno.txt > mutagen.txt
awk '{for(i=1;i<=NF;i++) if ($i=="DISULFID") print $1 "\t" $(i+1) "\t" $(i+2)}' uniprot_anno.txt > disulfide.txt
awk '{for(i=1;i<=NF;i++) if ($i=="MUTAGEN") print $1 "\t" $(i+1)}' uniprot_anno.txt > mutagen.txt
awk '{for(i=1;i<=NF;i++) if ($i=="CA_BIND") print $1 "\t" $(i+1) "\t" $(i+2)}' uniprot_anno.txt > ca_bind.txt
awk '{for(i=1;i<=NF;i++) if ($i=="DNA_BIND") print $1 "\t" $(i+1) "\t" $(i+2)}' uniprot_anno.txt > dna_bind.txt
awk '{for(i=1;i<=NF;i++) if ($i=="METAL") print $1 "\t" $(i+1)}' uniprot_anno.txt > metal.txt
awk '{for(i=1;i<=NF;i++) if ($i=="LIPID") print $1 "\t" $(i+1)}' uniprot_anno.txt > lipid.txt
awk '{for(i=1;i<=NF;i++) if ($i==" SITE") print $1 "\t" $(i+1)}' uniprot_anno.txt > site.txt






