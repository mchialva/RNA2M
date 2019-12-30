#!/bin/bash
#### Download nt.fa database, taxonomy and create taxoner64 index ####
# must have Taxoner and Taxoner64 installed
# download nt and unzip
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip nt.gz

# download taxonomy and unzip
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
unzip -d ./ taxdmp.zip nodes.dmp names.dmp

# parse names.dmp for compatibility with taxoner64
grep "scientific name" names.dmp | awk -F\\t\|\\t '{printf("%s\t%s\n",$1,$3)}' > names

# download gi to taxid conversion and unzip
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz 
gunzip -d ./ gi_taxid_nucl.dmp

# make bowtie2 index for taxoner64
mkdir bowtie2
cd bowtie2

# parse fasta file with dbcreator (taxoner)
dbcreator -n nt -g gi_taxid_nucl.dmp -o nodes.dmp

# Index with bowtie2 (large genome indexing)
for i in *.fasta;
do basename=$(basename "$i" .fasta);
bowtie2-build --large-index --threads 50 "$i" $basename;
done
