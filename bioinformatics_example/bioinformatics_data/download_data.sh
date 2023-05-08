#!/bin/bash

#SBATCH -J data_download
#SBATCH -N 1
#SBATCH -t 06:00:00
#SBATCH --mem=8GB
#SBATCH --partition=main

#Download eQTL data
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000013/QTD000110/QTD000110.cc.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000013/QTD000110/QTD000110.cc.tsv.gz.tbi
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000013/QTD000110/QTD000110.permuted.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000013/QTD000110/QTD000110.credible_sets.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000013/QTD000110/QTD000110.lbf_variable.txt.gz

#Downlaod pQTL data
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000035/QTD000584/QTD000584.cc.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000035/QTD000584/QTD000584.cc.tsv.gz.tbi
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000035/QTD000584/QTD000584.permuted.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000035/QTD000584/QTD000584.lbf_variable.txt.gz
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000035/QTD000584/QTD000584.credible_sets.tsv.gz

#Download gene and protein metadata
wget https://zenodo.org/record/7808390/files/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz
wget https://zenodo.org/record/7808390/files/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz

#Download GWAS data
wget https://zenodo.org/record/7901534/files/VitD.clpp_combined.tsv.gz
wget https://zenodo.org/record/7901534/files/VitD.coloc3_combined.tsv.gz
wget https://zenodo.org/record/7901534/files/VitD.coloc5_combined.tsv.gz