#!/bin/bash

#SBATCH -J coloc
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --mem=8GB
#SBATCH --partition=main

module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.7.3
module load squashfs/4.4

nextflow run main_gwas.nf \
    -profile tartu_hpc \
    --eqtl_tsv /path/to/eqtl_samples.tsv \
    --gwas_tsv /path/to/gwas_samples.tsv \
    --chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X \
    --n_chromosomes 23 \
    --chr_batches 1 \
    --metadata_file /path/to/metadata_file.tsv.gz \
    -w ./work_gwas \
    --outdir ./results_gwas \
    -resume