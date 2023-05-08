# coloc
**NB! Not a polished product.**
A workflow for running colocalisation analysis on eQTL vs pQTL or eQTL vs GWAS data. The analyses performed are CLPP (colocalisation posterior probability), coloc.abf (here referred to as coloc v3) and coloc.susie (here referred to as coloc v5).
To run the workflow, submit `run_pqtl.sh` or `run_gwas` to HPC.

## Usage
Input parameters:
*  **--eqtl_tsv** TSV file with 6 columns. Example in `sample_data/sample_eqtl.tsv`.
*  **eqtl_dataset_id** Dataset name or identifier.
*  **eqtl_sumstats** Summary statistics file containing the following columns: molecular_trait_id, chromosome, position, ref, alt, variant, ma_samples, maf, pvalue, beta, se, type, ac, an, r2, molecular_trait_object_id, gene_id, median_tpm, rsid
*  **eqtl_sumstats_index** Tabix-generated index file for the summary statistics file.
*  **eqtl_permuted** Permutation file containing the following columns: molecular_trait_object_id, molecular_trait_id, n_traits, n_variants, variant, chromosome, position, pvalue, beta, p_perm, p_beta
*  **eqtl_credible_sets** Credible sets file containing the following columns: molecular_trait_id, gene_id, cs_id, variant, rsid, cs_size, pip, pvalue, beta, se, z, cs_min_r2, region
*  **eqtl_lbf** LBF variable file containing the following columns: molecular_trait_id, region, variant, chromosome, position, lbf_variable1, lbf_variable2, lbf_variable3, lbf_variable4, lbf_variable5, lbf_variable6, lbf_variable7, lbf_variable8, lbf_variable9, lbf_variable10
*  **--pqtl_tsv** (only for `run_pqtl.sh). TSV file formatted exactly like eqtl_tsv, only the prefix changes from "eqtl" to "pqtl". Example in `sample_data/sample_pqtl.tsv`.
*  **--gwas_tsv** (only for `run_gwas.sh`). TSV file with 5 columns. Example in `sample_data/sample_gwas.tsv`.
*  **gwas_subset** Subset name or identifier.
*  **gwas_clpp** TSV file containing the following columns: molecular_trait_id, region, variant, chromosome, position, ref, alt, cs_id, cs_index, pip, z, cs_min_r2
*  **gwas_coloc3** TSV file containing the following columns: molecular_trait_id, region, variant, ref, alt, chromosome, position, maf, beta, se, an, ac, n, log10p, info
*  **gwas_coloc5** TSV file containing the following columns:
molecular_trait_id, region, variant, chromosome, position, lbf_variable1, lbf_variable2, lbf_variable3, lbf_variable4, lbf_variable5, lbf_variable6, lbf_variable7, lbf_variable8, lbf_variable9, lbf_variable10
*  **--chromosomes** comma-separated list of chromosomes.
*  **--n_chromosomes** number of chromosomes in previous list
*  **--chr_batches** number of splits to create for each chromosome (1 would analyse the entire chromosome, 3 would split each chromosome into 3 different jobs).
*  **--metadata_file** File to use for reference file creation. [Examples](https://zenodo.org/record/3366011#.ZBgpJrRBwl4)
*  **--outdir** Directory to write the results to.

## Bioinformatics example
1. Download necessary files by navigating to `bioinformatics_example/bioinformatics_data` and running `sbatch download_data.sh` in the terminal (**takes a while**).
2. Navigate to main folder, start workflow execution from terminal with the command `sbatch run_bioinformatics_gwas.sh` or `sbatch run_bioinformatics_pqtl.sh`
3. The results will appear in `results_bioinformatics_gwas` or `results_bioinformatics_pqtl`, respectively.
