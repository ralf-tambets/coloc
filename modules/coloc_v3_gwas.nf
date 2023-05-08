process run_coloc_v3{
    container = 'quay.io/eqtlcatalogue/coloc:v23.05.1'

    input:
    tuple val(eqtl_dataset_id), file(reference_file), file(eqtl_ss), file(eqtl_ss_index), file(eqtl_permuted), file(eqtl_credible_sets), file(eqtl_lbf), val(chr), val(batch), val(gwas_subset), file(gwas_clpp), file(gwas_coloc3), file(gwas_coloc5)

    output:
    tuple val(eqtl_dataset_id), val("${eqtl_dataset_id}_${gwas_subset}"), file("${eqtl_dataset_id}_${gwas_subset}_${chr}_${batch}_${params.chr_batches}.coloc.v3.tsv")

    script:
    """
    Rscript $baseDir/bin/coloc_v3_gwas.R \\
        --eqtl_file=$eqtl_ss \\
        --gwas_file=$gwas_coloc3 \\
        --reference_file=$reference_file \\
        --chromosome=$chr \\
        --chunk='${batch} ${params.n_batches}' \\
        --output_prefix=${eqtl_dataset_id}_${gwas_subset}_${chr}_${batch}_${params.chr_batches}.coloc.v3.tsv \\
        --outdir=.
    """
}

process merge_coloc_v3_results{
    publishDir "${params.outdir}/coloc_v3_results_merged/${eqtl_dataset_id}", mode: 'copy'
    container 'quay.io/eqtlcatalogue/colocalisation:v21.01.1'
    
    input:
    tuple val(eqtl_dataset_id), val(eqtl_gwas_id), file(eqtl_gwas_coloc_results_batch_files)

    output:
    tuple val(eqtl_gwas_id), file("${eqtl_gwas_id}.coloc.v3.txt.gz")

    script:
    """
    awk 'NR == 1 || FNR > 1{print}' ${eqtl_gwas_coloc_results_batch_files.join(' ')} | awk -F'\t' 'BEGIN{OFS=FS}NR==1{print \$0; next}{print \$0 | "sort -k9,9gr"}' | bgzip -c > ${eqtl_gwas_id}.coloc.v3.txt.gz
    """
}