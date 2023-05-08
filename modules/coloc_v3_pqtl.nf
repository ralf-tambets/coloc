process run_coloc_v3{
    container = 'quay.io/eqtlcatalogue/coloc:v23.05.1'

    input:
    tuple val(eqtl_dataset_id), file(reference_file), file(eqtl_ss), file(eqtl_ss_index), file(eqtl_permuted), file(eqtl_credible_sets), file(eqtl_lbf), val(chr), val(batch), val(pqtl_dataset_id), file(pqtl_ss), file(pqtl_ss_index), file(pqtl_permuted), file(pqtl_credible_sets), file(pqtl_lbf) 

    output:
    tuple val(eqtl_dataset_id), val("${eqtl_dataset_id}_${pqtl_dataset_id}"), file("${eqtl_dataset_id}_${pqtl_dataset_id}_${chr}_${batch}_${params.chr_batches}.coloc.v3.tsv")

    script:
    """
    Rscript $baseDir/bin/coloc_v3_pqtl.R \\
        --eqtl_file=$eqtl_ss \\
        --pqtl_file=$pqtl_ss \\
        --reference_file=$reference_file \\
        --chromosome=$chr \\
        --chunk='${batch} ${params.n_batches}' \\
        --output_prefix=${eqtl_dataset_id}_${pqtl_dataset_id}_${chr}_${batch}_${params.chr_batches}.coloc.v3.tsv \\
        --outdir=.
    """
}

process merge_coloc_v3_results{
    publishDir "${params.outdir}/coloc_v3_results_merged/${eqtl_dataset_id}", mode: 'copy'
    container 'quay.io/eqtlcatalogue/colocalisation:v21.01.1'

    input:
    tuple val(eqtl_dataset_id), val(eqtl_pqtl_id), file(eqtl_pqtl_coloc_results_batch_files)

    output:
    tuple val(eqtl_pqtl_id), file("${eqtl_pqtl_id}.coloc.v3.txt.gz")

    script:
    """
    awk 'NR == 1 || FNR > 1{print}' ${eqtl_pqtl_coloc_results_batch_files.join(' ')} | awk -F'\t' 'BEGIN{OFS=FS}NR==1{print \$0; next}{print \$0 | "sort -k8,8gr"}' | bgzip -c > ${eqtl_pqtl_id}.coloc.v3.txt.gz
    """
}