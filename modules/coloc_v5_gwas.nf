process tabix_eqtl_by_chr_position{
    container 'quay.io/eqtlcatalogue/colocalisation:v21.01.1'

    input:
    tuple val(eqtl_dataset_id), file(eqtl_ss), file(eqtl_ss_index), file(eqtl_permuted), file(eqtl_credible_sets), file(eqtl_lbf)

    output:
    tuple val(eqtl_dataset_id), file("sorted_${eqtl_lbf.simpleName}.txt.gz"), file("sorted_${eqtl_lbf.simpleName}.txt.gz.tbi")
    
    script:
    """
    zcat ${eqtl_lbf} | awk -F'\t' 'BEGIN{OFS=FS}NR==1{print \$0; next}{print \$0 | "sort -k4,4n -k5,5n"}' | bgzip > sorted_${eqtl_lbf.simpleName}.txt.gz
    tabix -s4 -b5 -e5 -S1 -f sorted_${eqtl_lbf.simpleName}.txt.gz
    """
}

process tabix_gwas_by_chr_position{
    container 'quay.io/eqtlcatalogue/colocalisation:v21.01.1'

    input:
    tuple val(gwas_subset), file(gwas_clpp), file(gwas_coloc3), file(gwas_coloc5)

    output:
    tuple val(gwas_subset), file("sorted_${gwas_coloc5.simpleName}.txt.gz"), file("sorted_${gwas_coloc5.simpleName}.txt.gz.tbi")
    
    script:
    """
    zcat ${gwas_coloc5} | awk -F'\t' 'BEGIN{OFS=FS}NR==1{print \$0; next}{print \$0 | "sort -k4,4n -k5,5n"}' | bgzip > sorted_${gwas_coloc5.simpleName}.txt.gz
    tabix -s4 -b5 -e5 -S1 -f sorted_${gwas_coloc5.simpleName}.txt.gz
    """
}

process run_coloc_v5{
    container = 'quay.io/eqtlcatalogue/coloc:v23.05.1'

    input:
    tuple val(chr), val(batch), val(eqtl_dataset_id), file(reference_file), file(eqtl_file), file(eqtl_index), val(gwas_subset), file(gwas_file), file(gwas_index)

    output:
    tuple val(eqtl_dataset_id), val("${eqtl_dataset_id}_${gwas_subset}"), file("${eqtl_dataset_id}_${gwas_subset}_${chr}_${batch}_${params.chr_batches}.coloc.v5.tsv")

    script:
    """
    Rscript $baseDir/bin/coloc_v5_gwas.R \\
        --eqtl_file=$eqtl_file \\
        --gwas_file=$gwas_file \\
        --reference=$reference_file \\
        --chromosome=${chr} \\
        --chunk='${batch} ${params.chr_batches}' \\
        --output_prefix=${eqtl_dataset_id}_${gwas_subset}_${chr}_${batch}_${params.chr_batches}.coloc.v5.tsv \\
        --outdir=.
    """
}

process merge_coloc_v5_results{
    publishDir "${params.outdir}/coloc_v5_results_merged/${eqtl_dataset_id}", mode: 'copy'
    container 'quay.io/eqtlcatalogue/colocalisation:v21.01.1'

    input:
    tuple val(eqtl_dataset_id), val(gwas_eqtl_id), file(gwas_eqtl_coloc_results_batch_files)

    output:
    tuple val(gwas_eqtl_id), file("${gwas_eqtl_id}.coloc.v5.txt.gz")

    script:
    """
    awk 'NR == 1 || FNR > 1{print}' ${gwas_eqtl_coloc_results_batch_files.join(' ')} | awk -F'\t' 'BEGIN{OFS=FS}NR==1{print \$0; next}{print \$0 | "sort -k11,11gr"}' | bgzip -c > ${gwas_eqtl_id}.coloc.v5.txt.gz
    """
}