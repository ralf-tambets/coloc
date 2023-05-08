nextflow.enable.dsl=2

Channel.fromPath(params.eqtl_tsv)
    .ifEmpty { error "Cannot find eQTL file: ${params.eqtl_tsv}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.eqtl_dataset_id, file(row.eqtl_sumstats), file(row.eqtl_sumstats_index), file(row.eqtl_permuted), file(row.eqtl_credible_sets), file(row.eqtl_lbf)]}
    .set{ eqtl_raw_ch }


Channel.fromPath(params.gwas_tsv)
    .ifEmpty { error "Cannot find GWAS_variant_list: ${params.gwas_tsv}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.gwas_subset, file(row.gwas_clpp), file(row.gwas_coloc3), file(row.gwas_coloc5)]}
    .set{ gwas_raw_ch }

Channel.fromList(params.chromosomes.split(",") as List)
    .set{ chromosome_ch }

process create_reference_file{
    container = 'quay.io/eqtlcatalogue/coloc:v23.05.1'
    
    input:
    tuple val(eqtl_dataset_id), file(eqtl_sumstats), file(eqtl_sumstats_index), file(eqtl_permuted), file(eqtl_credible_sets), file(eqtl_lbf)

    output:
    tuple val(eqtl_dataset_id), file("${eqtl_dataset_id}.gene_positions.tsv")

    script:
    """
    Rscript $baseDir/bin/gene_meta.R \\
        --permutation_file ${eqtl_permuted} \\
        --metadata_file ${params.metadata_file} \\
        --qtl_dataset ${eqtl_dataset_id};
    """
}

process run_clpp{
    publishDir "${params.outdir}/clpp_results/${eqtl_dataset_id}", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/coloc:v23.05.1'

    input:
    tuple val(eqtl_dataset_id), file(eqtl_sumstats), file(eqtl_sumstats_index), file(eqtl_permuted), file(eqtl_credible_sets), file(eqtl_lbf), val(gwas_subset), file(gwas_clpp), file(gwas_coloc3), file(gwas_coloc5)

    output:
    file("${eqtl_dataset_id}_${gwas_subset}.clpp.tsv.gz")

    script:
    """
    Rscript $baseDir/bin/clpp_gwas.R \\
        --eqtl_file ${eqtl_credible_sets} \\
        --gwas_file ${gwas_clpp} \\
        --output_prefix '${eqtl_dataset_id}_${gwas_subset}.clpp.tsv.gz' \\
        --outdir .
    """
}

include { tabix_eqtl_by_chr_position; tabix_gwas_by_chr_position; run_coloc_v5; merge_coloc_v5_results } from './modules/coloc_v5_gwas.nf'
include { run_coloc_v3; merge_coloc_v3_results } from './modules/coloc_v3_gwas.nf'


workflow coloc_v5{
    take:
        reference_file_ch
    
    main:
    tabix_eqtl_by_chr_position(eqtl_raw_ch)
    tabix_gwas_by_chr_position(gwas_raw_ch)
    run_coloc_v5(
        chromosome_ch
        .combine(Channel.of(1..params.chr_batches))
        .combine(reference_file_ch.combine(tabix_eqtl_by_chr_position.out, by: 0))
        .combine(tabix_gwas_by_chr_position.out)
        )
    merge_coloc_v5_results(run_coloc_v5.out.groupTuple(by: [1, 0], size: params.n_chromosomes*params.chr_batches, sort: true))

}

workflow coloc_v3{
    take:
        reference_file_ch
    
    main:
    run_coloc_v3(
        reference_file_ch.combine(eqtl_raw_ch.combine(chromosome_ch), by: 0)
        .combine(Channel.of(1..params.chr_batches))
        .combine(gwas_raw_ch)
        )
    merge_coloc_v3_results(run_coloc_v3.out.groupTuple(by: [1, 0], size: params.n_chromosomes*params.chr_batches, sort: true))
}

workflow{
    create_reference_file(eqtl_raw_ch)
    coloc_v5(create_reference_file.out)
    coloc_v3(create_reference_file.out)
    run_clpp(eqtl_raw_ch.combine(gwas_raw_ch))
}