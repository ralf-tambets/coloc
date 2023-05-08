nextflow.enable.dsl=2

Channel.fromPath(params.eqtl_tsv)
    .ifEmpty { error "Cannot find eQTL file: ${params.eqtl_tsv}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.eqtl_dataset_id, file(row.eqtl_sumstats), file(row.eqtl_sumstats_index), file(row.eqtl_permuted), file(row.eqtl_credible_sets), file(row.eqtl_lbf)]}
    .set{ eqtl_raw_ch }

Channel.fromPath(params.pqtl_tsv)
    .ifEmpty { error "Cannot find pQTL file: ${params.pqtl_tsv}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.pqtl_dataset_id, file(row.pqtl_sumstats), file(row.pqtl_sumstats_index), file(row.pqtl_permuted), file(row.pqtl_credible_sets), file(row.pqtl_lbf)]}
    .set{ pqtl_raw_ch }

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
    tuple val(eqtl_dataset_id), file(eqtl_sumstats), file(eqtl_sumstats_index), file(eqtl_permuted), file(eqtl_credible_sets), file(eqtl_lbf), val(pqtl_dataset_id), file(pqtl_sumstats), file(pqtl_sumstats_index), file(pqtl_permuted), file(pqtl_credible_sets), file(pqtl_lbf)

    output:
    file("${eqtl_dataset_id}_${pqtl_dataset_id}.clpp.tsv.gz")

    script:
    """
    Rscript $baseDir/bin/clpp_pqtl.R \\
        --eqtl_file ${eqtl_credible_sets} \\
        --pqtl_file ${pqtl_credible_sets} \\
        --output_prefix '${eqtl_dataset_id}_${pqtl_dataset_id}.clpp.tsv.gz' \\
        --outdir .
    """
}

include { tabix_by_chr_position as tabix_eqtl_by_chr_position; tabix_by_chr_position as tabix_pqtl_by_chr_position } from './modules/coloc_v5_pqtl.nf'
include { run_coloc_v5; merge_coloc_v5_results } from './modules/coloc_v5_pqtl.nf'
include { run_coloc_v3; merge_coloc_v3_results } from './modules/coloc_v3_pqtl.nf'


workflow coloc_v5{
    take:
        reference_file_ch
    
    main:
    tabix_eqtl_by_chr_position(eqtl_raw_ch)
    tabix_pqtl_by_chr_position(pqtl_raw_ch)
    run_coloc_v5(
        chromosome_ch
        .combine(Channel.of(1..params.chr_batches))
        .combine(reference_file_ch.combine(tabix_eqtl_by_chr_position.out, by: 0))
        .combine(tabix_pqtl_by_chr_position.out)
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
        .combine(pqtl_raw_ch)
        )
    merge_coloc_v3_results(run_coloc_v3.out.groupTuple(by: [1, 0], size: params.n_chromosomes*params.chr_batches, sort: true))
}

workflow{
    create_reference_file(eqtl_raw_ch)
    coloc_v5(create_reference_file.out)
    coloc_v3(create_reference_file.out)
    run_clpp(eqtl_raw_ch.combine(pqtl_raw_ch))
}