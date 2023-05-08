#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("coloc"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("seqminer"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--eqtl_file"), type="character", default=NULL,
              help="File indexed by chromosome name", metavar = "type"),
  make_option(c("--gwas_file"), type="character", default=NULL,
              help="File indexed by chromosome", metavar = "type"),
  make_option(c("--reference_file"), type="character", default=NULL,
              help="Gene metadata file", metavar = "type"),
  make_option(c("--chromosome"), type="character", default=NULL,
              help="Chromosome", metavar = "type"),
  make_option(c("--chunk"), type="character", default="1 1", 
              help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed. [default \"%default\"]", metavar = "type"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Prefix of the output files.", metavar = "type"),
  make_option(c("--outdir"), type="character", default="./coloc_results/",
              help="Path to the output directory. [default \"%default\"]", metavar = "type")
  )

opt <- optparse::parse_args(OptionParser(option_list=option_list))

eqtl_file = opt$eqtl_file
gwas_file = opt$gwas_file
reference_file = opt$reference_file
chromosome = opt$chromosome
chunk = opt$chunk
output_prefix = opt$output_prefix
outdir = opt$outdir

message("__________ OPTIONS __________")
message("eqtl_file: ", eqtl_file)
message("gwas_file: ", gwas_file)
message("reference_file: ", reference_file)
message("chromosome:", chromosome)
message("chunk: ", chunk)
message("output_prefix: ", output_prefix)
message("output_file_path: ", outdir)

splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

splitIntoChunks <- function(chunk_number, n_chunks, n_total){
  chunk_size = floor(n_total/(n_chunks))
  batches = splitIntoBatches(n_total,chunk_size)
  batches[batches > n_chunks] = n_chunks
  selected_batch = batches == chunk_number
  return(selected_batch)
}

coloc_in_region <- function(eqtl_gene_name, eqtl_gene_for_coloc, gwas_gene_table){
  results = tibble(
    nsnps = numeric(),
    hit1 = character(),
    hit2 = character(),
    PP.H0.abf = numeric(),
    PP.H1.abf = numeric(),
    PP.H2.abf = numeric(),
    PP.H3.abf = numeric(),
    PP.H4.abf = numeric(),
    idx1 = numeric(),
    idx2 = numeric(),
    gene = character(),
    gwas = character(),
    gwas_region = character()
  )
  if (nrow(gwas_gene_table) == 0){
    return(results)
  }
  for (selected_gene_region in unique(gwas_gene_table$molecular_trait_id)){
    selected_gene_table <- gwas_gene_table %>% filter(molecular_trait_id == selected_gene_region)
    if (nrow(selected_gene_table) < 300){
      next
    }
    selected_gene_for_coloc <- t(as.matrix(selected_gene_table %>% select(paste0("lbf_variable",1:10))))
    colnames(selected_gene_for_coloc) <- selected_gene_table$variant

    skip_to_next <- FALSE

    tryCatch(
    expr = {
        coloc_results <- coloc.bf_bf(eqtl_gene_for_coloc, selected_gene_for_coloc, p1=1e-4, p2=1e-4, p12=5e-6)$summary
    },
    error = function(e){ 
        message(e)
        skip_to_next <- TRUE
    }
    )

    if (skip_to_next) {
      next
    }

    if (!is.null(coloc_results)){
      coloc_results$gene <- rep(eqtl_gene_name, nrow(coloc_results))
      coloc_results$gwas <- rep(selected_gene_table$molecular_trait_id[1], nrow(coloc_results))
      coloc_results$gwas_region <- rep(selected_gene_region, nrow(coloc_results))
      results <- rbind(results, coloc_results) 
    }
  }
  return(results)
}

get_genes_in_region <- function(eqtl_gene_name, gwas_file, reference_table){
  gene_reference <- reference_table[reference_table$molecular_trait_id == eqtl_gene_name,]
  if (nrow(gene_reference) == 0){
    message(paste("No genes of name", eqtl_gene_name, "in reference."))
    return(data.frame())
  }
  tabixrange = paste0(gene_reference$chromosome,":",gene_reference$startpos,"-",gene_reference$endpos)
  genes_in_region <- tabix.read.table(gwas_file, tabixRange=tabixrange)
  if (nrow(genes_in_region) == 0){
    message(paste("No genes in region", tabixrange))
    return(data.frame())
  }
  colnames(genes_in_region) <- c("molecular_trait_id", "region", "variant", "chromosome", "position", paste0("lbf_variable",1:10))
  return(genes_in_region)
}

coloc_per_gene <- function(eqtl_gene_name, eqtl_gene_table, gwas_file, reference_table){
  selected_gene_table <- eqtl_gene_table %>% filter(molecular_trait_id == eqtl_gene_name)
  if (nrow(selected_gene_table)==0){
    message(paste("No genes of name", eqtl_gene_name, "in eqtl_gene_table."))
    return(tibble(
      nsnps = numeric(),
      hit1 = character(),
      hit2 = character(),
      PP.H0.abf = numeric(),
      PP.H1.abf = numeric(),
      PP.H2.abf = numeric(),
      PP.H3.abf = numeric(),
      PP.H4.abf = numeric(),
      idx1 = numeric(),
      idx2 = numeric(),
      gene = character(),
      gwas = character(),
      gwas_region = character()
    ))
  }
  eqtl_gene_for_coloc <- t(as.matrix(selected_gene_table %>% select(paste0("lbf_variable",1:10))))
  colnames(eqtl_gene_for_coloc) <- selected_gene_table$variant
  gwas_genes_in_region <- get_genes_in_region(eqtl_gene_name, gwas_file, reference_table)
  
  coloc_results <- coloc_in_region(eqtl_gene_name,eqtl_gene_for_coloc, gwas_genes_in_region)
  return(coloc_results)
}

eqtl_gene_table <- tabix.read.table(eqtl_file, tabixRange=paste0(chromosome, ":1-2147483647"))
coloc_results <- data.frame()
if (!is.null(eqtl_gene_table) && nrow(eqtl_gene_table)>0){
  colnames(eqtl_gene_table) <- c("molecular_trait_id", "region", "variant", "chromosome", "position", paste0("lbf_variable",1:10))
  
  reference_table <- read_tsv(reference_file)
  
  eqtl_gene_table <- eqtl_gene_table %>%
    filter(molecular_trait_id %in% reference_table$molecular_trait_id) 

  eqtl_gene_names <- unique(eqtl_gene_table$molecular_trait_id)

  chunk_vector = strsplit(chunk, split = " ") %>% unlist() %>% as.numeric()
  chunk_id = chunk_vector[1]
  n_chunks = chunk_vector[2]
  if (length(eqtl_gene_names) < n_chunks) {
    n_chunks = length(eqtl_gene_names)
  }
  if (chunk_id <= n_chunks){
    selected_chunk = splitIntoChunks(chunk_id, n_chunks, length(eqtl_gene_names))
    selected_eqtl_gene_names = data.frame(eqtl_gene_names[selected_chunk])[1]
      
    coloc_results <- do.call("rbind", apply(selected_eqtl_gene_names, 1, coloc_per_gene,
                                            eqtl_gene_table = eqtl_gene_table,
                                            gwas_file = gwas_file,
                                            reference_table = reference_table))
  }
} else {
  message(paste ("Empty table for tabixrange", paste0(chromosome, ":1-2147483647"), "!"))
}

if (!dir.exists(outdir)) dir.create(outdir)

if (!is.null(output_prefix)) {
  file_name = file.path(outdir, output_prefix)
} else {
  file_name = file.path(outdir, "no_output_prefix_given.tsv")
}

if(nrow(coloc_results) > 0){
  message(" ## write colocalisation results to ", file_name )
  coloc_results = coloc_results %>% dplyr::select(c("gene", "gwas", "gwas_region", "hit1", "hit2"), everything())
  write_tsv(coloc_results, file_name)
} else {
  file.create(file_name)
  message("Coloc results are empty or null!")
}
