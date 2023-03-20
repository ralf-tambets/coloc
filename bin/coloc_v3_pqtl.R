#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("coloc"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--eqtl_file"), type="character", default=NULL,
              help="eQTL file indexed by chromosome", metavar = "type"),
  make_option(c("--pqtl_file"), type="character", default=NULL,
              help="pQTL file indexed by chromosome", metavar = "type"),
  make_option(c("--reference_file"), type="character", default=NULL,
              help="Gene metadata file", metavar = "type"),
  make_option(c("--chromosome"), type="character", default=NULL,
              help="Chromosome number", metavar="type"),
  make_option(c("--chunk"), type="character", default="1 1",
              help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed. [default \"%default\"]", metavar = "type"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Prefix of the output files.", metavar = "type"),
  make_option(c("--outdir"), type="character", default="./coloc_results/",
              help="Path to the output directory. [default \"%default\"]", metavar = "type")
)

opt <- optparse::parse_args(OptionParser(option_list=option_list))

eqtl_file = opt$eqtl_file
pqtl_file = opt$pqtl_file
reference_file = opt$reference_file
selected_chromosome = opt$chromosome
chunk = opt$chunk
output_prefix = opt$output_prefix
outdir = opt$outdir

message("__________ OPTIONS __________")
message("eqtl_file: ", eqtl_file)
message("pqtl_file: ", pqtl_file)
message("reference_file: ", reference_file)
message("chromosome: ", selected_chromosome)
message("chunk: ", chunk)
message("output_prefix: ", output_prefix)
message("output_file_path: ", outdir)

splitIntoBatches <- function(n, batch_size){
  message(paste0("n",n))
  message(paste0("batch_size", batch_size))
  n_batches = ceiling(n/batch_size)
  message("n_batches", n_batches)
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


coloc_in_region <- function(eqtl_gene_name, selected_eqtl_gene, pqtl_genes_in_region){
  results = tibble(
    nsnps = numeric(),
    PP.H0.abf = numeric(),
    PP.H1.abf = numeric(),
    PP.H2.abf = numeric(),
    PP.H3.abf = numeric(),
    PP.H4.abf = numeric(),
    eqtl = character(),
    pqtl = character()
  )
  
  for (pqtl_gene in unique(pqtl_genes_in_region$molecular_trait_id)){
    message(pqtl_gene)
    selected_pqtl_gene <- pqtl_genes_in_region %>% filter(molecular_trait_id == pqtl_gene) %>%
      select(-rsid) %>%
      distinct() %>%
      mutate(id = paste(chromosome, position, sep = ":")) %>% 
      group_by(id) %>% 
      mutate(row_count = n()) %>% dplyr::ungroup() %>% 
      filter(row_count == 1) %>%
      filter(!is.nan(se)) %>%
      filter(!is.na(se)) %>%
      select(molecular_trait_id, variant, maf, beta, se, an) %>%
      mutate(maf = as.numeric(maf), beta = as.numeric(beta), se = as.numeric(se), an = as.numeric(an))
    
    n_rows_before = max(nrow(selected_eqtl_gene), nrow(selected_pqtl_gene))

    if (nrow(selected_pqtl_gene)<300){
      next
    }
    
    #TO CIRCUMNAVIGATE THE FOLLOWING ERROR: Error in coloc.abf(dataset1 = eqtl_dataset, dataset2 = pqtl_dataset, p1 = 1e-04, dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified
    selected_eqtl_gene_clone <- selected_eqtl_gene %>% filter(variant %in% selected_pqtl_gene$variant)
    selected_pqtl_gene <- selected_pqtl_gene %>% filter(variant %in% selected_eqtl_gene_clone$variant)

    n_rows_after = nrow(selected_eqtl_gene_clone)

    if (n_rows_after/n_rows_before < 0.1){
      message(paste0("Too few common variants!"))
      message(eqtl_gene_name)
      message(pqtl_gene)
      next
    }

    varbeta = selected_eqtl_gene_clone$se^2
    N = selected_eqtl_gene_clone$an/2
    MAF = selected_eqtl_gene_clone$maf

    #sdY.est checks
    is_oneover_na = all(is.na(1/varbeta))
    is_nvx_na = all(is.na(2*N*MAF*(1-MAF)))

    if (is_oneover_na | is_nvx_na) {
      message("One column is NA")
      next
    }
    
    eqtl_dataset <- list(varbeta = varbeta,
                         N = N,
                         MAF = MAF,
                         type="quant",
                         beta = selected_eqtl_gene_clone$beta,
                         snp=selected_eqtl_gene_clone$variant)
    
    varbeta = selected_pqtl_gene$se^2
    N = selected_pqtl_gene$an/2
    MAF = selected_pqtl_gene$maf

    #sdY.est checks
    is_oneover_na = all(is.na(1/varbeta))
    is_nvx_na = all(is.na(2*N*MAF*(1-MAF)))

    if (is_oneover_na | is_nvx_na) {
      message("One column is NA")
      next
    }
    
    pqtl_dataset <- list(varbeta = varbeta,
                         N = N,
                         MAF = MAF,
                         type="quant",
                         beta = selected_pqtl_gene$beta,
                         snp=selected_pqtl_gene$variant)
    

    coloc_results <- coloc.abf(dataset1 = eqtl_dataset, dataset2 = pqtl_dataset, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
    coloc_results <- as_tibble(t(as.data.frame(coloc_results$summary)))    
    


    if (is.null(coloc_results)){
      next
    }
    coloc_results$eqtl = rep(eqtl_gene_name, nrow(coloc_results))
    coloc_results$pqtl = rep(selected_pqtl_gene$molecular_trait_id[1], nrow(coloc_results))
    results <- rbind(results, coloc_results)
  }
  
  return(results)
}

coloc_per_gene <- function(eqtl_gene_name, eqtl_gene_table, pqtl_file, reference_table){
  print(eqtl_gene_name)
  coloc_results <- tibble(
    nsnps = numeric(),
    PP.H0.abf = numeric(),
    PP.H1.abf = numeric(),
    PP.H2.abf = numeric(),
    PP.H3.abf = numeric(),
    PP.H4.abf = numeric(),
    eqtl = character(),
    pqtl = character()
  )
  selected_eqtl_gene <- eqtl_gene_table %>% filter(molecular_trait_id == eqtl_gene_name)%>%
    select(-rsid) %>%
    distinct() %>%
    mutate(id = paste(chromosome, position, sep = ":")) %>% 
    group_by(id) %>% 
    mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    filter(row_count == 1) %>%
    filter(!is.nan(se)) %>%
    filter(!is.na(se)) %>%
    select(molecular_trait_id, variant, maf, beta, se, an)
  message(paste0("eqtl geen valitud ", nrow(selected_eqtl_gene)))
  if (nrow(selected_eqtl_gene)<300){
    message("Not enough rows.")
    return(coloc_results)
  }
  

  gene_reference <- reference_table[reference_table$molecular_trait_id == eqtl_gene_name,]
  pqtl_genes_in_region <- pqtl_gene_table %>%
    filter(chromosome == gene_reference$chromosome) %>%
    filter(position >= gene_reference$startpos) %>%
    filter(position <= gene_reference$endpos)
  if (nrow(pqtl_genes_in_region) == 0){
    return(coloc_results)
  }
  message(paste0("colociks valmis ", nrow(pqtl_genes_in_region)))
  coloc_results <- coloc_in_region(eqtl_gene_name, selected_eqtl_gene, pqtl_genes_in_region)
  return(coloc_results)
}

coloc_results <- data.frame()

eqtl_gene_table = read_tsv(eqtl_file) %>%
  filter(chromosome == selected_chromosome)
pqtl_gene_table = read_tsv(pqtl_file) %>%
  filter(chromosome == selected_chromosome)

if (!is.null(eqtl_gene_table)) {
  if (nrow(eqtl_gene_table)>0){
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
      selected_eqtl_gene_names = data.frame(eqtl_gene_names[selected_chunk])
        
      coloc_results <- do.call("rbind", apply(selected_eqtl_gene_names, 1, coloc_per_gene,
                                              eqtl_gene_table = eqtl_gene_table,
                                              pqtl_file = pqtl_file,
                                              reference_table = reference_table))
    }
  } else {
    message("eqtl file is empty!")
  }

} else {
  message("eqtl file is null!")
}



if (!dir.exists(outdir)) dir.create(outdir)

if (!is.null(output_prefix)) {
  file_name = file.path(outdir, output_prefix)
} else {
  file_name = file.path(outdir, "no_output_prefix_given.tsv")
}

if(nrow(coloc_results > 0)){
  message(" ## write colocalisation results to ", file_name )
  coloc_results = coloc_results %>% dplyr::select(c("eqtl", "pqtl"), everything())
  write_tsv(coloc_results, file_name)
} else {
  file.create(file_name)
  message("Coloc results are empty!")
}
