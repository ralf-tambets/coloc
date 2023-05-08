#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--eqtl_file"), type="character", default=NULL,
              help="eQTL file indexed by chromosome", metavar = "type"),
  make_option(c("--gwas_file"), type="character", default=NULL,
              help="GWAS file indexed by chromosome", metavar = "type"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Prefix of the output files.", metavar = "type"),
  make_option(c("--outdir"), type="character", default=".",
              help="Path to the output directory. Defaults to work folder.", metavar = "type")
)

opt <- optparse::parse_args(OptionParser(option_list=option_list))

eqtl_file = opt$eqtl_file
gwas_file = opt$gwas_file
output_prefix = opt$output_prefix
outdir = opt$outdir

message("__________ OPTIONS __________")
message("eqtl_file: ", eqtl_file)
message("gwas_file: ", gwas_file)
message("output_prefix: ", output_prefix)
message("outdir: ", outdir)

eqtl <- read_tsv(eqtl_file)
eqtl <- eqtl %>%
  select(molecular_trait_id, variant, pip, cs_id) %>% rename(eqtl_id = molecular_trait_id, eqtl_pip=pip, eqtl_cs_id = cs_id) %>%
  distinct()

gwas <- read_tsv(gwas_file)
gwas <- gwas %>%
  select(molecular_trait_id, region, variant, pip, cs_id) %>% rename(gwas_id = molecular_trait_id, gwas_region = region, gwas_pip=pip, gwas_cs_id = cs_id) %>%
  distinct()

joined <- inner_join(eqtl, gwas, by=c("variant")) %>%
  mutate(unique_id = paste0(eqtl_id, "_", gwas_id)) %>%
  mutate(clpp_variant = eqtl_pip*gwas_pip) %>%
  distinct()
  
grouped <- joined %>%
  group_by(unique_id, eqtl_cs_id) %>%
  reframe(unique_id, eqtl_cs_id, clpp_cs = sum(clpp_variant), n_snps = n()) %>%
  ungroup() %>%
  distinct()

clpp_results <- grouped %>%
  select(unique_id, eqtl_cs_id, clpp_cs, n_snps) %>%
  left_join(joined, by=c("unique_id", "eqtl_cs_id")) %>%
  ungroup() %>%
  select(eqtl_id, gwas_id, eqtl_cs_id, gwas_region, variant, clpp_variant, clpp_cs, n_snps) %>%
  distinct() %>%
  arrange(desc(clpp_cs), desc(clpp_variant))

if (!dir.exists(outdir)) dir.create(outdir)

if (!is.null(output_prefix)) {
  clpp_file_name = file.path(outdir, output_prefix)
} else {
  clpp_file_name = file.path(outdir, "no_clpp_output_prefix_given.tsv")
}

if(!is.null(clpp_results) && nrow(clpp_results) > 0){
  message(" ## write colocalisation results to ", clpp_file_name )
  write_tsv(clpp_results, clpp_file_name)
} else {
  file.create(clpp_file_name)
  message("Coloc results are empty or null!")
}