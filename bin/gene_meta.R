suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--permutation_file"), type="character", default=NULL,
              help="Permutation file name", metavar = "type"),
  make_option(c("--metadata_file"), type="character", default=NULL,
              help="Metadata file name", metavar = "type"),
  make_option(c("--qtl_dataset"), type="character", default=NULL,
              help="QTL dataset name", metavar = "type")
)

opt <- optparse::parse_args(OptionParser(option_list=option_list))

permutation_file = opt$permutation_file
metadata_file = opt$metadata_file
qtl_dataset = opt$qtl_dataset

message("__________ OPTIONS __________")
message("permutation_file:", permutation_file)
message("metadata_file: ", metadata_file)
message("qtl_dataset: ", qtl_dataset)

permutations <- read_tsv(permutation_file)

permutations <- permutations %>%
  mutate(FDR = p.adjust(p=p_beta, method='fdr')) %>%
  filter(FDR < 0.01) %>%
  select(molecular_trait_object_id, molecular_trait_id) %>%
  distinct()

metadata <- read_tsv(metadata_file)
  
joined <- inner_join(permutations, metadata, by=c("molecular_trait_id"="phenotype_id")) %>%
  select(molecular_trait_id, chromosome, phenotype_pos) %>%
  mutate(startpos = pmax(1, phenotype_pos-1000000), endpos = phenotype_pos+1000000) %>%
  select(-phenotype_pos) %>%
  distinct()

if(nrow(joined) == 0) {
  stop("Position file not created!")
}

write_tsv(joined, paste0(qtl_dataset,".gene_positions.tsv"))