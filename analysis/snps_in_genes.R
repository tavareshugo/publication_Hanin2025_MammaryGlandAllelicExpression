library(tidyverse)
library(vroom)
library(slider)
library(valr)

# import variants
variants <- vroom("data/external/strains/B6xCAST.snps_and_indels.merged.bed.gz")

# get gene annotation - last GRCm38 version
mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl",
                   host="nov2020.archive.ensembl.org")

# gene information
gene_annot <- getBM(attributes = c("ensembl_gene_id",
                                   # "chromosome_name","strand",
                                   # "start_position", "end_position",
                                   "external_gene_name",
                                   "percentage_gene_gc_content",
                                   "gene_biotype"),
                    mart = mart)

transcript_annot <- getBM(attributes = c("ensembl_gene_id",
                                         "ensembl_transcript_id",
                                         "ensembl_exon_id",
                                         "exon_chrom_start",
                                         "exon_chrom_end",
                                         "chromosome_name",
                                         "strand",
                                         "transcript_biotype"),
                          mart = mart)



# Variants across genome --------------------------------------------------

variants <- variants %>%
  group_by(chrom) %>%
  arrange(chrom, start) %>%
  mutate(roll_count = slide_index_dbl(ref, start, length, .before = 100e3, .after = 100e3)) %>%
  ungroup()

variants <- variants %>%
  mutate(type = ifelse(str_length(alt) != str_length(ref), "indel", "snp"))

variants_merged <- variants %>%
  group_by(chrom) %>%
  bed_merge()
