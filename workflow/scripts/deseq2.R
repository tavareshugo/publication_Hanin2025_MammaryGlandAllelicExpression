library(DESeq2)
library(biomaRt)
library(tximport)
library(vroom)
library(dplyr)
#library(tidyr)
library(tibble)


# Prepare tx2gene tables --------------------------------------------------

# transcript-to-gene file
tx2gene <- vroom("data/external/reference/transcript2gene.tsv",
                 col_types = c("cc")) %>%
  # ensure columns are in the correct order for tximport
  select(transcript, gene)

# make a diploid version of the table
tx2gene_diploid <- tx2gene %>%
  mutate(transcript_diploid = paste0("C.", transcript),
         gene_diploid = paste0("C.", gene))
tx2gene_diploid <- tx2gene %>%
  mutate(transcript_diploid = paste0("B.", transcript),
         gene_diploid = paste0("B.", gene)) %>%
  bind_rows(tx2gene_diploid)

# list of input files
salmon_files <- list.files("results/salmon_diploid/",
                           pattern = "quant.sf",
                           full.names = TRUE,
                           recursive = TRUE)
names(salmon_files) <- basename(dirname(salmon_files))



# Prepare colData ---------------------------------------------------------

samples <- vroom("read_info_test.csv")
samples <- samples %>%
  distinct(sample, cell_type, stage, timepoint, cross, animal_id)
samples <- samples %>% column_to_rownames("sample")

# ensure it's the same order as files
samples <- samples[names(salmon_files), ]


# Gene-level Allele Specific -----------------------------------------------------------

txi_gene_allele <- tximport(salmon_files, type = "salmon",
                            tx2gene = tx2gene_diploid %>%
                              select(transcript_diploid, gene_diploid) %>%
                              distinct(),
                            importer = function(x) vroom::vroom(x,
                                                                progress = FALSE,
                                                                col_types = c("ciddd")))

dds_gene_allele <- DESeqDataSetFromTximport(txi_gene_allele,
                                            colData = samples,
                                            design = ~ 1)


# Gene-level Regular -----------------------------------------------------------

txi_gene_regular <- tximport(salmon_files, type = "salmon",
                             tx2gene = tx2gene_diploid %>%
                               select(transcript_diploid, gene) %>%
                               distinct(),
                             importer = function(x) vroom::vroom(x,
                                                                 progress = FALSE,
                                                                 col_types = c("ciddd")))

dds_gene_regular <- DESeqDataSetFromTximport(txi_gene_regular,
                                             colData = samples,
                                             design = ~ 1)


# Isoform-level Allele Specific -----------------------------------------------------------

# in this case we turn off summarisation with txOut = TRUE
txi_isoform_allele <- tximport(salmon_files, type = "salmon",
                               tx2gene = tx2gene_diploid %>%
                                 select(transcript_diploid, gene_diploid) %>%
                                 distinct(),
                               importer = function(x) vroom::vroom(x,
                                                                   progress = FALSE,
                                                                   col_types = c("ciddd")),
                               txOut = TRUE)

dds_isoform_allele <- DESeqDataSetFromTximport(txi_isoform_allele,
                                               colData = samples,
                                               design = ~ 1)


# Isoform-level Regular -----------------------------------------------------------

txi_isoform_regular <- tximport(salmon_files, type = "salmon",
                             tx2gene = tx2gene_diploid %>%
                               select(transcript_diploid, transcript) %>%
                               distinct(),
                             importer = function(x) vroom::vroom(x,
                                                                 progress = FALSE,
                                                                 col_types = c("ciddd")))

dds_isoform_regular <- DESeqDataSetFromTximport(txi_isoform_regular,
                                             colData = samples,
                                             design = ~ 1)


# Save Objects ------------------------------------------------------------

saveRDS(dds_gene_allele,
        "data/processed/dds_gene_allele.rds")

saveRDS(dds_gene_regular,
        "data/processed/dds_gene_regular.rds")

saveRDS(dds_isoform_allele,
        "data/processed/dds_isoform_allele.rds")

saveRDS(dds_isoform_regular,
        "data/processed/dds_isoform_regular.rds")
