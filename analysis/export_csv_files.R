# this script exports the counts matrices to CSV files for collaborators
# output is sent to a temporary dir as we do not need to keep these long-term

library(DESeq2)
library(tidyverse)

dds <- readRDS("results/DESeqDataSet/dds_gene_regular.rds")

# get the gene names
gene_meta <- rowData(dds) |> 
  as_tibble(rownames = "gene_id") |> 
  select(gene_id, gene_name = external_gene_name)

# export TPM
assay(dds, "tpm") |> 
  as_tibble(rownames = "gene_id") |> 
  full_join(gene_meta, by = "gene_id") |> 
  select(gene_id, gene_name, everything()) |> 
  write_csv("temp/hybrid_samples_tpm.csv")

# export normalised logcounts
assay(dds, "logcounts") |> 
  as_tibble(rownames = "gene_id") |> 
  full_join(gene_meta, by = "gene_id") |> 
  select(gene_id, gene_name, everything()) |> 
  write_csv("temp/hybrid_samples_logcounts.csv")

# export raw counts
counts(dds) |> 
  as_tibble(rownames = "gene_id") |> 
  full_join(gene_meta, by = "gene_id") |> 
  select(gene_id, gene_name, everything()) |> 
  write_csv("temp/hybrid_samples_counts.csv")

