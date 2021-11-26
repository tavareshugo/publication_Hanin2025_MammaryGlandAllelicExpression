setwd("~/mount/hpc_uni/mammary_gland_transcriptomes/test_pipeline/")

library(tidyverse)
theme_set(theme_minimal(base_size = 16))

regular <- read_tsv("results/salmon/stromal_t1_bc_M00323192/quant.sf")
diploid <- read_tsv("results/salmon_diploid/stromal_t1_bc_M00323192/quant.sf")

regular <- regular %>% select(gene = Name, both = NumReads)

diploid <- diploid %>%
  select(Name, NumReads) %>%
  separate(Name, c("strain", "gene")) %>%
  pivot_wider(names_from = strain, values_from = NumReads)

# compare regular salmon quantification (against standard genome) with diploid quantification
inner_join(regular, diploid, by = "gene") %>%
  ggplot(aes(both + 1, B + C + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Regular genome (counts + 1)", y = "Diploid genome (counts + 1)",
       subtitle = "Transcript-level quantification",
       caption = "Both axis on a log-scale")

# compare two alleles
inner_join(regular, diploid, by = "gene") %>%
  ggplot(aes(B + 1, C + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()


# Summarise per Gene ------------------------------------------------------

# transcript-to-gene
tx2gene <- read_tsv("data/external/reference/transcript2gene.tsv") %>%
  # ensure columns are in the correct order for tximport
  select(transcript, gene)

tx2gene_diploid <- tx2gene %>%
  mutate(transcript = paste0("C_", transcript),
         gene = paste0("C_", gene))
tx2gene_diploid <- tx2gene %>%
  mutate(transcript = paste0("B_", transcript),
         gene = paste0("B_", gene)) %>%
  bind_rows(tx2gene_diploid)

tx2gene_diploid %>%
  mutate(allele = transcript) %>%
  mutate(transcript = str_remove(transcript, "^._"),
        gene = str_remove(gene, "^._")) %>%
  select(gene, transcript, allele) %>%
  write_tsv("data/external/strains/allele_to_gene_map.tsv", col_names = FALSE)

library(tximport)
regular_gene <- tximport("results/salmon/stromal_t1_bc_M00323192/quant.sf",
                 tx2gene = tx2gene, type = "salmon") %>%
  pluck("counts") %>% as_tibble(rownames = "gene") %>%
  rename(counts = V1)

diploid_gene <- tximport("results/salmon_diploid/stromal_t1_bc_M00323192/quant.sf",
                  tx2gene = tx2gene_diploid, type = "salmon") %>%
  pluck("counts") %>% as_tibble(rownames = "gene") %>%
  rename(counts = V1) %>%
  separate(gene, c("strain", "gene")) %>%
  pivot_wider(names_from = strain, values_from = counts)

# compare the two approaches at the gene-level - very good correlation
diploid_gene %>%
  inner_join(regular_gene, by = "gene") %>%
  ggplot(aes(counts + 1, B + C + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown", size = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Regular genome (counts + 1)", y = "Diploid genome (counts + 1)",
       subtitle = "Gene-level quantification",
       caption = "Both axis on a log-scale")

# compare two allele's quantification
diploid_gene %>%
  ggplot(aes(B + 1, C + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()


# Compare with SNPsplit ---------------------------------------------------

# get SNP-split counts (from Russell)
snpsplit1 <- read_tsv("results/SLX-18042.D701rna_D505rna.H3YKJDSXY.s_4.r_1_val_1_CAST_EiJ_N-masked_hisat2_srtd.genome1.bam_featureCounts_counts.txt",
                     comment = "#")
snpsplit2 <- read_tsv("results/SLX-18042.D701rna_D505rna.H3YKJDSXY.s_4.r_1_val_1_CAST_EiJ_N-masked_hisat2_srtd.genome2.bam_featureCounts_counts.txt",
                      comment = "#")
snpsplit <- tibble(gene = snpsplit1$Geneid, B = snpsplit1[[7]], C = snpsplit2[[7]])
rm(snpsplit1, snpsplit2)

# the two alleles with SNP-split
snpsplit %>%
  ggplot(aes(B + 1, C + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

# compare SNP-split to SALMON
snpsplit %>%
  inner_join(diploid_gene, by = "gene", suffix = c("_split", "_salmon")) %>%
  ggplot(aes(B_split + 1, B_salmon + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

# ratio of counts
snpsplit %>%
  inner_join(diploid_gene, by = "gene", suffix = c("_split", "_salmon")) %>%
  ggplot(aes(log2((C_split + 1) / (B_split + 1)),
             log2((C_salmon + 1) / (B_salmon + 1)))) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown")

# MA plots
snpsplit %>%
  inner_join(diploid_gene, by = "gene", suffix = c("_split", "_salmon")) %>%
  ggplot(aes(log2((C_split + 1 + B_split + 1)/2),
             log2((C_split + 1) / (B_split + 1)))) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "brown")

snpsplit %>%
  inner_join(diploid_gene, by = "gene", suffix = c("_salmon", "_salmon")) %>%
  ggplot(aes(log2((C_salmon + 1 + B_salmon + 1)/2),
             log2((C_salmon + 1) / (B_salmon + 1)))) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "brown")



# Compare with feature counts ---------------------------------------------

featurects <- read_tsv("../check_strandness/results/featureCounts/rev_includedup",
                       comment = "#")
featurects <- tibble(gene = featurects$Geneid, cts = featurects[[7]])

#
snpsplit %>%
  inner_join(featurects) %>%
  ggplot(aes(B + C + 1, cts)) +
  geom_point() +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

regular_gene %>%
  inner_join(featurects) %>%
  ggplot(aes(counts, cts)) +
  geom_point() +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

diploid_gene %>%
  inner_join(featurects) %>%
  ggplot(aes(B + C + 1, cts)) +
  geom_point() +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()



