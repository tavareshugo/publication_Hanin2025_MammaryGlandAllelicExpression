# I'm a bad person
setwd("~/mount/hpc_uni/mammary_gland_transcriptomes/test_pipeline/")

library(tximport)
library(tidyverse)
theme_set(theme_minimal(base_size = 16))

# transcript-to-gene
tx2gene <- read_tsv("data/external/reference/transcript2gene.tsv") %>%
  # ensure columns are in the correct order for tximport
  select(transcript, gene)


# Prepare tx2gene tables --------------------------------------------------

# make a diploid version of the table
tx2gene_diploid <- tx2gene %>%
  mutate(transcript = paste0("C_", transcript),
         gene = paste0("C_", gene))
tx2gene_diploid <- tx2gene %>%
  mutate(transcript = paste0("B_", transcript),
         gene = paste0("B_", gene)) %>%
  bind_rows(tx2gene_diploid)

# tx2gene_diploid %>%
#   mutate(allele = transcript) %>%
#   mutate(transcript = str_remove(transcript, "^._"),
#          gene = str_remove(gene, "^._")) %>%
#   select(gene, transcript, allele) %>%
#   write_tsv("data/external/strains/allele_to_gene_map.tsv", col_names = FALSE)


# Transcript quantification ----------------------------------

# ground truth
truth_isoform <- read_tsv("results/rsem-simulate-reads/stromal_t1_bc_M00323192.sim.isoforms.results")
truth_isoform <- truth_isoform %>% select(transcript = transcript_id, truth = count)

# regular salmon quantification
regular_isoform <- read_tsv("results/salmon_sim/stromal_t1_bc_M00323192/quant.sf")
regular_isoform <- regular_isoform %>% select(transcript = Name, regular = NumReads)

# quantification against diploid transcriptome
# we map each transcript allele to its transcript ID to get tximport summaries
diploid_isoform <- tximport("results/salmon_diploid_sim/stromal_t1_bc_M00323192/quant.sf",
                            tx2gene = tx2gene_diploid %>%
                              mutate(gene = str_remove(transcript, "^._")),
                            type = "salmon") %>%
  pluck("counts") %>% as_tibble(rownames = "transcript") %>%
  rename(diploid = V1)

isoform <- Reduce(full_join, list(truth_isoform, regular_isoform, diploid_isoform))
rm(regular_isoform, diploid_isoform, truth_isoform)

# compare regular salmon quantification (against standard genome) with diploid quantification
isoform %>%
  ggplot(aes(regular + 1, diploid + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Regular genome (counts + 1)", y = "Diploid genome (counts + 1)",
       subtitle = "Transcript-level quantification",
       caption = "Both axis on a log-scale")

# compare difference to ground truth
isoform %>%
  ggplot(aes(log2((regular+1)/(truth+1)), log2((diploid+1)/(truth+1)))) +
  geom_point(alpha = 0.1, aes(colour = truth == 0)) +
  geom_abline(colour = "brown")

# estimated vs truth
isoform %>%
  pivot_longer(regular:diploid) %>%
  mutate(lfc = log2((value+1)/(truth + 1)),
         class = case_when(truth == 0 & value == 0 ~ "truth & estimated == 0",
                           truth == 0 & value != 0 ~ "truth == 0",
                           truth != 0 & value == 0 ~ "estimated == 0",
                           TRUE ~ "count > 0")) %>%
  ggplot(aes(log2(truth + 1), lfc)) +
  geom_point(aes(colour = class), alpha = 0.1) +
  facet_grid(name ~ .) +
  scale_colour_manual(values = c("black", "brown", "steelblue", "orange2")) +
  labs(y = "log2(estimated/truth)")

# diploid approach vs ground truth
isoform %>%
  ggplot(aes(truth + 1, diploid + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

# regular quantification vs truth
isoform %>%
  ggplot(aes(truth + 1, regular + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

# The diploid genome seems to do better at zero counts
isoform %>%
  filter(truth == 0 & regular != 0 & diploid == 0) %>% nrow()
isoform %>%
  filter(truth == 0 & regular == 0 & diploid != 0) %>% nrow()

# The diploid genome seems to do better at zero counts
isoform %>%
  filter(truth != 0 & regular == 0 & diploid != 0) %>% nrow()
isoform %>%
  filter(truth != 0 & regular != 0 & diploid == 0) %>% nrow()


# Summarise per Gene ------------------------------------------------------

# ground truth
truth_gene <- read_tsv("results/rsem-simulate-reads/stromal_t1_bc_M00323192.sim.genes.results")
truth_gene <- truth_gene %>% select(gene = gene_id, truth = count)

# regular salmon quantification
regular_gene <- tximport("results/salmon_sim/stromal_t1_bc_M00323192/quant.sf",
                         tx2gene = tx2gene, type = "salmon") %>%
  pluck("counts") %>% as_tibble(rownames = "gene") %>%
  rename(regular = V1)

# quantification against diploid transcriptome
# we map each transcript allele to its transcript ID to get tximport summaries
diploid_gene <- tximport("results/salmon_diploid_sim/stromal_t1_bc_M00323192/quant.sf",
                            tx2gene = tx2gene_diploid %>%
                              mutate(gene = str_remove(gene, "^._")),
                            type = "salmon") %>%
  pluck("counts") %>% as_tibble(rownames = "gene") %>%
  rename(diploid = V1)

gene <- Reduce(full_join, list(truth_gene, regular_gene, diploid_gene))
rm(regular_gene, diploid_gene, truth_gene)

# compare regular salmon quantification (against standard genome) with diploid quantification
gene %>%
  ggplot(aes(regular + 1, diploid + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Regular genome (counts + 1)", y = "Diploid genome (counts + 1)",
       subtitle = "Gene-level quantification",
       caption = "Both axis on a log-scale")

# compare difference to ground truth
gene %>%
  ggplot(aes(log2((regular+1)/(truth+1)), log2((diploid+1)/(truth+1)))) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown")

# estimated vs truth
gene %>%
  pivot_longer(regular:diploid) %>%
  mutate(lfc = log2((value+1)/(truth + 1)),
         class = case_when(truth == 0 & value == 0 ~ "truth & estimated == 0",
                           truth == 0 & value != 0 ~ "truth == 0",
                           truth != 0 & value == 0 ~ "estimated == 0",
                           TRUE ~ "count > 0")) %>%
  ggplot(aes(log2(truth + 1), lfc)) +
  geom_point(aes(colour = class), alpha = 0.1) +
  facet_grid(name ~ .) +
  scale_colour_manual(values = c("black", "brown", "steelblue", "orange2")) +
  labs(y = "log2(estimated/truth)")

# diploid approach vs ground truth
gene %>%
  ggplot(aes(truth + 1, diploid + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

# regular quantification vs truth
gene %>%
  ggplot(aes(truth + 1, regular + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

# The diploid genome seems to do better at zero counts
gene %>%
  filter(truth == 0 & regular != 0 & diploid == 0) %>% nrow()
gene %>%
  filter(truth == 0 & regular == 0 & diploid != 0) %>% nrow()

# The diploid genome seems to do better at zero counts
gene %>%
  filter(truth != 0 & regular == 0 & diploid != 0) %>% nrow()
gene %>%
  filter(truth != 0 & regular != 0 & diploid == 0) %>% nrow()



# Transcript allele specific ----------------------------------------------

# ground truth
truth_transcript_allele <- read_tsv("results/rsem-simulate-reads/stromal_t1_bc_M00323192.sim.alleles.results")
truth_transcript_allele <- truth_transcript_allele %>% select(allele = allele_id,
                                                              truth = count)

# quantification against diploid transcriptome
# we map each transcript allele to its transcript ID to get tximport summaries
diploid_transcript_allele <- read_tsv("results/salmon_diploid_sim/stromal_t1_bc_M00323192/quant.sf") %>%
  select(allele = Name, diploid = NumReads)

transcript_allele <- full_join(truth_transcript_allele, diploid_transcript_allele)
rm(diploid_transcript_allele, truth_transcript_allele)

# diploid approach vs ground truth
transcript_allele %>%
  ggplot(aes(truth + 1, diploid + 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

transcript_allele %>%
  ggplot(aes(log2((diploid+truth+2)/2), log2((diploid+1)/(truth+1)))) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "brown")

transcript_allele %>%
  ggplot(aes(log2((diploid+1)/(truth+1)))) +
  geom_density()


# gene allele specific ----------------------------------------------

# ground truth
truth_gene_allele <- read_tsv("results/rsem-simulate-reads/stromal_t1_bc_M00323192.sim.alleles.results")
truth_gene_allele <- truth_gene_allele %>%
  separate(allele_id, c("strain", "transcript"), sep = "_") %>%
  mutate(gene = paste(strain, gene_id, sep = "_")) %>%
  select(gene, count)
truth_gene_allele <- truth_gene_allele %>%
  group_by(gene) %>%
  summarise(truth = sum(count)) %>%
  ungroup()

# quantification against diploid geneome
# we map each gene allele to its gene ID to get tximport summaries
diploid_gene_allele <- tximport("results/salmon_diploid_sim/stromal_t1_bc_M00323192/quant.sf",
                                tx2gene = tx2gene_diploid,
                                type = "salmon") %>%
  pluck("counts") %>% as_tibble(rownames = "gene") %>%
  rename(diploid = V1)

gene_allele <- full_join(truth_gene_allele, diploid_gene_allele)
rm(diploid_gene_allele, truth_gene_allele)

# diploid approach vs ground truth
gene_allele %>%
  ggplot(aes(truth, diploid)) +
  geom_point(alpha = 0.1) +
  geom_abline(colour = "brown") +
  scale_x_log10() + scale_y_log10()

gene_allele %>%
  ggplot(aes(log2(truth+1), log2((diploid+1)/(truth+1)))) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "brown")

gene_allele %>%
  ggplot(aes(log2((diploid+1)/(truth+1)))) +
  geom_density()

gene_allele %>%
  separate(gene, c("strain", "gene"), sep = "_") %>%
  pivot_wider(names_from = strain, values_from = c(truth, diploid)) %>%
  ggplot(aes(log2((truth_B+1)/(truth_C+1)), log2((diploid_B+1)/(diploid_C+1)))) +
  geom_point(alpha = 0.1)

# % of allele estimate
gene_allele %>%
  separate(gene, c("strain", "gene"), sep = "_") %>%
  group_by(gene) %>%
  summarise(truth_fraction = (truth[strain == "B"])/sum(truth),
            total_truth = sum(truth),
            diploid_fraction = (diploid[strain == "B"])/sum(diploid),
            total_diploid = sum(diploid)) %>%
  full_join(snps_in_genes, by = "gene") %>%
  ggplot(aes(truth_fraction, diploid_fraction)) +
  geom_point(aes(colour = n_variants)) +
  geom_abline(colour = "grey") +
  scale_colour_viridis_c(trans = "log10") +
  labs(subtitle = "Allele ratio (gene-level)", colour = "# Variants")

gene_allele %>%
  separate(gene, c("strain", "gene"), sep = "_") %>%
  group_by(gene) %>%
  summarise(truth_fraction = (truth[strain == "B"])/sum(truth),
            total_truth = sum(truth),
            diploid_fraction = (diploid[strain == "B"])/sum(diploid),
            total_diploid = sum(diploid)) %>%
  full_join(snps_in_genes, by = "gene") %>%
  filter(n_variants == 0 & (near(diploid_fraction, 0) | near(diploid_fraction, 1))) %>%
  filter(total_truth > 0)

# distribution % of allele estimate
gene_allele %>%
  separate(gene, c("strain", "gene"), sep = "_") %>%
  group_by(gene) %>%
  summarise(truth_fraction = (truth[strain == "B"])/sum(truth),
            total_truth = sum(truth),
            diploid_fraction = (diploid[strain == "B"])/sum(diploid),
            total_diploid = sum(diploid)) %>%
  full_join(snps_in_genes, by = "gene") %>%
  pivot_longer(c(truth_fraction, diploid_fraction)) %>%
  mutate(name = ifelse(name == "diploid_fraction", "Estimated", "Truth")) %>%
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.1) +
  labs(x = "Allelic Fraction", fill = "")

gene_allele %>%
  separate(gene, c("strain", "gene"), sep = "_") %>%
  group_by(gene) %>%
  summarise(truth_fraction = (truth[strain == "B"])/sum(truth),
            total_truth = sum(truth),
            diploid_fraction = (diploid[strain == "B"])/sum(diploid),
            total_diploid = sum(diploid)) %>%
  full_join(snps_in_genes, by = "gene") %>%
  filter(n_variants == 0) %>%
  pivot_longer(c(truth_fraction, diploid_fraction)) %>%
  ggplot(aes(value, fill = name)) +
  geom_density(alpha = 0.1) +
  labs(x = "Allelic Fraction")

# Get SNPs ----------------------------------------------------------------

library(valr)
# read SNPs and Indels
snps <- vroom::vroom("data/external/strains/B6xCAST.snps_and_indels.merged.bed.gz",
                     col_names = c("chrom", "start", "end", "name"),
                     col_types = "ciic")
snps <- snps %>%
  mutate(type = case_when(str_detect(name, "SNP") & !str_detect(name, "INDEL") ~ "SNP",
                          str_detect(name, "INDEL") ~ "INDEL",
                          TRUE ~ NA_character_))

# read annotation
transcript_annot <- vroom::vroom("data/external/reference/transcripts.bed.gz",
                                 col_names = c("chrom", "start", "end", "transcript"),
                                 col_types = "ciic")
transcript_annot <- transcript_annot %>%
  full_join(tx2gene, by = "transcript")

gene_annot <- transcript_annot %>%
  group_by(gene) %>%
  bed_merge()

# intersect genes and variants
snps_in_genes <- gene_annot %>%
  group_by(gene) %>%
  bed_intersect(snps, suffix = c("", "_snp")) %>%
  group_by(gene) %>%
  summarise(n_variants = n(),
            n_snps = sum(type_snp == "SNP"),
            n_indels = sum(type_snp == "INDEL")) %>%
  ungroup() %>%
  # join to complete gene table to add "zeros"
  full_join(tx2gene %>% distinct(gene), by = "gene") %>%
  replace_na(list(n_variants = 0, n_snps = 0, n_indels = 0))

# fraction of genes with variants
sum(snps_in_genes$n_variants > 0)/nrow(snps_in_genes)

# distribution of # variants per gene
snps_in_genes %>%
  count(n_variants, name = "count") %>%
  ggplot(aes(n_variants, count/sum(count))) +
  geom_segment(aes(xend = n_variants, y = 0, yend = count/sum(count)),
               colour = "grey") +
  geom_point(colour = "brown", size = 3) +
  scale_x_log10() +
  labs(x = "# Variants per Gene", y = "Fraction")

snps_in_genes %>%
  pivot_longer(c(n_snps, n_indels)) %>%
  count(name, value, name = "count") %>%
  filter(value > 0) %>%
  ggplot(aes(value, count/sum(count))) +
  geom_line(aes(group = name), colour = "grey") +
  geom_point(aes(colour = name), size = 2) +
  scale_x_log10() +
  labs(x = "# Variants per Gene", y = "Fraction")

# number of genes with each type of variant
sum(snps_in_genes$n_indels > 0 & snps_in_genes$n_snps > 0)  # with indel and SNP
sum(snps_in_genes$n_indels > 0 & snps_in_genes$n_snps == 0) # with indel but no SNP
sum(snps_in_genes$n_snps > 0 & snps_in_genes$n_indels == 0) # with SNP but no indel
sum(snps_in_genes$n_variants == 0) # no variants

# distribution size of indels
# actually this is relative to the reference, so this is the size of deletions in CAST
# and all the ones where difference == 1 are insertions in CAST
snps %>%
  filter(type == "INDEL") %>%
  ggplot(aes(end - start)) +
  geom_bar()
