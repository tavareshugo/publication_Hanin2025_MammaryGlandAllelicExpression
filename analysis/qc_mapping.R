# Analysis of QC metrics

library(tidyverse)
library(patchwork)
theme_set(theme_minimal() + theme(text = element_text(size = 16)))

cell_colours <- RColorBrewer::brewer.pal(6, "Dark2")
names(cell_colours) <- c("adipocytes", "basal", "endothelial", "luminal differentiated", "luminal progenitors", "stromal")



# Read data ---------------------------------------------------------------

sample_info <- read_csv("read_info.csv")

# summarise to include only categories of interest for our QC
sample_info <- sample_info %>%
  count(sample, cell_type, stage, cross, animal_id, rin_category, rin)

# multiqc data
multiqc <- jsonlite::read_json("results/qc/multiqc_report_data/multiqc_data.json")


# Percent Duplicates ---------------------------------------

# parse multiqc
pct_dup <- multiqc %>%
  pluck("report_general_stats_data") %>%
  pluck(2) %>%
  map_dbl(~ .x$PERCENT_DUPLICATION) %>%
  enframe(name = "sample", value = "pct_dup")

pct_dup %>%
  full_join(sample_info, by = "sample") %>%
  ggplot(aes(pct_dup*100)) +
  geom_density(aes(colour = cell_type)) +
  # geom_density(fill = "steelblue", alpha = 0.3) +
  geom_rug(alpha = 0.2) +
  scale_colour_manual(values = cell_colours) +
  labs(x = "% duplicates",
       caption = "from Picard software", colour = "Cell")


# percent mapped (salmon) ----------------------------------

# parse multiqc
pct_mapped <- multiqc %>%
  pluck("report_general_stats_data") %>%
  pluck(3) %>%
  map_dbl(~ .x$percent_mapped) %>%
  enframe(name = "sample", value = "pct_mapped")

pct_mapped %>%
  full_join(sample_info, by = "sample") %>%
  ggplot(aes(pct_mapped)) +
  # geom_density(fill = "steelblue", alpha = 0.3) +
  geom_rug(alpha = 0.2) +
  geom_density(aes(colour = cell_type)) +
  scale_colour_manual(values = cell_colours) +
  labs(x = "% mapped",
       caption = "Using salmon software against 'diploidised' transcriptome")

# also get pct mapped from STAR (just to get an idea of the difference)
pct_mapped_genome <- list.files("results/star/", pattern = "Log.final.out", recursive = TRUE, full.names = TRUE)
names(pct_mapped_genome) <- pct_mapped_genome %>%
  basename() %>% str_remove(".Log.final.out")
pct_mapped_genome <- pct_mapped_genome %>%
  map_dfr(function(i){
    read_table(i, skip = 9, n_max = 1, col_names = FALSE, col_types = cols()) %>%
      select(pct_mapped = X6) %>%
      mutate(pct_mapped = as.numeric(str_remove(pct_mapped, "%")))
  }, .id = "sample")

pct_mapped_genome %>%
  full_join(sample_info, by = "sample") %>%
  ggplot(aes(pct_mapped)) +
  # geom_density(fill = "steelblue", alpha = 0.3) +
  geom_rug(alpha = 0.2) +
  geom_density(aes(colour = cell_type)) +
  scale_colour_manual(values = cell_colours) +
  labs(x = "% mapped",
       caption = "Using STAR software against B6 reference genome")


# percent rRNA -------------------------

# parse multiqc
pct_rrna <- multiqc %>%
  pluck("report_general_stats_data") %>%
  pluck(5) %>%
  map_dbl(~ .x$rRNA_pct) %>%
  enframe(name = "sample", value = "pct_rrna")

pct_rrna %>%
  full_join(sample_info, by = "sample") %>%
  ggplot(aes(pct_rrna)) +
  # geom_density(fill = "steelblue", alpha = 0.3) +
  geom_rug(alpha = 0.2) +
  geom_density(aes(colour = cell_type)) +
  scale_colour_manual(values = cell_colours) +
  labs(x = "% rRNA",
       caption = "Using rRNAs from 'rnacentral.org' database")


# Gene body coverage ------------------------------------------------------

# tidy data from json
gene_body_cov <- multiqc %>%
  pluck("report_data_sources") %>%
  pluck("RSeQC") %>%
  pluck("gene_body_coverage") %>%
  map_dfr(~ read_tsv(paste0("results/qc/rseqc/", basename(.x)),
                     col_types = cols(.default = col_double(), Percentile = col_character()))) %>%
  pivot_longer(-Percentile, names_to = "percentile", values_to = "value") %>%
  rename(sample = Percentile) %>%
  mutate(percentile = as.numeric(percentile),
         sample = str_remove(sample, "\\..*")) %>%
  group_by(sample) %>%
  mutate(scaled = (value - min(value))/(max(value) - min(value))) %>%
  ungroup()

# plot of gene body coverage (partitioned by RIN categorisation)
p <- gene_body_cov %>%
  left_join(sample_info, by = "sample") %>%
  mutate(rin_category = factor(rin_category,
                               levels = c("lo", "med lo", "med hi", "hi"))) %>%
  ggplot(aes(percentile, scaled, group = sample)) +
  geom_line(aes(colour = cell_type)) +
  facet_grid(~ rin_category) +
  scale_colour_manual(values = cell_colours) +
  labs(x = "Gene body percentile (5' -> 3')",
       y = "Scaled fraction coverage",
       colour = "Cell",
       caption = "Split by RIN category.")
plotly::ggplotly(p)

# problematic samples
gene_cov_fail <- gene_body_cov %>%
  group_by(sample) %>%
  filter(scaled > .9) %>%
  filter(percentile < 15) %>%
  ungroup() %>%
  distinct(sample) %>%
  mutate(genecov_fail = TRUE)

# plot just those
gene_cov_fail %>%
  left_join(gene_body_cov) %>%
  left_join(sample_info, by = "sample") %>%
  ggplot(aes(percentile, scaled, group = sample)) +
  geom_line(alpha = 0.8)

# how many per group
gene_body_cov %>%
  group_by(sample) %>%
  filter(scaled > .9) %>%
  filter(percentile < 15) %>%
  ungroup() %>%
  distinct(sample) %>%
  left_join(sample_info, by = "sample") %>%
  count(cell_type, stage, cross) %>%
  ggplot(aes(cell_type, stage)) +
  geom_label(aes(label = n)) +
  facet_grid(~ cross) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# relationship with RIN
gene_body_cov %>%
  group_by(sample) %>%
  filter(scaled == max(scaled)) %>%
  summarise(percentile = min(percentile)) %>%
  left_join(sample_info, by = "sample") %>%
  mutate(rin = ifelse(is.na(rin), 0, rin)) %>%
  ggplot(aes(percentile, rin)) +
  geom_point(aes(colour = cell_type)) +
  scale_colour_manual(values = cell_colours)



# Read distribution -------------------------------------------------------

read_dist <- multiqc %>%
  pluck("report_data_sources") %>%
  pluck("RSeQC") %>%
  pluck("read_distribution") %>%
  map_dfr(~ read_table(paste0("results/qc/rseqc/", basename(.x)),
                     skip = 4, n_max = 10,
                     col_types = cols(.default = col_double(), Group = col_character())),
          .id = "sample") %>%
  select(-X5) %>%
  janitor::clean_names() %>%
  mutate(sample = str_remove(sample, "\\..*")) %>%
  group_by(sample) %>%
  mutate(tag_pct = tag_count/sum(tag_count)) %>%
  ungroup() %>%
  mutate(group = factor(group,
                        levels = c("TSS_up_10kb", "TSS_up_5kb", "TSS_up_1kb",
                                   "5'UTR_Exons", "CDS_Exons", "3'UTR_Exons",
                                   "Introns",
                                   "TES_down_1kb", "TES_down_5kb", "TES_down_10kb")))

read_dist %>%
  left_join(sample_info) %>%
  ggplot(aes(tag_pct, sample)) +
  geom_col(aes(fill = group)) +
  theme(axis.text.y = element_blank())

# exon pct
read_dist %>%
  left_join(sample_info) %>%
  filter(str_detect(group, "Exons")) %>%
  group_by(sample, cell_type) %>%
  summarise(tag_pct = sum(tag_pct)) %>%
  ungroup() %>%
  mutate(sample = fct_reorder(sample, tag_pct)) %>%
  ggplot(aes(tag_pct, sample)) +
  geom_col(aes(fill = cell_type)) +
  scale_fill_manual(values = cell_colours) +
  theme(axis.text.y = element_blank())

# intron pct
read_dist %>%
  filter(str_detect(group, "Intron")) %>%
  mutate(sample = fct_reorder(sample, tag_pct)) %>%
  ggplot(aes(tag_pct, sample)) +
  geom_col() +
  theme(axis.text.y = element_blank())

# relationship between mapped and intron reads
pct_mapped_genome %>%
  rename(genome = pct_mapped) %>%
  full_join(pct_mapped, by = "sample") %>%
  left_join(sample_info, by = "sample") %>%
  rename(transcriptome = pct_mapped) %>%
  mutate(dif = transcriptome/genome) %>%
  full_join(read_dist %>% filter(str_detect(group, "Intron")), by = "sample") %>%
  ggplot(aes(tag_pct, dif)) +
  geom_point(aes(colour = cell_type)) +
  scale_colour_manual(values = cell_colours) +
  labs(x = "% reads overlapping introns",
       y = "Aligned ratio (transcriptome / genome)",
       colour = "Cell")

# look at mapping rates between genome or transcriptome
pct_mapped_genome %>%
  rename(genome = pct_mapped) %>%
  full_join(pct_mapped, by = "sample") %>%
  rename(transcriptome = pct_mapped) %>%
  pivot_longer(c(genome, transcriptome)) %>%
  full_join(sample_info, by = "sample") %>%
  ggplot(aes(value)) +
  geom_density(aes(fill = name), alpha = 0.2) +
  geom_rug(alpha = 0.2, aes(colour = name)) +
  labs(x = "% mapped",
       colour = "Aligned to", fill = "Aligned to",
       caption = "Using STAR software against B6 reference genome and salmon software against 'diploidised' transcriptome.")



# Save summary QC table -------------------------------------------------

read_dist %>%
  filter(str_detect(group, "Intron")) %>%
  distinct(sample, pct_intron = tag_pct) %>%
  full_join(pct_dup, by = "sample") %>%
  full_join(pct_mapped, by = "sample") %>%
  full_join(pct_rrna, by = "sample") %>%
  full_join(gene_cov_fail, by = "sample") %>%
  mutate(genecov_fail = ifelse(is.na(genecov_fail), FALSE, genecov_fail)) %>%
  write_csv("results/qc/mapping_qc_summary.csv")
