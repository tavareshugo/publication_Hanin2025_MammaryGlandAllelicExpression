library(DESeq2)
library(tidyverse)
library(UpSetR)
library(patchwork)
theme_set(theme_minimal(base_size = 16))

# vector of colours for plotting developmental stages
stage_colours <- c("nulliparous" = "grey",
                   "gestation d5.5" = "#a6cee3",
                   "gestation d9.5" = "#9ecae1",
                   "gestation d14.5" = "#3182bd",
                   "lactation d5" = "#e5f5e0",
                   "lactation d10" = "#a1d99b",
                   "lactation d15" = "#31a354",
                   "involution d1" = "#fee6ce",
                   "involution d6" = "#fdae6b",
                   "involution d14" = "#e6550d")

# colours for cells
cell_colours <- RColorBrewer::brewer.pal(6, "Dark2")
names(cell_colours) <- c("adipocytes", "basal", "endothelial", "luminal differentiated", "luminal progenitors", "stromal")


# Read Data ---------------------------------------------------------------

# deseq dataset
dds <- readRDS("results/DESeqDataSet/dds_gene_regular.rds")

# fetch gene annotation
annot <- dds %>% rowData() %>% as_tibble(rownames = "gene") %>%
  mutate(external_gene_name = ifelse(is.na(external_gene_name), gene, external_gene_name))

# imprinted genes from Tucci et al
tucci <- readxl::read_excel("data/external/tucci_et_al_sup1_imprinted_genes.xlsx",
                            sheet = "Mouse_All_Known IG_NoStatus")

# imprinted genes from Xu et al
xu <- read_csv("data/external/xu_et_al_imprints.csv")

# other known imprints
imprints <- c("DLK1", "MEST", "IGF2R", "IGF2", "GRB10", "GNAS", "UBE3A", "CDKN1C", "SNRPN", "NNAT", "PEG3", "COPG")

# isolde data
isolde <- read_csv("results/ase_isolde/isolde_all.csv") %>%
  rename(gene = name) %>%
  mutate(cell_type = case_when(cell_type == "luminalp" ~ "luminal progenitors",
                               cell_type == "luminald" ~ "luminal differentiated",
                               TRUE ~ cell_type)) %>%
  left_join(annot, by = "gene") %>%
  left_join(dds %>% colData() %>% as_tibble() %>%
              count(cell_type, timepoint, stage) %>% rename(n_reps = n),
            by = c("cell_type", "timepoint")) %>%
  mutate(chromosome_name = factor(chromosome_name, levels = c(1:19, "X", "Y")))

# isolde strain bias analysis
strain <- read_csv("results/ase_isolde/isolde_all_strain.csv") %>%
  rename(gene = name) %>%
  mutate(cell_type = case_when(cell_type == "luminalp" ~ "luminal progenitors",
                               cell_type == "luminald" ~ "luminal differentiated",
                               TRUE ~ cell_type)) %>%
  left_join(annot, by = "gene") %>%
  left_join(dds %>% colData() %>% as_tibble() %>%
              count(cell_type, timepoint, stage) %>% rename(n_reps = n),
            by = c("cell_type", "timepoint")) %>%
  mutate(chromosome_name = factor(chromosome_name, levels = c(1:19, "X", "Y")))


# chromosome sizes (for plotting ideograms)
chrom_sizes <- read_tsv("data/external/chrom_sizes.tsv") %>%
  rename(chromosome_name = chrom) %>%
  mutate(start_position = 0) %>%
  filter(chromosome_name != "MT" & chromosome_name != "Y") %>%
  mutate(chromosome_name = factor(chromosome_name, levels = c(1:19, "X", "Y")))

# IGRs from Andergassen et al
igrs <- readxl::read_excel("data/external/elife-25125-supp1-v2.xlsx",
                           sheet = "K", range = "A5:A32",
                           col_names = "external_gene_name")

# get the coordinates
# PWS/AS is the Prader-Willi syndrome/Angelman syndrome region near Snrpn
igrs <- igrs %>%
  mutate(external_gene_name = ifelse(external_gene_name == "Pws/As", "Snrpn", external_gene_name)) %>%
  left_join(annot %>% distinct(external_gene_name, chromosome_name, start_position)) %>%
  mutate(chromosome_name = factor(chromosome_name, levels = c(1:19, "X", "Y")))

# get normalised gene expression - only for genes that isolde analysed
expr <- dds[unique(c(isolde$gene, strain$gene)), ] %>%
  assay("vst") %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(dds %>% colData() %>% as_tibble(rownames = "sample"), by = "sample")

summarised_expr <- expr %>%
  group_by(gene, cell_type, stage, timepoint) %>%
  summarise(mean_expr = mean(expr), median_expr = median(expr))

# allele-specific counts
dds_isolde <- readRDS("results/DESeqDataSet/dds_gene_allele.rds")

# estimate size factors for normalisation (ISoLDE recommends doing this,
# although I wonder if size factors should be estimated at the sample level (rather than allele-specific level))
dds_isolde <- estimateSizeFactors(dds_isolde)



# ISOlDE ------------------------------------------------------------------

# distribution of allele difference - all pooled together
wrap_plots(
  isolde %>%
    ggplot(aes(diff_prop)) +
    geom_histogram(binwidth = 0.1) +
    labs(x = "", subtitle = "All genes"),
  isolde %>% filter(status == "ASE") %>%
    ggplot(aes(diff_prop)) +
    geom_histogram(binwidth = 0.1) +
    labs(x = "Difference (M - P)", subtitle = "ASE genes"),
  ncol = 1
) &
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

# distribution per condition
isolde %>%
  ggplot(aes(diff_prop)) +
  geom_density(aes(colour = cell_type)) +
  facet_wrap(~ stage) +
  scale_colour_manual(values = cell_colours) +
  labs(x = "Difference (M - P)", colour = "Cell")


# a kind of MA plot
isolde %>%
  left_join(summarised_expr, by = c("gene", "cell_type", "stage", "timepoint")) %>%
  arrange(status == "ASE") %>%
  ggplot(aes(mean_expr, diff_prop)) +
  geom_point(aes(colour = status == "ASE"), size = 1) +
  facet_grid(cell_type ~ stage) +
  scale_colour_manual(values = c("FALSE" = "grey", "TRUE" = "black")) +
  labs(x = "Mean expression", y = "Difference (M - P)", colour = "ASE")

# overlap
isolde %>%
  filter(status == "ASE") %>%
  mutate(group = paste(cell_type, stage, sep = " | ")) %>%
  with(split(gene, group)) %>%
  fromList() %>%
  upset(order.by = "freq", nsets = 10)

# barplot
isolde %>%
  filter(status == "ASE") %>%
  count(gene, external_gene_name) %>%
  count(n) %>%
  ggplot(aes(n, nn)) +
  geom_segment(aes(xend = n, y = 0, yend = nn), colour = "grey48", size = 1) +
  geom_point(colour = "steelblue", size = 3) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  labs(x = "Number of conditions a gene has ASE",
       y = "count", caption = "All genes flagged as ASE by ISOlDE")

# Number of ASE genes
wrap_plots(
  isolde %>%
    filter(status == "ASE") %>%
    count(gene, external_gene_name) %>%
    count(n) %>%
    ggplot(aes(n, nn)) +
    geom_segment(aes(xend = n, y = 0, yend = nn), colour = "grey48", size = 1) +
    geom_point(colour = "steelblue", size = 3) +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    labs(x = "",
         y = "count", subtitle = "All genes flagged as ASE by ISOlDE"),
  isolde %>%
    filter(status == "ASE" & abs(diff_prop) > 0.7) %>%
    count(gene, external_gene_name) %>%
    count(n) %>%
    ggplot(aes(n, nn)) +
    geom_segment(aes(xend = n, y = 0, yend = nn), colour = "grey48", size = 1) +
    geom_point(colour = "steelblue", size = 3) +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    labs(x = "Number of conditions a gene has ASE",
         y = "count", subtitle = "Genes with > 0.7 difference"),
  ncol = 1
)


# Gene biotype ---------------

# Count ASE genes per biotype
biotype_counts <- isolde %>%
  group_by(gene, gene_biotype) %>%
  summarise(ase = any(status == "ASE" & abs(diff_prop) > 0.7)) %>%
  ungroup()

# Visualise counts
biotype_counts %>%
  group_by(gene_biotype) %>%
  summarise(n_ase = sum(ase),
            n = n(),
            pct_ase = mean(ase)*100) %>%
  mutate(gene_biotype = gene_biotype %>% str_replace_all("_", " ")) %>%
  ggplot(aes(pct_ase, gene_biotype)) +
  geom_col()

# permutation to get 95% expected intervals
permutation_test <- map_dfr(1:500, function(i){
  biotype_counts %>%
    # shuffle
    mutate(ase = sample(ase)) %>%
    group_by(gene_biotype) %>%
    summarise(n_ase = sum(ase),
              n = n(),
              pct_ase = mean(ase)*100) %>%
    ungroup() %>% mutate(rep = i)
})

permutation_summary <- permutation_test %>%
  group_by(gene_biotype) %>%
  summarise(lo = quantile(pct_ase, 0.025),
            hi = quantile(pct_ase, 0.975))

biotype_counts %>%
  group_by(gene_biotype) %>%
  summarise(n_ase = sum(ase),
            n = n(),
            pct_ase = mean(ase)*100) %>%
  filter(n_ase > 0) %>%
  mutate(gene_biotype = fct_reorder(gene_biotype, pct_ase)) %>%
  ggplot() +
  geom_point(aes(pct_ase, gene_biotype, colour = "Observed"), size = 2) +
  geom_linerange(data = permutation_summary,
                 aes(xmin = lo, xmax = hi, y = gene_biotype, colour = "95% \nexpected\nrange"),
                 size = 3, alpha = 0.2) +
  scale_colour_manual(values = c("steelblue", "black")) +
  labs(x = "% ASE genes", y = "Gene biotype", colour = "",
       caption = "Shaded area is the 95% expected interval from 500 permutations.")



# count ASE per cell/timepoint -------

# split between parents
isolde %>%
  filter(status == "ASE" & abs(diff_prop) >= 0.7) %>%
  count(status, origin, cell_type, stage, n_reps) %>%
  mutate(miss_reps = ifelse(n_reps == 8, NA, n)) %>%
  ggplot(aes(stage, n)) +
  geom_line(aes(group = interaction(origin, cell_type),
                colour = origin)) +
  geom_point(aes(y = miss_reps), shape = 8) +
  facet_grid(~ cell_type) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "", y = "Number of genes",
       subtitle = "ASE difference > 0.7",
       caption = "* = samples with fewer than 8 replicates") +
  scale_colour_manual(values = c("M" = "purple", "P" = "green4"))

# total per cell type
isolde %>%
  filter(status == "ASE" & abs(diff_prop) >= 0.7) %>%
  count(status, cell_type, stage, n_reps) %>%
  mutate(miss_reps = ifelse(n_reps == 8, NA, n)) %>%
  ggplot(aes(stage, n)) +
  geom_line(aes(group = interaction(cell_type),
                colour = cell_type)) +
  geom_point(aes(y = miss_reps), shape = 8) +
  # facet_grid(origin ~ .) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "", y = "Number of genes",
       subtitle = "ASE difference > 0.7",
       colour = "Cell",
       caption = "* = samples with fewer than 8 replicates") +
  scale_colour_manual(values = cell_colours)

# as proportion of genes tested by ISoLDE
isolde %>%
  group_by(cell_type, stage) %>%
  summarise(ase_hi = sum(abs(diff_prop) > 0.7 & status == "ASE", na.rm = TRUE),
            ase_all = sum(status == "ASE"),
            tested = sum(!is.na(diff_prop)),
            n_reps = unique(n_reps)) %>%
  ungroup() %>%
  mutate(ase_hi_prop = ase_hi/tested,
         miss_reps = ifelse(n_reps != 8, ase_hi_prop, NA)) %>%
  ggplot(aes(stage, ase_hi_prop*100)) +
  geom_line(aes(group = interaction(cell_type),
                colour = cell_type)) +
  geom_point(aes(y = miss_reps*100), shape = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "", y = "% of tested genes",
       subtitle = "ASE difference > 0.7",
       colour = "Cell",
       caption = "* = samples with fewer than 8 replicates") +
  scale_colour_manual(values = cell_colours)



# ideogram ------

# make list of genes with ASE
ase_genes <- isolde %>%
  filter(status == "ASE" & chromosome_name != "Y" & abs(diff_prop) > 0.7) %>%
  group_by(gene) %>% filter(n() >= 2) %>% ungroup() %>%
  distinct(gene) %>% pull(gene)

# ASE genes in 4 or more conditions
isolde %>%
  filter(status == "ASE" & chromosome_name != "Y") %>%
  group_by(gene) %>% filter(n() > 3) %>% ungroup() %>%
  mutate(external_gene_name = ifelse(is.na(external_gene_name), gene, external_gene_name)) %>%
  distinct(gene, external_gene_name, chromosome_name, start_position, end_position,
           origin) %>%
  ggplot(aes(x = chromosome_name, y = start_position/1e6)) +
  geom_segment(data = chrom_sizes,
               aes(xend = chromosome_name, yend = size/1e6),
               colour = "lightgrey", size = 4, lineend = "round") +
  geom_point(data = igrs, shape = "diamond", size = 3) +
  # geom_text(data = igrs, aes(label = external_gene_name), nudge_x = -0.7) +
  geom_point(alpha = 0.3, shape = "_", size = 8) +
  ggrepel::geom_text_repel(aes(label = external_gene_name, colour = origin),
                           max.overlaps = 20) +
  # geom_jitter(width = 0.05, height = 0, alpha = 0.2) +
  scale_y_continuous(breaks = seq(0, 500, by = 50)) +
  scale_colour_manual(values = c("M" = "purple", "P" = "green4")) +
  theme(panel.grid = element_blank()) +
  labs(x = "Chromosome", y = "Mb")

# all
isolde %>%
  filter(status == "ASE" & chromosome_name != "Y" & abs(diff_prop) > 0.7) %>%
  group_by(gene) %>% filter(n() >= 2) %>% ungroup() %>%
  mutate(external_gene_name = ifelse(is.na(external_gene_name), gene, external_gene_name)) %>%
  distinct(gene, external_gene_name, chromosome_name, start_position, end_position,
           origin) %>%
  ggplot(aes(x = chromosome_name, y = start_position/1e6)) +
  geom_segment(data = chrom_sizes,
               aes(xend = chromosome_name, yend = size/1e6),
               colour = "lightgrey", size = 4, lineend = "round") +
  geom_point(data = igrs, shape = "diamond", size = 3) +
  geom_point(alpha = 0.3, shape = "_", size = 8) +
  ggrepel::geom_text_repel(aes(label = external_gene_name, colour = origin),
                           max.overlaps = 20) +
  scale_y_continuous(breaks = seq(0, 500, by = 50)) +
  scale_colour_manual(values = c("M" = "purple", "P" = "green4")) +
  theme(panel.grid = element_blank()) +
  labs(x = "Chromosome", y = "Mb",
       subtitle = "ASE difference > 0.7 in at least 2 conditions")

# location of all these genes per condition
isolde %>%
  filter(status == "ASE" & chromosome_name != "Y" & abs(diff_prop) > 0.7) %>%
  group_by(gene) %>% filter(n() >= 2) %>% ungroup() %>%
  mutate(external_gene_name = ifelse(is.na(external_gene_name), gene, external_gene_name)) %>%
  ggplot(aes(x = start_position/1e6, y = timepoint)) +
  geom_segment(data = chrom_sizes,
               aes(xend = size/1e6, y = "t0", yend = "t0"),
               colour = NA) +
  # geom_point(data = igrs, shape = "diamond", size = 3) +
  # geom_text(data = igrs, aes(label = external_gene_name), nudge_x = -0.7) +
  geom_point(alpha = 0.3, size = 1) +
  # ggrepel::geom_text_repel(aes(label = external_gene_name, colour = origin),
  #                          max.overlaps = 20) +
  # geom_jitter(width = 0.05, height = 0, alpha = 0.2) +
  scale_x_continuous(breaks = seq(0, 500, by = 50)) +
  # scale_colour_manual(values = c("M" = "purple", "P" = "green4")) +
  theme(panel.grid = element_blank()) +
  labs(x = "Mb", y = "Stage") +
  facet_grid(cell_type ~ chromosome_name) +
  theme_bw()

isolde %>%
  filter(chromosome_name != "Y" & gene %in% ase_genes) %>%
  filter(status %in% c("ASE", "BA", "UN")) %>%
  # filter((status == "ASE" & abs(diff_prop) > 0.7) | (status == "BA")) %>%
  # group_by(gene) %>% filter(sum(status == "ASE") > 0) %>% ungroup() %>%
  group_by(gene) %>% filter(n_distinct(interaction(cell_type, stage)) >= 2) %>% ungroup() %>%
  mutate(external_gene_name = ifelse(is.na(external_gene_name), gene, external_gene_name)) %>%
  ggplot(aes(x = start_position/1e6, y = timepoint)) +
  geom_segment(data = chrom_sizes,
               aes(xend = size/1e6, y = "t0", yend = "t0"),
               colour = NA) +
  geom_point(aes(colour = status), size = 1) +
  # ggrepel::geom_text_repel(aes(label = external_gene_name, colour = origin),
  #                          max.overlaps = 20) +
  scale_x_continuous(breaks = seq(0, 500, by = 50)) +
  labs(x = "Mb", y = "Stage") +
  facet_grid(cell_type ~ chromosome_name) +
  theme_bw()



# Incoherent signal -----

# list of genes with both ASE and BA status
incoherent_genes <- isolde %>%
  filter(gene %in% ase_genes) %>%
  filter(status %in% c("ASE", "BA")) %>%
  group_by(gene) %>%
  filter(n_distinct(status) > 1) %>%
  ungroup() %>%
  distinct(gene) %>% pull()

# visualise them
isolde %>%
  filter(gene %in% incoherent_genes) %>%
  # plot
  ggplot(aes(stage, diff_prop)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = cell_type, group = cell_type)) +
  # geom_linerange(aes(ymin = diff_prop - variability, ymax = diff_prop + variability)) +
  geom_point(aes(colour = cell_type)) +
  facet_wrap(~ external_gene_name) +
  scale_colour_manual(values = cell_colours) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "", y = "Difference (M - P)", colour = "Cell")



# ASE trends -------------------------------------------------------

# Same genes shown in ideogram
isolde %>% filter(gene %in% ase_genes) %>%
  # plot
  ggplot(aes(stage, diff_prop)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = cell_type, group = cell_type)) +
  # geom_linerange(aes(ymin = diff_prop - variability, ymax = diff_prop + variability)) +
  geom_point(aes(colour = cell_type)) +
  facet_wrap(~ external_gene_name) +
  scale_colour_manual(values = cell_colours) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "", y = "Difference (M - P)")



# genes with DE in Zfp57 dataset
zfp57_degs <- c("Csdc2", "Chil1", "Peg12", "Peg3", "Plagl1", "Zdbf2",
                "H19", "Igf2", "Msi1", "Gm7072", "Rasgrf1", "Rcvrn",
                "Meg3", "Mirg", "Rian", "B830012L14Rik", "Gm37899")
isolde %>%
  filter(external_gene_name %in% zfp57_degs) %>%
  # plot
  ggplot(aes(stage, diff_prop)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = cell_type, group = cell_type)) +
  # geom_linerange(aes(ymin = diff_prop - variability, ymax = diff_prop + variability)) +
  geom_point(aes(colour = cell_type)) +
  facet_wrap(~ external_gene_name) +
  scale_colour_manual(values = cell_colours) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "", y = "Difference (M - P)")


# DEGs in Zfp57 dataset with GO mammary development
isolde %>%
  filter(external_gene_name %in% c("Tbx3", "Bmp4", "Robo1", "Sostdc1", "Csmd1")) %>%
  # plot
  ggplot(aes(stage, diff_prop)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = cell_type, group = cell_type)) +
  # geom_linerange(aes(ymin = diff_prop - variability, ymax = diff_prop + variability)) +
  geom_point(aes(colour = cell_type)) +
  facet_wrap(~ external_gene_name) +
  scale_colour_manual(values = cell_colours) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "", y = "Difference (M - P)")

# DEGs in Zfp57 dataset with GO vasculature development
isolde %>%
  filter(external_gene_name %in% c("Dcn", "Ccl2", "Cldn5", "Col4a3", "Cnmd", "Hhip", "Bmp4")) %>%
  # plot
  ggplot(aes(stage, diff_prop)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = cell_type, group = cell_type)) +
  # geom_linerange(aes(ymin = diff_prop - variability, ymax = diff_prop + variability)) +
  geom_point(aes(colour = cell_type)) +
  facet_wrap(~ external_gene_name) +
  scale_colour_manual(values = cell_colours) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "", y = "Difference (M - P)")


# Example genes ----

# Dlk1 - detected mostly in stromal cells
isolde %>% filter(gene %in% ase_genes) %>%
  filter(external_gene_name == "Dlk1") %>%
  count(cell_type, is.na(diff_prop))


# Dlk1 overall expression
expr %>%
  filter(gene == "ENSMUSG00000040856") %>%
  ggplot(aes(stage, expr)) +
  ggbeeswarm::geom_quasirandom(aes(colour = stage)) +
  geom_line(stat = "summary", fun = "median", aes(group = 1)) +
  facet_grid(~ cell_type) +
  scale_colour_manual(values = stage_colours) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Normalised expression",
       caption = "Line shows the median")

dds_isolde["ENSMUSG00000040856", ] %>%
  assay("counts") %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "sample", values_to = "expr") %>%
  full_join(colData(dds_isolde) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  filter(cell_type %in% c("stromal", "basal")) %>%
  ggplot(aes(allele_parent, log2(expr + 1))) +
  ggbeeswarm::geom_quasirandom(aes(colour = cell_type),
                               dodge.width = 0.7) +
  facet_wrap( ~ stage, ncol = 5) +
  scale_colour_manual(values = cell_colours[c("basal", "stromal")]) +
  labs(x = "Parent allele", colour = "Cell") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


counts(dds_isolde, normalized = TRUE)["ENSMUSG00000040856", , drop = FALSE] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "sample", values_to = "expr") %>%
  full_join(colData(dds_isolde) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  filter(cell_type %in% c("stromal", "basal")) %>%
  ggplot(aes(allele_parent, log2(expr + 1))) +
  ggbeeswarm::geom_quasirandom(aes(colour = cell_type, shape = cross), size = 2) +
  facet_wrap(~ stage, ncol = 5) +
  scale_colour_manual(values = cell_colours[c("stromal", "basal")]) +
  # geom_hline(yintercept = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(15, 17)) +
  labs(x = "Allele", colour = "Cell")



# Tpbgl gene
expr %>%
  filter(gene == "ENSMUSG00000096606") %>%
  ggplot(aes(stage, expr)) +
  ggbeeswarm::geom_quasirandom(aes(colour = stage)) +
  geom_line(stat = "summary", fun = "median", aes(group = 1)) +
  facet_grid(~ cell_type) +
  scale_colour_manual(values = stage_colours) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Normalised expression",
       caption = "Line shows the median")

counts(dds_isolde, normalized = TRUE)["ENSMUSG00000096606", , drop = FALSE] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "sample", values_to = "expr") %>%
  full_join(colData(dds_isolde) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  filter(cell_type %in% c("stromal", "basal")) %>%
  ggplot(aes(allele_parent, log2(expr + 1))) +
  ggbeeswarm::geom_quasirandom(aes(colour = cell_type, shape = cross), size = 2) +
  facet_wrap(~ stage, ncol = 5) +
  scale_colour_manual(values = cell_colours[c("stromal", "basal")]) +
  # geom_hline(yintercept = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(15, 17)) +
  labs(x = "Allele", colour = "Cell")


# Mcts2 gene
expr %>%
  filter(gene == "ENSMUSG00000042814") %>%
  ggplot(aes(stage, expr)) +
  ggbeeswarm::geom_quasirandom(aes(colour = stage)) +
  geom_line(stat = "summary", fun = "median", aes(group = 1)) +
  facet_grid(~ cell_type) +
  scale_colour_manual(values = stage_colours) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Normalised expression",
       caption = "Line shows the median")

counts(dds_isolde, normalized = TRUE)["ENSMUSG00000042814", , drop = FALSE] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "sample", values_to = "expr") %>%
  full_join(colData(dds_isolde) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  filter(cell_type %in% c("adipocytes", "luminal progenitors")) %>%
  ggplot(aes(allele_parent, log2(expr + 1))) +
  ggbeeswarm::geom_quasirandom(aes(colour = cell_type, shape = cross), size = 2) +
  facet_wrap(~ stage, ncol = 5) +
  scale_colour_manual(values = c("adipocytes" = "#1B9E77", "luminal progenitors" = "black")) +
  # geom_hline(yintercept = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(15, 17)) +
  labs(x = "Allele", colour = "Cell")


# Wap gene
expr %>%
  filter(gene == "ENSMUSG00000000381") %>%
  ggplot(aes(stage, expr)) +
  ggbeeswarm::geom_quasirandom(aes(colour = stage)) +
  geom_line(stat = "summary", fun = "median", aes(group = 1)) +
  facet_grid(~ cell_type) +
  scale_colour_manual(values = stage_colours) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Normalised expression",
       caption = "Line shows the median")

counts(dds_isolde, normalized = TRUE)["ENSMUSG00000000381", , drop = FALSE] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "sample", values_to = "expr") %>%
  full_join(colData(dds_isolde) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  filter(cell_type %in% c("luminal differentiated", "endothelial")) %>%
  ggplot(aes(allele_parent, log2(expr + 1))) +
  ggbeeswarm::geom_quasirandom(aes(colour = cell_type, shape = cross), size = 2) +
  facet_wrap(~ stage, ncol = 5) +
  scale_colour_manual(values = cell_colours[c("luminal differentiated", "endothelial")]) +
  # geom_hline(yintercept = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(15, 17)) +
  labs(x = "Allele", colour = "Cell")

isolde %>%
  filter(gene == "ENSMUSG00000000381") %>%
  ggplot(aes(stage, diff_prop)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(colour = cell_type)) +
  geom_line(aes(colour = cell_type, group = cell_type)) +
  scale_y_continuous(limits = c(-1, 1))


# Ramp3 gene
expr %>%
  filter(gene == "ENSMUSG00000041046") %>%
  ggplot(aes(stage, expr)) +
  ggbeeswarm::geom_quasirandom(aes(colour = stage)) +
  geom_line(stat = "summary", fun = "median", aes(group = 1)) +
  facet_grid(~ cell_type) +
  scale_colour_manual(values = stage_colours) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Normalised expression",
       caption = "Line shows the median")

counts(dds_isolde, normalized = TRUE)["ENSMUSG00000041046", , drop = FALSE] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "sample", values_to = "expr") %>%
  full_join(colData(dds_isolde) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  filter(cell_type %in% c("endothelial")) %>%
  ggplot(aes(allele_parent, log2(expr + 1))) +
  ggbeeswarm::geom_quasirandom(aes(colour = cell_type, shape = cross), size = 2) +
  facet_wrap(~ stage, ncol = 5) +
  scale_colour_manual(values = cell_colours[c("luminal differentiated", "endothelial")]) +
  # geom_hline(yintercept = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(15, 17)) +
  labs(x = "Allele", colour = "Cell")


# Show all these three genes
expr %>%
  left_join(annot %>% select(gene, external_gene_name), by = "gene") %>%
  filter(external_gene_name %in% c("Tbrg4", "Wap", "Ramp3")) %>%
  ggplot(aes(stage, expr)) +
  ggbeeswarm::geom_quasirandom(aes(colour = stage)) +
  geom_line(stat = "summary", fun = "median", aes(group = 1)) +
  facet_grid(external_gene_name ~ cell_type) +
  scale_colour_manual(values = stage_colours) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(colour = "black")) +
  labs(y = "Normalised expression",
       caption = "Line shows the median")
isolde %>%
  filter(external_gene_name %in% c("Tbrg4", "Wap", "Ramp3")) %>%
  # filter(gene == "ENSMUSG00000041046") %>%
  ggplot(aes(stage, diff_prop)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(colour = cell_type), size = 1.5) +
  geom_line(aes(colour = cell_type, group = cell_type)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_colour_manual(values = cell_colours) +
  facet_wrap(~ external_gene_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Strain-biases ---------------

# distribution of allele difference - all pooled together
wrap_plots(
  strain %>%
    ggplot(aes(diff_prop)) +
    geom_histogram(binwidth = 0.1) +
    labs(x = "", subtitle = "All genes"),
  strain %>% filter(status == "ASE") %>%
    ggplot(aes(diff_prop)) +
    geom_histogram(binwidth = 0.1) +
    labs(x = "Difference (M - P)", subtitle = "ASE genes"),
  ncol = 1
) &
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

strain %>%
  filter(status == "ASE") %>%
  distinct(gene, cell_type, stage) %>%
  inner_join(isolde %>% filter(gene %in% ase_genes),
             by = c("gene", "cell_type", "stage")) %>%
  distinct(gene, external_gene_name)



# DEGs compared to ASE genes -----

diffexp <- read_csv("results/diffexp/")

