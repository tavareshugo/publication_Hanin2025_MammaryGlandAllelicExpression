library(DESeq2)
library(tidyverse)
library(patchwork)
library(ggfortify)
library(ggbeeswarm)
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

cell_colours <- RColorBrewer::brewer.pal(6, "Dark2")
names(cell_colours) <- c("adipocytes", "basal", "endothelial", "luminal differentiated", "luminal progenitors", "stromal")

cross_colours <- c("bc" = "black", "cb" = "brown")


# Read Data ---------------------------------------------------------------

# DDS object with gene-level quantification
dds <- readRDS("results/DESeqDataSet/dds_gene_regular.rds")

# add QC metrics
dds$sum <- colSums(counts(dds))
dds$detected <- colSums(counts(dds) > 0)
rowData(dds)$sum <- rowSums(counts(dds))
rowData(dds)$detected <- rowSums(counts(dds) > 0)

# add mapping QC to colData - file produced from the qc_mapping.R script
map_qc <- read_csv("results/qc/mapping_qc_summary.csv") %>%
  column_to_rownames("sample")
colData(dds) <- cbind(colData(dds), map_qc[colnames(dds), ])


# check number of replicates in each group
colData(dds) %>%
  as_tibble(rownames = "sample") %>%
  distinct(cell_type, stage, timepoint, cross, animal_id) %>%
  mutate(stage = reorder(stage, rank(timepoint))) %>%
  count(cell_type, stage, cross) %>%
  group_by(cell_type, stage) %>%
  mutate(complete = sum(n) != 8) %>%
  ungroup() %>%
  ggplot(aes(cell_type, stage)) +
  geom_label(aes(label = n, fill = complete)) +
  facet_grid(~ cross) +
  scale_fill_manual(values = c("grey", "lightpink")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))



# Filtering criteria -------------------------------------------------------

# total counts
qplot(rowData(dds)$sum + 1) + scale_x_log10() + labs(x = "Total Counts Per Gene (+1)")
qplot(dds$sum/1e6) + scale_x_log10() + labs(x = "Total Counts Per Sample (Millions)")

# genes detected
qplot(rowData(dds)$detected) + labs(title = "Number of samples a gene is detected in")
qplot(dds$detected) + labs(title = "Number of genes detected per sample")

# correlation
colData(dds) %>% as_tibble() %>%
  ggplot(aes(sum/1e6, detected, colour = cell_type)) +
  geom_point() +
  stat_ellipse() +
  scale_x_log10() +
  labs(x = "Total Counts Per Sample (Millions)",
       y = "Number of genes detected") +
  ggthemes::scale_colour_tableau()

# samples with lower number of detected genes
# mostly adipocytes (maybe they have lower complexity transcriptome?)
# and the endothelial, which is the sample with less than 1M counts
which(colSums(counts(dds) > 0) < 12500)
which(colSums(counts(dds)) < 1e6)

# number of samples/genes retained with different filters
filtered_numbers <- crossing(min_counts = 1:5, min_samples = seq(0, 300, by = 10),
                             min_genes = c(0, 5000, 10000, 15000)) %>%
  rowwise() %>%
  mutate(ngenes = sum(rowSums(counts(dds) >= min_counts) >= min_samples),
         nsamples = sum(colSums(counts(dds) >= min_counts) >= min_genes))

filtered_numbers %>%
  distinct(min_samples, min_counts, ngenes) %>%
  ggplot(aes(min_samples, ngenes)) +
  geom_line(aes(colour = factor(min_counts))) +
  geom_vline(xintercept = 80) +
  labs(x = "Threshold for minimum samples detected",
       y = "Number of genes retained",
       colour = "Count\nthreshold")

filtered_numbers %>%
  distinct(min_genes, min_counts, nsamples) %>%
  ggplot(aes(min_genes, nsamples)) +
  geom_line(aes(colour = factor(min_counts))) +
  labs(x = "Threshold for minimum genes detected",
       y = "Number of samples retained",
       colour = "Count\nthreshold")

# keep genes detected in at least 80 samples
# remove sample with less than 1M reads
dds <- dds[rowSums(counts(dds) >= 1) > 80, colSums(counts(dds)) > 1e6]
# dds <- dds[rowMedians(counts(dds)) > 1, colSums(counts(dds)) > 1e6]


# Normalisation -----------------------------------------------------------

# VST seems to be struggling to stabilize the variance
qplot(rowMeans(assay(dds, "vst")),
      rowSds(assay(dds, "vst"))) +
  geom_smooth(se = FALSE)

# normalise accounting for design
dds$cell_type <- factor(dds$cell_type)
dds$stage <- factor(dds$stage)
dds$timepoint <- factor(dds$timepoint)
dds$cross <- factor(dds$cross)
design(dds) <- ~ 0 + cell_type:timepoint:cross
assay(dds, "vst_design") <- varianceStabilizingTransformation(counts(dds),
                                                              blind = FALSE)

qplot(rowMeans(assay(dds, "vst_design")),
      rowSds(assay(dds, "vst_design"))) +
  geom_smooth(se = FALSE)

# comparing these two
qplot(rowMeans(assay(dds, "vst")),
      rowMeans(assay(dds, "vst_design"))) +
  geom_abline()

# distribution of normalised counts
assay(dds, "vst") %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  full_join(dds %>% colData() %>% as_tibble(rownames = "sample")) %>%
  ggplot(aes(expr)) +
  geom_density(aes(group = sample,
                   colour = cell_type)) +
  scale_colour_manual(values = cell_colours) +
  facet_wrap(~ timepoint)

assay(dds, "vst") %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  full_join(dds %>% colData() %>% as_tibble(rownames = "sample")) %>%
  ggplot(aes(expr)) +
  geom_density(aes(group = sample,
                   colour = stage)) +
  geom_density(data = tibble(expr = rowMeans(assay(dds, "vst"))), size = 1) +
  scale_colour_manual(values = stage_colours) +
  facet_wrap(~ cell_type) +
  labs(x = "Normalised Expression (log2)", colour = "Stage",
       caption = "black line is the distribution of gene average expression across all samples")

assay(dds, "vst") %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  full_join(dds %>% colData() %>% as_tibble(rownames = "sample")) %>%
  ggplot(aes(expr, stage, height = ..density..)) +
  ggridges::geom_density_ridges(aes(group = sample, colour = stage),
                                stat = "density", alpha = 0.1) +
  scale_colour_manual(values = stage_colours) +
  facet_wrap(~ cell_type)

# PCA ---------------------------------------------------------------------

topn <- order(rowVars(assay(dds, "vst")), decreasing = TRUE)[1:1000]
pca <- dds[topn, ] %>%
  assay("vst") %>%
  t() %>% #scale() %>%
  prcomp()

pca %>% broom::tidy("eigenvalues") %>%
  filter(PC %in% 1:20) %>%
  ggplot(aes(factor(PC), cumulative*100)) +
  geom_col(aes(y = percent*100)) +
  geom_line(aes(group = 1))

wrap_plots(
  pca$x %>%
    as_tibble(rownames = "sample") %>%
    select(sample, PC1:PC5) %>%
    left_join(colData(dds) %>% as_tibble(rownames = "sample"),
              by = "sample") %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(colour = cell_type)) +
    scale_colour_manual(values = cell_colours),
  pca$x %>%
    as_tibble(rownames = "sample") %>%
    select(sample, PC1:PC5) %>%
    left_join(colData(dds) %>% as_tibble(rownames = "sample"),
              by = "sample") %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(colour = stage)) +
    scale_colour_manual(values = stage_colours),
  pca$x %>%
    as_tibble(rownames = "sample") %>%
    select(sample, PC1:PC5) %>%
    left_join(colData(dds) %>% as_tibble(rownames = "sample"),
              by = "sample") %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(colour = cross)) +
    scale_colour_manual(values = cross_colours)
) &
  coord_equal() &
  theme_classic(base_size = 16) &
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# first 3 PCs against each other
wrap_plots(
  autoplot(pca,
           data = colData(dds) %>% as_tibble(rownames = "sample"),
           x = 1, y = 2, colour = "cell_type"),
  autoplot(pca,
           data = colData(dds) %>% as_tibble(rownames = "sample"),
           x = 1, y = 3, colour = "cell_type"),
  autoplot(pca,
           data = colData(dds) %>% as_tibble(rownames = "sample"),
           x = 2, y = 3, colour = "cell_type")
) +
  plot_layout(guides = "collect") & ggthemes::scale_colour_tableau() & coord_equal()


# highlight samples flagged as having bad gene coverage (RNA degradation?)
pca$x %>%
  as_tibble(rownames = "sample") %>%
  select(sample, PC1:PC5) %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"),
            by = "sample") %>%
  mutate(highlight = ifelse(genecov_fail, PC2, NA)) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(colour = cell_type)) +
  geom_point(aes(y = highlight), pch = 1, size = 3) +
  scale_colour_manual(values = cell_colours) +
  labs(caption = "Highlighted points have bad gene coverage.") +
  coord_equal()

# look at % RRNA
pca$x %>%
  as_tibble(rownames = "sample") %>%
  select(sample, PC1:PC5) %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"),
            by = "sample") %>%
  mutate(highlight = ifelse(genecov_fail, PC2, NA)) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(colour = pct_rrna)) +
  geom_point(aes(y = highlight), pch = 1, size = 3) +
  scale_colour_viridis_c()

# highlight Involution D1 samples
# where a dip in ASE genes was observed in adipocytes, luminal (both) and endothelial
pca$x %>%
  as_tibble(rownames = "sample") %>%
  select(sample, PC1:PC5) %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"),
            by = "sample") %>%
  mutate(highlight = ifelse(stage == "involution d1", PC2, NA)) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(colour = cell_type)) +
  geom_point(aes(y = highlight), pch = 1, size = 3) +
  scale_colour_manual(values = cell_colours) +
  labs(caption = "Highlighted points are involution D1.",
       colour = "Cell") +
  coord_equal()



# UMAP --------------------------------------------------------------------

library(uwot)
topn <- order(rowVars(assay(dds, "vst")), decreasing = TRUE)[1:1000]
set.seed(20210401)
sample_umap <- dds[topn, ] %>%
  assay("vst") %>%
  t() %>%
  umap() %>%
  as_tibble() %>%
  bind_cols(colData(dds) %>% as_tibble())

wrap_plots(
  sample_umap %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(colour = cell_type)) +
    scale_colour_manual(values = cell_colours),
  sample_umap %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(colour = stage)) +
    scale_colour_manual(values = stage_colours),
  sample_umap %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(colour = cross)) +
    scale_colour_manual(values = cross_colours)
) &
  coord_equal() &
  theme_classic(base_size = 16) &
  theme(axis.text = element_blank(), axis.ticks = element_blank()) &
  labs(x = "UMAP1", y = "UMAP2")

sample_umap %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = cell_type)) +
  # scale_colour_viridis_d()
  ggthemes::scale_colour_tableau() +
  coord_equal() + theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2")

sample_umap %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = stage)) +
  scale_colour_manual(values = stage_colours) +
  coord_equal() + theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2")

set.seed(20210401)
sample_umap <- dds[topn, ] %>%
  assay("vst") %>%
  t() %>%
  umap(n_components = 3) %>%
  as_tibble() %>%
  bind_cols(colData(dds) %>% as_tibble())

library(rgl); options(rgl.printRglwidget = TRUE)
plot3d(sample_umap$V1, sample_umap$V2, sample_umap$V3,
       col = stage_colours[sample_umap$stage],
       type = "p", radius = .2, axes = FALSE, xlab = "", ylab = "", zlab = "")
play3d( spin3d( axis = c(0, 0, 1), rpm = 8), duration = 1/8*60)


# hclust ------------------------------------------------------------------

topn <- order(rowVars(assay(dds, "vst")), decreasing = TRUE)[1:1000]
sample_clust <- dds[topn, ] %>%
  assay("vst") %>%
  t() %>% #scale() %>%
  dist() %>% hclust(method = "ward.D2")

plot(sample_clust)

dend_data <- ggdendro::dendro_data(sample_clust)

# Add sample information to dendrogram labels
dend_data$labels <- dend_data$labels %>%
  full_join(colData(dds) %>% as_tibble(rownames = "label") %>%
              distinct(label, timepoint, stage, cell_type, cross))

# Make a plot
# Note: the invisible point layer is for use with coord_polar
ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dend_data$labels,
             aes(x = x, y = y, colour = cell_type),
             size = 2, fill = "white") +
  geom_point(data = tibble(x = 0, y = 0), aes(x, y), alpha = 0, colour = "red") +
  # geom_text(data = dend_data$labels,
  #           aes(x = x, y = y, colour = cell_type, label = paste(cross, animal)),
  #           hjust = 0, nudge_y = 0.003, show.legend = FALSE) +
  scale_colour_manual(values = cell_colours) +
  scale_shape_manual(values = c(21, 19)) +
  scale_y_reverse(expand = c(0, 0.1)) +
  coord_flip() +
  coord_polar() +
  theme_void(base_size = 16) +
  labs(colour = "Cell Type")

ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dend_data$labels,
             aes(x = x, y = y, colour = stage),
             size = 2, fill = "white") +
  geom_point(data = tibble(x = 0, y = 0), aes(x, y), alpha = 0, colour = "red") +
  # geom_text(data = dend_data$labels,
  #           aes(x = x, y = y, colour = cell_type, label = paste(cross, animal)),
  #           hjust = 0, nudge_y = 0.003, show.legend = FALSE) +
  scale_colour_manual(values = stage_colours) +
  scale_shape_manual(values = c(21, 19)) +
  scale_y_reverse(expand = c(0, 0.1)) +
  coord_flip() +
  coord_polar() +
  theme_void(base_size = 16) +
  labs(colour = "Stage")


# Calculate correlations
cor_dists <- cor(assay(dds, "vst"))
diag(cor_dists) <- NA

# Convert to distance
# https://arxiv.org/pdf/1208.3145.pdf
# https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering
cor_dists <- as.dist( sqrt(1/2*(1 - cor_dists)) )

# Apply clustering
cor_cluster <- hclust(cor_dists, method = "ward.D2")

dend_data <- ggdendro::dendro_data(cor_cluster)

# Add sample information to dendrogram labels
dend_data$labels <- dend_data$labels %>%
  full_join(colData(dds) %>% as_tibble(rownames = "label") %>%
              distinct(label, timepoint, stage, cell_type, cross))

# Make a plot
# Note: the invisible point layer is for use with coord_polar
ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dend_data$labels,
             aes(x = x, y = y, colour = cell_type),
             size = 2, fill = "white") +
  geom_point(data = tibble(x = 0, y = 0), aes(x, y), alpha = 0, colour = "red") +
  # geom_text(data = dend_data$labels,
  #           aes(x = x, y = y, colour = cell_type, label = paste(cross, animal)),
  #           hjust = 0, nudge_y = 0.003, show.legend = FALSE) +
  scale_colour_manual(values = cell_colours) +
  scale_shape_manual(values = c(21, 19)) +
  scale_y_reverse(expand = c(0, 0.1)) +
  coord_flip() +
  coord_polar() +
  theme_void(base_size = 16) +
  labs(colour = "Cell Type")

ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dend_data$labels,
             aes(x = x, y = y, colour = stage),
             size = 2, fill = "white") +
  geom_point(data = tibble(x = 0, y = 0), aes(x, y), alpha = 0, colour = "red") +
  # geom_text(data = dend_data$labels,
  #           aes(x = x, y = y, colour = cell_type, label = paste(cross, animal)),
  #           hjust = 0, nudge_y = 0.003, show.legend = FALSE) +
  scale_colour_manual(values = stage_colours) +
  scale_shape_manual(values = c(21, 19)) +
  scale_y_reverse(expand = c(0, 0.1)) +
  coord_flip() +
  coord_polar() +
  theme_void(base_size = 16) +
  labs(colour = "Stage")



# Sample-Sample Correlation -------------------------------------------------

# heatmap suggests some samples with low correlation - what are they?
cor_samples <- cor(assay(dds, "vst"))
# cor_samples[upper.tri(cor_samples, diag = TRUE)] <- NA
hist(cor_samples)

# heatmap
pheatmap::pheatmap(cor_samples, cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE, show_colnames = FALSE,
                   annotation_col = tibble(sample = colnames(dds),
                                           cell_type = dds$cell_type,
                                           stage = dds$stage) %>%
                     column_to_rownames("sample"),
                   annotation_row = tibble(sample = colnames(dds),
                                           cell_type = dds$cell_type,
                                           stage = dds$stage) %>%
                     column_to_rownames("sample"),
                   annotation_colors = list(cell_type = cell_colours,
                                            stage = stage_colours))

pheatmap::pheatmap(cor_samples, clustering_method = "ward.D2",
                   show_rownames = FALSE, show_colnames = FALSE,
                   annotation_col = tibble(sample = colnames(dds),
                                           cell_type = dds$cell_type,
                                           stage = dds$stage) %>%
                     column_to_rownames("sample"),
                   annotation_row = tibble(sample = colnames(dds),
                                           cell_type = dds$cell_type,
                                           stage = dds$stage) %>%
                     column_to_rownames("sample"),
                   annotation_colors = list(cell_type = cell_colours,
                                            stage = stage_colours))

# check samples with highest number of low-correlations
cor_samples %>%
  as_tibble(rownames = "sample1") %>%
  pivot_longer(-sample1, names_to = "sample2", values_to = "cor") %>%
  drop_na() %>%
  filter(cor < 0.75) %>%
  distinct(sample1, sample2) %>%
  gather() %>%
  count(value) %>%
  arrange(desc(n))

# these are some of them
low_cor_samples <- c("adipocytes_t7_bc_M00324968",
  "adipocytes_t7_cb_M00323782",
  "adipocytes_t6_bc_M00314944",
  "adipocytes_t4_bc_M00314943",
  "stromal_t8_cb_M00338410")

cor_samples %>%
  as_tibble(rownames = "sample1") %>%
  pivot_longer(-sample1, names_to = "sample2", values_to = "cor") %>%
  drop_na() %>%
  filter(cor < 1) %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"),
            by = c("sample1" = "sample")) %>%
  mutate(highlight = ifelse(
    sample1 %in% low_cor_samples,
    paste(cell_type, stage, cross, sep = "; "),
    NA)) %>%
  arrange(!is.na(highlight)) %>%
  ggplot(aes(cor)) +
  geom_line(stat = "density",
            aes(group = sample1, colour = highlight,
                alpha = highlight)) +
  scale_colour_brewer(palette = "Set2", na.value = "grey50") +
  scale_alpha_manual(values = rep(1, 5), na.value = 0.5) +
  labs(x = "Pearson's Correlation",
       title = "Correlation between a given sample and all others",
       colour = "", alpha = "")

# highlight them in PCA
pca$x %>%
  as_tibble(rownames = "sample") %>%
  select(sample, PC1:PC5) %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"),
            by = "sample") %>%
  mutate(highlight = ifelse(
    sample %in% low_cor_samples,
    PC2,
    NA)) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(colour = cell_type)) +
  geom_point(aes(y = highlight), shape = 1, size = 3) +
  scale_colour_manual(values = cell_colours)


# Investigating a few genes -----------------------------------------------

# known imprints
imprints <- c("DLK1", "MEST", "IGF2R", "IGF2", "GRB10", "GNAS", "UBE3A", "CDKN1C", "SNRPN", "NNAT", "PEG3", "COPG2")

# export data for Boshra
assay(dds, "vst")[toupper(rowData(dds)$external_gene_name) %in% c(imprints, "ZFP57"), ] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  select(sample, gene, external_gene_name, cell_type, stage, timepoint, cross, expr) %>%
  write_csv("~/temp/rnaseq_normalised_expression_imprinted_genes.csv")

# visualise
assay(dds, "vst")[toupper(rowData(dds)$external_gene_name) %in% imprints, ] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  ggplot(aes(cell_type, expr)) +
  # geom_jitter(aes(colour = timepoint)) +
  ggbeeswarm::geom_quasirandom(aes(colour = stage)) +
  facet_wrap(~ external_gene_name) +
  scale_colour_manual(values = stage_colours) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

assay(dds, "vst")[toupper(rowData(dds)$external_gene_name) %in% imprints, ] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  ggplot(aes(timepoint, expr)) +
  ggbeeswarm::geom_quasirandom(aes(colour = cell_type)) +
  facet_wrap(~ external_gene_name) +
  ggthemes::scale_colour_tableau() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

assay(dds, "vst")[toupper(rowData(dds)$external_gene_name) %in% imprints, ] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  ggplot(aes(timepoint, expr)) +
  geom_point(colour = "steelblue") +
  geom_line(stat = "summary", fun = "median", aes(group = 1)) +
  # ggbeeswarm::geom_quasirandom(aes(colour = cross)) +
  facet_grid(external_gene_name ~ cell_type) +
  ggthemes::scale_colour_tableau() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic(base_size = 16) +
  theme(panel.border = element_rect(colour = "grey48", fill=NA, size=1))


# ZFP57
assay(dds, "vst")[which(rowData(dds)$external_gene_name == "Zfp57"), , drop = FALSE] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  ggplot(aes(timepoint, expr)) +
  geom_point() +
  geom_line(stat = "summary", fun = "median", aes(group = 1)) +
  facet_grid(external_gene_name ~ cell_type) +
  ggthemes::scale_colour_tableau() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

assay(dds, "vst") %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  separate(sample, c("cell_type", "timepoint", "cross", "animal_id")) %>%
  ggplot(aes(expr, group = sample, colour = cell_type)) +
  geom_density(alpha = 0.1) +
  ggthemes::scale_colour_tableau()



# Gene Expression Distributions ------------------------------------------

# distribution of normalised counts
assay(dds, "vst") %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  full_join(dds %>% colData() %>% as_tibble(rownames = "sample")) %>%
  # filter(str_detect(cell_type, "luminal") & str_detect(stage, "lactation")) %>%
  ggplot(aes(expr)) +
  geom_density(aes(group = sample,
                   colour = cell_type)) +
  geom_density(data = tibble(expr = rowMeans(assay(dds, "vst")))) +
  scale_colour_manual(values = cell_colours) +
  facet_grid(cell_type ~ stage) #+ coord_cartesian(ylim = c(0, 0.25))

# average expression
qplot(rowMeans(assay(dds, "vst")), geom = "density") +
  scale_x_continuous(breaks = seq(2, 30, 3))

qplot(rowMeans(assay(dds, "vst")), geom = "density") +
  coord_cartesian(ylim = c(0, 0.1)) +
  scale_x_continuous(breaks = seq(2, 30, 3))



# Reference Genes ---------------------------------------------------------

# Taken from https://dx.doi.org/10.1038%2Fsrep35595
ref_genes_mammary <- c("Arpc3", "Clock", "Ctbp1", "Phf7", "Prdx1", "Sugp2", "Taf11", "Usp7")
# all commonly used reference genes
ref_genes <- c("Prdx1", "Phf7", "Ctbp1", "Tbp", "Rpl13a", "Hprt", "Sugp2", "Gapdh", "Clock", "Usp7", "Taf11", "18S", "Sdha", "Actb", "Hmbs", "Arpc3", "B2m", "Ywhaz", "Pgk1", "Gusb")
all(ref_genes %in% rowData(dds)$external_gene_name)
ref_genes[!ref_genes %in% rowData(dds)$external_gene_name]

assay(dds, "vst")[which(rowData(dds)$external_gene_name %in% ref_genes), ] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  select(sample, gene, external_gene_name, cell_type, stage, timepoint, cross, expr) %>%
  ggplot(aes(stage, expr)) +
  geom_quasirandom(aes(colour = cell_type),
                   dodge.width = 0.7) +
  facet_wrap(~ external_gene_name) +
  scale_colour_manual(values = cell_colours) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Stage", y = "Normalised Expression", colour = "Cell Type") +
  scale_y_continuous(breaks = seq(2, 30, 3)) +
  theme(axis.text.x = element_text(size = 8))

assay(dds, "vst")[which(rowData(dds)$external_gene_name %in% ref_genes_mammary), ] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  select(sample, gene, external_gene_name, cell_type, stage, timepoint, cross, expr) %>%
  ggplot(aes(stage, expr)) +
  geom_quasirandom(aes(colour = cell_type),
                   dodge.width = 0.7) +
  facet_wrap(~ external_gene_name) +
  scale_colour_manual(values = cell_colours) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Stage", y = "Normalised Expression", colour = "Cell Type") +
  scale_y_continuous(breaks = seq(2, 30, 3))


# some cyclins
assay(dds, "vst")[grep("Ccn[a-d]", rowData(dds)$external_gene_name), ] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  select(sample, gene, external_gene_name, cell_type, stage, timepoint, cross, expr) %>%
  ggplot(aes(stage, expr)) +
  geom_quasirandom(aes(colour = cell_type),
                   dodge.width = 0.7) +
  facet_wrap(~ external_gene_name, scales = "free_y") +
  scale_colour_manual(values = cell_colours) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Stage", y = "Normalised Expression", colour = "Cell Type") +
  scale_y_continuous(breaks = seq(2, 30, 3))

# and cyclin-dependent kinases
assay(dds, "vst")[grep("Cdk", rowData(dds)$external_gene_name), ] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(rowData(dds) %>% as_tibble(rownames = "gene"), by = "gene") %>%
  left_join(colData(dds) %>% as_tibble(rownames = "sample"), by = "sample") %>%
  select(sample, gene, external_gene_name, cell_type, stage, timepoint, cross, expr) %>%
  filter(str_detect(cell_type, "luminal")) %>%
  ggplot(aes(stage, expr)) +
  geom_quasirandom(aes(colour = cell_type),
                   dodge.width = 0.7) +
  facet_wrap(~ external_gene_name, scales = "free_y") +
  scale_colour_manual(values = cell_colours) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Stage", y = "Normalised Expression", colour = "Cell Type") +
  scale_y_continuous(breaks = seq(2, 30, 3))
