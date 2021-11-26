library(DESeq2)
library(tidyverse)
library(Matrix)
library(igraph)
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


# Read data ---------------------------------------------------------------

dds <- readRDS("results/diffexp/dds.rds")

lfc <- read_csv("results/diffexp/lfc.csv")

# add gene names
lfc <- rowData(dds) %>% as_tibble(rownames = "gene") %>%
  distinct(gene, external_gene_name) %>%
  rename(gene_name = external_gene_name) %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene, gene_name)) %>%
  full_join(lfc, by = "gene")

rm(dds) # only needed this to fetch rowData


# MA plots ----------------------------------------------------------------

lfc %>%
  filter(cell == "adipocytes") %>%
  ggplot(aes(log10(baseMean), log2FoldChange)) +
  geom_point(size = 0.8) +
  facet_wrap( ~ contrast)


# Cluster genes -----------------------------------------------------------

# create correlation matrix
adjmat <- lfc %>%
  filter(padj < 0.01) %>%
  filter(cell == "endothelial") %>%
  group_by(gene) %>% filter(any(abs(log2FoldChange) > 1)) %>%
  mutate(group = paste(cell, contrast, sep = ".")) %>%
  select(gene, group, log2FoldChange) %>%
  pivot_wider(names_from = "group", values_from = "log2FoldChange") %>%
  column_to_rownames("gene") %>%
  as.matrix() %>%
  t() %>%
  cor()

# convert to 0/1 based on threshold and store as a sparseMatrix
adjmat[adjmat < 0.9] <- 0
adjmat[adjmat >= 0.9] <- 1
adjmat <- Matrix(adjmat, sparse = TRUE)

# build graph
g <- graph_from_adjacency_matrix(adjmat, mode = "undirected")

# degree distribution
plot(1:length(degree.distribution(g)), degree.distribution(g), log = "xy", pch = 19,
     xlab = "Degree", ylab = "p", las = 1)

# clustering
clust <- cluster_louvain(g)

tibble(gene = clust$names, cluster = clust$membership) %>%
  group_by(cluster) %>%
  filter(n() > 50) %>%
  ungroup() %>%
  inner_join(lfc, by = "gene") %>%
  filter(cell == "endothelial") %>%
  ggplot(aes(contrast, log2FoldChange)) +
  geom_line(aes(colour = cell, group = interaction(cell, gene)), alpha = 0.5) +
  facet_wrap(~ cluster) +
  scale_colour_manual(values = cell_colours)



lfc %>%
  filter(cell == "adipocytes" & contrast == "t1_t0") %>% nrow()
