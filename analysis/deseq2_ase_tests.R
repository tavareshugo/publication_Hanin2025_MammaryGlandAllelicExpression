library(DESeq2)
library(UpSetR)
library(tidyverse)
theme_set(theme_minimal(base_size = 14))


# Read Data ---------------------------------------------------------------

# imprinted genes from Tucci et al
tucci <- readxl::read_excel("data/external/tucci_et_al_sup1_imprinted_genes.xlsx",
                            sheet = "Mouse_All_Known IG_NoStatus")

# other known imprints
imprints <- c("DLK1", "MEST", "IGF2R", "IGF2", "GRB10", "GNAS", "UBE3A", "CDKN1C", "SNRPN", "NNAT", "PEG3", "COPG")

# DDS object with allele-specific counts at the gene level
dds_gene_allele <- readRDS("data/processed/dds_gene_allele.rds")

# check number of replicates in each group
colData(dds_gene_allele) %>%
  as_tibble(rownames = "sample") %>%
  distinct(cell_type, stage, timepoint, cross, animal_id) %>%
  mutate(stage = reorder(stage, rank(timepoint))) %>%
  count(cell_type, stage, cross) %>%
  group_by(cell_type, stage) %>%
  mutate(complete = sum(n) == 8) %>%
  ungroup() %>%
  ggplot(aes(cell_type, stage)) +
  geom_label(aes(label = n, fill = complete)) +
  facet_grid(~ cross) +
  scale_fill_manual(values = c("grey", "lightpink")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# filter to retain a few samples only for testing
keep <- dds_gene_allele$cell_type == "stromal" &
  dds_gene_allele$timepoint %in% c("t0", "t1", "t2", "t4", "t5", "t6")
dds <- dds_gene_allele[, keep]
colData(dds) <- droplevels(colData(dds))

# create a condition variable to analyse these groups separately (easier for analysis)
dds$condition <- paste(dds$cell_type,
                       dds$timepoint,
                       sep = ".")
dds$condition <- gsub("luminal differentiated", "luminald", dds$condition)
dds$condition <- gsub("luminal progenitors", "luminalp", dds$condition)

# factorize all columns
dds$condition <- factor(dds$condition)
dds$animal_id <- factor(gsub(" .*|\\(AF.*", "", dds$animal_id))
dds$allele_strain <- factor(dds$allele_strain)
dds$allele_parent <- factor(dds$allele_parent)

# design used is based on:
# https://rstudio-pubs-static.s3.amazonaws.com/275642_e9d578fe1f7a404aad0553f52236c0a4.html
# design(dds) <- ~ condition + condition:id_nested_within_cross + condition:allele_strain

# https://support.bioconductor.org/p/81030/
# when doing all cell types maybe the design should be
# ~ 0 + cell_type:animal_id + condition:allele_parent
design <- model.matrix(~ 0 + animal_id + condition:allele_parent, colData(dds))
design <- design[, !grepl("allele_parentpaternal", colnames(design))] # to make it full rank
design(dds) <- design

# filter lowly abundant genes
dds <- dds[rowSums(counts(dds) > 0) > 16, ]

# fit model
sizeFactors(dds) <- 1
dds <- DESeq(dds)

# save object for later exploration
saveRDS(dds, "data/processed/dds_ase_test.rds")

#
# # WRONG!!!! ORDER OF NAMES NOT CORRECT
# names(imprints) <- gene_annot %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   filter(toupper(external_gene_name) %in% imprints) %>%
#   pull(ensembl_gene_id)
#
# imprints_correct <- gene_annot[toupper(gene_annot$external_gene_name) %in% imprints, ]



# Results -----------------------------------------------------------------

# get results for ASE in all timepoints
contrasts <- resultsNames(dds)
contrasts <- contrasts[!grepl("animal_id", contrasts)]

res <- lapply(contrasts, function(i){
  dds %>%
    results(contrast = list(i)) %>%
    as_tibble(rownames = "gene") %>%
    mutate(contrast = i)
})
res <- bind_rows(res)

# calculate allele bias (maternal allele)
res <- res %>%
  mutate(ase = 2^log2FoldChange/(1 + 2^log2FoldChange))

# nicer naming for the contrast
res <- res %>%
  mutate(contrast = str_replace(contrast, ".*\\.t", "t") %>% str_remove("\\..*"))



# Exploratory Analysis ----------------------------------------------------

# distribution of ASE
res %>%
  ggplot(aes(ase)) +
  geom_density(aes(colour = contrast))

# correlation across timepoints
res %>%
  group_by(gene) %>%
  filter(any(padj < 0.05)) %>%
  ungroup() %>%
  select(gene, contrast, ase) %>%
  pivot_wider(names_from = contrast, values_from = ase) %>%
  ggplot(aes(t0, t2)) +
  geom_point(alpha = 0.1)

# MA plots
res %>%
  ggplot(aes(log10(baseMean), ase)) +
  geom_point(size = 0.5) +
  geom_point(data = . %>% filter(padj < 0.05), colour = "orange", size = 0.5) +
  geom_point(data = . %>% filter(padj < 0.01), colour = "brown", size = 0.5) +
  facet_wrap(~ contrast) +
  labs(x = "log2(mean expression)", y = "Maternal Allele Bias")

# checking some known imprints
res %>%
  full_join(rowData(dds) %>% as_tibble(rownames = "gene") %>% select(gene, external_gene_name),
            by = "gene") %>%
  ggplot(aes(log10(baseMean), ase)) +
  geom_point(size = 0.5, colour = "lightgrey") +
  geom_point(data = . %>% filter(toupper(external_gene_name) %in% imprints),
             colour = "brown", size = 0.5) +
  ggrepel::geom_text_repel(data = . %>% filter(toupper(external_gene_name) %in% imprints),
            aes(label = external_gene_name)) +
  facet_wrap(~ contrast) +
  labs(x = "log2(mean expression)", y = "Maternal Allele Bias")

# checking Tucci et al full list
res %>%
  full_join(rowData(dds) %>% as_tibble(rownames = "gene") %>% select(gene, external_gene_name),
            by = "gene") %>%
  ggplot(aes(log10(baseMean), ase)) +
  geom_point(size = 0.5, colour = "lightgrey") +
  geom_point(data = . %>% filter(toupper(external_gene_name) %in% toupper(tucci$Gene)),
             colour = "brown", size = 0.5) +
  facet_wrap(~ contrast) +
  labs(x = "log2(mean expression)", y = "Maternal Allele Bias")


# intersection plot
res %>%
  filter(padj < 0.01 & (ase > 0.9 | ase < 0.1)) %>%
  with(split(gene, contrast)) %>%
  fromList() %>%
  upset(order.by = "freq")



# Effect of time ----------------------------------------------------------

dds_lrt <- DESeq(dds, test = "LRT",
                 reduced = model.matrix(~ 0 + animal_id, colData(dds)))
res_lrt <- dds_lrt %>%
  results() %>%
  as_tibble(rownames = "gene")

dif_genes <- res_lrt %>%
  filter(padj < 0.01) %>%
  pull(gene)

# highlight in MA plots
res %>%
  filter(baseMean > 0) %>%
  ggplot(aes(log10(baseMean), ase)) +
  geom_point(size = 0.5) +
  geom_point(data = . %>% filter(gene %in% dif_genes), colour = "brown", size = 0.7) +
  facet_wrap(~ contrast) +
  labs(x = "log2(mean expression)", y = "Maternal Allele Bias")

res %>%
  filter(gene %in% dif_genes) %>%
  ggplot(aes(contrast, ase)) +
  geom_line(aes(group = gene))


# Read snakemake-processed file ----

dds <- readRDS("results/ase_deseq2/dds_ase.rds")

# get results for ASE in all conditions
contrasts <- resultsNames(dds)
contrasts <- contrasts[grepl("maternal", contrasts)]

res <- mclapply(contrasts, function(i){
  cat(i)
  dds %>%
    results(contrast = list(i)) %>%
    as_tibble(rownames = "gene") %>%
    mutate(contrast = i)
}, mc.cores = 3)
res <- bind_rows(res)

# calculate allele bias (maternal allele)
res <- res %>%
  mutate(ase = 2^log2FoldChange/(1 + 2^log2FoldChange))

# nicer naming for the contrast
res <- res %>%
  mutate(contrast = str_replace(contrast, ".*\\.t", "t") %>% str_remove("\\..*"))
