library(DESeq2)
library(tidyverse)


# Read Data ---------------------------------------------------------------

# DDS object with allele-specific counts at the gene level
dds <- readRDS("results/DESeqDataSet/dds_gene_allele.rds")


# Prepare Data ------------------------------------------------------------

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


# Fit model ------------------------------------------------------------

# design used is based on:
# https://rstudio-pubs-static.s3.amazonaws.com/275642_e9d578fe1f7a404aad0553f52236c0a4.html
# design(dds) <- ~ condition + condition:id_nested_within_cross + condition:allele_strain

# https://support.bioconductor.org/p/81030/
# when doing all cell types maybe the design should be
# ~ 0 + cell_type:animal_id + condition:allele_parent
design <- model.matrix(~ 0 + animal_id + condition:allele_parent, colData(dds))
design <- model.matrix(~ 0 + animal_id + cell_type + condition:allele_parent, colData(dds))
design <- model.matrix(~ 0 + timepoint + cell_type + cell_type:timepoint:allele_parent, colData(dds))
design <- model.matrix(~ 0 + cell_type:animal_id + cell_type:timepoint:allele_parent, colData(dds))
design <- design[, !grepl("allele_parentpaternal", colnames(design))] # to make it full rank
design <- design[, which(colSums(design) > 0)] # remove terms for missing replicates
design(dds) <- design

# filter lowly abundant genes
dds <- dds[rowSums(counts(dds) > 0) > 16, ]

# fit model
sizeFactors(dds) <- 1
dds <- DESeq(dds)


# Results ------------------------------------------------------------

# get results for ASE in all conditions
contrasts <- resultsNames(dds)
contrasts <- contrasts[grepl("maternal", contrasts)]

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


# save objects
saveRDS(dds, "results/ase_deseq2/dds_ase.rds")
write_csv(res, "results/ase_deseq2/deseq2_ase.csv")
