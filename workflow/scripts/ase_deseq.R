library(DESeq2)


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


# design used is based on:
# https://rstudio-pubs-static.s3.amazonaws.com/275642_e9d578fe1f7a404aad0553f52236c0a4.html
# design(dds) <- ~ condition + condition:id_nested_within_cross + condition:allele_strain

# https://support.bioconductor.org/p/81030/
# when doing all cell types maybe the design should be
# ~ 0 + cell_type:animal_id + condition:allele_parent
design <- model.matrix(~ 0 + animal_id + condition:allele_parent, colData(dds))
design <- model.matrix(~ 0 + animal_id + cell_type + condition:allele_parent, colData(dds))
design <- model.matrix(~ 0 + timepoint + cell_type + cell_type:timepoint:allele_parent, colData(dds))
design <- design[, !grepl("allele_parentpaternal", colnames(design))] # to make it full rank
design(dds) <- design

# filter lowly abundant genes
dds <- dds[rowSums(counts(dds) > 0) > 16, ]

# fit model
sizeFactors(dds) <- 1
dds <- DESeq(dds)


# save object for later exploration
saveRDS(dds, "results/ASE_DEseq/dds_ase_test.rds")

