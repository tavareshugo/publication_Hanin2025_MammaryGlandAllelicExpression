#----------------------------------#
# Fit DESeq model for standard     #
# differential expression analysis #
#----------------------------------#

library(DESeq2)

# Read Data ---------------------------------------------------------------

# DDS object with gene-level quantification
dds <- readRDS("results/DESeqDataSet/dds_gene_regular.rds")

# filtering out lowly-expressed genes and one sample with very few reads
dds <- dds[rowSums(counts(dds) >= 1) > 80, colSums(counts(dds)) > 1e6]


# Fit model ------------------------------------------------------------

# Fit full model
dds$cell <- factor(gsub(" ", "_", dds$cell_type))
design(dds) <- ~ 0 + cell:timepoint:cross

dds <- DESeq(dds)

saveRDS(dds, "results/diffexp/dds.rds")


# Run DE test between timepoints ---------------------------------------

# obtain table of results with LFC between timepoints per cell
contrast_combinations <- combn(unique(as.character(dds$timepoint)), 2)
lfc <- vector("list", nlevels(dds$cell_type)*ncol(contrast_combinations))
mod_mat <- model.matrix(design(dds), colData(dds))
i <- 1
for (cell in levels(dds$cell_type)){
  for (contrast in 1:ncol(contrast_combinations)){
    # timepoints being compared
    t1 <- contrast_combinations[, contrast][1]
    t2 <- contrast_combinations[, contrast][2]

    # coefficient weights for each timepoint
    t1_coef <- colMeans(mod_mat[dds$cell_type == cell & dds$timepoint == t1, ])
    t2_coef <- colMeans(mod_mat[dds$cell_type == cell & dds$timepoint == t2, ])

    # results table (with log2FC shrinkage)
    res <- lfcShrink(dds, type = "ashr", contrast = t2_coef - t1_coef)
    res["cell"] <- cell
    res["contrast"] <- paste0(t2, "_", t1)

    # add to results list
    lfc[[i]] <- res

    # increment counter
    i <- i + 1
  }
}
# clean environment
rm(contrast_combinations, mod_mat, i, cell, contrast, t2, t1, t2_coef, t1_coef, res)

lfc <- do.call("rbind", lfc)
lfc$gene <- rownames(lfc)

write.csv(lfc, "results/diffexp/lfc_timepoints.csv", row.names = FALSE)


# Run DE test between cells -------------------------------------------

# obtain table of results with LFC between cells within each timepoint
contrast_combinations <- combn(unique(as.character(dds$cell_type)), 2)
lfc <- vector("list", nlevels(dds$timepoint)*ncol(contrast_combinations))
mod_mat <- model.matrix(design(dds), colData(dds))
i <- 1
for (t in levels(dds$timepoint)){
  for (contrast in 1:ncol(contrast_combinations)){
    # cells being compared
    cell1 <- contrast_combinations[, contrast][1]
    cell2 <- contrast_combinations[, contrast][2]

    # coefficient weights for each cell type
    cell1_coef <- colMeans(mod_mat[dds$timepoint == t & dds$cell_type == cell1, ])
    cell2_coef <- colMeans(mod_mat[dds$timepoint == t & dds$cell_type == cell2, ])

    # results table (with log2FC shrinkage)
    res <- lfcShrink(dds, type = "ashr", contrast = cell2_coef - cell1_coef)
    res["timepoint"] <- t
    res["contrast"] <- paste0(cell2, "_", cell1)

    # add to results list
    lfc[[i]] <- res

    # increment counter
    i <- i + 1
  }
}
# clean environment
rm(contrast_combinations, mod_mat, i, t, contrast, cell2_coef, cell1_coef, res, cell2, cell1)

lfc <- do.call("rbind", lfc)
lfc$gene <- rownames(lfc)

write.csv(lfc, "results/diffexp/lfc_cells.csv", row.names = FALSE)

