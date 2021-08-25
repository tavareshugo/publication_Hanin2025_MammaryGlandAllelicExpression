#----------------------------------#
# Fit DESeq model for standard     #
# differential expression analysis #
#----------------------------------#

library(DESeq2)

# Read Data ---------------------------------------------------------------

# DDS object with gene-level quantification
dds <- readRDS("results/DESeqDataSet/dds_gene_regular.rds")


# Fit model ------------------------------------------------------------

# Fit full model
dds$cell <- factor(gsub(" ", "_", dds$cell_type))
design(dds) <- ~ 0 + cell:timepoint:cross

dds <- DESeq(dds)

saveRDS(dds, "results/diffexp/dds.rds")


# Run DE test between timepoints ---------------------------------------

# obtain table of results with LFC across consecutive timepoints per cell
lfc <- vector("list", nlevels(dds$cell_type)*9)
mod_mat <- model.matrix(design(dds), colData(dds))
i <- 1
for (cell in levels(dds$cell_type)){
  for (t in 1:9){
    # coefficient weights for each timepoint
    t2_coef <- colMeans(mod_mat[dds$cell_type == cell & dds$timepoint == paste0("t", t), ])
    t1_coef <- colMeans(mod_mat[dds$cell_type == cell & dds$timepoint == paste0("t", t-1), ])

    # results table (with log2FC shrinkage)
    res <- lfcShrink(dds, type = "ashr", contrast = t2_coef - t1_coef)
    res["cell"] <- cell
    res["contrast"] <- paste0("t", t, "_", "t", t-1)

    # add to results list
    lfc[[i]] <- res

    # increment counter
    i <- i + 1
  }
}
rm(mod_mat, i, cell, t, t2_coef, t1_coef, res) # clean environment

lfc <- do.call("rbind", lfc)
lfc$gene <- rownames(lfc)

write.csv(lfc, "results/diffexp/lfc_timepoints.csv", row.names = FALSE)


# Run DE test between cells -------------------------------------------

# obtain table of results with LFC between cells within each timepoint
contrast_combinations <- combn(unique(as.character(dds$cell_type)), 2)
lfc <- vector("list", nlevels(dds$cell_type)*9)
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
rm(mod_mat, i, t, cell2_coef, cell1_coef, res, cell2, cell1) # clean environment

lfc <- do.call("rbind", lfc)
lfc$gene <- rownames(lfc)

write.csv(lfc, "results/diffexp/lfc_cells.csv", row.names = FALSE)

