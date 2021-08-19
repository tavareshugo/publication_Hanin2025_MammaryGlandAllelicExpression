#----------------------------------#
# Fit DESeq model for standard     #
# differential expression analysis #
#----------------------------------#

library(DESeq2)

# Read Data ---------------------------------------------------------------

# DDS object with gene-level quantification
dds <- readRDS("results/DESeqDataSet/dds_gene_regular.rds")


# Run DE tests ------------------------------------------------------------

# LRT test for interaction
dds$cell <- factor(gsub(" ", "_", dds$cell_type))
design(dds) <- ~ 0 + cell:timepoint:cross

dds <- DESeq(dds, test = "LRT",
             reduced = ~ cross)

saveRDS(dds, "results/diffexp/dds.rds")

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

write.csv(lfc, "results/diffexp/lfc.csv", row.names = FALSE)
