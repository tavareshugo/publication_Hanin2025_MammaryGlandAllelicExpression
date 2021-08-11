library(DESeq2)
library(ISoLDE)

# Read Data ---------------------------------------------------------------

# DDS object with allele-specific counts at the gene level
dds <- readRDS("results/DESeqDataSet/dds_gene_allele.rds")


# Prepare Data for ISolDE -------------------------------------------------

# ISoLDE has functions to read counts from files
# I have checked the format, and basically this corresponds to pulling the counts matrix from the DDS object
# and creating what they call "sample target", which is basically some form of colData()

# reorder object to have alternating parents for each replicate
dds <- dds[, order(dds$cell_type, dds$timepoint, dds$animal_id, dds$allele_parent, dds$allele_strain)]

# create the "sample target" data frame
sample_target <- data.frame(sample = colnames(dds),
                            parent = dds$allele_parent,
                            strain = dds$allele_strain,
                            stringsAsFactors = TRUE)

# rename cell types for convenience in output without spaces
dds$cell_type <- gsub("luminal differentiated", "luminald", dds$cell_type)
dds$cell_type <- gsub("luminal progenitors", "luminalp", dds$cell_type)


# Run Tests ---------------------------------------------------------------

# Tests are run on each cell_type:timepoint combination

for (i_cell in unique(dds$cell_type)){
  for (i_timepoint in unique(dds$timepoint)) {
    cat("Processing cell: ", i_cell, "\n")
    cat("Processing timepoint: ", i_timepoint, "\n")

    # samples to keep
    samples_keep <- which(dds$cell_type == i_cell & dds$timepoint == i_timepoint)

    # estimate size factors for normalisation (ISoLDE recommends doing this,
    # although I wonder if size factors should be estimated at the sample level (rather than allele-specific level))
    dds_subset <- estimateSizeFactors(dds[, samples_keep])

    # filter
    filtered_cts <- filterT(rawASRcounts = counts(dds_subset),
                            normASRcounts = counts(dds_subset, normalized = TRUE),
                            target = sample_target[samples_keep, ],
                            tol_filter = 50,
                            bias = "parental")

    # test
    test_results <- isolde_test(bias = "parental",
                                method = "default",
                                asr_counts = filtered_cts$filteredASRcounts,
                                target = sample_target[samples_keep, ],
                                nboot = 5000, pcore = 100,
                                graph = FALSE,
                                prefix = paste0(i_cell, "_", i_timepoint),
                                outdir = "results/ase_isolde/")

    # rename output files to remove date
    outfile <- list.files("results/ase_isolde/", pattern = paste0(i_cell, "_", i_timepoint, "_ALL"),
                          full.names = TRUE)
    if (length(outfile) != 1) stop("problem with output file")
    file.rename(outfile, paste0("results/ase_isolde/", i_cell, "_", i_timepoint, "_isolde.tsv"))
  }
}



# Compile results ---------------------------------------------------------


compiled_results <- list()
for (i_cell in unique(dds$cell_type)){
  for (i_timepoint in unique(dds$timepoint)) {
    # read isolde results table
    res <- read.table(paste0("results/ase_isolde/", i_cell, "_", i_timepoint, "_isolde.tsv"),
                      header = TRUE)

    # add cell and timepoint columns
    res$cell_type <- i_cell
    res$timepoint <- i_timepoint

    # add to list
    compiled_results[[paste0(i_cell, "_", i_timepoint)]] <- res

    # remove original file
    unlink(paste0("results/ase_isolde/", i_cell, "_", i_timepoint, "_isolde.tsv"))
  }
}

# bind all data frames
compiled_results <- do.call("rbind", compiled_results)

# export
write.csv(compiled_results,
          file = "results/ase_isolde/isolde_all.csv",
          row.names = FALSE)
