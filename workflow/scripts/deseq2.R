packages <- c("valr", "DESeq2", "biomaRt", "tximport", "vroom", "dplyr", "tibble")
install.packages(setdiff(packages, rownames(installed.packages())), repos = "https://cran.ma.imperial.ac.uk/")
library(valr)
library(DESeq2)
library(biomaRt)
library(tximport)
library(vroom)
library(dplyr)
#library(tidyr)
library(tibble)

# create output directory 
if (!dir.exists("results/DESeqDataSet/")) dir.create("results/DESeqDataSet/", recursive = TRUE)


# Prepare tx2gene tables --------------------------------------------------

# transcript-to-gene file
tx2gene <- vroom("data/external/reference/transcript2gene.tsv",
                 col_types = c("cc")) %>%
  # ensure columns are in the correct order for tximport
  select(transcript, gene)

# make a diploid version of the table
tx2gene_diploid <- tx2gene %>%
  mutate(transcript_diploid = paste0("C.", transcript),
         gene_diploid = paste0("C.", gene))
tx2gene_diploid <- tx2gene %>%
  mutate(transcript_diploid = paste0("B.", transcript),
         gene_diploid = paste0("B.", gene)) %>%
  bind_rows(tx2gene_diploid)

# list of input files
salmon_files <- list.files("results/salmon_diploid/",
                           pattern = "quant.sf",
                           full.names = TRUE,
                           recursive = TRUE)
names(salmon_files) <- basename(dirname(salmon_files))



# Prepare colData ---------------------------------------------------------

samples <- vroom("read_info.csv")
# samples <- samples %>% filter(sample %in% names(salmon_files))
samples <- samples %>%
  distinct(sample, cell_type, stage, timepoint, cross, animal_id)

# encode animal ID as a nested level
samples <- samples %>%
  group_by(cell_type, stage) %>%
  mutate(id_nested = factor(1:n())) %>%
  group_by(cell_type, stage, cross) %>%
  mutate(id_nested_within_cross = factor(1:n())) %>%
  ungroup()

# encode variables as factors (needed for DESeq design)
samples <- samples %>%
  mutate(cell_type = factor(cell_type),
         stage = reorder(stage, rank(timepoint)),
         timepoint = factor(timepoint),
         cross = factor(cross),
         animal_id = factor(animal_id)) %>%
  column_to_rownames("sample")

# ensure it's the same order as files
samples <- samples[names(salmon_files), ]


# Prepare rowData ---------------------------------------------------------

# get gene annotation - this was mapped to GRCm38 assembly
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",
                   host="nov2020.archive.ensembl.org")

# gene information
gene_annot <- getBM(attributes = c("ensembl_gene_id",
                                   "chromosome_name",#"strand",
                                   "start_position", "end_position",
                                   "external_gene_name",
                                   "percentage_gene_gc_content",
                                   "gene_biotype"),
                    mart = mart)

# transcript annotations
transcript_annot <- getBM(attributes = c("ensembl_gene_id",
                                         "ensembl_transcript_id",
                                         "ensembl_exon_id",
                                         "exon_chrom_start",
                                         "exon_chrom_end",
                                         "chromosome_name",
                                         "strand",
                                         "transcript_biotype"),
                          mart = mart)

# convert strand to +/- which is required for GRanges object
transcript_annot$strand <- case_when(transcript_annot$strand == -1 ~ "-",
                                     transcript_annot$strand == 1 ~ "+",
                                     TRUE ~ NA_character_)



# SNP information ---------------------------------------------------------

# read SNPs/INDELs
snps <- vroom("data/external/strains/B6xCAST.snps_and_indels.merged.bed.gz",
              col_names = c("chrom", "start", "end", "name"),
              col_types = c("ciic"))

snps <- snps %>% mutate(indel = grepl("INDEL", name)) %>% select(-name)

# create gene intervals table
gene_intervals <- transcript_annot %>%
  # rename for valr
  rename(chrom = chromosome_name, start = exon_chrom_start, end = exon_chrom_end) %>%
  select(chrom, start, end, everything()) %>%
  # turn start coordinates to zero-based
  mutate(start = start - 1) %>%
  arrange(chrom, start, end) %>%
  # merge exons
  group_by(ensembl_gene_id) %>%
  bed_merge() %>%
  ungroup()

# calculate coverage with SNPs
gene_variant_coverage <- gene_intervals %>%
  bed_coverage(snps %>% filter(!indel)) %>%
  group_by(ensembl_gene_id) %>%
  summarise(n_snps = sum(.ints))

# calculate coverage with INDELs
gene_variant_coverage <- gene_intervals %>%
  bed_coverage(snps %>% filter(indel)) %>%
  group_by(ensembl_gene_id) %>%
  summarise(n_indels = sum(.ints)) %>%
  full_join(gene_variant_coverage, by = "ensembl_gene_id")

# add the two together
gene_variant_coverage <- gene_variant_coverage %>% mutate(n_variants = n_snps + n_indels)

# add to gene annotation table
gene_annot <- gene_annot %>%
  full_join(gene_variant_coverage, by = "ensembl_gene_id") %>%
  column_to_rownames("ensembl_gene_id")


# Same thing for transcripts
# create transcript intervals table
transcript_intervals <- transcript_annot %>%
  # rename for valr
  rename(chrom = chromosome_name, start = exon_chrom_start, end = exon_chrom_end) %>%
  select(chrom, start, end, everything()) %>%
  # turn start coordinates to zero-based
  mutate(start = start - 1) %>%
  arrange(chrom, start, end) %>%
  # merge exons
  group_by(ensembl_transcript_id) %>%
  bed_merge() %>%
  ungroup()

# calculate coverage with SNPs
transcript_variant_coverage <- transcript_intervals %>%
  bed_coverage(snps %>% filter(!indel)) %>%
  group_by(ensembl_transcript_id) %>%
  summarise(n_snps = sum(.ints))

# calculate coverage with INDELs
transcript_variant_coverage <- transcript_intervals %>%
  bed_coverage(snps %>% filter(indel)) %>%
  group_by(ensembl_transcript_id) %>%
  summarise(n_indels = sum(.ints)) %>%
  full_join(transcript_variant_coverage, by = "ensembl_transcript_id")

# add the two together
transcript_variant_coverage <- transcript_variant_coverage %>% mutate(n_variants = n_snps + n_indels)

# add to transcript annotation table
transcript_annot <- transcript_annot %>%
  full_join(transcript_variant_coverage, by = "ensembl_transcript_id")

# create a data.frame to add as rowData
transcript_rowdata <- transcript_annot %>%
  group_by(ensembl_transcript_id, ensembl_gene_id, chromosome_name, strand, transcript_biotype) %>%
  summarise(n_indels = sum(n_indels, na.rm = TRUE),
            n_snps = sum(n_snps, na.rm = TRUE),
            n_variants = sum(n_variants, na.rm = TRUE),
            start = min(exon_chrom_start),
            end = max(exon_chrom_end)) %>%
  ungroup() %>%
  column_to_rownames("ensembl_transcript_id")

# clean environment
rm(snps, gene_intervals, gene_variant_coverage, transcript_intervals, transcript_variant_coverage,
   mart)


# Gene-level Allele Specific -----------------------------------------------------------

txi_gene_allele <- tximport(salmon_files, type = "salmon",
                            tx2gene = tx2gene_diploid %>%
                              select(transcript_diploid, gene_diploid) %>%
                              distinct(),
                            importer = function(x) vroom::vroom(x,
                                                                progress = FALSE,
                                                                col_types = c("ciddd")))

dds_gene_allele <- DESeqDataSetFromTximport(txi_gene_allele,
                                            colData = samples,
                                            design = ~ 1)

# prepare allele-specific dds object
dds_gene_allele_b <- dds_gene_allele[grep("B\\.", rownames(dds_gene_allele)), ]
dds_gene_allele_c <- dds_gene_allele[grep("C\\.", rownames(dds_gene_allele)), ]
# add suffix to colnames
colnames(dds_gene_allele_b) <- paste0(colnames(dds_gene_allele_b), "_B")
colnames(dds_gene_allele_c) <- paste0(colnames(dds_gene_allele_c), "_C")
# remove prefix from rownames
rownames(dds_gene_allele_b) <- gsub("B\\.", "", rownames(dds_gene_allele_b))
rownames(dds_gene_allele_c) <- gsub("C\\.", "", rownames(dds_gene_allele_c))
# check rownames are all matching and in the same order
if(!all(rownames(dds_gene_allele_b) %in% rownames(dds_gene_allele_c)) |
   !all(rownames(dds_gene_allele_c) %in% rownames(dds_gene_allele_b))){
  stop("rownames issue with DDS object!")
}
dds_gene_allele_b <- dds_gene_allele_b[rownames(dds_gene_allele_c), ]
# bind together
dds_gene_allele <- cbind(dds_gene_allele_b, dds_gene_allele_c)
# add allele columns to colData
colData(dds_gene_allele)$allele_strain <- gsub(".*_", "", colnames(dds_gene_allele))
colData(dds_gene_allele)$allele_parent <- ifelse(tolower(dds_gene_allele$allele_strain) == substring(dds_gene_allele$cross, 1, 1), "maternal", "paternal")
rm(dds_gene_allele_c, dds_gene_allele_b)

# add rowData
rowData(dds_gene_allele) <- DataFrame(gene_annot[rownames(dds_gene_allele), ])

# # normalised counts
# assay(dds_gene_allele, "vst_blind") <- varianceStabilizingTransformation(counts(dds_gene_allele),
#                                                                          blind = TRUE)

# save object
saveRDS(dds_gene_allele,
        "results/DESeqDataSet/dds_gene_allele.rds")

# clean environment
rm(dds_gene_allele, txi_gene_allele)


# Gene-level Regular -----------------------------------------------------------

txi_gene_regular <- tximport(salmon_files, type = "salmon",
                             tx2gene = tx2gene_diploid %>%
                               select(transcript_diploid, gene) %>%
                               distinct(),
                             importer = function(x) vroom::vroom(x,
                                                                 progress = FALSE,
                                                                 col_types = c("ciddd")))

dds_gene_regular <- DESeqDataSetFromTximport(txi_gene_regular,
                                             colData = samples,
                                             design = ~ 1)

# add rowData
rowData(dds_gene_regular) <- DataFrame(gene_annot[rownames(dds_gene_regular), ])

# normalised counts
assay(dds_gene_regular, "vst_blind") <- varianceStabilizingTransformation(counts(dds_gene_regular),
                                                                          blind = TRUE)

# assay(dds_gene_regular, "rlog_blind") <- rlog(counts(dds_gene_regular),
#                                               blind = TRUE)

# save object
saveRDS(dds_gene_regular,
        "results/DESeqDataSet/dds_gene_regular.rds")

# clean environment
rm(txi_gene_regular, dds_gene_regular)



# Isoform-level Allele Specific -----------------------------------------------------------

# in this case we turn off summarisation with txOut = TRUE
txi_isoform_allele <- tximport(salmon_files, type = "salmon",
                               tx2gene = tx2gene_diploid %>%
                                 select(transcript_diploid, gene_diploid) %>%
                                 distinct(),
                               importer = function(x) vroom::vroom(x,
                                                                   progress = FALSE,
                                                                   col_types = c("ciddd")),
                               txOut = TRUE)

dds_isoform_allele <- DESeqDataSetFromTximport(txi_isoform_allele,
                                               colData = samples,
                                               design = ~ 1)


# prepare allele-specific dds object
# separate into two sets
dds_isoform_allele_b <- dds_isoform_allele[grep("B\\.", rownames(dds_isoform_allele)), ]
dds_isoform_allele_c <- dds_isoform_allele[grep("C\\.", rownames(dds_isoform_allele)), ]
# add suffix to colnames
colnames(dds_isoform_allele_b) <- paste0(colnames(dds_isoform_allele_b), "_B")
colnames(dds_isoform_allele_c) <- paste0(colnames(dds_isoform_allele_c), "_C")
# remove prefix from rownames
rownames(dds_isoform_allele_b) <- gsub("B\\.", "", rownames(dds_isoform_allele_b))
rownames(dds_isoform_allele_c) <- gsub("C\\.", "", rownames(dds_isoform_allele_c))
# check rownames are all matching and in the same order
if(!all(rownames(dds_isoform_allele_b) %in% rownames(dds_isoform_allele_c)) |
   !all(rownames(dds_isoform_allele_c) %in% rownames(dds_isoform_allele_b))){
  stop("rownames issue with DDS object!")
}
dds_isoform_allele_b <- dds_isoform_allele_b[rownames(dds_isoform_allele_c), ]
# bind together
dds_isoform_allele <- cbind(dds_isoform_allele_b, dds_isoform_allele_c)
# add allele columns to colData
colData(dds_isoform_allele)$allele_strain <- gsub(".*_", "", colnames(dds_isoform_allele))
colData(dds_isoform_allele)$allele_parent <- ifelse(tolower(dds_isoform_allele$allele_strain) == substring(dds_isoform_allele$cross, 1, 1), "maternal", "paternal")
rm(dds_isoform_allele_c, dds_isoform_allele_b)

# add rowData
rowData(dds_isoform_allele) <- DataFrame(transcript_rowdata[rownames(dds_isoform_allele), ])


# save object
saveRDS(dds_isoform_allele,
        "results/DESeqDataSet/dds_isoform_allele.rds")

# clean environment
rm(txi_isoform_allele, dds_isoform_allele)



# Isoform-level Regular -----------------------------------------------------------

txi_isoform_regular <- tximport(salmon_files, type = "salmon",
                             tx2gene = tx2gene_diploid %>%
                               select(transcript_diploid, transcript) %>%
                               distinct(),
                             importer = function(x) vroom::vroom(x,
                                                                 progress = FALSE,
                                                                 col_types = c("ciddd")))

dds_isoform_regular <- DESeqDataSetFromTximport(txi_isoform_regular,
                                             colData = samples,
                                             design = ~ 1)


# add rowData
rowData(dds_isoform_regular) <- DataFrame(transcript_rowdata[rownames(dds_isoform_regular), ])

# normalised counts
assay(dds_isoform_regular, "vst_blind") <- varianceStabilizingTransformation(counts(dds_isoform_regular),
                                                                          blind = TRUE)

# assay(dds_isoform_regular, "rlog_blind") <- rlog(counts(dds_isoform_regular),
#                                                  blind = TRUE)

# save object
saveRDS(dds_isoform_regular,
        "results/DESeqDataSet/dds_isoform_regular.rds")

