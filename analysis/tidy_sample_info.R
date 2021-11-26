# Tidy sample info table

library(tidyverse)
library(readxl)
library(janitor)
theme_set(theme_minimal() + theme(text = element_text(size = 16)))


# Read sample info --------------------------------------------------------

# master table with information about each sample
sample_info <- read_excel("data/raw/RNA seq - lane and index plan - final.xlsx",
                          sheet = "All samples in the experiment") %>%
  clean_names()

# tidy
sample_info <- sample_info %>%
  # a rogue column that was in the excel sheet
  select(-column1) %>%
  # lowercase some character columns to avoid case typos
  mutate(across(c(stage, cell_type, cross, rin_category), tolower)) %>%
  # convert variables that should be numeric
  mutate(across(c(rin, rna_qubit_concentration_ng_ul,
                  library_concentration_from_qubit_ng_ul,
                  average_library_size_according_to_the_bioanalyzer_bp),
                as.numeric)) %>%
  # convert date variable
  mutate(library_prep_date = as.Date(library_prep_date, "%d.%m.%y")) %>%
  # lane numbers are not all correct as these were adjusted after
  select(-lane_number)



# Read lane info ----------------------------------------------------------

# lane information in several separate sheets of the excel
lanes <- list(
  pool_1_4 = read_excel("data/raw/RNA seq - lane and index plan - final.xlsx",
                        sheet = "pooling lanes 1-4",
                        na = c("", "-", "didn't work")),
  pool_5_8 = read_excel("data/raw/RNA seq - lane and index plan - final.xlsx",
                        sheet = "pooling lanes 5-8",
                        na = c("", "-", "didn't work")),
  pool_9_12 = read_excel("data/raw/RNA seq - lane and index plan - final.xlsx",
                         sheet = "pooling lanes 9-12",
                         na = c("", "-", "didn't work")),
  pool_13_16 = read_excel("data/raw/RNA seq - lane and index plan - final.xlsx",
                          sheet = "pooling lanes 13-16",
                          na = c("", "-", "didn't work")),
  pool_17_20 = read_excel("data/raw/RNA seq - lane and index plan - final.xlsx",
                          sheet = "pooling lanes 17-20",
                          na = c("", "-", "didn't work"))
)

# tidy
lanes <- lanes %>%
  map(clean_names) %>%
  map(~ select(.x,
               sample_number:animal_id, lane_number,
               index_true_seq, slx_identifier)) %>%
  bind_rows(.id = "pool") %>%
  # lowercase some character columns to avoid case typos
  mutate(across(c(stage, cell_type, cross), tolower)) %>%
  # change index name to match fastq file names
  mutate(index_true_seq = str_replace(index_true_seq, "-", "_"))

# check everything is in common between the two tables
anti_join(
  sample_info,
  lanes,
  by = c("sample_number", "stage", "cell_type", "cross", "animal_id")
)

# remove samples that were not sequenced
lanes <- lanes %>%
  drop_na(slx_identifier, index_true_seq)


# FASTQ file information --------------------------------------------

# Read file (obtained with)
# find -maxdepth 3 -mindepth 3 -type f -name "SLX-*r_1.fq.gz" | grep "Nova"
fqs <- read_csv("data/raw/fastq_list.csv",
                col_names = "fq1")

fqs <- fqs %>%
  # add information about read 2
  mutate(fq2 = str_replace(fq1, "r_1.fq.gz$", "r_2.fq.gz")) %>%
  # get run and slx_identifier
  mutate(basedir = fq1 %>% dirname() %>% str_remove("^\\./")) %>%
  separate(basedir, c("run", "slx_identifier"), sep = "/") %>%
  # get index_true_seq
  mutate(index_true_seq = fq1 %>% basename %>%
           str_split("\\.") %>% map(pluck, 2) %>% unlist()) %>%
  # reorder columns
  select(slx_identifier, index_true_seq, run, fq1, fq2)


# Join tables --------------------------------------------------------

info_tidy <- fqs %>%
  # there was one sample for which we have no fastq files (SLX-19101 index D711rna-D505rna)
  # check with anti_join(lanes, fqs)
  left_join(lanes, by = c("slx_identifier", "index_true_seq")) %>%
  # there were several samples which were not made into libraries
  # check with anti_join(sample_info, lanes)
  left_join(sample_info,
            by = c("sample_number", "stage",
                   "cell_type", "cross", "animal_id"))

# tidy up a bit more
# add timepoint information (for easier naming of IDs)
info_tidy <- tribble(
  ~stage, ~timepoint,
  "nulliparous", "t0",
  "gestation d5.5", "t1",
  "gestation d9.5", "t2",
  "gestation d14.5", "t3",
  "lactation d5", "t4",
  "lactation d10", "t5",
  "lactation d15", "t6",
  "involution d1", "t7",
  "involution d6", "t8",
  "involution d14", "t9"
) %>%
  full_join(info_tidy, by = "stage") %>%
  # create a descriptive sample name
  mutate(sample_name = paste(str_remove(cell_type, " "),
                             timepoint,
                             cross,
                             str_remove(animal_id, ".\\(.*"), sep = "_")) %>%
  # reorder columns for convenience
  select(sample_name, sample_number, cell_type, stage, timepoint, cross, animal_id,
         slx_identifier, index_true_seq, run, fq1, fq2, pool, everything())

# add unique id and library name (which is the sample_number in this case)
info_tidy <- info_tidy %>%
  mutate(id = paste(slx_identifier, index_true_seq, sep = "."),
         library = paste0("LIB", sample_number)) %>%
  select(id, sample_name, library,
         fq1, fq2, everything())

# information about the sequencing samples (reads)
write_csv(info_tidy, "read_info.csv")

# information about the biological samples only
info_tidy %>%
  distinct(sample_name, cell_type, stage,
           timepoint, cross, animal_id, rin) %>%
  write_csv("sample_info.csv")


# Exploratory plots -------------------------------------------------------

# how many samples per combination
info_tidy %>%
  group_by(stage, cell_type, cross) %>%
  summarise(nanimal = n_distinct(animal_id),
            nlane = n_distinct(lane_number),
            nsample = n_distinct(sample_number)) %>%
  ggplot(aes(nanimal)) + geom_bar() +
  labs(caption = "number of replicates per factor combination (stage x cell x cross)")

# sort by replicate number
info_tidy %>%
  distinct(sample_name, cell_type, stage, cross) %>%
  count(stage, cell_type, cross) %>%
  arrange(n, cell_type, stage, cross)

info_tidy %>%
  distinct(sample_name, cell_type, stage, cross, sample_number) %>%
  group_by(cell_type, stage, cross) %>%
  filter(n() < 4)

# check the samples missing sequence data
anti_join(sample_info, info_tidy) %>%
  select(sample_number, stage, cell_type, cross, animal_id)

# samples re-sequenced
info_tidy %>%
  group_by(stage, cell_type, cross, animal_id) %>%
  summarise(nlanes = n_distinct(lane_number),
            nfqs = n_distinct(fq1),
            sample_number = unique(sample_number),
            nindex = n_distinct(index_true_seq)) %>%
  filter(nfqs > 1) %>%
  ggplot(aes(nfqs, paste(cell_type, stage, cross, animal_id))) +
  geom_col(aes(fill = factor(nindex))) +
  labs(x = "number fastq files", y = "", caption = "samples that were resequenced")

# do the resequenced samples have the same index?
info_tidy %>%
  group_by(stage, cell_type, cross, animal_id) %>%
  summarise(nlanes = n_distinct(lane_number),
            nfqs = n_distinct(fq1),
            sample_number = unique(sample_number),
            nindex = n_distinct(index_true_seq)) %>%
  filter(nfqs > 1) %>%
  ungroup() %>%
  count(nindex)

# how pools are distributed across runs
info_tidy %>%
  count(pool, run) %>%
  ggplot(aes(pool, run)) +
  geom_text(aes(label = n)) +
  theme_minimal()


