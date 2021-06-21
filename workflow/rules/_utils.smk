
def get_trimmed_fq1(wildcards):
  fqs = sample_info.loc[sample_info["sample"] == wildcards.sample, "id"].tolist()
  return expand("results/trimmed/{fqs}.R1.fq.gz", fqs = fqs)

def get_trimmed_fq2(wildcards):
  fqs = sample_info.loc[sample_info["sample"] == wildcards.sample, "id"].tolist()
  return expand("results/trimmed/{fqs}.R2.fq.gz", fqs = fqs)

# from: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastqc.html
def basename_without_ext(file_path):
    """Returns basename of file path, without the file extension."""

    base = os.path.basename(file_path)
    # Remove file extension(s) (similar to the internal fastqc approach)
    base = re.sub("\\.gz$", "", base)
    base = re.sub("\\.bz2$", "", base)
    base = re.sub("\\.txt$", "", base)
    base = re.sub("\\.fastq$", "", base)
    base = re.sub("\\.fq$", "", base)
    base = re.sub("\\.sam$", "", base)
    base = re.sub("\\.bam$", "", base)

    return base