rule DESeqDataSet:
  input:
    expand(
      "results/salmon_diploid/{sample}",
      sample = sample_info["sample"].unique()
    )
  output:
    "results/DESeqDataSet/dds_gene_allele.rds",
    "results/DESeqDataSet/dds_gene_regular.rds",
    "results/DESeqDataSet/dds_isoform_allele.rds",
    "results/DESeqDataSet/dds_isoform_regular.rds"
  params:
    runtime = "12:00:00",
  threads: 10
  log:
    "logs/deseq2/DESeqDataSet.log"
  conda:
    "../envs/deseq2.yaml"
  shell:
    "Rscript workflow/scripts/deseq2.R > {log} 2>&1"


rule ase_deseq2:
  input:
    "results/DESeqDataSet/dds_gene_allele.rds"
  output:
    "results/ase_deseq2/dds_ase.rds"
  params:
    runtime = "36:00:00",
  threads: 10
  log:
    "logs/deseq2/ase_deseq2.log"
  conda:
    "../envs/deseq2.yaml"
  shell:
    "Rscript workflow/scripts/ase_deseq.R > {log} 2>&1"


rule diffexp:
  input:
    "results/DESeqDataSet/dds_gene_regular.rds"
  output:
    "results/diffexp/lfc_timepoints.csv",
    "results/diffexp/lfc_cells.csv",
    "results/diffexp/dds.rds"
  params:
    runtime = "36:00:00"
  threads: 10
  log:
    "logs/deseq2/diffexp.log"
  conda:
    "../envs/deseq2.yaml"
  shell:
    "Rscript workflow/scripts/diffexp.R > {log} 2>&1"
