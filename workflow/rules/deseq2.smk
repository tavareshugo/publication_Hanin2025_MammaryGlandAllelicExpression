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
