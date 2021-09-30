
rule ase_isolde:
  input:
    "results/DESeqDataSet/dds_gene_allele.rds"
  output:
    "results/ase_isolde/isolde_all.csv"
  params:
    runtime = "12:00:00",
  threads: 10
  log:
    "logs/isolde/ase_isolde.log"
  conda:
    "../envs/isolde.yaml"
  shell:
    "Rscript workflow/scripts/ase_isolde.R > {log} 2>&1"


rule ase_isolde_strain:
  input:
    "results/DESeqDataSet/dds_gene_allele.rds"
  output:
    "results/ase_isolde/isolde_all_strain.csv"
  params:
    runtime = "12:00:00",
  threads: 10
  log:
    "logs/isolde/ase_isolde_strain.log"
  conda:
    "../envs/isolde.yaml"
  shell:
    "Rscript workflow/scripts/ase_isolde_strain.R > {log} 2>&1"

rule ase_isolde_isoform:
  input:
    "results/DESeqDataSet/dds_isoform_allele.rds"
  output:
    "results/ase_isolde_isoform/isolde_isoform_all.csv"
  params:
    runtime = "12:00:00",
  threads: 10
  log:
    "logs/isolde/ase_isolde_isoform.log"
  conda:
    "../envs/isolde.yaml"
  shell:
    "Rscript workflow/scripts/ase_isolde_isoform.R > {log} 2>&1"
