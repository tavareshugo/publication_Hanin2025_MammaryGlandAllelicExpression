#---- Read QC ----#

def get_fastqc_input(wildcards):
  filenames = sample_info["fq1"].tolist() + sample_info["fq2"].tolist()
  
  out = [x for x in filenames if basename_without_ext(x).endswith(wildcards.readfile)]

  if len(out) != 1:
    raise ValueError("Bug.")
  else:
    return out
  

rule fastqc:
  input:
    get_fastqc_input
  output:
    "results/qc/fastqc/{readfile}_fastqc.html",
    "results/qc/fastqc/{readfile}_fastqc.zip"
  params:
    runtime = "00:30:00",
    outdir = "results/qc/fastqc/"
  threads: 4
  priority: 1
  log:
    "logs/fastqc/{readfile}.log"
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "fastqc -o {params.outdir} -t {threads} {input} > {log} 2>&1"


#---- Alignment QC ----#
# These rules are based on https://github.com/snakemake-workflows/rna-seq-star-deseq2

# RSeQC scripts supported by multiqc:
# bam_stat
# gene_body_coverage
# infer_experiment
# inner_distance
# junction_annotation
# junction_saturation
# read_distribution
# read_duplication
# read_gc


rule rseqc_gtf2bed:
  input:
    "results/reference/annotation.gtf",
  output:
    bed="results/qc/rseqc/annotation.bed",
    db=temp("results/qc/rseqc/annotation.db")
  log:
    "logs/rseqc/gtf2bed.log",
  params:
    runtime = "01:00:00"
  conda:
    "../envs/rnaseq.yaml"
  script:
    "../scripts/gtf2bed.py"


# rule rseqc_gene_body_coverage:
#   input:
#     bams = expand("results/markdup/{sample}.sortedByCoord.markdup.bam", sample = sample_info["sample"].unique()),
#     bais = expand("results/markdup/{sample}.sortedByCoord.markdup.bam.bai", sample = sample_info["sample"].unique()),
#     bed = "results/qc/rseqc/annotation.bed"
#   output:
#     "results/qc/rseqc/all_samples.geneBodyCoverage.txt"
#   params:
#     runtime = "06:00:00",
#     bamdir = "results/markdup",
#     prefix = "results/qc/rseqc/all_samples"
#   priority: 1
#   log:
#     "logs/rseqc/geneBody_coverage.log"
#   conda:
#     "../envs/rnaseq.yaml"
#   shell:
#     '''
#     geneBody_coverage.py -i {params.bamdir} -r {input.bed} -o {params.prefix} \
#       > {log} 2>&1
#     '''

rule rseqc_gene_body_coverage:
  input:
    bam = "results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai = "results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
    bed = "results/qc/rseqc/annotation.bed"
  output:
    "results/qc/rseqc/{sample}.geneBodyCoverage.txt"
  params:
    runtime = "06:00:00",
    bamdir = "results/markdup",
    prefix = "results/qc/rseqc/{sample}"
  threads: 4  # for the HPC to make sure there's enough memory allocated
  priority: 1
  log:
    "logs/rseqc/geneBody_coverage/{sample}.log"
  conda:
    "../envs/rnaseq.yaml"
  shell:
    '''
    geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params.prefix} \
      > {log} 2>&1
    '''


rule rseqc_junction_annotation:
  input:
    bam="results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai="results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
    bed="results/qc/rseqc/annotation.bed"
  output:
    "results/qc/rseqc/{sample}.junctionanno.junction.bed"
  params:
    runtime = "01:00:00",
    extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
    prefix="results/qc/rseqc/{sample}.junctionanno"
  log:
    "logs/rseqc/junction_annotation/{sample}.log"
  threads: 4  # for the HPC to make sure there's enough memory allocated
  priority: 1
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
    "> {log} 2>&1"


rule rseqc_junction_saturation:
  input:
    bam="results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai="results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
    bed="results/qc/rseqc/annotation.bed",
  output:
    "results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
  params:
    runtime = "01:00:00",
    extra=r"-q 255",
    prefix="results/qc/rseqc/{sample}.junctionsat"
  threads: 4  # for the HPC to make sure there's enough memory allocated
  priority: 1
  log:
    "logs/rseqc/junction_saturation/{sample}.log",
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
    "> {log} 2>&1"


rule rseqc_stat:
  input:
    bam="results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai="results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
  output:
    "results/qc/rseqc/{sample}.stats.txt",
  params:
    runtime = "01:00:00",
  threads: 4  # for the HPC to make sure there's enough memory allocated
  priority: 1
  log:
    "logs/rseqc/bam_stat/{sample}.log",
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "bam_stat.py -i {input.bam} > {output} 2> {log}"


rule rseqc_infer:
  input:
    bam="results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai="results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
    bed="results/qc/rseqc/annotation.bed",
  output:
    "results/qc/rseqc/{sample}.infer_experiment.txt",
  priority: 1
  params:
    runtime = "01:00:00",
  threads: 4  # for the HPC to make sure there's enough memory allocated
  log:
    "logs/rseqc/infer_experiment/{sample}.log",
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "infer_experiment.py -r {input.bed} -i {input.bam} -s 1000000 > {output} 2> {log}"


rule rseqc_innerdis:
  input:
    bam="results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai="results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
    bed="results/qc/rseqc/annotation.bed",
  output:
    "results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
  params:
    runtime = "01:00:00",
    prefix="results/qc/rseqc/{sample}.inner_distance_freq"
  threads: 4  # for the HPC to make sure there's enough memory allocated
  priority: 1
  log:
    "logs/rseqc/inner_distance/{sample}.log",
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
  input:
    bam="results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai="results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
    bed="results/qc/rseqc/annotation.bed",
  output:
    "results/qc/rseqc/{sample}.readdistribution.txt",
  params:
    runtime = "01:00:00",
  threads: 4  # for the HPC to make sure there's enough memory allocated
  priority: 1
  log:
    "logs/rseqc/read_distribution/{sample}.log",
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
  input:
    bam="results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai="results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
  output:
    "results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
  params:
    runtime = "01:00:00",
    prefix="results/qc/rseqc/{sample}.readdup"
  threads: 4  # for the HPC to make sure there's enough memory allocated
  priority: 1
  log:
    "logs/rseqc/read_duplication/{sample}.log",
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "read_duplication.py -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
  input:
    bam="results/markdup/{sample}.sortedByCoord.markdup.bam",
    bai="results/markdup/{sample}.sortedByCoord.markdup.bam.bai",
  output:
    "results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
  params:
    runtime = "01:00:00",
    prefix="results/qc/rseqc/{sample}.readgc"
  threads: 4  # for the HPC to make sure there's enough memory allocated
  priority: 1
  log:
    "logs/rseqc/read_GC/{sample}.log",
  conda:
    "../envs/rnaseq.yaml"
  shell:
    "read_GC.py -i {input.bam} -o {params.prefix} > {log} 2>&1"



#---- MultiQC ----#

rule multiqc:
  input:
    expand(
      "results/sortmerna/{id}", 
      id = sample_info["sample"].unique()
    ),
    expand(
      "results/qc/fastqc/{readfile}_fastqc.zip", 
      readfile = [basename_without_ext(i) for i in sample_info["fq1"].tolist() + sample_info["fq2"].tolist()]
    ),
    expand(
      "logs/cutadapt/{id}.log", 
      id = sample_info["id"]
    ), 
    expand(
      "results/salmon_diploid/{sample}", 
      sample = sample_info["sample"].unique()
    ), 
    expand(
      "results/markdup/{sample}.metrics.txt",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.geneBodyCoverage.txt",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.junctionanno.junction.bed",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.infer_experiment.txt",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.stats.txt",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.readdistribution.txt",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
      sample = sample_info["sample"].unique(),
    ),
    expand(
      "logs/rseqc/junction_annotation/{sample}.log",
      sample = sample_info["sample"].unique(),
    ),
  output:
    "results/qc/multiqc_report.html",
  params:
    runtime = "01:00:00",
    qcdirs = ["results/qc/fastqc", "logs/cutadapt/", "results/salmon_diploid/", "results/sortmerna/", "results/markdup/", "results/qc/rseqc/"]
  conda:
    "../envs/rnaseq.yaml"
  log:
    "logs/multiqc.log",
  shell:
    '''
    multiqc -o $(dirname {output}) -n multiqc_report.html {params.qcdirs} \
      > {log} 2>&1
    '''
