
# decompress files if zipped
localrules: prepare_reference_files # we execute this rule locally if running on HPC
rule prepare_reference_files:
  input:
    genome = config["reference"]["genome"],
    gtf = config["reference"]["annotation"]
  output:
    genome = temp("results/reference/genome.fa"),
    gtf = temp("results/reference/annotation.gtf")
  params:
    runtime = "00:30:00"
  log:
    "logs/prepare_reference_files/out.log"
  threads: 1
  run:
    if input.genome.endswith('.gz'):
      shell("zcat {input.genome} > {output.genome} 2> {log}")
    else:
      shell("cp {input.genome} > {output.genome} 2> {log}")
    
    if input.gtf.endswith('.gz'):
      shell("zcat {input.gtf} > {output.gtf} 2> {log}")
    else:
      shell("cp {input.gtf} > {output.gtf} 2> {log}")


# make diploid genome
rule diploidome:
  input:
    genome = "results/reference/genome.fa",
    gtf = "results/reference/annotation.gtf",
    vcf = config["variants"]["vcf"]
  output:
    idx = directory("results/reference/diploidome/salmon_index"), 
    genome = "results/reference/diploidome/diploid_genome.fa.gz",
    transcritps = "results/reference/diploidome/diploid_transcripts.fa.gz"
  params:
    runtime="02:00:00",
    sample = config["variants"]["sample"],
    hap1 = config["variants"]["hap1"],
    hap2 = config["variants"]["hap2"]
  log:
    "logs/reference/diploidome.log"
  threads: 10
  conda:
    "../envs/rnaseq.yaml"
  shell:
    '''
    bash workflow/scripts/diploidomes.sh \
      --ref {input.genome} \
      --gtf {input.gtf} \
      --vcf {input.vcf} \
      --sample {params.sample} \
      --hap1 {params.hap1} \
      --hap2 {params.hap2} \
      --index \
      --cpus {threads} \
      --outdir $(dirname {output.idx}) \
      > {log} 2>&1
    '''

rule star_index:
  input:
    genome = "results/reference/genome.fa",
    gtf = "results/reference/annotation.gtf"
  output:
    directory("results/reference/star_index")
  params:
    runtime="01:00:00",
    extra = config["params"]["star-index"]
  log:
    "logs/reference/star_index.log"
  threads: 10
  conda:
    "../envs/rnaseq.yaml"
  shell:
    '''
    STAR \
      --runMode genomeGenerate \
      --runThreadN {threads} \
      --genomeDir {output} \
      --genomeFastaFiles {input.genome} \
      --sjdbGTFfile {input.gtf} \
      {params.extra} \
      > {log} 2>&1
    '''
