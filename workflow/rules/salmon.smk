
# SALMON quantification on the diploid genome
rule salmon_diploid:
  input:
    idx = "results/reference/diploidome/salmon_index",
    fq1 = get_trimmed_fq1,
    fq2 = get_trimmed_fq2
  output:
    directory("results/salmon_diploid/{sample}")
  params:
    runtime = "02:00:00",
    extras = config["params"]["salmon"]
  log:
    "logs/salmon_diploid/{sample}.log"
  threads: 10
  conda:
    "../envs/rnaseq.yaml"
  shell:
    '''
    salmon quant \
      {params.extras} \
      -i {input.idx} \
      -1 {input.fq1} \
      -2 {input.fq2} \
      --output {output} \
      --seqBias --gcBias --posBias --validateMappings \
      --threads {threads} \
      >{log} 2>&1
    '''
  

# SALMON quantification on the regular genome
rule salmon_regular:
  input:
    idx = "results/reference/salmon_index",
    fq1 = get_trimmed_fq1,
    fq2 = get_trimmed_fq2
  output:
    directory("results/salmon_regular/{sample}")
  params:
    runtime = "02:00:00",
    extras = config["params"]["salmon"]
  log:
    "logs/salmon_regular/{sample}.log"
  threads: 10
  conda:
    "../envs/salmon.yaml"
  shell:
    '''
    salmon quant \
      {params.extras} \
      -i {input.idx} \
      -1 {input.fq1} \
      -2 {input.fq2} \
      --output {output} \
      --seqBias --gcBias --posBias --validateMappings \
      --threads {threads} \
      >{log} 2>&1
    '''
  
