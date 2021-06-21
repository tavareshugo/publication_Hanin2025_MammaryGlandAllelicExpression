
rule cutadapt:
  input:
    fq1 = lambda wildcards: sample_info.loc[sample_info["id"] == wildcards.id, "fq1"],
    fq2 = lambda wildcards: sample_info.loc[sample_info["id"] == wildcards.id, "fq2"]
  output:
    fq1 = temp("results/trimmed/{id}.R1.fq.gz"),
    fq2 = temp("results/trimmed/{id}.R2.fq.gz")
  params:
    runtime="02:00:00",
    # note: had to use an anonymous function here due to an obscure error from snakemake
    extras = lambda wildcards: config["params"]["cutadapt-pe"]
  log:
    "logs/cutadapt/{id}.log"
  threads: 10
  conda:
    "../envs/rnaseq.yaml"
  shell:
    '''
    cutadapt {params.extras} \
      --cores {threads} \
      -o {output.fq1} -p {output.fq2} \
      {input.fq1} {input.fq2} \
      >{log} 2>&1
    '''

