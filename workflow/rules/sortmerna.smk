
rule sortmerna:
  input:
    fq1 = get_trimmed_fq1,
    fq2 = get_trimmed_fq2,
    rrnadb = config["rRNA"]
  output:
    directory("results/sortmerna/{sample}"),
  params:
    runtime="36:00:00"
  log:
    "logs/sortmerna/{sample}.log"
  threads: 20
  conda:
    "../envs/rnaseq.yaml"
  shell:
    '''
    # concatenate in case there's multiple files
    # cat {input.fq1} > results/sortmerna/{wildcards.sample}.R1.tmp
    # cat {input.fq2} > results/sortmerna/{wildcards.sample}.R2.tmp
    
    # sample 10M reads
    mkdir -p results/sortmerna/{wildcards.sample}/temp/
    seqtk sample -s 20210607 <(cat {input.fq1}) 10000000 > results/sortmerna/{wildcards.sample}/temp/{wildcards.sample}.R1.fq
    seqtk sample -s 20210607 <(cat {input.fq2}) 10000000 > results/sortmerna/{wildcards.sample}/temp/{wildcards.sample}.R2.fq

    # merge into interleaved file with basename only (this makes multiqc output neater)
    seqtk mergepe results/sortmerna/{wildcards.sample}/temp/{wildcards.sample}.R1.fq results/sortmerna/{wildcards.sample}/temp/{wildcards.sample}.R2.fq > results/sortmerna/{wildcards.sample}/temp/{wildcards.sample}

    # identify rRNA reads
    sortmerna \
      --ref {input.rrnadb} \
      -L 10 \
      --workdir results/sortmerna/{wildcards.sample}/ \
      --reads results/sortmerna/{wildcards.sample}/temp/{wildcards.sample} \
      --paired \
      --threads {threads} \
      >{log} 2>&1

    # --reads results/sortmerna/{wildcards.sample}.R1.tmp \
    # --reads results/sortmerna/{wildcards.sample}.R2.tmp \
    # --paired_in \
    # --fastx \
    # --aligned results/sortmerna/{wildcards.sample}/out/{wildcards.sample}_aligned \
    # --other results/sortmerna/{wildcards.sample}/out/{wildcards.sample}_other \
    # --out2 \

    # clean all output files apart from log, to save space
    rm results/sortmerna/{wildcards.sample}/out/aligned.blast
    rm -r results/sortmerna/{wildcards.sample}/temp
    rm -r results/sortmerna/{wildcards.sample}/idx
    rm -r results/sortmerna/{wildcards.sample}/kvdb
    rm -r results/sortmerna/{wildcards.sample}/readb
    '''
