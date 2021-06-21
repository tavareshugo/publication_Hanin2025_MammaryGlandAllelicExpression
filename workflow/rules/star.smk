def get_read_group(wildcards):
  '''
  Return a string with read group information for output alignment.
  This is done so that MarkDuplicates correctly deduplicates at the library level.
  See section 3.2 of STAR manual for how read group information is defined.
  '''
  sample_id = sample_info.loc[sample_info["sample"] == wildcards.sample, "id"]
  sample_name = sample_info.loc[sample_info["sample"] == wildcards.sample, "sample"]
  library_id = sample_info.loc[sample_info["sample"] == wildcards.sample, "library"]
  
  read_group = []
  for i in zip(sample_id, sample_name, library_id):
    # note these are separated by TAB character (not space)
    read_group.append("ID:{}	SM:{}	LB:{}	PL:illumina".format(i[0], i[1], i[2]))
    
  return " , ".join(read_group)


rule star:
  input:
    idx = "results/reference/star_index",
    fq1 = get_trimmed_fq1,
    fq2 = get_trimmed_fq2
  output:
    bam = temp("results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"),
    bai = temp("results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai")
  params:
    runtime = "03:00:00",
    rg = get_read_group,
    # format input for STAR
    fq1 = lambda wildcards: ",".join(get_trimmed_fq1(wildcards)),
    fq2 = lambda wildcards: ",".join(get_trimmed_fq2(wildcards))
  log:
    "logs/star/{sample}.log"
  threads: 10
  conda:
    "../envs/rnaseq.yaml"
  shell:
    '''
    STAR \
      --runMode alignReads \
      --quantMode GeneCounts \
      --outSAMtype BAM SortedByCoordinate \
      --runThreadN {threads} \
      --genomeDir {input.idx} \
      --readFilesCommand zcat \
      --readFilesIn {params.fq1} {params.fq2} \
      --outSAMattrRGline {params.rg} \
      --outFileNamePrefix results/star/{wildcards.sample}/{wildcards.sample}. \
      >{log} 2>&1

    samtools index {output.bam} >>{log} 2>&1
    '''


rule markdup:
  input:
    bam = "results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    bai = "results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
  output:
    bam = temp("results/markdup/{sample}.sortedByCoord.markdup.bam"),
    bai = temp("results/markdup/{sample}.sortedByCoord.markdup.bam.bai"),
    metrics = "results/markdup/{sample}.metrics.txt"
  params:
    runtime = "06:00:00"
  log:
    "logs/markdup/{sample}.log"
  threads: 10
  conda:
    "../envs/rnaseq.yaml"
  shell:
    '''
    picard -Xms1g -Xmx20g MarkDuplicates \
      INPUT={input.bam} \
      OUTPUT={output.bam} \
      METRICS_FILE={output.metrics} \
      MAX_RECORDS_IN_RAM=100000 \
      >{log} 2>&1

    samtools index {output.bam} >>{log} 2>&1
    '''


