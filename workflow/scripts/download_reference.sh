#!/bin/bash
#SBATCH -A LEYSER-SL2-CPU
#SBATCH -J prepare_reference
#SBATCH -D /rds/user/hm533/hpc-work/mammary_gland_transcriptomes/test_pipeline/  # your working directory
#SBATCH -o logs/slurm/00-prepare_reference.log
#SBATCH -p skylake             # or `skylake-himem`
#SBATCH -c 1                   # max 32 CPUs; default 1 CPU
#SBATCH --mem-per-cpu=5980MB   # max 5980MB (or 12030MB for skilake-himem)
#SBATCH -t 03:00:00            # HH:MM:SS with maximum 12:00:00 for SL1 or 36:00:00 for SL2; default 00:10::00

#### Activate Software Environment ####

# conda create -n rnaseq -c bioconda gffread=0.12.1 cutadapt=3.4 salmon=1.4 fastqc=0.11.9 multiqc=1.10.1  sortmerna=4.3.2
source activate rnaseq


#### Input Variables ####

CPUS="$SLURM_CPUS_PER_TASK"    # number of CPUs to use
OUTDIR="data/external/reference/"  # output directory name


#### prepare directories ####

# create directory
mkdir -p "$OUTDIR"

# work from this directory
cd "$OUTDIR"


#### Download reference genome ####

# Using "primary assembly" files (advised in STAR manual)
# Using release M25, which is the last GRCm38 release

# Genome 
wget -O genome.fa.gz http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz

# GFF3/GTF
wget -O annotation.gtf.gz ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
wget -O annotation.gff3.gz ftp://ftp.ensembl.org/pub/release-102/gff3/mus_musculus/Mus_musculus.GRCm38.102.gff3.gz

# Extract transcript sequences from genome
# the reason we extract them ourselves is that ENSEMBL includes transcripts in assembly patches and excludes TECs (to be experimentally confirmed)
gzip -dc genome.fa.gz > genome.fa # gffread needs the files to be decompressed
zcat annotation.gtf.gz | gffread - -w - -g genome.fa | gzip > transcripts.fa.gz
rm genome.fa # remove decompressed file

# create a table with transcript-to-gene IDs
printf "transcript\tgene\n" > transcript2gene.tsv
zcat annotation.gtf.gz | \
  awk -F "\t" '$3 == "transcript" {print $9}' | \
  awk -F ";" '{print $3,$1}' | \
  sed 's/"//g' | sed 's/ *transcript_id *//' | sed 's/ *gene_id */\t/' \
  >> transcript2gene.tsv

# create BED file (useful for downstream analysis)
zcat annotation.gtf.gz | \
  awk 'OFS="\t" {if ($3 == "exon" || $3 == "five_prime_utr" || $3 == "three_prime_utr") { print $1,$4-1,$5,$14}}' | \
  tr -d '";' | \
  sort -k1,1 -k2n,2n | \
  gzip > transcripts.bed.gz


#### salmon indexing ####

# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
mkdir salmon_index 

# prepare "gentrome"
zcat genome.fa.gz | grep "^>" | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
cat transcripts.fa.gz genome.fa.gz > gentrome.fa.gz

# index
salmon index -t gentrome.fa.gz -d decoys.txt -p ${CPUS} -i salmon_index 

# remove unnecessary files
rm gentrome.fa.gz
rm decoys.txt


#### STAR indexing ####

mkdir -p star_index

# STAR wants these uncompressed - could I use process substitution? `<(zcat genome.fa.gz)`
zcat genome.fa.gz > genome.fa 
zcat annotation.gtf.gz > annotation.gtf

STAR \
  --runMode genomeGenerate \
  --runThreadN ${CPUS} \
  --genomeDir star_index \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile annotation.gtf \
  --sjdbOverhang 97

rm genome.fa; rm annotation.gtf


