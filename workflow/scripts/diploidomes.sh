#!/bin/bash

#----- Parse Input Arguments -----#

# from https://stackoverflow.com/a/39376824
# usage function
function usage()
{
   cat << HEREDOC

   Usage: $progname [options] [--ref genome.fa.gz] [--gtf annotation.gtf] [--vcf variants.vcf]

   required arguments:
     -r, --ref STR            reference genome in fasta format (can be compressed)
     -a, --gtf STR            annotation GTF file (can NOT be compressed - see notes)
     -v, --vcf STR            variants VCF file (can be compressed)
     -s, --sample STR         sample name from which to extract variants in VCF file 
                              (see 'bcftools consensus' documentation for details)
    
    optional arguments:
     --hap1 STR               name for haplotype 1. Default: HAP1
     --hap2 STR               name for haplotype 2. Default: HAP2
     -o, --outdir STR         output directory. Default: diploidome
     -i, --index              create Salmon index. Default: off
     -c, --cpus NUM           number of CPUs (used for salmon indexing). Default: 1
     -k, --keep-intermediate  keep intermediate haplotype-specific files. Default: off

     -h, --help               print this message

     Notes:
       The GTF file should not be compressed. However, you can pass a compressed 
       file by using the following syntax (where "annotation.gtf.gz" is your file):
         --gtf <(zcat annotation.gtf.gz)

HEREDOC
}  

# initialize default variables
progname=$(basename $0) # for getopt
HAP1="HAP1"
HAP2="HAP2"
OUTDIR="diploidome"
CPUS="1"
KEEP="0"
INDEX="0"

# REF="data/genome.fa.gz"
# GTF="data/annotation.gtf.gz"
# VCF="data/B6xCAST.snps_and_indels.vcf.gz"
# SAMPLE="B6xCAST"

# use getopt and store the output into $OPTS
# note the use of -o for the short options, --long for the long name options
# and a : for any option that takes a parameter
OPTS=$(getopt -o "hr:a:v:s:o:ic:k" --long "help,ref:,gtf:,vcf:,sample:,outdir:,index,cpus:,keep-intermediate,hap1:,hap2:" -n "$progname" -- "$@")
if [ $? != 0 ] ; then echo "Error in command line arguments." >&2 ; usage; exit 1 ; fi
eval set -- "$OPTS"

while true; do
  # uncomment the next line to see how shift is working
  # echo "\$1:\"$1\" \$2:\"$2\""
  case "$1" in
    -h | --help ) usage; exit; ;;
    -r | --ref ) REF="$2"; shift 2 ;;
    -a | --gtf ) GTF="$2"; shift 2 ;;
    -v | --vcf ) VCF="$2"; shift 2 ;;
    -s | --sample ) SAMPLE="$2"; shift 2 ;;
    -o | --outdir ) OUTDIR="$2"; shift 2 ;;
    --hap1 ) HAP1="$2"; shift 2 ;;
    --hap2 ) HAP2="$2"; shift 2 ;;
    -i | --index ) INDEX=1; shift ;;
    -c | --cpus ) CPUS="$2"; shift 2 ;;
    -k | --keep-intermediate ) KEEP=1; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

# if (( $verbose > 0 )); then

#    # print out all the parameters we read in
#    cat <<EOM
#    num=$num_str
#    time=$time_str
#    verbose=$verbose
#    dryrun=$dryrun
# EOM
# fi

#----- Check Inputs -----#

if [ ! -f "$REF" ]; then
  printf "\n ERROR: reference file does not exist $REF\n"
  usage; exit 1
fi

if [ ! -f "$GTF" ]; then
  printf "\n ERROR: annotation file does not exist $GTF"
  usage; exit 1
fi

if [ ! -f "$VCF" ]; then
  printf "\n ERROR: variant file does not exist $VCF\n"
  usage; exit 1
fi

if bcftools view --header-only --samples $SAMPLE $VCF > /dev/null 2>&1 ; then
  echo ""
else
  printf "\n ERROR: sample could not be found in VCF file: $SAMPLE\n"
  usage; exit 1
fi


# print out all the parameters we read in
cat <<EOM 
Running with options:
  --ref $REF
  --gtf $GTF
  --vcf $VCF
  --sample $SAMPLE
  --outdir $OUTDIR
  --hap1 $HAP1
  --hap2 $HAP2
  --index $INDEX
  --cpus $CPUS
  --keep-intermediate $KEEP

EOM


#----- Check Software -----#
# https://stackoverflow.com/a/677212/5023162

command -v bcftools >/dev/null 2>&1 || { echo >&2 "Please install bcftools: 'conda install -c bioconda bcftools'."; exit 1; }
command -v STAR >/dev/null 2>&1 || { echo >&2 "Please install STAR: 'conda install -c bioconda star'."; exit 1; }
command -v gffread >/dev/null 2>&1 || { echo >&2 "Please install gffread: 'conda install -c bioconda gffread'."; exit 1; }
command -v salmon >/dev/null 2>&1 || { echo >&2 "Please install salmon: 'conda install -c bioconda salmon'."; exit 1; }



#----- Prepare Output -----#

# current working directory
WD=$(pwd)

# output directory
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/${HAP1}
mkdir -p ${OUTDIR}/${HAP2}

# prefixes for haplotype-specific files
HAP1PREFIX="${OUTDIR}/${HAP1}/${HAP1}"
HAP2PREFIX="${OUTDIR}/${HAP2}/${HAP2}"



#----- Prepare Genome 1 -----#

echo "Preparing genome for $HAP1..."

# Create HAP1 "genome"
bcftools consensus \
  --chain "${HAP1PREFIX}.chain" \
  --fasta-ref "${REF}" \
  --haplotype 1 \
  --sample "${SAMPLE}" \
  --prefix "${HAP1}." "${VCF}" |\
  gzip > "${HAP1PREFIX}.genome.fa.gz"

echo "  $HAP1 reference genome fasta: DONE"

# liftover variants to create HAP1 GTF
if [ "${GTF: -3}" == ".gz" ]; then
  # if compressed need to zcat it first
  echo "  Decompressing GTF file"
  zcat ${GTF} > ${HAP1PREFIX}.gtf
  STAR --runMode liftOver \
    --genomeChainFiles "${HAP1PREFIX}.chain" \
    --sjdbGTFfile ${HAP1PREFIX}.gtf \
    --outFileNamePrefix "${HAP1PREFIX}."
  rm ${HAP1PREFIX}.gtf
else
  STAR --runMode liftOver \
    --genomeChainFiles "${HAP1PREFIX}.chain" \
    --sjdbGTFfile ${GTF} \
    --outFileNamePrefix "${HAP1PREFIX}."
fi

echo "  $HAP1 GTF liftOver: DONE"

# get transcript fasta sequences for HAP1
gzip -dc ${HAP1PREFIX}.genome.fa.gz | sed "s/^>${HAP1}\./>/" > ${HAP1PREFIX}.fa
cat ${HAP1PREFIX}.GTFliftOver_1.gtf | \
  gffread - -w - -g ${HAP1PREFIX}.fa | \
  sed "s/^>/>${HAP1}./" | \
  gzip > ${HAP1PREFIX}.transcripts.fa.gz
rm ${HAP1PREFIX}.fa # remove temporary decompressed file

echo "  $HAP1 transcript fasta: DONE"


#----- Prepare Genome 2 -----#

echo "Preparing genome for $HAP2..."

# Create HAP2 "genome"
bcftools consensus \
  --chain "${HAP2PREFIX}.chain" \
  --fasta-ref "${REF}" \
  --haplotype 2 \
  --sample "${SAMPLE}" \
  --prefix "${HAP2}." "${VCF}" |\
  gzip > "${HAP2PREFIX}.genome.fa.gz"

echo "  $HAP2 reference genome fasta: DONE"

# liftover variants to create HAP2 GTF
if [ "${GTF: -3}" == ".gz" ]; then
  # if compressed need to zcat it first
  STAR --runMode liftOver \
    --genomeChainFiles "${HAP2PREFIX}.chain" \
    --sjdbGTFfile <(zcat ${GTF}) \
    --outFileNamePrefix "${HAP2PREFIX}."
else
  STAR --runMode liftOver \
    --genomeChainFiles "${HAP2PREFIX}.chain" \
    --sjdbGTFfile ${GTF} \
    --outFileNamePrefix "${HAP2PREFIX}."
fi

echo "  $HAP2 GTF liftOver: DONE"

# get transcript fasta sequences for HAP2
gzip -dc ${HAP2PREFIX}.genome.fa.gz | sed "s/^>${HAP2}\./>/" > ${HAP2PREFIX}.fa
cat ${HAP2PREFIX}.GTFliftOver_1.gtf | \
  gffread - -w - -g ${HAP2PREFIX}.fa | \
  sed "s/^>/>${HAP2}./" | \
  gzip > ${HAP2PREFIX}.transcripts.fa.gz
rm ${HAP2PREFIX}.fa # remove temporary decompressed file

echo "  $HAP2 transcript fasta: DONE"




#----- Concatenate Diploid Genomes -----#

echo "Creating diploid genome/transcriptome..."

# concatenate to create diploid genome
cat "${HAP1PREFIX}.genome.fa.gz" "${HAP2PREFIX}.genome.fa.gz" > ${OUTDIR}/diploid_genome.fa.gz

# create diploid transcript FASTA
cat ${HAP1PREFIX}.transcripts.fa.gz ${HAP2PREFIX}.transcripts.fa.gz > ${OUTDIR}/diploid_transcripts.fa.gz



#----- Salmon Indexing -----#
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

if (( $INDEX > 0 )); then
  echo "Creating Salmon index..."

  cd ${OUTDIR}
  mkdir salmon_index 

  # prepare "gentrome"
  zcat diploid_genome.fa.gz | grep "^>" | cut -d " " -f 1 | sed 's/>//g' > diploid_decoys.txt
  cat diploid_transcripts.fa.gz diploid_genome.fa.gz > diploid_gentrome.fa.gz

  # index
  # --keepDuplicates is used so that identical transcripts (alleles with no SNPs) get quantified separately
  # see: https://combine-lab.github.io/salmon/faq/ (search for --keepDuplicates)
  salmon index \
    --keepDuplicates \
    -t diploid_gentrome.fa.gz \
    -d diploid_decoys.txt \
    -p ${CPUS} \
    -i salmon_index \
    > salmon_index.log

  # remove unnecessary files
  rm diploid_gentrome.fa.gz
  rm diploid_decoys.txt
fi


#----- Clean -----#

if (( $KEEP > 0 )); then
  echo "Intermediate files for each haplotype kept in:"
  echo "${WD}/${OUTDIR}/${HAP1}"
  echo "${WD}/${OUTDIR}/${HAP2}"
else 
  echo "Removing intermediate files for each haplotype from:"
  echo "${WD}/${OUTDIR}/${HAP1}"
  echo "${WD}/${OUTDIR}/${HAP2}"
  rm -r "${WD}/${OUTDIR}/${HAP1}" "${WD}/${OUTDIR}/${HAP2}"
fi