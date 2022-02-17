#!/bin/bash

# Mapping (version 1.1)
## Tom van Schaik, 180110

# - Functionality:
# A simple mapping module using bwa mem, including sorting, marking of duplicated
# reads and indexing. 

# - Input
# * Required input:
#  * Read input file (.fastq(.gz), or stdin)
#  * 
#  * Output directory
# * Optional output:
#  * Paired-end file (default False)
#  * Basename for all output files (default basename input file)
#  * Cores for bowtie (default 8)
#  * Quiet (default False)
#  * Mark duplicates (default False)
#  * Local mapping (default False)

# - Output
# * Standard output:
#  * Bam file
#  * Mapping statistics file (without any real format..)
#  * Bam index file

# - Method
# bwa mem mapping, followed by samtools-based conversions.

# - Version
# 1.0 - First version
# 1.1 - BWA mem mapping, updated samtools version (and script)

# - To do
# * Create a structured statistics file, instead of this huge list.
# * Create some basic statistics figures.

version="1.1"

# External tools:
BWA="bwa"
SAMTOOLS="samtools"         # Note: samtools version 1.6 instead of 0.1something
BAMUTILS="bam"

###############################################################################
## Part 1 - read input arguments ##############################################
###############################################################################

# the usage and getopts section is based upon:
# http://rsalveti.wordpress.com/2007/04/03/bash-parsing-arguments-with-getopts/
usage()
{
  cat << EOF
  usage: $0 options

  This script takes a read file(s) and a genome index, returning an alignment file 
  together with a statistics file. 

  OPTIONS:
  -h      Show this message
  -r      The reads input file, can be fasta or fastq (or - for stdin) (required)
  -i      Index bwa base (required)
  -o      Output directory (required)
  -p      Paired-end read file (default: None)
  -b      Basename (default: Read file basename)
  -c      Cores to be used (default: 8)
  -q      Quiet (default: False)
  -d      Mark duplicates using bamutils (not remove) (default: False)
  -f      Log file (instead of stderr) (default: False)
EOF
}

#  -l      Local mapping (default: False)


CORES=8
while getopts “hr:i:o:p:b:c:qdlf:” OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    r)
      READS_FNAME="$OPTARG"
      ;;
    i)
      INDEX_BASE="$OPTARG"
      ;;
    o)
      OUTPUT_DIR="$OPTARG"
      ;;
    p)
      PAIRED_FNAME="$OPTARG"
      ;;
    b)
      BASENAME="$OPTARG"
      ;;
    c)
      CORES="$OPTARG"
      ;;
    q)
      QUIET=1
      ;;
    d)
      DEDUP=1
      ;;
    # l)
    #   LOCAL=1
    #   ;;
    f)
      LOG_FILE="$OPTARG"
      ;;
    ?)
      usage
      exit
      ;;
  esac
done

if [[ -z $READS_FNAME ]] || [[ -z $INDEX_BASE ]] || [[ -z $OUTPUT_DIR ]]; then
  usage
  exit 1
fi

if [[ $READS_FNAME == "-" ]]; then
  if [[ -z $BASENAME ]]; then
    echo "With stdin, a basename is required."
    usage
    exit 1
  fi
fi

if [[ ! -a $INDEX_BASE.pac ]]; then
  echo "Cannot find the index."
  usage
  exit 1
fi

###############################################################################
## Part 2 - Common functions ##################################################
###############################################################################

# Stop script upon errors
set -e 

# Logging (to stderr)
log () {
  if [[ -z $QUIET ]]; then
    if [[ -z $LOG_FILE ]]; then
      (>&2 echo -e $*)
    else
      echo -e $* >> $LOG_FILE
    fi
  fi
}

# Check exit code
CheckExit() {
  if [ $1 -ne 0 ]; then
    echo "$2" 1>&2
    exit 1
  fi
}

###############################################################################
## Part 3 - Set-up input / output #############################################
###############################################################################

if [[ -z $BASENAME ]]; then
  BASENAME=$(basename $READS_FNAME | cut -d"." -f1)
fi

OUT_BASE=$OUTPUT_DIR/$BASENAME
OUT=$OUTPUT_DIR/$BASENAME.bam
OUT_STAT=$OUTPUT_DIR/$BASENAME.mapping.statistics.txt
OUT_ERR=$OUTPUT_DIR/$BASENAME.mapping.stderr.txt

# Create bowtie call
BWA_CALL="$BWA mem -v 2 -t $CORES $INDEX_BASE"

# if [[ -n $LOCAL ]]; then
#   BWA_CALL+=" --very-sensitive-local"
# else
#   BWA_CALL+=" --very-sensitive"
# fi

if [[ -z $PAIRED_FNAME ]]; then
  BWA_CALL+=" $READS_FNAME"
else
  BWA_CALL+=" $READS_FNAME $PAIRED_FNAME"
fi

# Create dam dedup call if required
if [[ -n $DEDUP ]]; then
  DEDUP_TMP="$OUTPUT_DIR/$BASENAME.tmp.bam"
  DEDUP_CALL="$BAMUTILS dedup --in $OUT --log - --out $DEDUP_TMP"
fi

BWA_VERSION=`$BWA 2>&1 | awk '/Version/{ print $2 }'`
SAMTOOLS_VERSION=`$SAMTOOLS 2>&1 | awk '/Version/{ print $2 }'`
BAM_VERSION=`$BAMUTILS 2>&1 | awk '/Version/{ print $2 }' | sed "s/;//"`
  
###############################################################################
## Part 4 - Execute commands ##################################################
###############################################################################

log "----- mapping: version $version -----"

mkdir -p $OUTPUT_DIR

log "Using bwa mem version $BWA_VERSION"
log "Using samtools version $SAMTOOLS_VERSION"
log "Call: `echo "( $BWA_CALL | $SAMTOOLS view -Sb - | $SAMTOOLS sort -T /tmp/$BASENAME - -o $OUT_BASE.bam )"`\n"

if [[ -z $QUIET ]]; then
  if [[ -z $LOG_FILE ]]; then
    ( $BWA_CALL | $SAMTOOLS view -Sb - | $SAMTOOLS sort -T /tmp/$BASENAME - -o $OUT_BASE.bam )
  else
    ( $BWA_CALL | $SAMTOOLS view -Sb - | $SAMTOOLS sort -T /tmp/$BASENAME - -o $OUT_BASE.bam) 2>> $LOG_FILE
  fi
else
  ( $BWA_CALL | $SAMTOOLS view -Sb - | $SAMTOOLS sort -T /tmp/$BASENAME - -o $OUT_BASE.bam ) 2>> /dev/null
fi

# Put stderr to file if log file is given
#( $BOWTIE_CALL | $SAMTOOLS view -Sb - | $SAMTOOLS sort - $OUT_BASE )
CheckExit $? "BWA mapping failed"


log "------------------------------------------------\n"

if [[ -n $DEDUP ]]; then
  log "Using bamutils dedup version $BAM_VERSION"
  log "Call: `echo "$DEDUP_CALL"`"
  
  $DEDUP_CALL 2>> /dev/null
  CheckExit $? "Dedup failed."
  
  log "Call: mv $DEDUP_TMP $OUT"
  mv $DEDUP_TMP $OUT
  
  log "------------------------------------------------\n"
fi

log "Call: $SAMTOOLS index $OUT"
$SAMTOOLS index $OUT
log "Call: $SAMTOOLS flagstat $OUT > $OUT_STAT"
$SAMTOOLS flagstat $OUT > $OUT_STAT

echo "## Versions ##" >> $OUT_STAT
echo "# mapping.sh: $version" >> $OUT_STAT
echo "# bwa: $BWA_VERSION" >> $OUT_STAT
echo "# samtools: $SAMTOOLS_VERSION" >> $OUT_STAT
echo "# bam dedup: $BAM_VERSION" >> $OUT_STAT

log "------------------------------------------------\n"

log "Done.\n"
