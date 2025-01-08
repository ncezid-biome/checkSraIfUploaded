#!/bin/bash

set -e
set -o pipefail
set -u

if [ $# -ne 1 ]; then
    echo "Usage: $0 list_of_genomes.tsv.gz"
    echo "  where list of genomes is three columns: SRS_acc, R1.fastq.gz, R2.fastq.gz"
    exit 1
fi

# Check for dependencies
which efetch esearch fastq-dump prefetch 

mkdir -pv checkSraUploaded
zcat $1 | perl -lane '$line=$_; for(1..800){$line.="\t".<>; chomp($line); } print $line' > checkSraUploaded/list_of_genomes.tsv
num_jobs=$(cat checkSraUploaded/list_of_genomes.tsv | wc -l)
qsub -N checkSraUploaded -q edlb.q -e checkSraUploaded -o checkSraUploaded -V -l h_vmem=4G -l h_rt=1:00:00 -cwd -t 1-$num_jobs -tc 8 \
  -v "num_jobs=$num_jobs"  -v "CTRL_FILE=checkSraUploaded/list_of_genomes.tsv"  <<- "END_OF_SCRIPT"
    # This is a "here document."  It gets submitted as though it were a 
    # separate file. The here document ends right before END_OF_SCRIPT
    set -e
    set -o pipefail
    set -u

    line=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
    echo $line | xargs -P 1 -n 3 bash -c '
      SRS=$0
      R1=$1
      R2=$2
      if [ ! -f $R1 ]; then
          echo "Missing $R1"
          exit 0
      fi
      if [ ! -f $R2 ]; then
          echo "Missing $R2"
          exit 0
      fi
      perl ../scripts/checkLocalVsNcbiWithMash.pl --R1 $R1 --R2 $R2 --SRS $SRS
    ' 
END_OF_SCRIPT