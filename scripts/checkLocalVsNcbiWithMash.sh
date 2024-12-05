#!/bin/bash
set -e
set -o pipefail
set -u

if [ $# -lt 1 ]; then
  echo "Usage: $0 <localR1> <localR2> <SRS>"
  exit 0
fi

mashThreshold=0.95
k=32
stackSize=10000

tempdir=$(mktemp -d -t checkSra.XXXXXX)
trap 'rm -rf $tempdir' EXIT

localR1=$1
localR2=$2
SRS=$3

sample=$(basename $localR1 .fastq.gz)

echo "NOTE: testing $sample ($SRS)" >&2
SRR=$(esearch -db sra -query $SRS | efetch -format xml | xtract -pattern EXPERIMENT_PACKAGE -element RUN@accession)
echo "NOTE: $sample => $SRR ($SRS)" >&2
fasterq-dump $SRR --threads 1 --outdir $tempdir --split-files --skip-technical
ncbiR1=$tempdir/${SRR}_1.fastq
ncbiR2=${ncbiR1/_1/_2}
distR1=$(mash dist -s $stackSize -k $k $ncbiR1 $localR1 | cut -f 3)
distR2=$(mash dist -s $stackSize -k $k $ncbiR2 $localR2 | cut -f 3)

echo "NOTE: $sample($SRR) mash dists: $distR1 $distR2" >&2

if (( $(echo "$distR1 > $mashThreshold" | bc -l) )) && (( $(echo "$distR2 > $mashThreshold" | bc -l) )); then
  echo "FOUND ON NCBI: $SRR $sample" >&2
  exit 0
else
  echo "NOT FOUND ON NCBI: $SRR $sample" >&2
  exit 1
fi
