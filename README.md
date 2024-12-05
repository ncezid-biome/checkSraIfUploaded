# Check NCBI if you already uploaded a read set

I have a large number of fastq files locally and most of them have been uploaded to sra. Some might not have been. Can I check against SRA if any particular fastq file has been uploaded already? I feel like it might involve downloading a large spreadsheet of every sra entry for a given pathogen, checking the number of spots against the number of reads I have, and then if it matches, verify with a hashsum.

## Get the list of possible accessions from NCBI

```bash
# Get all metadata from the Listeria NCBI Pathogens
cat ~/GWA/projects/SneakerNet/ncbiPathogens/data/ftp.ncbi.nlm.nih.gov/pathogen/Results/Listeria/latest_snps/Metadata/PDG000000001.3974.metadata.tsv | sort | gzip -c9 > PDG000000001.3974.metadata.tsv.gz
# extract any SRA run ID
zcat PDG000000001.3974.metadata.tsv | perl -lane 'while(/([DES]RR\d+)/g){print $1;}' > SRA_run.txt
zcat PDG000000001.3974.metadata.tsv | perl -lane 'while(/([DES]RS\d+)/g){print $1;}' > SRS_run.txt
# Get the number of spots per SRA ID
nohup efetch -input SRA_run.txt -db SRA > Listeria.SRR.csv 2> Listeria.SRR.csv.log &
nohup efetch -input SRS_run.txt -db SRA > Listeria.SRS.csv 2> Listeria.SRS.csv.log &

# Make a two column file with SRR identifier and number of spots
# Create gz file to help with disk I/O times when searching later
cut -f4,25 -d, Listeria.SRS.csv Listeria.SRR.csv | sort -nr | gzip -9c > spots.csv.gz
```

## Check my number of spots vs NCBI

```bash
zcat testdata/PNUSAL001024_[12].fastq.gz | fasten_metrics | sed 's/^/#/' | column -t
#totalLength  numReads  avgReadLength  avgQual
#305569340    2094942   145.86052      37.695873

# Take the number of reads and divide by 2 to get spots
zgrep $((2094942/2)) spots.csv.gz | sed 's/^/#/'
#1047471,SRS716008
```

Automate it

```bash
# Loop through all R1 files
dir=testdata
for R1 in $dir/*_R1_*.fastq.gz $dir/*_1.fastq.gz; do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    R2=${R2/_R1_/_R2_}

    echo "Looking at $R1 and $R2 ..."

    # get the number of spots in the local fastq 
    numReads=$(zcat $R1 $R2 | fasten_metrics | tail -n 1 | cut -f 2)
    SRS=$(zgrep $(($numReads/2)) spots.csv.gz | cut -f2 -d, | sort -n | uniq)

    # Even though only one SRS is expected, more could be found and so loop through them
    # If the SRS is found, then mark for later verification
    (
        for srs_acc in $SRS; do
            echo -e "$srs_acc\t$R1\t$R2"
        done
    ) | gzip -9c >> spots_found_on_ncbi.txt.gz
    # If the SRS is not found, then it is at least not found in NCBI Pathogens and so
    # mark it for safe keeping
    (
        if [ ! "$SRS" ]; then
            echo $R1 $R2
        fi
    ) >> not_on_ncbi.txt
done

# Some results seem to be duplicated a ton and so let's mostly ignore them
zcat spots_found_on_ncbi.txt.gz | perl -lane '
  next if($F[1]++ > 10); 
  print;
' > spots_found_on_ncbi.txt && \
  gzip -fv9 spots_found_on_ncbi.txt

# TODO download fastq from NCBI and compare against local fastq file
# To recap: not_on_ncbi.txt represents reads not already on NCBI for sure, but
# spots_found_on_ncbi.txt represents reads that _might_ be on NCBI.
# We need to verify the reads from spots_found_on_ncbi.txt and either:
#   * (not verified) add these reads to not_on_ncbi.txt 
#   * or, (verified) add these reads to a new file found_on_ncbi.txt
mashThreshold=0.95
k=32
stackSize=10000
parentTempdir="fasterq-dump"
zcat spots_found_on_ncbi.txt.gz | shuf | xargs -P 4 -L 1 bash -c '
  SRS=$0
  localR1=$1
  localR2=$2
  if [[ $(stat -c%s "$localR1") -lt 1000 || $(stat -c%s "$localR2") -lt 1000 ]]; then
    echo "Error: One or both files are smaller than 1000 bytes." >&2
    ls -lhL $localR1 $localR2 >&2
    exit 1
  fi

  perl scripts/checkLocalVsNcbiWithMash.pl --R1 $localR1 --R2 $localR2 --SRS $SRS
' > checkNcbiWithMash.tsv
```

## Notices

### Public Domain Notice

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

### License Standard Notice

### Privacy Notice

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

### Contributing Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Notice

This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).
