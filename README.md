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
    R2=${R1/_1/_2}
    R2=${R2/_R1_/_R2_}

    echo "Looking at $R1 and $R2 ..."

    # get the number of spots in the local fastq 
    numReads=$(zcat $R1 $R2 | fasten_metrics | tail -n 1 | cut -f 2)
    SRS=$(zgrep $(($numReads/2)) spots.csv.gz | cut -f2 -d,)

    # Even though only one SRS is expected, more could be found and so loop through them
    # If the SRS is found, then mark for later verification
    for srs_acc in $SRS; do
        echo $R1 $R2 >> spots_found_on_ncbi.txt
    done
    # If the SRS is not found, then it is at least not found in NCBI Pathogens and so
    # mark it for safe keeping
    if [ ! "$SRS" ]; then
        echo $R1 $R2 >> not_on_ncbi.txt
    fi
done

# TODO download fastq from NCBI and compare against local fastq file
```
