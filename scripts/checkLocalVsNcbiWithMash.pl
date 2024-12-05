#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use Getopt::Long qw/ GetOptions /;

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help SRS=s R1=s R2=s)) or die $!;
  die usage() if($$settings{help});
  
  $$settings{mashThreshold}=0.05;
  $$settings{k}=32;
  $$settings{stackSize}=10000;

  $$settings{tempdir} //= tempdir("checkSra.XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $was_found = checkSra($$settings{SRS}, $$settings{R1}, $$settings{R2}, $settings);
  print join("\t",$was_found, $$settings{SRS}, $$settings{R1}, $$settings{R2})."\n";

  return 0;
}

# Return 1 if found on SRA, 0 if not found on SRA
sub checkSra{
  my($SRS, $localR1, $localR2, $settings)=@_;

  my $xtractElement = 'RUN@accession'; # avoid interpretation of @ in the command below
  my $SRR = `esearch -db sra -query $SRS | efetch -format xml | xtract -pattern EXPERIMENT_PACKAGE -element $xtractElement`;
  chomp($SRR);
  
  system("fasterq-dump $SRR --threads 1 --outdir $$settings{tempdir} --split-files --skip-technical");
  if($?){
    die "ERROR with fasterq-dump";
  }
  system("gzip -1 $$settings{tempdir}/${SRR}_[12].fastq");
  if($?){
    die "ERROR with gzip";
  }
  my $ncbiR1="$$settings{tempdir}/${SRR}_1.fastq.gz";
  my $ncbiR2="$$settings{tempdir}/${SRR}_2.fastq.gz";

  my $mashR1 = mashDist($ncbiR1, $localR1, $settings); 
  my $mashR2 = mashDist($ncbiR2, $localR2, $settings);
  my $distR1 = $$mashR1{dist};
  my $distR2 = $$mashR2{dist};

  if($distR1 < $$settings{mashThreshold} && $distR2 < $$settings{mashThreshold}){
    logmsg "FOUND ON NCBI: $SRR $SRS" if($$settings{verbose});
    return 1;
  } else {
    logmsg "NOT FOUND ON NCBI: $SRR $SRS" if($$settings{verbose});
    return 0;
  }
}

sub mashDist{
  my($file1, $file2, $settings)=@_;

  # Set a default return hash
  my $mashHash = {
    dist       =>  1,
    pval       =>  1,
    fracHashes => "-1/1",
  };

  for my $path($file1, $file2){
    # get the file size of the path
    my $size=-s $path;
    if($size < 1000){
      logmsg "WARNING: $path is less than 1kb. Skipping mash for this file.";
      return $mashHash;
    }
    # if file magic shows that it isn't a fastq.gz file, skip it.
    my $magic=`file $path`;
    if($magic !~ /gzip compressed data/){
      logmsg "WARNING: $path is not a gzipped file. Skipping mash for this file.";
      return $mashHash;
    }
    my $f=basename($path);
    system("mash sketch -s $$settings{stackSize} -k $$settings{k} -o $$settings{tempdir}/$f.msh $path");
    die "ERROR with mash sketch $path" if $?;
  }
  my $output=`mash dist $$settings{tempdir}/*.msh`;
  die "ERROR with mash dist $file1 / $file2" if $?;
  chomp($output);
  my($tempPath1,$tempPath2,$dist,$pval,$fracHashes)=split(/\t/,$output);
  #logmsg join("\t",$file1,$file2,$dist,$pval,$fracHashes);
  if(!defined($dist)){
    logmsg "NOTE: dist not defined for $file1 $file2: setting to 1";
    $dist=1;
  }
  $$mashHash{dist} = $dist;
  $$mashHash{pval} = $pval;
  $$mashHash{fracHashes} = $fracHashes;
  return $mashHash;
}