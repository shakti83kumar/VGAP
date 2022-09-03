#!/usr/bin/perl -w
use strict;
##
## raw reads is processed by fastp program
##
## Illumina FASTQ files use the following naming scheme:
## <sample name>_<barcode sequence>_L<lane>_R<read number>_<set number>.fastq.gz
## example: SARSCOV2_S1_L001_R1_001.fastq.gz
##
##
## ************************ user input ************************
##
##
my $inFolder       = '/media/shakti/SHAKTI4/SARSCoV2_based_analysis/20220902_R51'; #folder path of raw reads;
my $rawreadfmt = 'fastq.gz';
my $pear  = '/media/shakti/SHAKTI4/Softwares/pear-0.9.11-linux-x86_64/bin/pear'; #path of fastp program:
my $p = 2; #number of processors; 
my $qual = 30;
my $len = 10; ## Minimum read size after trimming;
my $basewordsnum = 2; #number of starting words to inculde to make base word
my $lane = 'L001'; # L<lane number>
my $read1 = 'R1'; # R<read number>
my $read2 = 'R2'; # R<read number>
my $setnum = '001';  # <set number>
##
## ************************************************************
my $outFolder     = $inFolder.'_01mergedreads'; #folder path of clean raw reads
mkdir($outFolder, 0755);
my $inFolderNamePath = &Strip($inFolder);
my $outFolderNamePath = &Strip($outFolder);
#my (@fileArray, %baseNameHash);
opendir(INDIR, $inFolderNamePath)||die("Could not open the $inFolderNamePath\n");
opendir(INDIR, $inFolderNamePath)||die("Could not open the $inFolderNamePath\n");
while (my $file = readdir INDIR) 
 {
   $file = &Strip($file);
   next if($file =~ /^\.+/);
   if($file =~ /$rawreadfmt$/)
    {
       my @fileNameArray = split(/\_/, $file);
       #my $baseName = join("_", $fileNameArray[0]);
       my $baseName = join("_", @fileNameArray[0..$basewordsnum-1]);
       #$baseNameHash{$baseName} += 1;
       #my $R2 = $inFolderNamePath.'/'.$baseName.'_L001_R1_001.'.$rawreadfmt;
       #my $R2 = $inFolderNamePath.'/'.$baseName.'_L001_R2_001.'.$rawreadfmt;
       my $R1 = $inFolderNamePath.'/'.$baseName;
       my $R2 = $inFolderNamePath.'/'.$baseName;
       if($lane ne '')
       	{
       		$R1 = $R1.'_'.$lane;
       		$R2 = $R2.'_'.$lane;
       	}
       if(($read1 ne '')&&($read2 ne ''))
       	{
       		$R1 = $R1.'_'.$read1;
       		$R2 = $R2.'_'.$read2;
       	}
       if($setnum ne '')
       	{
       		$R1 = $R1.'_'.$setnum;
        		$R2 = $R2.'_'.$setnum;      		
       	}
       $R1 = $R1.'.'.$rawreadfmt;
       $R2 = $R2.'.'.$rawreadfmt;
       my $R12 = $outFolderNamePath.'/'.$baseName; 
       print $R1,"\n";
       print $R2,"\n";
       print $R12,"\n";
       system("$pear -j $p -q $qual -t $len -f $R1 -r $R2 -o $R12");
       #push(@fileArray, $file);
    }
 }
closedir(INDIR);
#undef(@fileArray);
#undef(%baseNameHash);
###........Removing either spaces or newline characters........................
sub Strip
 {
    my $stripVal = shift;
    $stripVal =~ s/^\s+|\s+$//;
    return($stripVal);
 }
