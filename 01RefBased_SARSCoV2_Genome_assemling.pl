#!/usr/bin/perl -w
use strict;


## ************************ user input ************************
##
my $dateNrun       = '20220902R51'; #Date and Run of sequencing
my $inFolder       = '/media/shakti/SHAKTI4/SARSCoV2_based_analysis/20220902_R51_01mergedreads'; #folder path of raw reads
my $outfolder      = '/media/shakti/SHAKTI4/SARSCoV2_based_analysis/20220902_R51_02mergedassembly'; #folder path of raw reads
my $RefSeqidxHisat = '/media/shakti/SHAKTI4/SARSCoV2_based_analysis/SARSCoV2ref/NC_045512.2'; #path of reference sequence indexed by hisat2 
my $RefSeqidxSamtl = '/media/shakti/SHAKTI4/SARSCoV2_based_analysis/SARSCoV2ref/NC_045512.2.fa'; #path of reference sequence indexed by samtools
my $searchfileExt = 'fastq';
my $AgvReadlength  = 70; # find from FASTQC program
##
## *********************** user input *************************
mkdir($outfolder, 0755);
&AssemblingRefBased($RefSeqidxHisat, $RefSeqidxSamtl, $inFolder);
sub AssemblingRefBased
	{
		my($RefSeqidxH, $RefSeqS, $in) = @_;
		my %totalSamples;
		opendir(DIR, $in) or die("could not open the folder $in");
		while(my $file = readdir(DIR))
			{	
				if($file =~ m/$searchfileExt$/)
					{
						print $file,"\n";
						my @filaA = split(/\./, $file);
						my $baseName =  $filaA[0];
						$totalSamples{$baseName} += 1;
					}
			}
		foreach my $eachSample (keys %totalSamples)
			{
				print $eachSample,"\n";
				my $assemblyfolderName = $outfolder.'/'.$eachSample;
				mkdir($assemblyfolderName, 0755);
				my $assembledfqfile = $in.'/'.$eachSample.'.assembled.'.$searchfileExt;
				my $unassembledfwdfqfile = $in.'/'.$eachSample.'.unassembled.forward.'.$searchfileExt;
				my $unassembledrwdfqfile = $in.'/'.$eachSample.'.unassembled.reverse.'.$searchfileExt;
				my $hisat2Sam = $assemblyfolderName.'/'.$eachSample.'.sam';
				my $sam2bam   = $assemblyfolderName.'/'.$eachSample.'.bam';
				my $mapq20    = $assemblyfolderName.'/'.$eachSample.'_mapped_q20.bam';
				my $mapidx    = $assemblyfolderName.'/'.$eachSample.'_mapped_q20.bai';
				my $mapfq     = $assemblyfolderName.'/'.$eachSample.'_mapped_q20.fq';
				my $idxstats  = $assemblyfolderName.'/'.$eachSample.'_idxstats';
				my $mapconsfq = $assemblyfolderName.'/'.$eachSample.'_consensus.fq';
				my $genomefas = $assemblyfolderName.'/'.$eachSample.'_consensus.fa';
				system("hisat2 -x $RefSeqidxH -q -U $assembledfqfile -U $unassembledfwdfqfile -U $unassembledrwdfqfile -S $hisat2Sam");
				system("samtools sort $hisat2Sam -o $sam2bam -O BAM");
				system("samtools view -F 0x4 -q 20 -b $sam2bam >$mapq20");
				system("samtools index $mapq20 $mapidx");
				system("samtools idxstats $mapq20  >$idxstats");
				system("bamToFastq -i $mapq20 -fq $mapfq");
				#system("bcftools mpileup -f $refseqS $mapq20 -o -b $mpilupbcf");
				system("bcftools mpileup -f $RefSeqS $mapq20 | bcftools call -c | vcfutils.pl vcf2fq >$mapconsfq");
				system("seqtk seq -A -U $mapconsfq >$genomefas"); 
				my $coverge =  &CalculateCoverage($idxstats);
				print $coverge, "\n";
				&RenamingAssembledSeq($genomefas, $assemblyfolderName, $eachSample, $coverge);
			}
	}
### ........................... Genome Coverage calculation ........................
sub CalculateCoverage
	{
		my $idxfile = shift;
		print $idxfile, "\n";
		my $cov = 0;
		my $check = 0;
		my $NoMapRead;
		my $RefSeqlen;
		open(IDX, $idxfile) or die "Could not open $idxfile $!";
		while(my $ln = <IDX>)
			{
				if( $check == 0)
					{
						$ln =~ s/^\s+|\s+$//;
						my @lnA = split(/\s+/, $ln);
						$RefSeqlen = $lnA[1];
						$NoMapRead = $lnA[2];
					}
				$check++;
			}
		$cov = ($NoMapRead*$AgvReadlength)/$RefSeqlen;
		$cov = int($cov);
		return($cov);
		close(IDX);
				
	}

### 
### ..........................Renaming the Assembled sequnece ....................
sub RenamingAssembledSeq
	{
		my($fasG, $folderPath, $sampleName, $coverge) = @_;
		#my $outputfile = $folderPath.'/'.$sampleName.'/'.$dateNrun.'_'.$sampleName.'_'.$coverge.'x'.'_final.fa';
		my $outputfile = $folderPath.'/'.$dateNrun.'_'.$sampleName.'_'.$coverge.'x'.'_final.fa';
		#my $outputfile = $folderPath.$sampleName.'/'.$sampleName.'_'.$coverge.'x'.'_final.fa';
		open(GENOME, $fasG) or die "Could not open $fasG $!";
		open(OUT, ">$outputfile") or die "Could not open $outputfile $!";
		while(my $Gln = <GENOME>)
			{
				if($Gln =~ /^>/)
					{
						print OUT '>'.$dateNrun.'_'.$sampleName.'_'.$coverge.'x'.'_final.fa',"\n";
						#print OUT '>'.$sampleName.'_'.$coverge.'x'.'_final.fa',"\n";
					}
				else
					{
						print OUT $Gln;
					}
			}
		close(GENOME);
		close(OUT);
	}

