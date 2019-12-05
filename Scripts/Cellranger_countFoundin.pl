#!/usr/bin/perl -w


##=============================================================================================
## USAGE: either run it using nohup or generate the commands to run later one by on or modify it according to your job submission system
## perl Cellranger_countFoundin.pl --sampleFile Sample_test.txt --refDirName /media/root/dataE/FOUNDIN_vikas/Database/GRCh38_foundinGTF  --resultDir /media/root/dataE/FOUNDIN_vikas/Counts_cellRanger/ --localMem 50
## perl Cellranger_countFoundin.pl --sampleFile ../FASTQ_all/Sample_names.txt --refDirName /media/root/dataE/FOUNDIN_vikas/Database/GRCh38_foundinGTF  --resultDir /media/root/dataE/FOUNDIN_vikas/Counts_cellRanger/ --localMem 50 > Commands_run27.sh
## 
##
## sampleFile should be the tab sep 3 column file. First column contains the sample ID (not the full fastq file name but the first part) example D17-8753_S1_L001_I1_001.fastq.gz should be only written D17-8753 (notice part before first underscore) and second column is the directory with fastq files for that sample and the third column is the ID for output.
## refDirName is the directory contaning genome reference build by cellranger mkref. Example command
##
## cellranger mkref --genome=GRCh38_foundinGTF --fasta=./fasta/genome.fa --genes=./GTF_foundin/gencode_v29.lncipedia_v5_2_hc.annotation.gtf
##
##
## Author:   Vikas Bansal
## Date:     20.11.2019
##============================================================================================
#





use warnings;
use strict;

use Getopt::Long;


my $sample_file;
my $result_dir;
my $ref_dir;
my $local_mem;

GetOptions(
    "sampleFile=s"   => \$sample_file,
    "resultDir=s"    => \$result_dir,
    "refDirName=s"    => \$ref_dir,
    "localMem=s"    => \$local_mem,	
);

open FILE, "$sample_file";

while(<FILE>){
	if(!/^\#/){
	   chomp;
	   (my $sampleID, my $sample_folder, my $out_file)=split();
	    

		 #`cd $result_dir && nohup /media/root/dataE/FOUNDIN_vikas/Tools/cellranger-3.1.0/cellranger count --id=$out_file --transcriptome=$ref_dir --fastqs=$sample_folder --sample=$sampleID --localmem=$local_mem > $out_file.err &`;
		  print "cd $result_dir && /media/root/dataE/FOUNDIN_vikas/Tools/cellranger-3.1.0/cellranger count --id=$out_file --transcriptome=$ref_dir --fastqs=$sample_folder --sample=$sampleID --localmem=$local_mem > $out_file.err\n" ;

	} 
}

