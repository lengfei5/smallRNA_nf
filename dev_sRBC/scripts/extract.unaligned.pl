#!/usr/bin/perl

###########################
# Thomas R Burkard        #
# IMP/IMBA Bioinformatics #
###########################

use strict;
use Getopt::Long;

my $bam;
my $fastq;

GetOptions("b=s"=>\$bam,
    "f=s"=>\$fastq); #rm strand errors                                                                                                                                                                                               

print STDERR "Reading BAM ...";
my $bid = {};
my $cntUniqChr = {};
my $cntMultiChr = {};

open(BAM, "samtools view $bam |"); 
while (my $l=<BAM>) { 

    my @cell=split "\t", $l;
    $bid->{$cell[0]}++;

    if ($l =~ /NH:i:(\d+)/)
    {
	if ($1 == 1)
	{
	    $cntUniqChr->{$cell[2]}++;
	} else {
	    $cntMultiChr->{$cell[2]}+=1/$1;
	}
    }
}
close(BAM);
print STDERR "done\n";

my $cntUniq = 0;
my $cntMulti = 0;

foreach my $id (keys %$bid)
{
    if ($bid->{$id} == 1)
    {
	$cntUniq++;
    } else {
	$cntMulti++;
    }
}


print STDERR "Reading FASTQ ...";
open(FASTQ, "$fastq");
my $cntUnalign = 0;
while (my $line = <FASTQ>)
{
    my $id = $line;
    chomp $id;
    $id=~s/^@//;
    $id=~s/^(\S+).*/$1/; #SRR repository

    if (! (exists $bid->{$id})) #unaligned
    {
	print $line;
	my $tmp = <FASTQ>;
	print $tmp;
	$tmp = <FASTQ>;
	print $tmp;
	$tmp = <FASTQ>;
	print $tmp;


	$cntUnalign++;
    } else {
	my $tmp = <FASTQ>;
	$tmp = <FASTQ>;
	$tmp = <FASTQ>;

    }
}
close(FASTQ);
print STDERR "done\n";

open (OUT, ">$bam.virusStat.txt");
print OUT "Total count (uniq):\t$cntUniq\n";
print OUT "Total count (multimap):\t$cntMulti\n";
print OUT "Total count (unaligned):\t$cntUnalign\n";
close OUT;
