#!/usr/bin/perl

use strict;
use Getopt::Long;
Getopt::Long::Configure ("bundling"); #allow case sensitive!

my $minLength=0;
my $maxLength=-1;

#min/max sequence length
GetOptions('m:f' => \$minLength,
	   'M:f' => \$maxLength);

print STDERR "Min. length: $minLength\n";
print STDERR "Max. length: $maxLength\n";

# trim random 4mer at 5' and 3'
sub trim53 {
    my ($fastq) = @_;
    my @fqEntry = split "\t", $fastq;

    my $fqTrimmed = "";
    if (length($fqEntry[1]) >= 8) #else shorter than random 2x4-mers
    {
	#extract UMI & sequence
	my $umi=$fqEntry[1]; 
	$umi=~s/(....).*(....)/\1\2/; 
	$fqEntry[1]=~s/....(.*)..../\1/; #sequence
	$fqEntry[3]=~s/....(.*)..../\1/; #quality
	$fqEntry[0].="_UMI:$umi"; #id
	$fqTrimmed=join "\n", @fqEntry; 
	$fqTrimmed.="\n";
    } 

    return($fqTrimmed);
}

# check length constrains
sub checkLength {

    my ($fqTrimmed) = @_;
    my @fqEntry = split "\n", $fqTrimmed;
    
    my $seq = $fqEntry[1];
    if (length($seq) >= $minLength)
    {
	if ($maxLength == -1 || length($seq) <= $maxLength)
	{
	    print "$fqTrimmed";
	}
    }   
}
    
#main
while (<>)
{
    chomp;
    my $fqTrimmed = trim53($_);
    checkLength($fqTrimmed) if ($fqTrimmed ne "");
}
