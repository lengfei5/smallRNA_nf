#!/usr/bin/perl

use strict;
use Getopt::Long;
Getopt::Long::Configure ("bundling"); #allow case sensitive!

my $minLength=0;
my $maxLength=-1;
my $trim5=0; #random 5' length
my $trim3=0; #random 3' length

#min/max sequence length
GetOptions('m:f' => \$minLength,
	   'M:f' => \$maxLength,
           '5:f' => \$trim5,
           '3:f' => \$trim3);

print STDERR "Min. length: $minLength\n";
print STDERR "Max. length: $maxLength\n";
print STDERR "Random 5' length: $trim5\n";
print STDERR "Random 3' length: $trim3\n";

# trim random 4mer at 5' and 3'
sub trim53 {
    my ($fastq) = @_;
    my @fqEntry = split "\t", $fastq;

    my $fqTrimmed = "";
    if (length($fqEntry[1]) >= ($trim5+$trim3)) #else shorter than random regions
    {
	#extract UMI & sequence
	my $umi=$fqEntry[1];
	if ($trim3 == 0)
        {
            $umi=substr($umi, 0, $trim5);
            $fqEntry[1]=substr($fqEntry[1], $trim5); #sequence
            $fqEntry[3]=substr($fqEntry[3], $trim5); #quality
        } else {
            $umi=substr($umi, 0, $trim5) . substr($umi, -$trim3);
            $fqEntry[1]=substr($fqEntry[1], $trim5, -$trim3); #sequence
            $fqEntry[3]=substr($fqEntry[3], $trim5, -$trim3); #quality
        }
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
