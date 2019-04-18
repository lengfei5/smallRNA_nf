#!/usr/bin/perl

###########################
# Thomas R Burkard        #
# IMP/IMBA Bioinformatics #
###########################

#replace NH:i: by FR:f: (fraction). For use with tailor alignments.
#Attention multihits without NH:i: will get counted several times with HTSeq-count!!!

use strict;
use Getopt::Long;

my $maxNH=-1;

GetOptions('m=s'=>\$maxNH) || die <<END;
  Usage: $0 [options] -
  Options:
          -m        Maximum number of multihits
END

my $cnt = 0;
while (my $l = <>)
{
    $cnt++;
    if ($l=~/NH:i:(\d+)/)
    {
	my $f=1/$1;
	next if (($maxNH > 0) && ($1 > $maxNH));
	$l=~s/NH:i:\d+/FR:f:$f/;
#	$total += $f;
    }
    print $l;
}

#print STDERR "$total\n";
