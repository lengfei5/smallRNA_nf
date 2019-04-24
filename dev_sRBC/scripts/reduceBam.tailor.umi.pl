#!/usr/bin/perl -w

###########################
# Thomas R Burkard        #
# IMP/IMBA Bioinformatics #
###########################

use List::Util qw(sum);
use strict;

use Getopt::Long;

my $gtf="";
my $sam="";
GetOptions(
    'g=s'       => \$gtf,
    's=s'       => \$sam,
    ) || die <<END;
Usage: $0 [options] -
  Options:
          -g        GTF for annotation
          -s        SAM input
END


my $gname = {};
my $gtype = {};
my $chromosome = {};
open(GTF, "$gtf") || die "No GTF file provided";
while (<GTF>)
{
    chomp;
    
    my @field = split "\t";
    $field[0] =~ s/\s//g;

    $chromosome->{$field[0]} = 1;

    if ($field[-1] =~ /gene_id "(.*?)"; transcript_id ".*?"; gene_name "(.*?)";/)
    {
	$gname -> {$1} = $2;	    
	$gtype -> {$1} = $field[1];
    } elsif ($field[-1] =~ /gene_id "(.*?)";/)
    {
	$gname -> {$1} = $1;	    
	$gtype -> {$1} = $field[1];
    }
}
close(GTF);


print "chr\tstart\tend\tstrand\tfbgn\tgname\ttype\tsequence\trepeat\ttail\ttailLen\tcount\tUMInum\tUMIfr\n";

my $cntFeature = {};
foreach my $chrom (sort keys %$chromosome)
{
    my $seq = {};
    open(SAM, "$sam")  || die "No SAM (" . $sam . ") file provided";
    while (<SAM>)
    {
	chomp;

	my @field = split "\t";
	next if ($field[2] ne $chrom); #chromsome by chromosome to save memory but slower...


	my $umi = $field[0];
	$umi =~ s/.*_UMI://;
	$umi = "N" if ($umi =~ /N/); #collapse all N containing UMIs
	
	my $fbgn = "none";
	if (/XF:Z:mRNA-(\S+)/) {
	    $fbgn = $1;
	} elsif (/XF:Z:(\S+)/) {
	    $fbgn = $1;
	} 

	my $fr = 0;
	if (/NH:i:(\d+)/)
	{
	    $fr = 1/$1;
	} elsif (/FR:f:(\S+)/) {  #fraction see bam.NH2fraction.pl
	    $fr = $1;
	} 

	my $repeat="";
	if ($fr == 1)
	{
	    $repeat="U";
	} else {
	    $repeat="R";
	}

	my $tail = "";
	if (/TL:Z:(\S+)/)
	{
	    $tail = $1;
	}

	$field[3] -= 1; #1based to 0based

	$field[5] =~ s/(\d+I)//g; #remove Inserts which don't contribute to the end mapping position
	$field[5] =~ s/(\d+S)//g; #remove Softclip which don't contribute to the end mapping position
	my @nt = split "[A-Z]", $field[5];
	my $ntSum = sum @nt;
	my $end = $field[3] + $ntSum;

	my $strand = "";
	if ($field[1] & 0x0010) #reverse complement
	{
	    $strand = "-";
	    $field[9] =~ tr/gatcGATC/ctagCTAG/;
	    $field[9] = reverse($field[9]);

	} else {
	    $strand = "+";
	}
    
	$seq -> {"$field[2]"} -> {"$field[3]"} -> {$end} -> {$strand} -> {$fbgn} -> {$repeat} -> {$tail} -> {"$field[9]"}->{"cnt"} += $fr;
	$seq -> {"$field[2]"} -> {"$field[3]"} -> {$end} -> {$strand} -> {$fbgn} -> {$repeat} -> {$tail} -> {"$field[9]"}->{"umi"}->{$umi} = $fr;
    
	$fbgn =~ s/ambiguous.*/ambiguous/;
	$cntFeature -> {"$fbgn"} += $fr;
    }
    close(SAM);

    foreach my $chr (sort keys %$seq)
    {
	foreach my $pos (sort keys %{$seq->{$chr}})
	{
	    foreach my $end (sort keys %{$seq->{$chr}->{$pos}})
	    {
		foreach my $strand (sort keys %{$seq->{$chr}->{$pos}->{$end}})
		{
		    
		    foreach my $fbgn (sort keys %{$seq->{$chr}->{$pos}->{$end}->{$strand}})
		    {
			foreach my $repeat (sort keys %{$seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}})
			{
			    foreach my $tail (sort keys %{$seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}})
			    {
				foreach my $s (sort keys %{$seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}})
				{
				    if ( exists $gname->{$fbgn} )
				    {
					my $umiNum = keys %{$seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}->{$s}->{"umi"}};
					my $UMIfr = 0;
					foreach my $umi (keys %{$seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}->{$s}->{"umi"}})
					{
					    $UMIfr += $seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}->{$s}->{"umi"}->{$umi};
					}
					print "$chr\t$pos\t$end\t$strand\t$fbgn\t" . $gname->{$fbgn} . "\t" . $gtype->{$fbgn} . ";\t$s\t$repeat\t$tail\t" . length($tail) . "\t"  . $seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}->{$s}->{"cnt"} . "\t$umiNum\t$UMIfr\n";
				    } else {
					my $type = "";
					my $gn = "";
					if ($fbgn =~ /ambiguous\[(.*)\]/)
					{
					    my @amb = split '\\+', $1;
					    my $t = {};
					    my $g = {};
					    foreach my $id (@amb)
					    {
						$g->{$gname->{$id}} ++;
						$t->{$gtype->{$id}} ++;
					    }
					
					    foreach my $key (sort keys %$g)
					    {
						$gn .= "$key|";
					    }
					    $gn =~s/\|$//;

					    #$gn =~ join "|", sort keys %$g;
					
					    foreach my $key (sort keys %$t)
					    {
						$type .= "$key;" . $t->{$key} . "|";
					    }
					    $type =~s/\|$//;

					} else {
					    $type=$fbgn;
					    $gn = $fbgn;
					}
					
					my $umiNum = keys %{$seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}->{$s}->{"umi"}};
					my $UMIfr = 0;
					foreach my $umi (keys %{$seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}->{$s}->{"umi"}})
					{
					    $UMIfr += $seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}->{$s}->{"umi"}->{$umi};
					}
					
					print "$chr\t$pos\t$end\t$strand\t$fbgn\t$gn\t$type\t$s\t$repeat\t$tail\t" . length($tail) . "\t"  . $seq->{$chr}->{$pos}->{$end}->{$strand}->{$fbgn}->{$repeat}->{$tail}->{$s}->{"cnt"} . "\t$umiNum\t$UMIfr\n";
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}


foreach my $fbgn (sort keys %$cntFeature)
{
    print STDERR "$fbgn\t" . $cntFeature->{$fbgn} . "\n";
}
