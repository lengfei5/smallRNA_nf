#!/usr/bin/perl

use POSIX qw/floor/;

while(<>)
{
    @cells = split "\t";

    $mid = floor(($cells[4]+$cells[3])/2);

    @out2=@cells;

    if ($cells[6] eq "+")
    {
	$cells[8] =~ s/";/-5p";/g;
	$out2[8] =~ s/";/-3p";/g;
    } else {
	$cells[8] =~ s/";/-3p";/g;
	$out2[8] =~ s/";/-5p";/g;
    }

    $cells[4] = $mid;
    $out2[3] = $mid+1;

    print join "\t", @cells;
    print join "\t", @out2;

}

    
	
