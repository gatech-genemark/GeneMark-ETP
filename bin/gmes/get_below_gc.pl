#!/usr/bin/env perl

use strict;
use warnings;

my %h;

while(my $line = <>)
{
	$h{'A'} = 0;
	$h{'C'} = 0;
	$h{'G'} = 0;
	$h{'T'} = 0;
	my $size = length($line);

	for( my $i = 2; $i < $size; $i += 3)
	{
		$h{ substr( $line, $i, 1) } += 1;
	}

	my $count = $h{'C'}+$h{'G'} + $h{'A'}+$h{'T'};

	if ( $count == 0 )
	{
#		print $line;
		next;
	}

	my $gc = ($h{'C'}+$h{'G'})/($h{'C'}+$h{'G'} + $h{'A'}+$h{'T'});

	if ( $gc < 0.65 )
	{
		print $line;
#		print  $gc ."\n";
	}
	else
	{
#		print  $gc ."\n";
	}
}

