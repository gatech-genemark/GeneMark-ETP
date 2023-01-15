#!/usr/bin/env perl
# -------------
# AL
# -------------

use strict;
use warnings;

my $gff = shift;
my $trace = shift;

my %h = LoadTracs($trace);

open(my $IN, $gff) or die;
while(<$IN>)
{
	next if /^#/;
	next if /^\s*$/;

        my @arr = split('\t');

	my $id = $arr[0];

	my $new_id = '';
	if ( exists $h{$id}{"sid"} )
	{
		$new_id = $h{$id}{"sid"};
	}
	else
	{
		die "error, not found in trace: $id\n";
	}

	$arr[0] = $h{$id}{"sid"};
        $arr[3] = $arr[3] + $h{$id}{"L"} - 1;
        $arr[4] = $arr[4] + $h{$id}{"L"} - 1;

	my $str = join( "\t", @arr);

	print $str;
}


sub LoadTracs
{
	my $fname = shift;

	my %h = ();

	open(my $IN, $fname) or die;
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		if ( /(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/ )
		{
			my $current = $1;
			my $orignal = $2;
			my $L = $3;
			my $R = $4;

			$h{$current}{"sid"} = $orignal;
			$h{$current}{"L"} = $L;
			$h{$current}{"R"} = $R;
		}
		else {die "error $_";}	
	}
	close $IN;

	return %h;
}

