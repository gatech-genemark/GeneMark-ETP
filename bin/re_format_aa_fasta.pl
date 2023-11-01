#!/usr/bin/env perl
# -------------
# Alex Lomsadze
# 2023
# 
# Re-format FASTA file with protein sequences
# The script takes the first word in the defline as a unique identifier.
# Sequence is uppercased.
# Only "A-Z", "*", and "-" letters are allowed in sequence.
# There is no strict checking for the allowed AA alphabet
# -------------

use strict;
use warnings;

my $in_name = shift;
my $out_name = shift;

if (( ! defined $in_name )or( ! defined $out_name))
{
	print "Usage: [in file] [out file]\n";
	exit;
}

die "error, input and output file names are identical\n" if ( $in_name eq $out_name );

my %h = ();
my @arr = ();
my $id = '';

open(my $IN, $in_name) or die "error on open file $in_name: $!\n";
while(<$IN>)
{
	next if /^\s*$/;

	if (/^>\s*(\S+)/)
	{
		$id = $1;

		if ( exists $h{$id} )
		{
			die "error, duplicated protein ID $id was found in the FASTA file $in_name\n";
		}
		else
		{
			$h{$id} = '';
			push @arr, $id;
		}
	}
	else
	{
		chomp;
		my $seq = uc $_;
		$seq =~ s/\s//g;

		if ( $seq =~ /^[A-Z\*\-]+$/ )
		{
			$h{$id} .= $seq;
		}
		else
		{
			die "error, an unexpected symbol was found in the file $in_name in: $seq\n";
		}
	}
}
close $IN;

foreach my $key (@arr)
{
	$h{$key} =~ s/\*$//;

	if ( $h{$key} =~ /\*/ )
	{
		print "warning, stop codon found in protein sequence, record $key\n";
	}
}

open(my $OUT, ">", $out_name) or die "error on open file $out_name\n";
foreach my $key (@arr)
{
	print $OUT (">". $key ."\n");
	print $OUT ($h{$key} ."\n");
}
close $OUT;

exit 0;
