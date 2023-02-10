#!/usr/bin/env perl
# =============
# Alex Lomsadze
# GaTech
# =============

use strict;
use warnings;

use Data::Dumper;

my $stat_fasta = shift;
my $L = shift;
my $R = shift;

my %stat = LoadStat($stat_fasta);

SafeIDs( "low.ids",    \%stat, 0, $L - 1 );
SafeIDs( "medium.ids", \%stat, $L, $R );
SafeIDs( "high.ids",   \%stat, $R + 1 , 100 );

# -------------
sub SafeIDs
{
	my $fname = shift;
	my $ref = shift;
	my $L = shift;
	my $R = shift;

	open( my $OUT, ">", $fname ) or die "error on open file $fname: $!\n";
	foreach my $id ( keys %{$ref} )
	{
		if (( $ref->{$id} >= $L )and( $ref->{$id} <= $R ))
		{
			print $OUT $id ."\n";
		}
	}
	close $OUT;
}
# -------------
sub LoadStat
{
	my $fname = shift;

	my %h = ();

	open(my $IN, $fname) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		my @arr = split('\s');

		my $id = $arr[0];
		$id =~ s/>//;

		$h{$id} = $arr[2];
	}
	close($IN);

#	print Dumper(\%h);

	return %h
}
# -------------

