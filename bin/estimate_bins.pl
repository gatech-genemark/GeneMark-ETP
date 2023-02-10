#!/usr/bin/env perl
# =============
# Alex Lomsadze
# GaTech
# =============

use strict;
use warnings;

use Data::Dumper;

my $hc_region = shift;
my $hc_trace = shift;
my $stat_fasta = shift;
my $verbose = shift;

$verbose = 0 if ( ! defined $verbose );

my %gids = LoadRegions($hc_region);
my %trace = LoadTrace($hc_trace);
my %stat = LoadStat($stat_fasta);
my %r2gid = RegionToGid(\%trace, \%gids);
my %gid2gc = Gid2gc( \%r2gid, \%stat );

my @data = CollectGC(\%stat);
my @bins = SplitByBins(\@data);

SafeIDs( "low.ids",    \%gid2gc, 0, $bins[0]        );   # 0,    $bins[0] - 1
SafeIDs( "medium.ids", \%gid2gc, $bins[0], $bins[1] );   # $bins[0], $bins[1]
SafeIDs( "high.ids",   \%gid2gc, $bins[1]     , 100 );   # $bins[1] + 1 , 100

print join( " ", @bins ) ."\n";

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
sub SplitByBins
{
	my $ref = shift;

	my $total = SumArr($ref, 0, 100);

	my $bin_size = 9;

	my $bin_pos_min = 4;
	my $bin_pos_max = 96;
	my $bin_pos_best = 0;
	my $bin_sum_best = 0;

	print "# bin-center   number of genes in bin $bin_size\n" if $verbose;

	for( my $pos = $bin_pos_min; $pos <= $bin_pos_max; $pos += 1)
	{
		my $sum = 0;

		for( my $i = $pos - 4; $i <= $pos + 4; $i += 1 )
		{
			$sum += $ref->[$i];
		}

		print $pos ."\t". $sum ."\n" if $verbose;

		if ( $sum >= $bin_sum_best )
		{
			$bin_sum_best = $sum;
			$bin_pos_best = $pos;
		}
	}

	my $L_bin = 0;
	my $M_bin = $bin_sum_best;
	my $H_bin = 0;

#	if ( $bin_sum_best/$total < 0.8 )
	{
		$L_bin = SumArr($ref, 0, $bin_pos_best - 5 );
		$H_bin = SumArr($ref, $bin_pos_best + 5, 100 );
	}

	print "# total: low medium high : - - -\n" if $verbose;
	print "$total : $L_bin $M_bin $H_bin : ". ($bin_pos_best - 4) ." ". $bin_pos_best ." ". ($bin_pos_best + 4) ."\n" if $verbose;

	if (( $L_bin >= 0 )and( $H_bin >= 0 )and( $M_bin >= 0 ))
	{
		return ($bin_pos_best - 4), ($bin_pos_best + 4); 
	}
	else
	{
		die "error, add missin section\n";
	}
}
# -------------
sub SumArr
{
	my $ref = shift;
	my $L = shift;
	my $R = shift;

	my $total = 0;

	for( my $i = $L; $i <= $R; $i += 1 )
	{
		$total += $ref->[$i];
	}

	return $total;
}
# -------------
sub CollectGC
{
	my $ref = shift;

	my @arr = ();

	for my $i (0..100)
	{
		$arr[$i] = 0;
	}

	foreach my $key (keys %{$ref})
	{
		$arr[$ref->{$key}] += 1;
	}

	if ( $verbose )
	{
		print "# GC  number of genes\n";
		for my $i (0..100)
		{
			print $i ."\t". $arr[$i] ."\n";
		}
	}

#	print Dumper(\@arr);

	return @arr;
}
# -------------
sub Gid2gc
{
	my $ref_gid = shift;
	my $ref_gc = shift;;

	my %h = ();

	foreach my $key (keys %{$ref_gid} )
	{
		$h{ $ref_gid->{$key} } = $ref_gc->{$key};
	}

#	print Dumper(\%h);

	return %h;
}
# -------------
sub RegionToGid
{
	my $reg = shift;
	my $gid = shift;

	my %h = ();

	foreach my $key (keys %{$reg})
	{
		if ( exists $gid->{$key} )
		{
			$h{ $reg->{$key} } = $gid->{$key};
		}
		else { die "error, gid not found for region: $key $reg->{$key}\n"; }
	}

#	print Dumper(\%h);

	return %h;
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
sub LoadTrace
{
	my $fname = shift;

	my %h = ();

	open(my $IN, $fname) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		my @arr = split('\s');

		my $key = $arr[1] ." ". $arr[2] ." ". $arr[3];

		$h{$key} = $arr[0];
	}
	close($IN);

#	print Dumper(\%h);

	return %h;
}
# -------------
sub LoadRegions
{
	my $fname = shift;

	my %h = ();

	open(my $IN, $fname) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		my @arr = split('\t');

		my $key = $arr[0] ." ". $arr[3] ." ". $arr[4];

		my $id = '';
		if ( /gene_id \"(\S+?)\"/ )
		{
			$id = $1;
		}

		$h{$key} = $id;
	}
	close($IN);

#	print Dumper(\%h);

	return %h;
}
# -------------

