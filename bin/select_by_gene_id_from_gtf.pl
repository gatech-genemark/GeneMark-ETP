#!/usr/bin/env perl
# ---------------------------
# Alex Lomsadze
# GaTech
#
# Select from GTF file lines with transcript ID
# ---------------------------

use strict;
use warnings;

my $v = 1;
my %h;

my $ids_fname = shift;
my $gtf_in = shift;
my $gtf_out = shift;

LoadToHash( $ids_fname, \%h );
SelectFromGFF( $gtf_in, $gtf_out, \%h );
CheckForFound(\%h ) if $v;

print "# done\n" if $v;

# -----------------------------
sub CheckForFound
{
	my $ref = shift;

	my $count_not_found = 0;

	foreach my $key (keys %{$ref})
	{
		if ( $ref->{$key} !=  10 )
		{
			$count_not_found += 1;
#			print "# warning, not found:  $key\n"
		}
	}

	print "# not found in input: $count_not_found\n" if $v;
}
# -----------------------------
sub SelectFromGFF
{
	my $name_in = shift;
	my $name_out = shift;
	my $ref = shift;

	open( my $IN, $name_in ) or die "error on open file $name_in: $!\n";
	open( my $OUT, ">", $name_out ) or die "error on open file $name_out: $!\n";
	while(<$IN>)
	{
		if ( /gene_id \"(\S+?)\";/ )
		{
			my $id = $1;

			if ( exists $ref->{$id} )
			{
				print $OUT $_;

				$ref->{$id} = 10;
			}
		}
		
	}
	close $OUT;
	close $IN;
}
# -----------------------------
sub LoadToHash
{
	my $name = shift;
	my $ref = shift;

	open( my $IN, $name ) or die "error on open file $name: $!\n";
	while(<$IN>)
	{
		next if /^\s*$/;

		if (/^(\S+)\s*/)
		{
			if ( ! exists $ref->{$1} )
			{
				$ref->{$1} = 1;
			}
			else { die "error, unexpected duplication of transcrip ID was detected: $1"; }
		}
		else
			{ die "error, unexpected line format found: $_\n"; }
	}
	close $IN;

	if ($v)
	{
		print "# from file $name parsed IDs: ". (scalar keys %{$ref} ) ."\n"
	}
}
# -----------------------------

