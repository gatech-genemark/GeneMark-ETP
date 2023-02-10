#!/usr/bin/env perl
# -----------------------
# Alex Lomsadse
# 2022
# GaTech
#
# Filter out overlap with HC genes
# -----------------------

use strict;
use warnings;

use Getopt::Long qw(GetOptions);

# -------------
my $hc = '';
my $hints = '';
my $out = '';

my $verbose = 0;
my $debug = 0;
my $warnings = 0;
# -------------
Usage() if ( @ARGV < 1 );
ParseCMD();

my %tr = ();

LoadHC($hc, \%tr);

my %gn = ();

MoveToGenes(\%tr, \%gn);

my %seq = MoveToSeqid(\%gn);

CreateOutput($hints, \%seq, $out);

# ==============
sub CreateOutput
{
	my $in = shift;
	my $ref = shift;
	my $out = shift;

	open( my $IN, $in ) or die "error on open file $in: $!\n";
	open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
	while(my $line = <$IN>)
	{
		my @arr = split('\t', $line);

		if ( exists $ref->{$arr[0]} )
		{
			if (( $ref->{$arr[0]}{'L'} < $arr[3] )and( $arr[4] < $ref->{$arr[0]}{'R'} ))
			{
				print $OUT $line;
			}
			else
			{
				;
			}
		}
		else
		{
			print $OUT $line;
		}
	}
	close $OUT;
	close $IN;
}

# --------------
sub MoveToSeqid
{
	my $ref = shift;

	my %h = ();

	foreach my $gid (keys %{$ref})
	{
		my $seqid = $ref->{$gid}{"seqid"};

		$h{$seqid}{"L"} = 0;
		$h{$seqid}{"R"} = 1000000000;

		my $strand = $ref->{$gid}{"strand"};

		if ( $strand eq "+" )
		{
			$h{$seqid}{"R"} = $ref->{$gid}{"L"} + 3;
		}
		else
		{
			$h{$seqid}{"L"} = $ref->{$gid}{"R"} - 3;
		}
	}

	return %h;
}
# --------------
sub MoveToGenes
{
	my $ref_tr = shift;
	my $ref_gn = shift;

	foreach my $tid (keys %{$ref_tr} )
	{
		my $current_gene = $ref_tr->{$tid}{"gid"};

		if ( ! exists $ref_gn->{$current_gene} )
		{
			$ref_gn->{$current_gene}{"gid"}    = $ref_tr->{$tid}{"gid"};
			$ref_gn->{$current_gene}{"tid"}    = $ref_tr->{$tid}{"tid"};
			$ref_gn->{$current_gene}{"seqid"}  = $ref_tr->{$tid}{"seqid"};
			$ref_gn->{$current_gene}{"L"}      = $ref_tr->{$tid}{"L"};
			$ref_gn->{$current_gene}{"R"}      = $ref_tr->{$tid}{"R"};
			$ref_gn->{$current_gene}{"strand"} = $ref_tr->{$tid}{"strand"};
		}
		else
		{
			if (( $ref_tr->{$tid}{"L"} < $ref_gn->{$current_gene}{"L"} )or( $ref_gn->{$current_gene}{"R"} < $ref_tr->{$tid}{"R"} ))
			{
				$ref_gn->{$current_gene}{"tid"} = $ref_tr->{$tid}{"tid"};
				$ref_gn->{$current_gene}{"L"} = $ref_tr->{$tid}{"L"};
				$ref_gn->{$current_gene}{"R"} = $ref_tr->{$tid}{"R"};
			}
			elsif (( $ref_tr->{$tid}{"L"} == $ref_gn->{$current_gene}{"L"} )and( $ref_gn->{$current_gene}{"R"} == $ref_tr->{$tid}{"R"} ))
			{
				my $f = 1;
				if ( $ref_gn->{$current_gene}{"tid"} =~ /\.(\d+)$/ )
				{
					$f = $1;
				}

				my $s = 1;
				if ( $ref_tr->{$tid}{"tid"} =~ /\.(\d+)$/ )
				{
					$s = $1;
				}

				if ( $f > $s )
				{
					$ref_gn->{$current_gene}{"tid"} = $ref_tr->{$tid}{"tid"};
				}
			}
		}
	}

	if ($verbose)
	{
		my $count = scalar (keys %{$ref_gn});
		print "# number of genes in set: $count\n";
	}
}
# --------------
sub LoadHC
{
	my $fname = shift;
	my $ref = shift;

	die "error, file name is empty: LoadHC\n" if ! $fname;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		my @arr = split('\t');

		my $tid = '';
		my $gid = '';
		
		if ( $arr[8] =~ /transcript_id \"(\S+)\"/ )
		{
			$tid = $1;
		}
		else { die "error, no transcript_id: $_"; }

		if ( $arr[8] =~ /gene_id \"(\S+)\"/ )
		{
			$gid = $1;
		}
		else { die "error, no gene_id: $_"; }

		if ( ! exists $ref->{$tid} )
		{
			$ref->{$tid}{"gid"} = $gid;
			$ref->{$tid}{"tid"} = $tid;
			$ref->{$tid}{"seqid"} = $arr[0];
			$ref->{$tid}{"L"} = $arr[3];
			$ref->{$tid}{"R"} = $arr[4];
			$ref->{$tid}{"strand"} = $arr[6]
		}
		else
		{
			if ( $ref->{$tid}{"strand"} ne $arr[6] )
				{ die "error, mismatch in strand: $tid $gid\n"; }
			if ( $ref->{$tid}{"seqid"} ne $arr[0] )
			{
				print "error, mismatch in seqid: $tid $gid\n" if $verbose; 
				next;
			}

			$ref->{$tid}{"L"} = $arr[3] if ( $arr[3] < $ref->{$tid}{"L"} );
			$ref->{$tid}{"R"} = $arr[4] if ( $ref->{$tid}{"R"} < $arr[4] );
		}
	}

	if ($verbose)
	{
		my $count = scalar (keys %{$ref});
		print "# number of transcripts in file: $count $fname\n";
	}
}
# --------------
sub Usage
{
	print "Usage: $0 --hc [] --hints [] --out []\n";
	print "  --hc []\n";
	print "  --hints []\n";
	print "  --out []\n";
	print "Optional:\n";
	print "  --verbose\n";
	print "  --warnings\n";
	print "  --debug\n";

	exit 1;
}
# --------------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'hc=s'       => \$hc,
		'hints=s'    => \$hints,
		'out=s'      => \$out,
		'verbose'    => \$verbose,
		'warnings'   => \$warnings,
		'debug'      => \$debug,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line: $ARGV[0]\n" if( @ARGV > 0 );

	$verbose = 1 if $debug;
}
# ==============

