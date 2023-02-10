#!/usr/bin/env perl
# -----------------------
# Alex Lomsadse
# 2022
# GaTech
#
# Create regions from coordiantes of high confidence genes
# -----------------------

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Data::Dumper;

# -------------
my $hcc = '';
my $hcp = '';
my $out = '';

my $margin = 0;

my $verbose = 0;
my $debug = 0;
my $warnings = 0;
# -------------
Usage() if ( @ARGV < 1 );
ParseCMD();

my %trc = ();
my %trp = ();

LoadHC($hcc, \%trc);
LoadHC($hcp, \%trp);

my %gnc = ();
my %gnp = ();

MoveToGenes(\%trc, \%gnc);
MoveToGenes(\%trp, \%gnp);

PartialToPoint(\%gnp, \%gnc);

CreateOutput($out, \%gnc);

# ==============
sub CreateOutput
{
	my $fname = shift;
	my $ref = shift;

	my $count = 0;

	my @gid_sorted = sort {
		( $ref->{$a}{"seqid"} cmp $ref->{$b}{"seqid"} )or
		( $ref->{$a}{"L"} <=> $ref->{$b}{"L"} )or
		( $ref->{$a}{"R"} <=> $ref->{$b}{"R"} )
	} keys %{$ref};

	open( my $OUT, ">", $fname ) or die "error on open file $fname: $!\n";
	foreach my $gid (@gid_sorted)
	{
		my $L = $ref->{$gid}{"L"};
		my $R = $ref->{$gid}{"R"};

		print $OUT $ref->{$gid}{"seqid"} ."\thc_gene\tmrna\t". $L ."\t". $R ."\t0\t". $ref->{$gid}{"strand"} ."\t.\tgene_id \"". $ref->{$gid}{"gid"} ."\"; transcript_id \"". $ref->{$gid}{"tid"} ."\";\n";
		$count += 1;
	}
	close $OUT;

	if ($verbose)
	{
		print "# number regions in the output file: $count $fname\n";
	}
}
# --------------
sub PartialToPoint
{
	my $ref = shift;
	my $target = shift;

	my $count = 0;

	foreach my $gid (keys %{$ref})
	{
		if ( ! exists $target->{$gid} )
		{
			$target->{$gid}{"gid"}    = $ref->{$gid}{"gid"};
                        $target->{$gid}{"tid"}    = $ref->{$gid}{"tid"};
                        $target->{$gid}{"seqid"}  = $ref->{$gid}{"seqid"};
			$target->{$gid}{"strand"} = $ref->{$gid}{"strand"};
			if ( $target->{$gid}{"strand"} eq "+" )
			{
				$target->{$gid}{"L"} = $ref->{$gid}{"R"};
				$target->{$gid}{"R"} = $ref->{$gid}{"R"};
			}
			else
			{
				$target->{$gid}{"L"} = $ref->{$gid}{"L"};
				$target->{$gid}{"R"} = $ref->{$gid}{"L"};
			}

			$count += 1;
		}
		else
			{ die "error, not supported yet\n"; }
	}

	if ($verbose)
	{
		print "# number of incomplete regions added: $count\n";
	}
}
# --------------
sub MoveToGenes
{
	my $ref_tr = shift;
	my $ref_gn = shift;

	my $by_rigions = 1;

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
			$ref_gn->{$current_gene}{"cds"}    = $ref_tr->{$tid}{"cds"};
		}
		else
		{
if($by_rigions) # by gene region
{
			if ($ref_gn->{$current_gene}{"L"} > $ref_tr->{$tid}{"L"})
			{
				$ref_gn->{$current_gene}{"L"} = $ref_tr->{$tid}{"L"};
			}

			if ($ref_gn->{$current_gene}{"R"} < $ref_tr->{$tid}{"R"})
			{
				$ref_gn->{$current_gene}{"R"} = $ref_tr->{$tid}{"R"};
			}

			if ( $ref_gn->{$current_gene}{"cds"} < $ref_tr->{$tid}{"cds"} )
			{
				$ref_gn->{$current_gene}{"tid"} = $ref_tr->{$tid}{"tid"};
				$ref_gn->{$current_gene}{"cds"} = $ref_tr->{$tid}{"cds"};
			}

			if ( $ref_gn->{$current_gene}{"cds"} == $ref_tr->{$tid}{"cds"} )
			{
				my $f = 0;
				if ( $ref_gn->{$current_gene}{"tid"} =~ /\.(\d+)$/ )
				{
					$f = $1;
				}

				my $s = 0;
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

if(!$by_rigions) # by cds length
{
			if ( $ref_gn->{$current_gene}{"cds"} < $ref_tr->{$tid}{"cds"} )
			{
				$ref_gn->{$current_gene}{"tid"} = $ref_tr->{$tid}{"tid"};
				$ref_gn->{$current_gene}{"L"}   = $ref_tr->{$tid}{"L"};
				$ref_gn->{$current_gene}{"R"}   = $ref_tr->{$tid}{"R"};
				$ref_gn->{$current_gene}{"cds"} = $ref_tr->{$tid}{"cds"};
			}
			elsif ( $ref_gn->{$current_gene}{"cds"} == $ref_tr->{$tid}{"cds"} )
			{
				my $f = 0;
				if ( $ref_gn->{$current_gene}{"tid"} =~ /\.(\d+)$/ )
				{
					$f = $1;
				}

				my $s = 0;
				if ( $ref_tr->{$tid}{"tid"} =~ /\.(\d+)$/ )
				{
					$s = $1;
				}

				if ( $f > $s )
				{
					$ref_gn->{$current_gene}{"tid"} = $ref_tr->{$tid}{"tid"};
					$ref_gn->{$current_gene}{"L"}   = $ref_tr->{$tid}{"L"};
					$ref_gn->{$current_gene}{"R"}   = $ref_tr->{$tid}{"R"};
				}
			}
}
		}
	}

	foreach my $gid (keys %{$ref_gn})
	{
		$ref_gn->{$gid}{"L"} -= $margin;
		$ref_gn->{$gid}{"L"} = 1 if ( $ref_gn->{$gid}{"L"} < 1 );

		$ref_gn->{$gid}{"R"} += $margin;
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
			$ref->{$tid}{"strand"} = $arr[6];
			$ref->{$tid}{"cds"} = 0;

			if ( $arr[2] eq "CDS" )
			{
				$ref->{$tid}{"cds"} = $arr[4] - $arr[3] + 1;
			}
		}
		else
		{
			if ( $ref->{$tid}{"strand"} ne $arr[6] )
				{ die "error, mismatch in strand: $tid $gid\n"; }
			if ( $ref->{$tid}{"seqid"} ne $arr[0] )
				{ die "error, mismatch in seqid: $tid $gid\n"; }

			$ref->{$tid}{"L"} = $arr[3] if ( $arr[3] < $ref->{$tid}{"L"} );
			$ref->{$tid}{"R"} = $arr[4] if ( $ref->{$tid}{"R"} < $arr[4] );

			if ( $arr[2] eq "CDS" )
			{
				$ref->{$tid}{"cds"} += $arr[4] - $arr[3] + 1;
			}
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
	print "Usage: $0 --hcc [] --hcp [] --out []\n";
	print "  --hcc []\n";
	print "  --hcp []\n";
	print "  --out []\n";
	print "Optional:\n";
	print "  --margin [$margin] margin into intergenic region\n";
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
		'hcc=s'      => \$hcc,
		'hcp=s'      => \$hcp,
		'out=s'      => \$out,
		'margin=i'   => \$margin,
		'verbose'    => \$verbose,
		'warnings'   => \$warnings,
		'debug'      => \$debug,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line: $ARGV[0]\n" if( @ARGV > 0 );

	$verbose = 1 if $debug;
}
# ==============

