#!/usr/bin/env perl
# -----------------------
# Alex Lomsadse
# 2022
# GaTech
#
# Select genes for training
# Input: GTF with some annotation
# Output: GTF with one isoform per gene
# 	- longest CDS
#	- name
#	- complete gene
#	- supported
# -----------------------

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Data::Dumper;

# -------------
my $in = '';
my $out = '';
my $gc = '';
my $hc = '';

my $verbose = 0;
my $debug = 0;
my $warnings = 0;
# -------------
Usage() if ( @ARGV < 1 );
ParseCMD();

my %trans = ();
my %gn = ();

my $def_status = 1;
$def_status = 0 if ($gc or $hc);

LoadTranscripts($in, \%trans);
SelectTranscript(\%trans, \%gn);

if ($gc)
{
	FilterByGClabel($gc, \%gn);
}

if ($hc)
{
	FilterByID($hc, \%gn);
}

CreateOutput($in, \%gn, $out);

# ==============
sub FilterByID
{
	my $name = shift;
	my $ref = shift;

	my %list = LoadList($name);

	my $count = 0;

	foreach my $key (keys %{$ref})
	{
		if ( exists $list{$key} )
		{
			$ref->{$key}{"status"} = 1;
			$count += 1;
		}
	}

	print "# genes found in list $name: $count\n" if $verbose;
}
# --------------
sub LoadList
{
	my $fname = shift;

	my %h = ();

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		if (/(\S+)/)
		{
			$h{$1} = 1;
		}
	}
	close $IN;

	return %h;
}
# --------------
sub FilterByGClabel
{
	my $label = shift;
	my $ref = shift;

	my $count = 0;

	foreach my $key (keys %{$ref})
	{
		my $bin = '';

		if ( $ref->{$key}{"tid"} =~ /^\d+_t([LMH])$/ )
		{
			$bin = $1;

			if ( $bin eq $label )
			{
				$ref->{$key}{"status"} = 1;
				$count += 1;
			}
		}
	}

	print "# genes found with label $label: $count\n" if $verbose;
}
# --------------
sub CreateOutput
{
	my $fname_in = shift;
	my $ref = shift;
	my $fname_out = shift;

	my $count = 0;

	my %h = ();
	foreach my $key (keys %{$ref})
	{
		if ( $ref->{$key}{"status"} )
		{
			$h{ $ref->{$key}{"tid"} } = $key;

			$count += 1;
		}
	}

	print "# genes found for training: $count\n" if $verbose;

	open( my $IN, $fname_in ) or die "error on open file $fname_in: $!\n";
	open( my $OUT, ">", $fname_out ) or die "error on open file $fname_out: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		my @arr = split('\t');

		my $tid = '';

		if ( $arr[8] =~ /transcript_id \"(\S+)\"/ )
		{
			$tid = $1;
		}
		else { die "error, no transcript_id: $_"; }

		if ( exists $h{$tid} )
		{
			print $OUT $_;
		}
	}
	close $OUT;
	close $IN;

	if ($verbose)
	{
		;
	}
}
# --------------
sub SelectTranscript
{
	my $ref_tr = shift;
	my $ref_gn = shift;

	foreach my $tid (keys %{$ref_tr} )
	{
		my $current_gene = $ref_tr->{$tid}{"gid"};

		if ( ! exists $ref_gn->{$current_gene} )
		{
			$ref_gn->{$current_gene}{"gid"}   = $ref_tr->{$tid}{"gid"};
			$ref_gn->{$current_gene}{"tid"}   = $ref_tr->{$tid}{"tid"};
			$ref_gn->{$current_gene}{"cds"}   = $ref_tr->{$tid}{"cds"};
			$ref_gn->{$current_gene}{"start"} = $ref_tr->{$tid}{"start"};
			$ref_gn->{$current_gene}{"stop"}  = $ref_tr->{$tid}{"stop"};
			$ref_gn->{$current_gene}{"status"} = $def_status;
		}
		else
		{
			next if !$ref_tr->{$tid}{"start"};
			next if !$ref_tr->{$tid}{"stop"};

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
	}

	my $partial_removed = 0;

	foreach my $gid (keys %{$ref_gn} )
	{
		if ( !$ref_gn->{$gid}{"start"} or !$ref_gn->{$gid}{"stop"} )
		{
			delete($ref_gn->{$gid});
			$partial_removed += 1;
		}
	}

	if ($verbose)
	{
		my $count = scalar (keys %{$ref_gn});
		print "# number of genes in set: $count\n";
		print "# removed partial: $partial_removed\n";
	}
}
# --------------
sub LoadTranscripts
{
	my $fname = shift;
	my $ref = shift;

	die "error, file name is empty: LoadTranscripts\n" if ! $fname;

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
			$ref->{$tid}{"cds"} = 0;
			$ref->{$tid}{"start"} = 0;
			$ref->{$tid}{"stop"} = 0;

			if ( $arr[2] eq "CDS" )
			{
				$ref->{$tid}{"cds"} = $arr[4] - $arr[3] + 1;
			}
			$ref->{$tid}{"start"} = 1 if ( $arr[2] eq "start_codon" );
			$ref->{$tid}{"stop"} = 1 if ( $arr[2] eq "stop_codon" );
		}
		else
		{
			if ( $arr[2] eq "CDS" )
			{
				$ref->{$tid}{"cds"} += $arr[4] - $arr[3] + 1;
			}
			$ref->{$tid}{"start"} = 1 if ( $arr[2] eq "start_codon" );
			$ref->{$tid}{"stop"} = 1 if ( $arr[2] eq "stop_codon" );
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
	print "Usage: $0 --in [] --out []\n";
	print "  --in []\n";
	print "  --out []\n";
	print "Optional:\n";
	print "  --gc [label]\n";
	print "  --hc [name]\n";
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
		'in=s'       => \$in,
		'out=s'      => \$out,
		'gc=s'       => \$gc,
		'hc=s'       => \$hc,
		'verbose'    => \$verbose,
		'warnings'   => \$warnings,
		'debug'      => \$debug,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line: $ARGV[0]\n" if( @ARGV > 0 );

	$verbose = 1 if $debug;
}
# ==============

