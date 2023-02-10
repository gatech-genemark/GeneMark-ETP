#!/usr/bin/env perl
# ------------------------------------------------------------
# Alex Lomsadze
# GeorgiaTech
# 2018
# Run GeneMark-ES on PBS 
# ------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use FindBin qw($Bin);
use Cwd qw( abs_path cwd );
use File::Temp qw( tempfile );
use Data::Dumper;

# ------------------------------------------------------------
my $in = '';
my $out_dir = ".";
my $pbs = "Y";
my $cores = 8;
my $verbose = '';
my $debug = '';
# ------------------------------------------------------------

ParseCMD();

my $es = "$Bin/gmes_petap.pl";

if ( ! -e $es ) { die "error, $es executable not found\n"; }
$es = abs_path($es);

CheckBeforeRun();

my $com = "cd $out_dir\n";

# edit this

$com .= "$es --ES --seq $in --max_intergenic 50000 --cores $cores --v > loginfo ";

# run

if ( $pbs )
{
	RunOnPBS( $com );
}
else
{
	system( $com ) and  die "error in code execution\n";
}

exit 0;

# ------------------------------------------------------------
sub RunOnPBS
{
	my $com = shift;

	my $pbs = UniqTmpFile();
	my $pbs_log = $pbs.".log";

my $text = "#!/bin/bash
#PBS -N es
#PBS -o $pbs_log
#PBS -j oe
#PBS -l nodes=1:ppn=$cores
#PBS -l walltime=24:00:00

$com
";

	open( my $OUT, ">$pbs" ) or die "error on create script for PBS: $pbs\n";
	print $OUT $text;
	close $OUT;
	chmod  0755, $pbs;

	my $id = `qsub $pbs`;

	chomp $id;
	print "qsub ID: $id\n" if $verbose;
}
# ------------------------------------------------
sub UniqTmpFile
{
	my ( $fh, $tmp_name ) = tempfile( "pbs_XXXXX" );
	if ( !fileno($fh) ) { die "error, can't open temporally file: $!\n"; }
	close $fh;
	chmod 0755, $fh;
	return $tmp_name;
}
# ------------------------------------------------------------
sub CheckBeforeRun
{
	if ( ! -e $in ) { die "error, file not found: $in\n"; }

	$pbs = uc $pbs;
	if ( $pbs ne "Y" and $pbs ne "N" )
	{
		die "error on --pbs option value: $pbs\n";
	} 

	if ( ! -e $out_dir )
	{
		mkdir $out_dir or die "error on create directory: $out_dir\n";;
	}
	
	$in = abs_path($in);
	$out_dir = abs_path($out_dir);
}
# ------------------------------------------------------------
sub ParseCMD
{
	if( @ARGV == 0 ) { Usage(); exit 1; }

	my $opt_results = GetOptions
	(
		'in=s'      => \$in,
		'out_dir=s' => \$out_dir,
		'pbs=s'     => \$pbs,
		'cores=i'   => \$cores,
		'verbose'   => \$verbose,
		'debug'     => \$debug,
	);

	if( !$opt_results ) { die "error on command line\n"; }
	if( @ARGV > 0 ) { die "error, unexpected argument found on command line: @ARGV\n"; }
	if ( !$in ) { die "error, required option not found: --in\n"; }

	$verbose = 1 if $debug;
};
# ------------------------------------------------------------
sub Usage
{
	my $pbs_status = '';
	if ( $pbs eq "Y" )
	{
		$pbs_status = "Y";
	}
	else
	{
		$pbs_status = "N";
	}

	my $txt =
"# -----------
Usage: $0  --in  [file_name]
  run GeneMark-ES

  --out_dir    [$out_dir]; output to this folder
  --pbs        [$pbs_status]; run job using PBS: Y or N
  --core       [$cores]; number of cores or threads per job
  --verbose
  --debug
# -----------
";
	print $txt;
	exit 1;
};
# ------------------------------------------------------------

