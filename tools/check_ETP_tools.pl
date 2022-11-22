#!/usr/bin/env perl
# -----------------
# Alex Lomsadze
# 2022
# Georgia Institute of Technology, Atlanta, Georgia, USA
# Copyright: GaTech
#
# This script performs installation test (sort of).
# Script checks if executable files return expected value on execution.
# -----------------

use strict;
use warnings;

print "# checking for presence of GeneMark-ETP dependencies\n";

my $status = 1;

# ----------------

CheckForVerison( "./bedtools", "--version" );
CheckForVerison( "./diamond", "version" );
CheckForVerison( "./fastq-dump", "-V" );
CheckForVerison( "./gffread", "2>&1" );
CheckForVerison( "./hisat2-align-l", "--version" );
CheckForVerison( "./hisat2-align-s", "--version" );
CheckForVerison( "./hisat2-build-l", "--version" );
CheckForVerison( "./hisat2-build-s", "--version" );
CheckForVerison( "./prefetch", "--version" );
CheckForVerison( "./samtools", "--version" );
CheckForVerison( "./stringtie", "2>&1" );
CheckForVerison( "./vdb-config", "--version" );

# ----------------
if ($status)
{
	print "# tools are OK\n";
}
else
{
	print "# WARNING, some issues were detected with required tools!\n";
}

# ===================
sub CheckForVerison
{
	my $cmd = shift;
	my $option = shift;

	my $info = `$cmd $option`;
	$info = uc $info;

	$cmd =~ s/^\.\///;
	$cmd = uc $cmd;

	if ($info =~ /$cmd / )
	{
		print "## $cmd OK\n";
	}
	else
	{
		print "## $cmd - WARNING! check it!\n";
		$status = 0;
	}
}
# ===================

