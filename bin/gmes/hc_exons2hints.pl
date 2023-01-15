#!/usr/bin/perl
# ----------------------
# AL
#
# reformat from ProtHint to GeneMark.hmm hints
# ---------------------

use strict;
use warnings;

my $in = shift;

open(my $IN, $in) or die "error on open file $in: $!\n";
while(<$IN>)
{
	if ( /\tCDS\t/ )
	{
		if ( /cds_type=Initial/ )
		{
			s/\tCDS\t/\tInitial\t/;
			print $_;
		}
		elsif ( /cds_type=Internal/ )
		{
			s/\tCDS\t/\tInternal\t/;
			print $_;
		}
		elsif ( /cds_type=Terminal/ )
		{
			s/\tCDS\t/\tTerminal\t/;
			print $_;
		}
		elsif ( /cds_type=Single/ )
		{
			s/\tCDS\t/\tSingle\t/;
			print $_;
		}
	}
}
close $IN;

