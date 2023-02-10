#!/usr/bin/env perl
# ========================
# Alex Lomsadze
# 2019
# Change path to Perl in GeneMark-ES suite
# ========================

use warnings;
use strict;

my $path = shift;

if ( ! $path )
{
	print "usage: $0 new_path_to_perl\n";
	print "   example:  $0 \"/usr/bin/perl\"\n";
	print "   example:  $0 \"/usr/bin/env perl\"\n";
	exit 1;
}

my @files =
(
"bed_to_gff.pl",
"bp_seq_select.pl",
"build_mod.pl",
"calc_introns_from_gtf.pl",
"gc_distr.pl",
"get_sequence_from_GTF.pl",
"gmes_petap.pl",
"histogram.pl",
# "hmm_to_gtf.pl",
"make_nt_freq_mat.pl",
"parse_ET.pl",
"parse_by_introns.pl",
"parse_gibbs.pl",
"parse_set.pl",
"predict_genes.pl",
# "reformat_fasta.pl",  # script was replaced by probuild
"reformat_gff.pl",
"rescale_gff.pl",
"rnaseq_introns_to_gff.pl",
"run_hmm_pbs.pl",
"scan_for_bp.pl",
"star_to_gff.pl",
"verify_evidence_gmhmm.pl",
"run_es.pl",
"get_below_gc.pl",
"compare_intervals_exact.pl",
);

foreach my $file (@files)
{
	if ( -e $file )
	{
		my $txt = "";
		open( my $IN, $file) or die "error on open file $file: $!";
		while(<$IN>)
		{
			$txt .= $_;
		}
		close $IN;

		if ( $txt !~ /^#!/ )
		{
			die "error, unexpected first line format found in $file\n";
		}

		if ( $txt =~ s/^#!.+\n?/#!$path\n/ )
		{
			open( my $OUT, '>' , $file) or die "error on open file $file: $!";
			print $OUT $txt;
			close $OUT;
		}
		else { die "error, unexpected error\n" };
	}
	else
	{
		die "error, required file not found: $file\n";
	}
}

print "done\n";
exit 0;


