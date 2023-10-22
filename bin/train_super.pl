#!/usr/bin/env perl
# =============
# Alex Lomsadze
# GaTech 2022
#
# build GeneMark.hmm model file from HC genes
# =============
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use Cwd qw(abs_path);
use YAML;
use Data::Dumper;

my $bin = $Bin;
# -------------
my $hc = '';
my $dna = '';
my $out = "output.mod";
my $bp = "";

my $verbose = 0;
my $warnings = 0;
my $debug = 0;
# -------------

Usage() if ( @ARGV < 1 );
ParseCMD();

symlink $hc, "hc.gff";
symlink $dna, "genome.fasta";

system("$bin/super.pl --gtf hc.gff --gseq genome.fasta --out set.out > loginfo");

if (!$bp)
{
	system("cp  $bin/super.cfg  run.cfg");
	system("cp  $bin/super.mod  ref.mod");
	system("cp  $bin/super.dur  intergenic.dur");
}
else
{
	system("cp  $bin/super_bp.cfg  run.cfg");
	system("cp  $bin/super_bp.mod  ref.mod");
	system("cp  $bin/super_bp.dur  intergenic.dur");
}

system("$bin/gmes/parse_set.pl  --set_in set.out --min_cds 30 --section ES_C --cfg run.cfg");

system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg  --section start_ATG  --format INI");
system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg  --section stop_TAA   --format TERM_TAA");
system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg  --section stop_TAG   --format TERM_TAG");
system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg  --section stop_TGA   --format TERM_TGA");

system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section donor_GT    --format DONOR_0    --phase 0  --quiet");
system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section acceptor_AG --format ACCEPTOR_0 --phase 0  --quiet");
system("cat  GT.mat > donor.mat");
system("cat  AG.mat > acceptor.mat");

system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section donor_GT    --format DONOR_1    --phase 1  --quiet");
system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section acceptor_AG --format ACCEPTOR_1 --phase 1  --quiet");
system("cat  GT.mat >> donor.mat");
system("cat  AG.mat >> acceptor.mat");

system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section donor_GT    --format DONOR_2    --phase 2  --quiet");
system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section acceptor_AG --format ACCEPTOR_2 --phase 2  --quiet");
system("cat  GT.mat >> donor.mat");
system("cat  AG.mat >> acceptor.mat");

system("mv donor.mat     GT.mat");
system("mv acceptor.mat  AG.mat");

system("$bin/gmes/make_nt_freq_mat.pl  --cfg run.cfg  --section donor_GC  --format DONOR_GC_0  --phase 0  --force");
system("cat  GC.mat > donor_GC.mat");
system("$bin/gmes/make_nt_freq_mat.pl  --cfg run.cfg  --section donor_GC  --format DONOR_GC_1  --phase 1  --force");
system("cat  GC.mat >> donor_GC.mat");
system("$bin/gmes/make_nt_freq_mat.pl  --cfg run.cfg  --section donor_GC  --format DONOR_GC_2  --phase 2  --force");
system("cat  GC.mat >> donor_GC.mat");

system("mv donor_GC.mat GC.mat");

my $str = " --MKCHAIN_L_MARGING 0 --MKCHAIN_R_MARGING 0 --MKCHAIN_PSEUDOCOUNTS 1 --MKCHAIN_PRECISION 8 --revcomp_non --ORDM 5 ";
system("$bin/gmes/probuild --cod cod.seq --non non.seq --mkmod_euk mkch  $str");

if (!$bp)
{
	my ($intron_min, $intron_max) = EstimateIntronMinMax( "intron.len" );
	print "intron $intron_min, $intron_max\n" if $verbose;
	UpdateCFG( "run.cfg", "intron_DUR", "max", $intron_max );
}

system("$bin/gmes/histogram.pl --cfg run.cfg --section intron_DUR");
system("$bin/gmes/histogram.pl --cfg run.cfg --section initial_DUR");
system("$bin/gmes/histogram.pl --cfg run.cfg --section internal_DUR");
system("$bin/gmes/histogram.pl --cfg run.cfg --section terminal_DUR");
system("$bin/gmes/histogram.pl --cfg run.cfg --section single_DUR");

if ($bp)
{
	my $cfg = YAML::LoadFile("run.cfg");

	system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section acceptor_short_AG --format ACC_BP_0 --phase 0");
	system("cat AG_SHORT.mat > acceptor_short.mat");
	system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section acceptor_short_AG --format ACC_BP_1 --phase 1" );
	system("cat AG_SHORT.mat >> acceptor_short.mat");
	system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section acceptor_short_AG --format ACC_BP_2 --phase 2" );
	system("cat AG_SHORT.mat >> acceptor_short.mat");
	system("mv acceptor_short.mat AG_SHORT.mat");

	system("$bin/gmes/bp_seq_select.pl --seq_in $cfg->{'ES_C'}->{'bp_region'} --seq_out $cfg->{'ES_C'}->{'gibbs_seq'}  --max_seq $cfg->{'ES_C'}->{'gibbs_seq_max'}  --bp_region_length  $cfg->{'ES_C'}->{'bp_region_length'}  --min_bp_region_length $cfg->{'ES_C'}->{'min_bp_region_length'}");

                # 9     motif length
                # -n    Use nucleic acid alphabet
                # -r    turn off reverse complements with DNA
                # -nopt Don't print Near Optimal output
                # -m    Do not maximize after near optimal sampling
                # -w    pseduocount weight
                # -Z    Don't write progress info
                # -s    random number generator seed

	system("$bin/gmes/Gibbs3 gibbs.seq 9 -n -r -o gibbs.out -nopt -m -w 0.001 -Z -s 1 -P $bin/gmes/prior.bp -F");

	system("$bin/gmes/parse_gibbs.pl --seq gibbs.seq --gibbs gibbs.out --motif_seq $cfg->{'branch_point'}->{'infile'} --spacer_len $cfg->{'spacer_DUR'}->{'in'} --spacer_seq spacer.seq --tr $cfg->{'ES_C'}->{'tr_bp'}");
	system("$bin/gmes/make_nt_freq_mat.pl --cfg run.cfg --section branch_point --format BRANCH");
	system("$bin/gmes/histogram.pl --cfg run.cfg --section spacer_DUR");
	system("$bin/gmes/scan_for_bp.pl --seq_in $cfg->{'ES_C'}->{'bp_region'}  --gibbs_in gibbs.out  --pos_out $cfg->{'prespacer_DUR'}->{'in'}");
	system("$bin/gmes/histogram.pl --cfg run.cfg --section prespacer_DUR");
	my $str = " --MKCHAIN_L_MARGING 0  --MKCHAIN_R_MARGING 0  --MKCHAIN_PSEUDOCOUNTS 1  --MKCHAIN_PRECISION 6 --ORDM 1 ";
	system("$bin/gmes/probuild --non spacer.seq --mkmod_non spacer.mkch  $str");
}

system("$bin/gmes/build_mod.pl --cfg run.cfg --section ES_C --def ref.mod  --out output.mod");

system("sed 's/INTERGENIC_TYPE \"INTERGENIC_DISTR\"/INTERGENIC_TYPE \"CONSTANT\"/' output.mod > tmp.mod");
system("mv tmp.mod $out");

# ===========
sub UpdateCFG
{
	my $fname = shift;
	my $section = shift;
	my $key = shift;
	my $value = shift;

	my $cfg = YAML::LoadFile( $fname );

	$cfg->{$section}->{$key} = $value;

	open( my $OUT, ">", $fname ) or die "error on open file $fname: $!\n";
        print $OUT Dump($cfg);
        close $OUT;
}
# -----------
sub EstimateIntronMinMax
{
	my $fname = shift;

	my $total = 0;
	my $R_10000 = 0;
	my $R_20000 = 0;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^\s*$/;

		if ( /^(\d+)\s*$/ )
		{
			my $len = $1;

			$total += 1;

			if ($len > 10000) { $R_10000 += 1; }
			if ($len > 20000) { $R_20000 += 1; }
		}
	}
	close $IN;

	my $min = 0;
	my $max = 10000;

	if ( $R_10000/$total < 0.05 )
	{
		$max = 10000;
	}
	else
	{
		$max = 30000;
	}

	print "# estimate: total R_10000 R_20000 ". $total ." ". $R_10000 ." ". $R_20000 ." ". ($R_10000/$total) ." ". ($R_20000/$total)  ."\n" if $debug;

	return ($min, $max);
}
# -----------
sub Usage
{
	print "Usage: $0 --hc [$hc] --dna [$dna] --out [$out]\n";
	print "  --hc [$hc]  high confidence genes GTF\n";
	print "  --dna [$dna] genome FASTA\n";
	print "  --out [$out] model file\n";
	print "  --bp  train fungi model\n";
	print "Developer:\n";
	print "  --verbose\n";
	print "  --debug\n";

	exit 1;
}
# -----------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
                'hc=s'    => \$hc,
                'dna=s'   => \$dna,
		'out=s'   => \$out,
		'bp'      => \$bp,
                'verbose' => \$verbose,
                'debug'   => \$debug,
        );

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line: $ARGV[0]\n" if( @ARGV > 0 );

	$verbose = 1 if $debug;

	die "error, option --hc is not set\n" if !$hc;
	die "error, file with training set is not found: $hc\n" if ( ! -e $hc);
	$hc = abs_path($hc);

	die "error, option --dna is not set\n" if !$dna;
	die "error, file with genome sequence is not found: $dna\n" if ( ! -e $dna);
	$dna = abs_path($dna);
}
# -----------

