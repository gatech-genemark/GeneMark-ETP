#!/usr/bin/env perl
# -------------
# Alex Lomsadze
# 2022
# 
# "GeneMark-ETP: Automatic Gene Finding in Eukaryotic Genomes in Consistence with Extrinsic Data"
# Tomas Bruna, Alexandre Lomsadze and Mark Borodovsky
# Publication in preparation 2022
# Georgia Institute of Technology, Atlanta, Georgia, USA
#
# Questions?
# Please contact Alex Lomsadze alexl@gatech.edu
#
# Copyright
# Georgia Institute of Technology, Atlanta, Georgia, USA
# -------------

use strict;
use warnings;

use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use Cwd qw(abs_path);
use YAML::XS qw(LoadFile);
use Data::Dumper;

my $bin = $Bin;
my $version = "1.08";
# -------------
my $cfg = '';
my $workdir = '.';
my $protdb_name = ''; # is initialized from $protdb_path

my $gc = 0;
my $bp = 0;
my $extend = 0;

my $penalty; # undefined = auto detection of penalty
my $local = '';
my $bam = '';
my $softmask = '';

my $cores = 64;
my $force = 0;

my $clean = 0;
my $clean_deep = 0;

my $paper = '0';

my $verbose = 0;
my $warnings = 0;
my $debug = 0;
# -------------
my $penalty_low;
my $penalty_medium;
my $penalty_high;
# -------------
# Species specific settings from --cfg YAML file
my $species = '';
my $genome_path = '';
my $RepeatMasker_path = '';
my @rnaseq_sets = ();
my $protdb_path = '';
my $annot_path = ''; # not used in this release
# -------------
# configuration for ab initio model training
my %train_cfg = ();
$train_cfg{"skip_extended"} = 0;         # skip extended training
$train_cfg{"extend_all"} = 0;            # run extended training for all
$train_cfg{"extend_low"} = 0;            # run extended training in low GC bin
$train_cfg{"extend_medium"} = 0;         # run extended training in medium GC bin 
$train_cfg{"extend_high"} = 0;           # run extended training in high GC bin
$train_cfg{"min_genes_to_skip"} = 30000; # threshold = number of genes in training = to skip extended training
$train_cfg{"max_iterations"} = 3;        # maximum number of iterations; 0 - initial prediction
$train_cfg{"iteration"} = 0;             # current iteration
$train_cfg{"model_0"} = '';              # model for initial prediction $workdir/$protdb_name/model/output.mod
$train_cfg{"model_0_low"} = '';          # model for initial prediction $workdir/$protdb_name/low/output.mod
$train_cfg{"model_0_medium"} = '';       # model for initial prediction $workdir/$protdb_name/medium/output.mod
$train_cfg{"model_0_high"} = '';         # model for initial prediction $workdir/$protdb_name/high/output.mod

# =============
Usage() if ( @ARGV < 1 );
ParseCMD();
LoadConfig($cfg);
Start() if $verbose;

PrepareGenome($genome_path) if 1;
if (!$softmask)
{
	MaskGenome($RepeatMasker_path) if 1;
}
$protdb_name = DownloadProteinDB($protdb_path) if 1;

# must run protdb_name initialization
$protdb_name = PrepareFolderForTraining( $protdb_path, $protdb_name );
print "# ProteinDB name: $protdb_name\n" if $verbose;

if (!$bam)
{
	PrepareRNASeq() if 1;
}
MapRNASeq()         if 1;
AssembleTrans()     if 1;
CreateIntronHints() if 1;
GeneMarkST()        if 1;
RunProtHintOnGMST($protdb_name)          if 1;
FilterGMST( $protdb_name, $protdb_name ) if 1;

my $hcc_genes      = "$workdir/rnaseq/hints/$protdb_name/complete.gtf";
my $hcp_genes      = "$workdir/rnaseq/hints/$protdb_name/incomplete.gtf";
my $rnaseq_hints   = "$workdir/rnaseq/hints/hintsfile_merged.gff";
my $prothint_evi   = "$workdir/rnaseq/hints/$protdb_name/prothint/evidence.gff";
my $prothint_hints = "$workdir/rnaseq/hints/$protdb_name/prothint/prothint.gff";

PrepareGenomeTraining($protdb_name)     if 1;
PrepareHCregions($protdb_name)          if 1;

# isGChetero($protdb_name); # just information in this version
if ($gc)
{
	TrainModelGC($protdb_name)      if 1;
}
else
{
	TrainModel($protdb_name)        if 1;
}

$train_cfg{"model_0"}        = "$workdir/$protdb_name/model/output.mod";
$train_cfg{"model_0_low"}    = "$workdir/$protdb_name/low/output.mod";
$train_cfg{"model_0_medium"} = "$workdir/$protdb_name/medium/output.mod";
$train_cfg{"model_0_high"}   = "$workdir/$protdb_name/high/output.mod";

if ( ! defined $penalty )
{
	EstinateMaskingPenalty($protdb_name)    if 1;
}
LoadPenalty($protdb_name);
PredictForProtHint($protdb_name)        if 1;
RunProtHint($protdb_name, $protdb_name) if 1;

PredictNonhc($protdb_name)              if 1;
ExtendEvidence($protdb_name)            if 1;

ITERATION:

PredictNonhc($protdb_name)              if 1;
Combine($protdb_name)                   if 1;
SelectSupported($protdb_name)           if 1;

EstimateAccPaper() if ($paper);

if ( !$gc and ExtendTraining($protdb_name) )
{
	print "# iteration ". $train_cfg{"iteration"} ."\n" if $verbose;

	ChDir($workdir);

	my $old_name = "$protdb_name/genemark.gtf";
	my $new_name = "$protdb_name/genemark.gtf_". $train_cfg{"iteration"};
	system ("cp $old_name $new_name" );
	$old_name = "$protdb_name/genemark_supported.gtf";
	$new_name = "$protdb_name/genemark_supported.gtf_". $train_cfg{"iteration"};
	system ("cp $old_name $new_name" );

	if ( $train_cfg{"iteration"} < $train_cfg{"max_iterations"} )
	{
		$force = 1;

		UpdateModel($protdb_name, $train_cfg{"iteration"}, "$workdir/$protdb_name/genemark_supported.gtf");

		$train_cfg{"model_etr"}        = "$workdir/$protdb_name/etr/". $train_cfg{"iteration"} ."/model/output.mod";
		$train_cfg{"model_etr_low"}    = "$workdir/$protdb_name/etr/". $train_cfg{"iteration"} ."/low/output.mod";
		$train_cfg{"model_etr_medium"} = "$workdir/$protdb_name/etr/". $train_cfg{"iteration"} ."/medium/output.mod";
		$train_cfg{"model_etr_high"}   = "$workdir/$protdb_name/etr/". $train_cfg{"iteration"} ."/high/output.mod";

		$train_cfg{"iteration"} += 1;

		if ( $train_cfg{"iteration"} > 0 )
		{
			ChDir($workdir);
			my $prev_i = $train_cfg{"iteration"} - 1;
			my $curr_i = $train_cfg{"iteration"};
			my $out_i  = "cds_". $prev_i ."_". $curr_i;
			system( "$bin/compare_intervals_exact.pl --f1 $workdir/$protdb_name/etr/$prev_i/model/training.gtf --f2 $workdir/$protdb_name/etr/$curr_i/model/training.gtf > $out_i" );

			my $F1 = F1_from_acc($out_i);

			if ($F1 >= 0)
			{
				goto LAST;
			}
		}

#		goto ITERATION;
	}
}

LAST:

Done() if $verbose;

# =============
sub F1_from_acc
{
	my $fname = shift;

	my $F1 = 0;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		;
	}
	close $IN;

	return $F1;
}
# -------------
sub EstimateAccPaper
{
	ChDir($workdir);

	if ( -e "annot/reliable.gtf" and -e "annot/union.gtf" and -e "annot/pseudo.gff3" )
	{
		system( "$bin/compare_intervals_exact.pl --f2 annot/reliable.gtf --f1 $protdb_name/genemark_supported.gtf --pseudo annot/pseudo.gff3 --gene > $protdb_name/sn_". $train_cfg{"iteration"} );
		system( "$bin/compare_intervals_exact.pl --f1 annot/union.gtf    --f2 $protdb_name/genemark_supported.gtf --pseudo annot/pseudo.gff3 --gene > $protdb_name/sp_". $train_cfg{"iteration"} );
	}
	elsif (  -e "annot/annot.gtf" and -e "annot/pseudo.gff3" )
	{
		system( "$bin/compare_intervals_exact.pl --f1 annot/annot.gtf    --f2 $protdb_name/genemark.gtf           --pseudo annot/pseudo.gff3 --gene > $protdb_name/snsp_". $train_cfg{"iteration"} );
	}
}
# -------------
sub UpdateModel
{
	my $wd = shift;
	my $itr = shift;
	my $prev = shift;

	print "\n### update training set: $itr\n" if $verbose;

	ChDir("$workdir/$wd");
	MkDir("$workdir/$wd/etr");
	ChDir("$workdir/$wd/etr");
	MkDir("$workdir/$wd/etr/$itr");
	ChDir("$workdir/$wd/etr/$itr");

	if (!$gc)
	{
		MkDir("$workdir/$wd/etr/$itr/model");
		ChDir("$workdir/$wd/etr/$itr/model");

		if ( $force or ( ! -e "output.mod" ))
		{
			system( "$bin/select_for_training.pl --in $prev --out training.gtf --verbose" );

			if(!$bp)
			{
				system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
			}
			else
			{
				system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta --bp" );
			}
		}
	}
	else
	{
		if ( $train_cfg{"extend_low"} )
		{
			MkDir("$workdir/$wd/etr/$itr/low");
			ChDir("$workdir/$wd/etr/$itr/low");
			if ( $force or ( ! -e "output.mod" ))
			{
				system( "$bin/select_for_training.pl --in $prev --out training.gtf --verbose --gc L --hc $workdir/$wd/hc/low.ids" );
				system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
			}
		}

		if ( $train_cfg{"extend_medium"} )
		{
			MkDir("$workdir/$wd/etr/$itr/medium");
			ChDir("$workdir/$wd/etr/$itr/medium");
			if ( $force or ( ! -e "output.mod" ))
			{
				system( "$bin/select_for_training.pl --in $prev --out training.gtf --verbose --gc M --hc $workdir/$wd/hc/medium.ids" );
				system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
			}
		}

		if ( $train_cfg{"extend_high"} )
		{
			MkDir("$workdir/$wd/etr/$itr/high");
			ChDir("$workdir/$wd/etr/$itr/high");
			if ( $force or ( ! -e "output.mod" ))
			{
				system( "$bin/select_for_training.pl --in $prev --out training.gtf --verbose --gc H --hc $workdir/$wd/hc/high.ids" );
				system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
			}
		}
	}

	print "### update training set: $itr ... done\n" if $verbose;
}
# -------------
sub ExtendTraining
{
	my $wd = shift;

	print "\n### check - run or not - extended training\n" if $verbose;
	
	if ( $train_cfg{"skip_extended"} )
	{
		print "# skipping extended training based on algo setup\n" if $verbose;
	}
	else
	{
		my $min = $train_cfg{"min_genes_to_skip"};

		ChDir("$workdir/$wd");

		my $count_low = 0;
		my $count_medium = 0;
		my $count_high = 0;
		my $count_all = 0;

		if (!$gc)
		{
			$count_all = `cut -f4 -d'"' ./model/training.gtf | uniq | sort | uniq | wc -l`;
			chomp $count_all;
		}
		else
		{
			$count_low    = CountLines("./hc/low.ids");
			$count_medium = CountLines("./hc/medium.ids");
			$count_high   = CountLines("./hc/high.ids");
			$count_all = $count_low + $count_medium + $count_high;
		}

		if (!$gc)
		{
			$train_cfg{"extend_all"} = 1 if $count_all < $min;
		}
		else
		{
			$train_cfg{"extend_low"}    = 1 if $count_low < $min;
			$train_cfg{"extend_medium"} = 1 if $count_medium < $min;
			$train_cfg{"extend_high"}   = 1 if $count_high < $min;
		}

		$train_cfg{"skip_extended"} = 0 if ( $train_cfg{"extend_all"} or $train_cfg{"extend_low"} or $train_cfg{"extend_medium"} or $train_cfg{"extend_high"} );

		if ( $verbose and $train_cfg{"skip_extended"} )
		{
			print "# skip extended training: $count_all = $count_low + $count_medium + $count_high\n";
		}
		elsif ( $verbose and !$train_cfg{"skip_extended"} )
		{
			print "# run extended training: $count_all = $count_low + $count_medium + $count_high\n";
		}
	}

	print "### check - run or not - extended training ... done\n" if $verbose;

	return ! $train_cfg{"skip_extended"};
}
# -------------
sub CountLines
{
	my $fname = shift;

	my $count = 0;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^\s*$/;
		$count += 1;
	}

	return $count;
}
# -------------
sub PrepareFolderForTraining
{
	my $fpath = shift;
	my $name = shift;

	die "error, file path to protein db is empty\n" if !$fpath;

	$name = GetFileNameFromPathName($fpath) if !$name;

	ChDir($workdir);
	MkDir("$workdir/$name");

	print "# folder for training: $name\n" if $verbose;

	return $name;
}
# -------------
sub LoadPenalty
{
	my $wd = shift;

	print "\n### load penalty values\n" if $debug;

	ChDir("$workdir/$wd");
	MkDir("penalty");
	ChDir("penalty");

	if ( defined $penalty )
	{
		if (!$gc)
		{
			system( "echo $penalty > penalty.value" );
		}
		else
		{
			system( "echo $penalty > penalty.value.low" );
			system( "echo $penalty > penalty.value.medium" );
			system( "echo $penalty > penalty.value.high" );
		}
	}

	if (!$gc)
	{
		$penalty = LoadPenaltyFromFile("penalty.value");
	}
	else
	{
		$penalty_low    = LoadPenaltyFromFile("penalty.value.low");
		$penalty_medium = LoadPenaltyFromFile("penalty.value.medium");
		$penalty_high   = LoadPenaltyFromFile("penalty.value.high");
	}

	print "\n### load penalty values ... done\n" if $debug;
}
# -------------
sub LoadPenaltyFromFile
{
	my $fname = shift;
	
	my $value = 0;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^\s*$/;
		next if /^#/;

		if( /^\s*(\S+)/ )
		{
			$value = $1;
		}
		else
			{ die "error, unexpected file format found in file penalty.value $_\n"; }
	}
	close $IN;

	print "# $fname = $value\n" if $debug;

	return $value;
}
# -------------
sub RunProtHintOnGMST
{
	my $name = shift;

	print "\n### generate ProtHint predictions with GeneMarkS-T seeds\n" if $verbose;

	ChDir("$workdir/rnaseq/hints");
	MkDir("$workdir/rnaseq/hints/$name");
	ChDir("$workdir/rnaseq/hints/$name");
	MkDir("$workdir/rnaseq/hints/$name/prothint");
	ChDir("$workdir/rnaseq/hints/$name/prothint");

	if ( CreateThis("evidence.gff") or CreateThis("prothint.gff") )
	{
		system( "$bin/gmes/ProtHint/bin/prothint.py  --geneSeeds $workdir/rnaseq/gmst/genome_gmst.gtf  $workdir/data/genome.softmasked.fasta  $workdir/data/$name --longGene 30000000 --longProtein 15000000 2> $workdir/prothint_gmst.log" );
	}

	StopIfNotFound("evidence.gff");
	StopIfNotFound("prothint.gff");

	print "### generate ProtHint predictions with GeneMarkS-T seeds ... done\n\n" if $verbose;
}
# =============
sub LoadConfig
{
	my $fname = shift;

	die "error, file name is empty in LoadConfig\n" if !$fname;

	my $config = LoadFile($fname);

	$species           = $config->{"species"}            if $config->{"species"};
	$genome_path       = $config->{"genome_path"}        if $config->{"genome_path"};
	$RepeatMasker_path = $config->{"RepeatMasker_path"}  if $config->{"RepeatMasker_path"};
	@rnaseq_sets       = ();
	@rnaseq_sets       = @{$config->{"rnaseq_sets"}}     if $config->{"rnaseq_sets"};
	$protdb_path       = $config->{"protdb_path"}        if $config->{"protdb_path"};
	$annot_path        = $config->{"annot_path"}         if $config->{"annot_path"};

	if ($debug)
	{
		print "\n";
		print "### from YAML species configuration file and command line\n";
		print "# cfg: $fname\n";
		print "# workdir: ". $workdir ."\n";
		print "# gc: ". $gc ."\n";
		print "# fungus: ". $bp ."\n";
		print "# extended ". $extend ."\n";
		if (defined $penalty) { print "# penalty: ". $penalty ."\n"; } else { print "# penalty: auto\n"; }
		print "# softmask: ". $softmask ."\n";
		print "# local: ". $local ."\n";
		print "# bam: ". $bam ."\n";
		print "# force: ". $force ."\n";
        print "# cores: ". $cores ."\n";
        print "# clean: ". $clean ."\n";
        print "# verbose: ". $verbose ."\n";
        print "# warnings: ". $warnings ."\n";
        print "# debug: ". $debug ."\n";
		print "# settings in cfg file:\n";
		print "#   species: ". $species ."\n";
		print "#   genome_path: ". $genome_path ."\n";
		print "#   RepeatMasker_path: ". $RepeatMasker_path ."\n";
		print "#   rnaseq_sets:";
		foreach my $name (@rnaseq_sets)
		{
			print " ". $name;
		}
		print "\n";
		print "#   protdb_path: ". $protdb_path ."\n";
		print "#   annot_path: ". $annot_path ."\n";
		print "### from YAML species configuration file and command line ... done\n\n";
	}
}
# -------------
sub Usage
{
	print "\n";
	print "ETP version $version\n";
	print "Usage: $0 --cfg [$cfg] --workdir [$workdir]\n";
	print "  --cfg [$cfg]      species configuration file in YAML format\n";
	print "  --workdir [$workdir] working directory; default is current directory\n";
	print "Optional parameters\n";
	print "Algo:\n";
	print "  --gc [$gc]      set algo into gc homogeneous [0] or heterogeneous [1] mode\n";
	print "  --fungus      run algorithm in fungus mode\n";
	print "  --extend [$extend]  run extended training\n";
	print "Algo modifiers:\n";
	print "  --softmask    input genome is in repeat soft-mask FASTA format\n";
	print "  --local [$local]    folder with RNA-Seq reads in FASTQ format; by default RNA-Seq files are loaded from SRA\n";
	print "  --bam [$bam]      file with reads to genome alignment in BAM format; by default BAM file is created by ETP\n";
	print "                BAM file must be sorted; BAM file must be compatible with StringTie2\n";
	print "  --penalty []  masking penalty; default mode is auto detection when penalty value is not provided\n";
	print "  --force       on algo re-run overwrite intermediate results\n";
	print "Run:\n";
	print "  --cores [$cores]  number of CPU/cores to use on multi CPU/core computer\n";
	print "  --clean       delete some temporary files\n";
	print "  --clean_deep  delete most of temporary files\n";
	print "Developer:\n";
	print "  --verbose\n";
	print "  --warnings\n";
	print "  --debug\n";
	print "  --paper [$paper]   parse genome FASTA ID's as in ETP publication\n";
	print "Output:\n";
	print "  Predicted genes are located in folder with protein database name\n";
	print "  GTF: genemark.gtf\n";
	print "\n";
	
	exit 1;
}
# -------------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'cfg=s'       => \$cfg,
		'workdir=s'   => \$workdir,
		'fungus'      => \$bp,
		'softmask'    => \$softmask,
		'local=s'     => \$local,
		'bam=s'       => \$bam,
		'cores=i'     => \$cores,
		'clean'       => \$clean,
		'clean_deep'  => \$clean_deep,
		'force'       => \$force,
		'verbose'     => \$verbose,
		'warnings'    => \$warnings,
		'debug'       => \$debug,
		'penalty=f'   => \$penalty,
		'gc=i'        => \$gc,
		'extend=i'    => \$extend,
		'paper=i'     => \$paper,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line: $ARGV[0]\n" if( @ARGV > 0 );

	$verbose = 1 if $debug;
	$warnings = 1 if $debug;
	$clean = 1 if $clean_deep;

	die "error, --cfg is not set\n" if !$cfg;
	StopIfNotFound($cfg);
	$cfg = abs_path($cfg);

	die "error, --workdir is not set\n" if !$workdir;
	StopIfNotFound($workdir);
	$workdir = abs_path($workdir);
	die "error, working directory must differ from installation directory: $workdir\n" if ( $bin eq $workdir);

	if ($local)
	{
		StopIfNotFound($local);
		$local = abs_path($local);
	}

	if ($bam)
	{
		StopIfNotFound($bam);
		$bam = abs_path($bam);
	}
}
# -------------
sub FilterGMST
{
	my $prot_db_name = shift;
	my $out_folder = shift;

	print "\n### filter gmst predictions\n" if $verbose;

	ChDir("$workdir/rnaseq/hints");
	MkDir("$workdir/rnaseq/hints/$out_folder");
	ChDir("$workdir/rnaseq/hints/$out_folder");

	my $in_gmst_gff = "$workdir/rnaseq/gmst/transcripts_merged.fasta.gff";
	my $in_tr_fasta = "$workdir/rnaseq/stringtie/transcripts_merged.fasta";
	my $in_tr_gff   = "$workdir/rnaseq/stringtie/transcripts_merged.gff";
	my $prot_db     = "$workdir/data/$prot_db_name";
	my $in_genome   = "$workdir/data/genome.softmasked.fasta";

	if ( CreateThis("complete.gtf") )
	{
		system( "$bin/GeneMarkSTFiltering/filter.py --minUnsupportedLogOdds 0 --PROTHINT_PATH $bin/gmes/ProtHint/bin --outputFolder . $in_gmst_gff $in_tr_fasta $in_tr_gff $prot_db $in_genome 2> $workdir/filter_gmst.log " );

		system( "$bin/compare_intervals_exact.pl --f1 complete.gtf --f2 complete.gtf --trans --shared12 --out complete.id --original 1" );
		system( "$bin/select_by_transcript_id_from_gtf.pl complete.id complete.gtf complete_uniq.gtf" );
		system( "mv complete_uniq.gtf complete.gtf" );
	}

	print "### filter gmst predictions ... done\n\n" if $verbose;
}
# -------------
sub GetFileNameFromPathName
{
	my $fpath = shift;

	die "error, file path is empty\n" if (!defined $fpath);

	my $fname = '';

	if (( $fpath =~ /.+\/(.+?)\.gz$/ )or( $fpath =~ /.+\/(.+?)$/ ))
	{
		$fname = $1;
	}
	else
		{ die "error on fname parsing from fpath: $fpath\n"; }

	return $fname;
}
# -------------
sub DownloadProteinDB
{
	my $fpath = shift;

	print "\n### download protein db\n" if $verbose;

	ChDir("$workdir/data");

	my $fname = FileFromPath($fpath);

	print "### download protein db ... done\n\n" if $verbose;

	return $fname;
}
# -------------
sub GeneMarkST
{
	print "\n### predict genes gmst\n" if $verbose;

	ChDir("$workdir/rnaseq");
	MkDir("$workdir/rnaseq/gmst");
	ChDir("$workdir/rnaseq/gmst");

	my $tseq = "$workdir/rnaseq/stringtie/transcripts_merged.fasta";
	StopIfNotFound($tseq);

	if( $force or ( ! -e "transcripts_merged.fasta.gff" ))
	{
		system("$bin/gmst/gmst.pl --format GFF $tseq" );
	}
	StopIfNotFound("transcripts_merged.fasta.gff");

	if ( $force or ( ! -e "genome_gmst.gtf" ))
	{
		system( "$bin/gms2hints.pl --tseq $tseq --tgff transcripts_merged.fasta.gff --ggtf ../stringtie/transcripts_merged.gff --out genome_gmst.gtf --long --oneiso" );
	}
	StopIfNotFound("genome_gmst.gtf");

	if ( $force or ( ! -e "genome_gmst_for_HC.gtf" ))
	{
		system( "$bin/gms2hints.pl --tseq $tseq --tgff transcripts_merged.fasta.gff --ggtf ../stringtie/transcripts_merged.gff --out genome_gmst_for_HC.gtf --nodup  --oneiso --min_cds 300");
	}
	StopIfNotFound("genome_gmst_for_HC.gtf");

	print "### predict genes gmst ... done\n\n" if $verbose;
}
# -------------
sub SelectHCset
{
	my $pcf = shift;
	my $pif = shift;
	my $gmstf = shift;
	my $out = shift;

	my %tid_pcf = ();
	my %gid_pcf = ();

	my %tid_pif = ();
	my %gid_pif = ();

	my %tid_gmstf = ();
	my %gid_gmstf = ();

	my $count_ext = 0;
	my $IN;

	open( $IN, $pcf ) or die "error on open file $pcf: $!\n";
	while(<$IN>)
	{
		my $tid = '';
		my $gid = '';
		if ( /transcript_id \"(\S+)\"/ ) { $tid = $1; }
		if ( /gene_id \"(\S+)\"/ )       { $gid = $1; }

		$tid_pcf{$tid} = $gid;
		$gid_pcf{$gid} = 1;
	}
	close $IN;

	open( $IN, $pif ) or die "error on open file $pif: $!\n";
	while(<$IN>)
	{
		my $tid = '';
		my $gid = '';
		if ( /transcript_id \"(\S+)\"/ ) { $tid = $1; }
		if ( /gene_id \"(\S+)\"/ )       { $gid = $1; }

		$tid_pif{$tid} = $gid;
		$gid_pif{$gid} = 1;
	}
	close $IN;

	open( $IN, $gmstf ) or die "error on open file $gmstf: $!\n";
	while(<$IN>)
	{
		if ( /complete/ and /LORF_UPSTOP/ )
		{
			my $tid = '';
			my $gid = '';

			if ( /transcript_id \"(\S+)\"/ ) { $tid = $1; }
			if ( /gene_id \"(\S+)\"/ )       { $gid = $1; }

			if ( ! exists $tid_pcf{$tid} and ! exists $gid_pif{$gid} )
			{
				$count_ext += 1 if ( ! exists $tid_gmstf{$tid} );

				$tid_gmstf{$tid} = $gid;
				$gid_gmstf{$gid} = 1;
			}
		}
	}
	close $IN;

	print "# extended by : $count_ext\n" if $debug;

	system( "cp $pcf $out" );
	open( my $OUT, ">>", $out ) or die "error on open file $out: $!\n";
	open( $IN, $gmstf ) or die "error on open file $gmstf: $!\n";
	while(<$IN>)
	{
		my $tid = '';
		if ( /transcript_id \"(\S+)\"/ ) { $tid = $1; }

		if ( exists $tid_gmstf{$tid} )
		{
			print $OUT $_;
		}
	}
	close $IN;
	close $OUT;
}
# -------------
sub ExtendHCset
{
	my $wd = shift;

	print "\n### extend HC set\n" if $verbose;

	ChDir("$workdir/rnaseq/hints/$wd");

	SelectHCset( "complete.gtf", "incomplete.gtf", "../../gmst/genome_gmst_for_HC.gtf", "complete_ext.gtf" );

	print "### extend HC set ... done\n" if $verbose;
}
# -------------
sub CreateIntronHints
{
	print "\n### create RNA-Seq hints\n" if $verbose;

	ChDir("$workdir/rnaseq");
	MkDir("$workdir/rnaseq/hints");
	ChDir("$workdir/rnaseq/hints");

	foreach my $set ( @rnaseq_sets )
	{
		print "# working on $set\n" if $verbose;

		my $in_bam = "$workdir/rnaseq/hisat2/mapping_". $set .".bam";
		my $out_hints = "bam2hints_". $set .".gff";

		if ( CreateThis($out_hints) )
		{
			StopIfNotFound($in_bam);
			system("$bin/bam2hints --intronsonly --in=$in_bam --out=$out_hints");
		}
	}

	if ( CreateThis("bam2hints_merged.gff") )
	{
		my @files = ();
		foreach my $set (@rnaseq_sets)
		{
			push @files, "bam2hints_". $set .".gff";
		}
		
		JoinIntronsGFF( "bam2hints_merged.gff", \@files );
		StopIfNotFound("$workdir/data/genome.fasta");
		system( "$bin/filterIntronsFindStrand.pl $workdir/data/genome.fasta bam2hints_merged.gff --score > hintsfile_merged.gff" );
	}

	print "### create RNA-Seq hints ... done\n\n" if $verbose;
}
# -------------
sub JoinIntronsGFF
{
	my $out = shift;
	my $ref = shift;

	my %h = ();

	foreach my $file (@{$ref})
	{
		open( my $IN, $file ) or die "error on open file $file: $!\n";
		while(<$IN>)
		{
			my @arr = split('\t');

			my $mult = 1;
			if ( $arr[8] =~ /mult=(\d+);/ )
			{
				$mult = $1;
			}
			
			my $key = $arr[0] ."\t". $arr[3] ."\t". $arr[4];

			$h{$key} += $mult;
		}
		close $IN;
	}

	open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
	foreach my $key (keys %h)
	{
		my @arr = split('\t', $key);

		print $OUT $arr[0] ."\tb2h\tintron\t". $arr[1] ."\t". $arr[2] ."\t0\t.\t.\tmult=". $h{$key} .";pri=4;src=E\n";
	}
	close $OUT;
}
# -------------
sub AssembleTrans
{
	print "\n### assemble transcripts\n" if $verbose;

	ChDir( "$workdir/rnaseq" );
	MkDir( "$workdir/rnaseq/stringtie" );
	ChDir( "$workdir/rnaseq/stringtie" );

	my $gff_list = '';

	foreach my $set ( @rnaseq_sets )
	{
		my $in_bam = "$workdir/rnaseq/hisat2/mapping_". $set .".bam";
		my $out_gff = "transcripts_". $set .".gff";
		$gff_list .= $out_gff ." ";
		StopIfNotFound($in_bam);

		print "# working on $set\n" if $verbose;
		if ( CreateThis($out_gff) )
		{
			system("stringtie -p $cores -o $out_gff $in_bam");
		}
	}

	if ( CreateThis("transcripts_merged.gff") )
	{
		system("stringtie --merge -o transcripts_merged.gff $gff_list");
	}

	if ( CreateThis("transcripts_merged.fasta") )
	{
		system("gffread -w transcripts_merged.fasta -g $workdir/data/genome.softmasked.fasta transcripts_merged.gff");
	}

	print "### assemble transcripts ... done\n\n" if $verbose;
}
# -------------
sub MapRNASeq
{
	print "\n### map RNA-Seq\n" if $verbose;

	MkDir("$workdir/rnaseq");
	ChDir("$workdir/rnaseq");
	MkDir("$workdir/rnaseq/hisat2");
	ChDir("$workdir/rnaseq/hisat2");

	if ( !$bam and CreateThis("genome.1.ht2") )
	{
		print "# building genome index\n" if $verbose;
		system("hisat2-build --quiet -p $cores $workdir/data/genome.softmasked.fasta genome");
		print "# building genome index ... done\n" if $verbose;
	}

	my $R1;
	my $R2;
	my $R;

	foreach my $set ( @rnaseq_sets )
	{
		print "# working on $set\n" if $verbose;

		$R  = $workdir ."/rnaseq/reads/". $set .".fastq";
		$R1 = $workdir ."/rnaseq/reads/". $set ."_1.fastq";
		$R2 = $workdir ."/rnaseq/reads/". $set ."_2.fastq";

		my $out_sam = "mapping_". $set .".sam";
		my $out_bam = "mapping_". $set .".bam";
		my $out_introns = "introns_". $set .".bed";

		if ( CreateThis($out_bam) )
		{
			if ($bam)
			{
				system( "ln -s $bam/$set.bam $out_bam" );
			}
			else
			{
				if ( CreateThis($out_sam) )
				{
					print "# mapping reads from $set\n" if $verbose;

					if ( -e $R )
					{
						system("hisat2 -x genome -U $R          --dta -p $cores -S $out_sam  --novel-splicesite-outfile $out_introns");
					}
					elsif ( -e $R1 and -e $R2 )
					{
						StopIfNotFound($R1);
						StopIfNotFound($R2);
				
						system("hisat2 -x genome -1 $R1 -2 $R2  --dta -p $cores -S $out_sam  --novel-splicesite-outfile $out_introns");
					}
					else
						{ die "error, reads are not found mapping step: $set\n"; }

					print "# mapping reads from $set ... done\n" if $verbose;
				}

				print "# from sam to sorted bam for $set\n" if $verbose;
				system("samtools sort -o $out_bam -@ $cores $out_sam");
				print "# from sam to sorted bam for $set ... done\n" if $verbose;
			}
		}

		if ( $clean and -e $out_sam ) { unlink $out_sam or die "error on delete file: $out_sam $!\n"; }
	}

	# not used in algo
	if ( $debug and CreateThis("introns_merged.bed") )
	{
		my $fnames_introns = '';
		foreach my $set ( @rnaseq_sets )
		{
			my $out_introns = "introns_". $set .".bed";
			$fnames_introns .= (" ". $out_introns);
		}
		system( "cat $fnames_introns | sort | uniq  > introns_merged.bed" );
	}

	print "### map RNA-Seq ... done\n\n" if $verbose;
}
# -------------
sub PrepareRNASeq
{
	print "\n### prepare RNA-Seq\n" if $verbose;

	ChDir($workdir);
	MkDir("$workdir/rnaseq");
	ChDir("$workdir/rnaseq");
	MkDir("$workdir/rnaseq/reads");
	ChDir("$workdir/rnaseq/reads");

	foreach my $set ( @rnaseq_sets )
	{
		print "# working on RNA-Seq set: $set\n" if $verbose;

		my $sra = $set ."/". $set .".sra";
		my $R1 = $set ."_1.fastq";
		my $R2 = $set ."_2.fastq";
		my $R  = $set .".fastq";

		if ((CreateThis($R1) or CreateThis($R2)) and  CreateThis($R))
		{
			if ($local)
			{
				if ( -e "$local/$R1" and -e "$local/$R2" )
				{
					system( "ln -s $local/$R1" );
					system( "ln -s $local/$R2" );
				}
				elsif ( -e "$local/$R" )
				{
					system( "ln -s $local/$R" );
				}
				else
					{ die "error, RNA-seq file for $set is not found at $local\n"; }
			}
			else
			{
				if ( CreateThis($sra) )
				{
					system( "prefetch --max-size  35G $set" );
				}

				StopIfNotFound("$workdir/rnaseq/reads/$set/$set.sra");

				system( "fastq-dump --split-3 $workdir/rnaseq/reads/$set/$set.sra" );
			}
		}

		if ( $clean and -e $sra ) { unlink $sra or die "error on delete file $sra: $!\n"; }
	}

	print "### prepare RNA-Seq ... done\n\n" if $verbose;
}
# -------------
sub MaskGenome
{
	my $fpath = shift;

	print "\n### mask genome\n" if $verbose;

	ChDir("$workdir/arx");

	# if no repeat masking information in file
	if (!$fpath)
	{
		system("touch ../data/repeats.gff");
		$penalty = 0;
	}
	else
	{
		my $fname = FileFromPath($fpath);;

		if ( CreateThis( "../data/repeats.gff" ))
		{
			RepeatMaskerOutputToGFF( $fname, "../data/repeats.gff" );
		}

		if ( $clean_deep) { unlink $fname or die "error on delete file $fname: $!\n"; }
	}

	ChDir("$workdir/data");

	if ( CreateThis( "genome.softmasked.fasta" ))
	{
		StopIfNotFound( "genome.fasta" );
		system("bedtools maskfasta -fi genome.fasta -bed repeats.gff -fo genome.softmasked.fasta -soft") and die "error on bedtools maskfasta\n";
	}

	StopIfNotFound( "genome.softmasked.fasta" );

	print "### mask genome ... done\n\n" if $verbose;
}
# -------------
sub RepeatMaskerOutputToGFF
{
	my $in = shift;
	my $out = shift;

	die "error, input file name is missing in RepeatMaskerOutputToGFF\n" if !$in;

	open( my $IN, $in ) or die "error on open file $in: $!\n";
	open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
	while(<$IN>)
	{
		next if ! /^\s*\d+/;
	
		my @arr = split(' ');

		my $seqid = $arr[4];
		my $label = "RM";
		my $type = "repeat";
		my $start = $arr[5];
		my $end = $arr[6];
		my $score = $arr[0];
		my $strand = $arr[8];
		my $phase = '.';
		my $attr = $arr[10];

		if ($strand eq '+')
		{
			;
		}
		elsif ($strand eq 'C')
		{
			$strand = '-';
		}
		else
			{ die "error, unexpected line format in file $in: $_\n"; }

		print $OUT $seqid ."\t". $label ."\t". $type ."\t". $start ."\t". $end ."\t". $score ."\t". $strand ."\t". $phase  ."\t". $attr ."\n"
	}
	close $OUT;
	close $IN;
}
# -------------
sub PrepareAnnotation
{
	my $fpath = shift;
	my $change_id = shift;
	
	print "\n### prepare annotation\n" if $verbose;

	ChDir( $workdir );
	MkDir( "arx" );
	MkDir( "annot" );
	ChDir( "arx" );

	my $fname = FileFromPath($fpath);

	system( "$bin/gff_to_gff_subset.pl --in $fname --out tmp_in.gff3 --list chr_annot.names --swap" ) and die "error on gff_to_gff_subset.pl\n";
exit;
	system( "$bin/probuild --stat_fasta stat_fasta  --seq $workdir/data/genome.softmasked.fasta" );
	PrepareGff3Header( "stat_fasta", "in_gencode.gff3" );
	system( "cat tmp_in.gff3 | grep -v '^#' >> in_gencode.gff3" );
	system( "gt gff3 -force -tidy -sort -retainids -checkids -o nice_in_gencode.gff3 -addintrons -v in_gencode.gff3");
exit;
	system( "$bin/enrich_gff.pl --in nice_in_gencode.gff3 --out annot.gff3 --gseq_in $workdir/data/genome.softmasked.fasta --v --warnings > stat.enrich");
	system( "$bin/gff3_to_gtf.pl annot.gff3 annot.gtf" ) and die "error on gff3_to_gtf.pl\n";

	move("annot.gff3", "../annot/annot.gff3") or die "error on move file: $!\n";
	move("annot.gtf", "../annot/annot.gtf") or die "error on move file: $!\n";

	print "### prepare annotation ... done\n\n" if $verbose;
}
# -------------
sub PrepareGff3Header
{
	my $in = shift;
	my $out = shift;

	open(my $IN, $in) or die "error on open file $in: $!\n";
	open(my $OUT, ">", $out) or die "error on open file $out: $!\n";
	print $OUT "##gff-version 3\n";
	while(<$IN>)
	{
		next if /^\s*$/;
		next if /^#/;

		if ( /^>(\S+)\s+(\d+)\s+/ )
		{
			my $seqid = $1;
			my $size = $2;
			print $OUT "##sequence-region ". $seqid ." 1 ". $size ."\n";
		}
		else
			{ die "error, unexpected line format found: $_"; }
	}
	close $OUT;
	close $IN;
}
# -------------
sub PrepareGenome
{
	my $fpath = shift;

	print "\n### prepare genome\n" if $verbose;

	die "error, file path is empty in PrepareGenome\n" if !$fpath;

	ChDir($workdir);
	MkDir("$workdir/data");
	MkDir("$workdir/arx");
	ChDir("$workdir/arx");

	my $fname = FileFromPath($fpath);

	SetChrNames( $fname, "chr.names" ) if ( CreateThis( "chr.names" ));

	if ( CreateThis( "../data/genome.fasta" ))
	{
		system( "$bin/probuild --reformat_fasta --in $fname --out ../data/genome.fasta --uppercase 1 --letters_per_line 60 --include_sid chr.names --first_w --swap_sid" ) and die "error on probuild in PrepareGenome\n";
	}

	if ($softmask and CreateThis( "../data/genome.softmasked.fasta" ))
	{
		system( "$bin/probuild --reformat_fasta --in $fname --out ../data/genome.softmasked.fasta --uppercase 0 --letters_per_line 60 --include_sid chr.names --first_w --swap_sid" ) and die "error on probuild in PrepareGenome\n";;
	}

	StopIfNotFound("../data/genome.fasta");
	StopIfNotFound("../data/genome.softmasked.fasta") if $softmask;

#	if ($clean_deep) { unlink $fname or die "error on delete file $fname: $!\n"; }

	print "### prepare genome ... done\n\n" if $verbose;
}
# -------------
sub ValueDuplicationInHash
{
	my $ref = shift;

	my $found = 0;
	my %h = ();

	foreach my $key (keys %{$ref})
	{
		if ( ! exists $h{ $ref->{$key} } )
		{
			$h{ $ref->{$key} } = $key;
		}
		else
		{
			$found = "$key $ref->{$key}\n";
			last;
		}
	}

	return $found;
}
# -------------
sub SetChrNames
{
	my $in = shift;
	my $out = shift;

	my %h = ();
	my @arr = ();

	# collect seqid's from FASTA deflines
	open( my $IN, $in ) or die "error on open file $in: $!\n";
	while(<$IN>)
	{
		if ( /^>/ )
		{
			if ( $paper and /mitochondrion|chloroplast|apicoplast|^>NW_/ )
			{
				print "warning, excluded FASTA record with defline: $_" if $warnings;
				next;
			}

			if ($paper and (( /^>(N[CT]_\S+) .* chromosome\s+(\S+)/ and !/Mus musculus/ )or( /^>(NC_\S+) .* chromosome\s+(\S+)/ and /Mus musculus/ )))
			{
				my $seqid = $1;
				my $chrid = $2;

				$chrid =~ s/,$//;

				# @arr is used to keep track the order of FASTA records
				push @arr, $seqid;

				if ( ! exists $h{$seqid} )
				{
					$h{$seqid} = $chrid;
				}
				else
					{ die "error, FASTA seqid duplication was detected: $seqid $_"; }

				print "# old-new seqid: ". $seqid ." ". $chrid ."\n" if $verbose;
			}
			elsif ( /^>\s*(\S+)\s*/ )
			{
				my $seqid = $1;
				my $chrid = $seqid;

				# @arr is used to keep track the order of FASTA records
				push @arr, $seqid;

				if ( ! exists $h{$seqid} )
				{
					$h{$seqid} = $chrid;
				}
				else
					{ die "error, FASTA seqid duplication was detected: $seqid $_"; }

				print "# old-new seqid: ". $seqid ." ". $chrid ."\n" if $verbose;
			}
			else
				{ print "warning, excluded FASTA record with defline: $_" if $warnings; }
		}
	}
	close $IN;

	if ( ! scalar @arr )
		{ die "error, FASTA selection was unsuccessful: $in\n"; }

	if ( my $result = ValueDuplicationInHash(\%h) )
		{ die "error, duplication was found in new ID's: $result\n"; }

	# create output
	open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
	foreach my $id (@arr)
	{
		print $OUT ($id ."\t". $h{$id} ."\n");
	}
	close $OUT;
}
# -------------
sub FileFromPath
{
	my $fpath = shift;

	print "# working on: $fpath\n" if $verbose;

	# '.gz' is excluded from file name
	my $fname = GetFileNameFromPathName($fpath);

	if ( $force or (! -e "$fname"))
	{
		if ( -e $fpath ) # if file is on local system
		{
			if ( $force or ( ! -e "$fname.gz" and ! -e $fname ))
			{
				system( "cp", $fpath, "./" );
			}
		}
		elsif ( $fpath =~ /^https:\/\// or $fpath =~ /^http:\/\// or $fpath =~ /^ftp:\/\// )
		{
			if ( $force or ( ! -e "$fname.gz" and ! -e $fname ))
			{
				system( "wget", "-q", $fpath ) and die "error on wget\n";
			}
		}
		elsif ( $fpath =~ /^(\S+)@/ )
		{
			if ( $force or ( ! -e "$fname.gz" and ! -e $fname ))
			{
				system( "scp", $fpath, "./" );
			}
	        }
		else
			{ die "error, unsupported method for file download: $fpath\n"; }
	}

	if (( $force or (! -e "$fname")) and -e "$fname.gz" )
	{
		system( "gunzip", "$fname.gz" ) and die "error on gunzip\n";
	}

	die "error, file not found: $fname\n" if ( ! -e $fname );

	print "# working on: $fpath ... done\n" if $verbose;

	return $fname;
}
# -------------
sub Start
{
	print "### ETP running\n";
	my $date = localtime();
	print "# ". $date ."\n";
}
# -------------
sub Done
{
	my $date = localtime();
	print "# ". $date ."\n";
	print "### ETP running ... done\n";
}
# -------------
sub StopIfNotFound
{
	my $fname = shift;

	die "error, file/forder name is empty\n" if !$fname;
    die "error, file/folder not found: ". $fname ."\n" if ( ! -e $fname );
}
# -------------
sub MkDir
{
	my $fname = shift;

	if ( ! $fname )
	{
		die "error, folder name is empty\n";
	}

	if ( ! -d $fname )
	{
		mkdir $fname or die "error on mkdir: ". $fname ."\n";
	}
}
# -------------
sub ChDir
{
	my $fname = shift;

	if ( ! $fname )
	{
		die "error, folder name is empty\n";
	}

	chdir $fname or die "error on chdir to: ". $fname ."\n";
}
# -------------
sub CreateThis
{
	my $fname = shift;

	if ( ! $fname )
	{
		die "error, name is empty\n";
	}

	if ( $force or ( ! -e $fname ))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
# -------------
sub ExtendEvidence
{
	my $wd = shift;

	print "\n### extend hints\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	my $threshold = 4;

	if (!$gc)
	{
		if ( $force or ! -e "rna_conflicts.gff" )
		{
			system( "$bin/printRnaAlternatives.py pred_m/genemark.gtf r_hints_nonhc.gtf --otherIntrons prothint/prothint.gff --minIntronScore $threshold > rna_conflicts.gff" );
		}
	}
	else
	{
		if ( $force or ! -e "rna_conflicts_low.gff" )
		{
			system( "$bin/printRnaAlternatives.py pred_m_low/genemark.gtf r_hints_nonhc.gtf --otherIntrons prothint/prothint.gff --minIntronScore $threshold > rna_conflicts_low.gff" );
		}
		if ( $force or ! -e "rna_conflicts_medium.gff" )
		{
			system( "$bin/printRnaAlternatives.py pred_m_medium/genemark.gtf r_hints_nonhc.gtf --otherIntrons prothint/prothint.gff --minIntronScore $threshold > rna_conflicts_medium.gff" );
		}
		if ( $force or ! -e "rna_conflicts_high.gff" )
		{
			system( "$bin/printRnaAlternatives.py pred_m_high/genemark.gtf r_hints_nonhc.gtf --otherIntrons prothint/prothint.gff --minIntronScore $threshold > rna_conflicts_high.gff" );
		}
	}
	
	print "### extend hints ... done\n" if $verbose;
}
# -------------
sub SelectSupported
{
	my $wd = shift;

	print "\n### select supported genes\n" if $verbose;

	ChDir( "$workdir/$wd/nonhc" );

	# prothint.gff           from ProtHint and GMS-T seeds
	# hintsfile_merged.gff   from RNA-Seq mapping
	# prothint.gff           from ProtHint on nonhc regions with GM.hmm seeds
	# HC complete            from HC module
	# HC partial             from HC module

	system( "$bin/format_back.pl prothint/prothint.gff nonhc.trace > prothint_nonhc.gff" );
	system( "cat $hcc_genes $hcp_genes $prothint_hints $rnaseq_hints prothint_nonhc.gff > allHints.gff" );
	system( "$bin/selectSupportedSubsets.py ../genemark.gtf allHints.gff --fullSupport /dev/null --anySupport ../genemark_supported.gtf --noSupport /dev/null" );

	print "### select supported genes ... done\n" if $verbose;
}
# -------------
sub EstinateMaskingPenalty
{
	my $wd = shift;

	print "\n### estimate masking penalty\n" if $verbose;

	ChDir("$workdir/$wd");
	MkDir("penalty");
	ChDir("penalty");

	my $estimated_penalty = 0;

	if (!$gc)
	{
		if ( $force or ! -e "penalty.value" )
		{	
			$estimated_penalty = `$bin/estimateMaskingPenalty.py  $hcc_genes  $workdir/data/genome.softmasked.fasta  ../model/output.mod  --GMES_PATH $bin/gmes  --threads $cores`;

			open( my $OUT, ">", "penalty.value" ) or die "error on open file penalty.value: $!\n";
			print $OUT "$estimated_penalty\n";
			close $OUT;
		}
	}
	else
	{
		if ( $force or ! -e "penalty.value.low" )
		{
			$estimated_penalty = `$bin/estimateMaskingPenalty.py  ../low/hc.gff  $workdir/data/genome.softmasked.fasta  ../low/output.mod  --GMES_PATH $bin/gmes  --threads $cores`;

			open( my $OUT, ">", "penalty.value.low" ) or die "error on open file penalty.value: $!\n";
			print $OUT "$estimated_penalty\n";
			close $OUT;
		}

		if ( $force or ! -e "penalty.value.medium" )
		{
			# medium GC
			$estimated_penalty = `$bin/estimateMaskingPenalty.py  ../medium/hc.gff  $workdir/data/genome.softmasked.fasta  ../medium/output.mod  --GMES_PATH $bin/gmes  --threads $cores`;

			open( my $OUT, ">", "penalty.value.medium" ) or die "error on open file penalty.value: $!\n";
			print $OUT "$estimated_penalty\n";
			close $OUT;
		}

		if ( $force or ! -e "penalty.value.high" )
		{
			$estimated_penalty = `$bin/estimateMaskingPenalty.py  ../high/hc.gff  $workdir/data/genome.softmasked.fasta  ../high/output.mod  --GMES_PATH $bin/gmes  --threads $cores`;

			open( my $OUT, ">", "penalty.value.high" ) or die "error on open file penalty.value: $!\n";
			print $OUT "$estimated_penalty\n";
			close $OUT;
		}
	}

	print "### estimate masking penalty ... done\n\n" if $verbose;
}
# -------------
sub isGChetero
{
	my $wd = shift;

	print "\n### check GC status\n" if $verbose;

	my $status = 0;

	ChDir("$workdir/$wd/hc");

	system( "$bin/probuild --stat_fasta stat.fasta --seq hc.fasta" );
	my $str = `$bin/estimate_bins.pl ../hc_regions.gtf  hc.trace  stat.fasta`;
	my $L = 0;
	my $R = 0;
	if ( $str =~ /(\d+)\s+(\d+)/ )
	{
		$L = $1;
		$R = $2;
	}
	print "# GC: $L $R\n" if $verbose;

	print "### check GC status ... done\n" if $verbose;

	return $status;
}
# -------------
sub Combine
{
	my $wd = shift;

	print "\n### combine hc and nonhc\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	if (!$gc)
	{
		system("$bin/format_back.pl  pred_m/genemark.gtf  nonhc.trace  > nonhc.gtf");
	}
	else
	{
		system("$bin/format_back.pl pred_m_low/genemark.gtf     nonhc.trace | sed 's/_g/_gL/' | sed 's/_t/_tL/'  > low.gtf");
		system("$bin/format_back.pl pred_m_medium/genemark.gtf  nonhc.trace | sed 's/_g/_gM/' | sed 's/_t/_tM/'  > medium.gtf");
		system("$bin/format_back.pl pred_m_high/genemark.gtf    nonhc.trace | sed 's/_g/_gH/' | sed 's/_t/_tH/'  > high.gtf");

		system("cat low.gtf medium.gtf high.gtf > nonhc.gtf");
	}

	ChDir("$workdir/$wd");

	system( "cat $hcc_genes nonhc/nonhc.gtf > genemark.gtf" );

	print "### combine hc and nonhc ... done\n\n" if $verbose;
}
# -------------
sub PredictNonhc
{
	my $wd = shift;

	my $mod;
	my $mod_low;
	my $mod_medium;
	my $mod_high;

	if (!$train_cfg{"iteration"})
	{
		$mod        = $train_cfg{"model_0"};
		$mod_low    = $train_cfg{"model_0_low"};
		$mod_medium = $train_cfg{"model_0_medium"};
		$mod_high   = $train_cfg{"model_0_high"};
	}
	else
	{
		$mod        = $train_cfg{"model_etr"};
		$mod_low    = $train_cfg{"model_etr_low"};
		$mod_medium = $train_cfg{"model_etr_medium"};
		$mod_high   = $train_cfg{"model_etr_high"};
	}

	print "\n### predict genes nonhc\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	# move to prothint section
	if($gc)
	{
		MkDir("prothint");
		system( "cat prothint_low/prothint.gff prothint_medium/prothint.gff prothint_high/prothint.gff > prothint/prothint.gff" );
		system( "cat prothint_low/evidence.gff prothint_medium/evidence.gff prothint_high/evidence.gff > prothint/evidence.gff" );
	}

	# introns - overlap rnaseq and prothint
	system( "$bin/compare_intervals_exact.pl --f1 prothint/prothint.gff --f2 r_hints_nonhc.gtf --out evi_all.gff --shared12 --original 1 --intron --no" );

	system( "grep -P \"\tstart_codon\t\" prothint/evidence.gff >> evi_all.gff" );
	system( "grep -P \"\tstop_codon\t\"  prothint/evidence.gff >> evi_all.gff" );

	system( "grep -P \"\tstart_codon\t\" p_evi_nonhc.gtf >> evi_all.gff" );
	system( "grep -P \"\tstop_codon\t\"  p_evi_nonhc.gtf >> evi_all.gff" );

	system( "$bin/filter_hints_partial.pl --hc hcp_hints_nonhc.gtf --hints evi_all.gff  --out evi.gff" );
	system( "cat hcp_hints_nonhc.gtf >> evi.gff" );

	if ( -e "rna_conflicts.gff" )
	{
		system( "cat rna_conflicts.gff >> evi.gff" );
	}

	if (!$gc)
	{
		MkDir("pred_m");
		ChDir("pred_m");

		if ( $force or ! -e "genemark.gtf" or $train_cfg{"extend_all"} )
		{
			print "#ITERATION $mod\n" if ( $verbose and $train_cfg{"extend_all"} );

			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty --predict_with $mod --seq ../nonhc.fasta --cores $cores --verbose --evi ../evi.gff --max_gap 40000 --max_mask 40000 > loginfo");
		}
	}
	else
	{
		MkDir("pred_m_low");
		MkDir("pred_m_medium");
		MkDir("pred_m_high");

		ChDir("$workdir/$wd/nonhc/pred_m_low");
		if ( $force or ! -e "genemark.gtf" or $train_cfg{"extend_low"} )
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_low --predict_with $mod_low --seq ../low.fasta --cores $cores --verbose --evi ../evi.gff  --max_gap 40000 --max_mask 40000 > loginfo");
		}

		ChDir("$workdir/$wd/nonhc/pred_m_medium");
		if ( $force or ! -e "genemark.gtf" or $train_cfg{"extend_medium"} )
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_medium --predict_with $mod_medium --seq ../medium.fasta --cores $cores --verbose --evi ../evi.gff  --max_gap 40000 --max_mask 40000 > loginfo");
		}

		ChDir("$workdir/$wd/nonhc/pred_m_high");
		if ( $force or ! -e "genemark.gtf" or $train_cfg{"extend_high"})
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_high --predict_with $mod_high --seq ../high.fasta --cores $cores --verbose --evi ../evi.gff  --max_gap 40000 --max_mask 40000 > loginfo");
		}
	}

	print "### predict genes nonhc ... done\n\n" if $verbose;
}
# -------------
sub RunProtHint
{
	my $wd = shift;
	my $prot = shift;

	print "\n### ProtHint\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	if (!$gc)
	{
		MkDir("prothint");
		ChDir("prothint");

		if ( $force or ! -e "evidence.gff" )
		{
			system("$bin/gmes/ProtHint/bin/prothint.py  --geneMarkGtf ../for_prothint/genemark.gtf ../nonhc.fasta  $workdir/data/$prot 2> log");
		}
	}
	else
	{
		MkDir("prothint_low");
		MkDir("prothint_medium");
		MkDir("prothint_high");

		ChDir("$workdir/$wd/nonhc/prothint_low");
		if ( $force or ! -e "evidence.gff" )
		{
			system("$bin/gmes/ProtHint/bin/prothint.py  --geneMarkGtf ../for_prothint_low/genemark.gtf ../low.fasta  $workdir/data/$prot 2> log");
		}

		ChDir("$workdir/$wd/nonhc/prothint_medium");
		if ( $force or ! -e "evidence.gff" )
		{
			system("$bin/gmes/ProtHint/bin/prothint.py  --geneMarkGtf ../for_prothint_medium/genemark.gtf ../medium.fasta  $workdir/data/$prot 2> log");
		}

		ChDir("$workdir/$wd/nonhc/prothint_high");
		if ( $force or ! -e "evidence.gff" )
		{
			system("$bin/gmes/ProtHint/bin/prothint.py  --geneMarkGtf ../for_prothint_high/genemark.gtf ../high.fasta  $workdir/data/$prot 2> log");

		}
	}

	print "### ProtHint ... done\n\n" if $verbose;
}
# -------------
sub PredictForProtHint
{
	my $wd = shift;

	print "\n### predict genes for ProtHint\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	my $mod;
	my $mod_low;
	my $mod_medium;
	my $mod_high;

#	if (!$train_cfg{"iteration"})
#	{
                $mod        = $train_cfg{"model_0"};
                $mod_low    = $train_cfg{"model_0_low"};
                $mod_medium = $train_cfg{"model_0_medium"};
                $mod_high   = $train_cfg{"model_0_high"};
#	}
#	else
#	{
#                $mod        = $train_cfg{"model_etr"};
#                $mod_low    = $train_cfg{"model_etr_low"};
#                $mod_medium = $train_cfg{"model_etr_medium"};
#                $mod_high   = $train_cfg{"model_etr_high"};
#	}

	# introns - overlap rnaseq and prothint
	system( "$bin/compare_intervals_exact.pl --f1 p_hints_nonhc.gtf --f2 r_hints_nonhc.gtf --out evi_ini.gff --shared12 --original 1 --intron --no" );
	system( "cat  p_evi_nonhc.gtf >> evi_ini.gff" );
	system( "cat  hcp_hints_nonhc.gtf >> evi_ini.gff" );

	if (!$gc)
	{
		MkDir("for_prothint");
		ChDir("for_prothint");

		if ( $force or ! -e "genemark.gtf" )
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty --predict_with $mod --seq ../nonhc.fasta --cores $cores --verbose --evi ../evi_ini.gff  --max_gap 40000 --max_mask 40000 > loginfo");

		}
	}
	else
	{
		MkDir("for_prothint_low");
		MkDir("for_prothint_medium");
		MkDir("for_prothint_high");

		ChDir("$workdir/$wd/nonhc/for_prothint_low");
		if ( $force or ! -e "genemark.gtf" )
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_low --predict_with $mod_low --seq ../low.fasta --cores $cores --verbose --evi ../evi_ini.gff --max_gap 40000 --max_mask 40000 > loginfo");
		}

		ChDir("$workdir/$wd/nonhc/for_prothint_medium");
		if ( $force or ! -e "genemark.gtf" )
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_medium --predict_with $mod_medium --seq ../medium.fasta --cores $cores --verbose --evi ../evi_ini.gff --max_gap 40000 --max_mask 40000 > loginfo");
		}

		ChDir("$workdir/$wd/nonhc/for_prothint_high");
		if ( $force or ! -e "genemark.gtf" )
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_high --predict_with $mod_high --seq ../high.fasta --cores $cores --verbose --evi ../evi_ini.gff --max_gap 40000 --max_mask 40000 > loginfo");
		}
	}

	print "### predict genes for ProtHint ... done\n\n" if $verbose;
}
# -------------
sub TranscriptIdFromRegions
{
	my $in = shift;
	my $out = shift;

	open( my $IN, $in ) or die "error on open file $in: $!\n";
	open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
	while(<$IN>)
	{
		my @arr = split('\t');

		# incomplete HC genes - one boundary 
		next if ( $arr[3] == $arr[4]);

		if( $arr[8] =~ /transcript_id \"(\S+)\";/ )
		{
			my $tid = $1;
			print $OUT $tid ."\n";
		}
		else { die "error in format\n"; }
	}
	close $OUT;
	close $IN;
}
# -------------
sub TrainModel
{
	my $wd = shift;

	print "\n### training abinitio parameters\n" if $verbose;

	ChDir("$workdir/$wd");
	MkDir("$workdir/$wd/model");
	ChDir("$workdir/$wd/model");

	if ( $force or ( ! -e "output.mod" ))
	{
		TranscriptIdFromRegions( "../hc_regions.gtf", "training.list" );
		system( "$bin/select_by_transcript_id_from_gtf.pl training.list $hcc_genes training.gtf" );

		if(!$bp)
		{
			system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.softmasked.fasta" );
		}
		else
		{
			system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.softmasked.fasta --bp" );
		}
	}

	print "### training abinitio parameters ... done\n\n" if $verbose;
}
# -------------
sub TrainModelGC
{
	my $wd = shift;

	print "\n### training GC abinitio parametwers\n" if $verbose;

	ChDir("$workdir/$wd/hc");

	system("$bin/probuild --stat_fasta stat.fasta --seq hc.fasta");
	my $str = `$bin/estimate_bins.pl ../hc_regions.gtf  hc.trace  stat.fasta`;
	my $L = 0;
	my $R = 0;
	if ( $str =~ /(\d+)\s+(\d+)/ )
	{
		$L = $1;
		$R = $2;
	}
	print "# GC: $L $R\n";

	# split nonhc regions by GC
	ChDir("$workdir/$wd/nonhc");

	system("$bin/probuild --stat_fasta stat.fasta --seq nonhc.fasta");
	system("$bin/split_by_bins.pl stat.fasta $L $R");

	JoinFasta( "regions", "low.fasta",    "low.ids" );
	JoinFasta( "regions", "medium.fasta", "medium.ids" );
	JoinFasta( "regions", "high.fasta",   "high.ids" );

	# train by GC
	ChDir("$workdir/$wd");
	MkDir("low");
	MkDir("medium");
	MkDir("high");

	ChDir("$workdir/$wd/low");
	if ( $force or ! -e "output.mod" )
	{
		system( "$bin/select_by_gene_id_from_gtf.pl ../hc/low.ids    $hcc_genes genes.gtf" );
		TranscriptIdFromRegions( "../hc_regions.gtf", "training.list" );
		system( "$bin/select_by_transcript_id_from_gtf.pl training.list genes.gtf training.gtf" );
		system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.softmasked.fasta" );
	}

	ChDir("$workdir/$wd/medium");
	if ( $force or ! -e "output.mod" )
	{
		system( "$bin/select_by_gene_id_from_gtf.pl ../hc/medium.ids $hcc_genes genes.gtf" );
		TranscriptIdFromRegions( "../hc_regions.gtf", "training.list" );
		system( "$bin/select_by_transcript_id_from_gtf.pl training.list genes.gtf training.gtf" );
		system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.softmasked.fasta" );
	}

	ChDir("$workdir/$wd/high");
	if ( $force or ! -e "output.mod" )
	{
		system( "$bin/select_by_gene_id_from_gtf.pl ../hc/high.ids   $hcc_genes genes.gtf" );
		TranscriptIdFromRegions( "../hc_regions.gtf", "training.list" );
		system( "$bin/select_by_transcript_id_from_gtf.pl training.list genes.gtf training.gtf" );
		system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.softmasked.fasta" );
	}

	print "### training GC abinitio parametwers ... done\n\n" if $verbose;
}
# -------------
sub PrepareGenomeTraining
{
	my $wd = shift;

	# from complete and partial high confidence genes (GFF) and sequence (FASTA)

	my $margin = 50; # for non-hc regions

	print "\n### prepare genome sequence for training\n" if $verbose;

	ChDir("$workdir/$wd");

	if ( $force or ( ! -e "hc_regions.gtf" ))
	{
		system( "$bin/create_regions.pl --hcc $hcc_genes --hcp $hcp_genes --out hc_regions.gtf --margin $margin" ) and die "error on create_regions.pl";
	}

	MkDir("$workdir/$wd/nonhc");
	ChDir("$workdir/$wd/nonhc");
	MkDir("$workdir/$wd/nonhc/regions");
	ChDir("$workdir/$wd/nonhc/regions");

	if ( $force or ( ! -e "../nonhc.fasta" ))
	{
		system( "$bin/probuild --cut nonhc --seq $workdir/data/genome.softmasked.fasta --regions ../../hc_regions.gtf --reverse --allow_x --trace ../nonhc.trace --soft --min_contig 120" ) and die "error on probuild in PrepareGenomeTraining";

		JoinNonhs( ".", "../nonhc.fasta" );
	}

	ChDir("$workdir/$wd/nonhc");

	if ( $force or ( ! -e "hcp_hints_nonhc.gtf" ))
	{
		# create hints from incomplete genes
		# remove start codons - when exist
		system("$bin/gmes/rescale_gff.pl --out hcp_hints_nonhc.gtf --in $hcp_genes      --trace nonhc.trace");
		system("grep -v start_codon hcp_hints_nonhc.gtf > hcp.tmp ; mv hcp.tmp hcp_hints_nonhc.gtf");
	}

	if ( $force or ( ! -e "r_hints_nonhc.gtf" ))
	{
		system("$bin/gmes/rescale_gff.pl --out r_hints_nonhc.gtf   --in $rnaseq_hints   --trace nonhc.trace");
	}

	if ( $force or ( ! -e "p_evi_nonhc.gtf" ))
	{
		system("$bin/gmes/rescale_gff.pl --out p_evi_nonhc.gtf     --in $prothint_evi   --trace nonhc.trace");
	}

	if ( $force or ( ! -e "p_hints_nonhc.gtf" ))
	{
		system("$bin/gmes/rescale_gff.pl --out p_hints_nonhc.gtf   --in $prothint_hints --trace nonhc.trace");
	}

	print "### prepare genome sequence for training ... done\n\n" if $verbose;
}
# -------------
sub PrepareHCregions
{
	my $wd = shift;

	print "\n### prepare HC regions training\n" if $verbose;

	ChDir("$workdir/$wd");
	MkDir("hc");
	ChDir("hc");
	MkDir("regions");
	ChDir("regions");

	if ( $force or ! -e "../hc.fasta" )
	{
		system( "$bin/probuild --cut hc --seq $workdir/data/genome.softmasked.fasta --regions ../../hc_regions.gtf  --allow_x --trace ../hc.trace --soft --min_contig 120" ) and die "error on parce hc regions";
		JoinFasta( ".", "../hc.fasta" );
	}

	print "### prepare HC regions training ... \n\n" if $verbose;
}
# -------------
sub JoinNonhs
{
	my $path = shift;
	my $out = shift;

	opendir( my $DIR, "$path" ) or  die "error on open dir $path: $!\n";
	my @list = readdir($DIR);
	closedir($DIR);

	open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
	foreach my $f (@list)
	{
		if ( $f =~ /^nonhc_(\d+)$/ )
		{
			open( my $IN, $f ) or die "error on open file $f: $!\n";
			while(<$IN>)
			{
				if ( /^>(\S+)/ )
				{
					my $id = $1;
					print $OUT ">$id\n";
				}
				else
				{
					print $OUT $_;
				}
			}
			close $IN;
		}
	}
	close $OUT;
}
# -------------
sub JoinFasta
{
	my $path = shift;
	my $out = shift;
	my $use_only_these = shift;

	my %h;

	if (defined $use_only_these)
	{
		open( my $IN, $use_only_these ) or die "error on open file $use_only_these: $!\n";
		while(<$IN>)
		{
			next if /^#/;
			next if /^\s*$/;

			if ( /^(\S+)/ )
			{
				$h{$1} = 1;
			}
		}
		close $IN;
	}

	opendir( my $DIR, "$path" ) or  die "error on open dir $path: $!\n";
	my @list = readdir($DIR);
	closedir($DIR);

	open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
	foreach my $f (@list)
	{
		if (( $f =~ /^nonhc_(\d+)$/ )or( $f =~ /^hc_(\d+)$/ ))
		{
			next if ( defined $use_only_these and !exists $h{$f} );
			
			open( my $IN, "$path/$f" ) or die "error on open file $path/$f: $!\n";
			while(<$IN>)
			{
				if ( /^>(\S+)/ )
				{
					my $id = $1;
					print $OUT ">$id\n";
				}
				else
				{
					print $OUT $_;
				}
			}
			close $IN;
		}
	}
	close $OUT;
}
# -------------

