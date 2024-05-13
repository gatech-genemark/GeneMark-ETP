#!/usr/bin/env perl
# -------------
# Alex Lomsadze
# 2023
# 
# "GeneMark-ETP: Automatic Gene Finding in Eukaryotic Genomes in Consistence with Extrinsic Data"
# Tomas Bruna, Alexandre Lomsadze and Mark Borodovsky
# Georgia Institute of Technology, Atlanta, Georgia, USA
#
# Publication bioRxiv
#    https://www.biorxiv.org/content/10.1101/2023.01.13.524024v1
#
# Questions?
#    Contact Alex Lomsadze alexl@gatech.edu
#
# Copyright
#    Georgia Institute of Technology, Atlanta, Georgia, USA
# -------------

use strict;
use warnings;

use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use Cwd qw(abs_path);
use YAML::XS qw(LoadFile);
use Parallel::ForkManager;
use Data::Dumper;

my $version = "1.02";

# -------------
my $bin = $Bin;          # folder with ETP code

my $cfg = '';            # required input configuration file with links to various types of input data
my $workdir = '.';       # working folder for whole ETP run

my $proc = '';           # name of working folders for protein DB related steps
                         # some of the ETP testing was done using the same RNA-Seq input and different protein DB inputs
                         # to accommodate such testing mode, the shared RNA-Seq part was reused,
                         # while all protein DB related steps were executed in "$proc" folders
                         # by default, "$proc" is initialized from protein DB name "$protdb_name"
my $protdb_name = '';    # is initialized from $protdb_path

my $gc = 0;              # "1" - run ETP in GC heterogeneous mode; "0" - homogeneous mode
my $bp = 0;              # "1" - run ETP in fungi mode; "0" - default general eukaryotic mode

my $extend = 1;          # "1" - run ETP with extended training (iterations); "0" - no iterations
my $paper = 0;           # "1" - parse input data following rules described in the ETP publication

my $penalty;             # if the variable is "undefined", then run ETP in masking penalty estimation mode

my $local = '';          # input data processing modifier
my $bam = '';            # input data processing modifier
my $softmask = '';       # input data processing modifier

my $cores = 64;
my $force = 0;
my $clean = 0;

my $verbose = 0;
my $warnings = 0;
my $debug = 0;

my $SAMTOOLS_M = 0;     # pass this to samtools
# -------------
# input data information from "$cfg" YAML file
my $species = '';
my $genome_path = '';
my $RepeatMasker_path = '';
my @rnaseq_sets = ();
my $protdb_path = '';
my $annot_path = '';     # this parameter is not used in this release
# -------------
# penalty values for individual GC bins in GC heterogeneous mode
my $penalty_low;
my $penalty_medium;
my $penalty_high;
# -------------
# configuration for ab initio model training
my %train_cfg = ();      # training configuration; this structure is updated during ETP run
InitTrainCfg();          # this part of the code should be improved
# =============

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckInstallation() if $debug;
LoadConfig($cfg);
ReportConfig($cfg) if $debug;
Start() if $verbose;

PrepareGenome($genome_path) if 1;
if (!$softmask)
{
	MaskGenome($RepeatMasker_path) if 1;
}
PrepareProteinDBFile($protdb_path) if 1;
# must run
($proc, $protdb_name) = SetTrainingPredictionFolder( $protdb_path, $proc );

if (!$bam)
{
	PrepareRNASeq() if 1;
}
MapRNASeq()         if 1;
AssembleTrans()     if 1;
CreateIntronHints() if 1;
GeneMarkST()        if 1;
RunProtHintOnGMST($proc,$protdb_name) if 1;
FilterGMST( $protdb_name, $proc ) if 1;

my $hcc_genes      = "$workdir/rnaseq/hints/$proc/complete.gtf";
my $hcp_genes      = "$workdir/rnaseq/hints/$proc/incomplete.gtf";
my $rnaseq_hints   = "$workdir/rnaseq/hints/hintsfile_merged.gff";
my $prothint_evi   = "$workdir/rnaseq/hints/$proc/prothint/evidence.gff";
my $prothint_hints = "$workdir/rnaseq/hints/$proc/prothint/prothint.gff";

PrepareGenomeTraining($proc) if 1;
PrepareHCregions($proc)      if 1;

# isGChetero($protdb_name) # which gc type to run
if ($gc)
{
	TrainModelGC($proc) if 1;
}
else
{
	TrainModel($proc)   if 1;
}

if (!defined $penalty)
{
	# minimum number of masked letter to run repeat penalty estimation
	my $MIN_FRACTION_MASKED = 0.01;
	$penalty = CheckMasked($MIN_FRACTION_MASKED) if 1;
}
if ( ! defined $penalty )
{
	EstinateMaskingPenalty($proc) if 1;
}

LoadPenalty($proc); # must run
PredictForProtHint($proc)        if 1;
RunProtHint($proc, $protdb_name) if 1;

unlink "$workdir/$proc/nonhc/rna_conflicts.gff"        if (-e "$workdir/$proc/nonhc/rna_conflicts.gff");
unlink "$workdir/$proc/nonhc/rna_conflicts_low.gff"    if (-e "$workdir/$proc/nonhc/rna_conflicts_low.gff");
unlink "$workdir/$proc/nonhc/rna_conflicts_medium.gff" if (-e "$workdir/$proc/nonhc/rna_conflicts_medium.gff");
unlink "$workdir/$proc/nonhc/rna_conflicts_high.gff"   if (-e "$workdir/$proc/nonhc/rna_conflicts_high.gff");

PredictNonhc($proc)    if 1;
ExtendEvidence($proc)  if 1;

ITERATION:

PredictNonhc($proc)    if 1;
Combine($proc)         if 1;
SelectSupported($proc) if 1;

EstimateAccPaper() if $paper;

if ($extend)
{
	RenameCurrentOutput( $train_cfg{"iteration"} );

	if( ! StopIterations($train_cfg{"iteration"}) )
	{
		if ( $train_cfg{"iteration"} < $train_cfg{"max_iterations"} )
		{
			$train_cfg{"iteration"} += 1;

			print "# iteration ". $train_cfg{"iteration"} ."\n" if $verbose;

			ChDir($workdir);
			UpdateModel($proc, $train_cfg{"iteration"}, "$workdir/$proc/genemark_supported.gtf");

			goto ITERATION;
		}
	}
}

Output($proc);
Done() if $verbose;

# =============
sub Output
{
	my $wd = shift;

	ChDir($workdir);

	if ( -e "$wd/genemark_supported.gtf" )
	{
		system( "cp $wd/genemark_supported.gtf genemark_supported.gtf" );
	}
	else
		{ die "error, output file $wd/genemark_supported.gtf not found\n"; }

	if ( -e "$wd/genemark.gtf" )
	{
		system( "cp $wd/genemark.gtf genemark.gtf" );
	}
	else
		{ die "error, output file $wd/genemark.gtf not found\n"; }
}
# -------------
sub RenameCurrentOutput
{
	my $id = shift;

	ChDir($workdir);

	my $old_name = "$proc/genemark.gtf";
	my $new_name = "$proc/genemark.gtf_". $id;
	system( "cp $old_name $new_name" );

	$old_name = "$proc/genemark_supported.gtf";
	$new_name = "$proc/genemark_supported.gtf_". $id;
	system( "cp $old_name $new_name" );
}
# -------------
sub InitTrainCfg
{
	$train_cfg{"run"} = 1;                   # run prediction step
	$train_cfg{"max_iterations"} = 3;        # maximum number of iterations; 0 - initial prediction
	$train_cfg{"iteration"} = 0;             # current iteration
	$train_cfg{"model"} = '';
	$train_cfg{"model_low"} = '';
	$train_cfg{"model_medium"} = '';
	$train_cfg{"model_high"} = '';
}
# -------------
sub StopIterations
{
	my $id = shift;

	return 0 if $id == 0;

	ChDir("$workdir/$proc");

	my $curr = "genemark_supported.gtf_". $id;
	my $prev = "genemark_supported.gtf_". ($id - 1);

	my $fname = "cds_". ($id - 1) ."_". $id;

	system( "$bin/compare_intervals_exact.pl --f1 $prev --f2 $curr > $fname" );

	my $f = 0;
	my $s = 0;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^\s*$/;

		if ( /^\d+\s+\d+\s+\d+\s+(\S+)/ )
		{
			my $value = $1;

			$f = $value if !$f;
			$s = $value if !$s;
		}
	}
	close $IN;

	my $F1 = 2*$f*$s/($f+$s);

	print "# $F1 for $id\n" if $verbose;

	return 1 if ( $F1 > 99 );
	return 0;
}
# -------------
sub EstimateAccPaper
{
	ChDir($workdir);

	print "### estimate accuracy\n" if $verbose;

	if ( -e "annot/reliable.gtf" and -e "annot/union.gtf" and -e "annot/pseudo.gff3" )
	{
		system( "$bin/compare_intervals_exact.pl --f2 annot/reliable.gtf --f1 $proc/genemark_supported.gtf --pseudo annot/pseudo.gff3 --gene > $proc/sn_". $train_cfg{"iteration"} );
		system( "$bin/compare_intervals_exact.pl --f1 annot/union.gtf    --f2 $proc/genemark_supported.gtf --pseudo annot/pseudo.gff3 --gene > $proc/sp_". $train_cfg{"iteration"} );
	}
	elsif (  -e "annot/annot.gtf" and -e "annot/pseudo.gff3" )
	{
		system( "$bin/compare_intervals_exact.pl --f1 annot/annot.gtf    --f2 $proc/genemark.gtf           --pseudo annot/pseudo.gff3 --gene > $proc/snsp_". $train_cfg{"iteration"} );
	}

	print "### estimate accuracy ... done\n" if $verbose;
}
# -------------
sub UpdateModel
{
	my $wd = shift;
	my $itr = shift;
	my $prev = shift;

	print "### update model: $itr\n" if $verbose;

	ChDir("$workdir/$wd");
	MkDir("$workdir/$wd/etr");
	ChDir("$workdir/$wd/etr");
	MkDir("$workdir/$wd/etr/$itr");
	ChDir("$workdir/$wd/etr/$itr");

	if (!$gc)
	{
		MkDir("$workdir/$wd/etr/$itr/model");
		ChDir("$workdir/$wd/etr/$itr/model");
		if ( CreateThis("output.mod"))
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
		MkDir("$workdir/$wd/etr/$itr/low");
		ChDir("$workdir/$wd/etr/$itr/low");
		if ( CreateThis("output.mod"))
		{
			system( "$bin/select_for_training.pl --in $prev --out training.gtf --verbose --gc L --hc $workdir/$wd/hc/low.ids" );
			system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
		}

		MkDir("$workdir/$wd/etr/$itr/medium");
		ChDir("$workdir/$wd/etr/$itr/medium");
		if ( CreateThis("output.mod"))
		{
			system( "$bin/select_for_training.pl --in $prev --out training.gtf --verbose --gc M --hc $workdir/$wd/hc/medium.ids" );
			system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
		}

		MkDir("$workdir/$wd/etr/$itr/high");
		ChDir("$workdir/$wd/etr/$itr/high");
		if ( CreateThis("output.mod"))
		{
			system( "$bin/select_for_training.pl --in $prev --out training.gtf --verbose --gc H --hc $workdir/$wd/hc/high.ids" );
			system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
		}
	}

	$train_cfg{"model"}        = "$workdir/$proc/etr/". $train_cfg{"iteration"} ."/model/output.mod";
	$train_cfg{"model_low"}    = "$workdir/$proc/etr/". $train_cfg{"iteration"} ."/low/output.mod";
	$train_cfg{"model_medium"} = "$workdir/$proc/etr/". $train_cfg{"iteration"} ."/medium/output.mod";
	$train_cfg{"model_high"}   = "$workdir/$proc/etr/". $train_cfg{"iteration"} ."/high/output.mod";

	print "### update model $itr ... done\n" if $verbose;
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
sub SetTrainingPredictionFolder
{
	my $fpath = shift;
	my $folder = shift;

	print "### set training-prediction folder\n" if $verbose;

	die "error, file path to protein db is empty\n" if !$fpath;

	my $fname = GetFileNameFromPathName($fpath);

	$folder = $fname if !$folder;

	ChDir($workdir);
	MkDir("$workdir/$folder");

	print "# training-prediction folder: $folder\n" if $verbose;
	print "### set training-prediction folder ... done\n\n" if $verbose;

	return ($folder, $fname);
}
# -------------
sub LoadPenalty
{
	my $wd = shift;

	print "### load penalty values\n" if $verbose;

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

	print "### load penalty values ... done\n\n" if $verbose;
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

	print "# penalty value loaded from $fname: $value\n" if $verbose;

	return $value;
}
# -------------
sub RunProtHintOnGMST
{
	my $name = shift;
	my $pfile = shift;

	print "### generate ProtHint predictions with GeneMarkS-T seeds\n" if $verbose;

	ChDir("$workdir/rnaseq/hints");
	MkDir("$workdir/rnaseq/hints/$name");
	ChDir("$workdir/rnaseq/hints/$name");
	MkDir("$workdir/rnaseq/hints/$name/prothint");
	ChDir("$workdir/rnaseq/hints/$name/prothint");

	StopIfNotFound("$workdir/data/$pfile");

	if ( CreateThis("evidence.gff") or CreateThis("prothint.gff") )
	{
		system( "$bin/gmes/ProtHint/bin/prothint.py --geneSeeds $workdir/rnaseq/gmst/genome_gmst.gtf $workdir/data/genome.fasta $workdir/data/$pfile --longGene 30000000 --longProtein 15000000 2> $workdir/prothint_gmst.log" );
		StopIfNotFound("evidence.gff");
		StopIfNotFound("prothint.gff");
	}
	else
	{
		print "# reusing ProtHint results\n" if $verbose;
	}

	print "### generate ProtHint predictions with GeneMarkS-T seeds ... done\n\n" if $verbose;
}
# -------------
sub LoadConfig
{
	my $fname = shift;

	die "error, file name is empty in LoadConfig\n" if !$fname;

	my $config = LoadFile($fname);

	$species           = $config->{"species"}            if $config->{"species"};
	$genome_path       = $config->{"genome_path"}        if $config->{"genome_path"};
	$RepeatMasker_path = $config->{"RepeatMasker_path"}  if $config->{"RepeatMasker_path"};
	@rnaseq_sets       = @{$config->{"rnaseq_sets"}}     if $config->{"rnaseq_sets"};
	$protdb_path       = $config->{"protdb_path"}        if $config->{"protdb_path"};
	$annot_path        = $config->{"annot_path"}         if $config->{"annot_path"};
}
# -------------
sub ReportConfig
{
	my $fname = shift;

	print "### Report config from YAML species configuration file and command line\n";
	print "# cfg: $fname\n";
	print "# workdir: ". $workdir ."\n";
	print "# proc: ". $proc ."\n";
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
	print "# paper: ". $paper ."\n";
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
	print "### Report config from YAML species configuration file and command line ... done\n\n";
}
# -------------
sub Usage
{
	print "\n";
	print "ETP version $version\n\n";
	print "Usage: $0 --cfg [$cfg]\n\n";
	print "  --cfg [$cfg]   configuration file in YAML format with links to input data\n";
	print "             see instructions about file format in the README file\n";
	print "             required parameter\n";
	print "  --workdir [$workdir]   working folder; default is current folder\n";
	print "                  output and large number of temporary files and folders will be created in working folder\n";
	print "\nAlgorithm modifiers:\n\n";
	print "  --gc [$gc]   run ETP in GC homogeneous [0] or heterogeneous [1] mode; default is GC homogeneous mode [0]\n";
	print "  --fungus   run ETP in fungus mode; default is general eukaryotic mode\n";
	print "  --extend [$extend]   run ETP in extended training mode [1] or just one training round [0]; default is extended training [1]\n";
	print "  --penalty []   masking penalty; default mode is auto estimation of masking parameter\n";
	print "\nAlgorithm input modifiers:\n\n";
	print "  --softmask   input genome is in repeat softmasked FASTA format\n";
	print "               default repeat setting in --cfg file is ignored and softmask genome information is used instead\n";
	print "  --local [$local]   folder where RNA-Seq reads in FASTQ format are located; by default RNA-Seq files are loaded from SRA\n";
	print "               RNA-Seq ID's are taken from --cfg file\n";
	print "  --bam [$bam]   folder where reads-to-genome alignments in BAM format are located; by default BAM files are created by ETP\n";
	print "             BAM files must be sorted; BAM file must be compatible with StringTie2 software; BAM file ID's are taken from --cfg file\n";
	print "  --paper [$paper]   select subset of input FASTA records as it is described in ETP publication [1]\n";
	print "                by default [0] all the input FASTA records are used in ETP\n";
	print "\nRun:\n\n";
	print "  --cores [$cores]  number of CPU/cores to use on multi CPU/core computer\n";
	print "  --clean       delete some temporary files during ETP run\n";
	print "  --verbose\n";
	print "  --warnings\n";
	print "  --debug\n";
	print "  --proc [$proc]   name of working folders for protein DB related steps\n";
	print "              by default, '--proc' is initialized from protein DB name in '--cfg' YAML file\n";
	print "  --samtools_m []  Set maximum memory per thread in samtools sort\n";
	print "\nOutput:\n\n";
	print "    predicted genes are located in working folder '--workdir'\n";
	print "    output files are named 'genemark.gtf' and 'genemark_supported.gtf'\n";
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
		'force'       => \$force,
		'verbose'     => \$verbose,
		'warnings'    => \$warnings,
		'debug'       => \$debug,
		'penalty=f'   => \$penalty,
		'gc=i'        => \$gc,
		'extend=i'    => \$extend,
		'paper=i'     => \$paper,
		'proc=s'      => \$proc,
		'samtools_m=s' => \$SAMTOOLS_M,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on the command line: $ARGV[0]\n" if( @ARGV > 0 );

	$verbose = 1 if $debug;
	$warnings = 1 if $debug;
	$verbose = 1 if $warnings;

	die "error, --cfg is not set\n" if !$cfg;
	StopIfNotFound($cfg);
	$cfg = abs_path($cfg);

	die "error, --workdir is not set\n" if !$workdir;
	StopIfNotFound($workdir);
	$workdir = abs_path($workdir);
	die "error, working directory must differ from ETP installation directory: $workdir\n" if ( $bin eq $workdir);

	if ($proc)
	{
		die "error, --proc option must be a folder name, not a path: $proc\n" if ( $proc =~ /\// );
		die "error, white space is not allowed in --proc option: $proc\n" if ( $proc =~ /\s/ );

		my %reserved = ( "data"=>"1", "arx"=>"1", "rnaseq"=>"1", "prothint_gmst.log"=>"1", "filter_gmst.log"=>"1", "genemark.gtf"=>"1", "genemark_supported.gtf"=>"1" );
		if ( exists $reserved{$proc} )
			{ die "error, option --proc is matching the reserved word\n"; }
	}

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

	print "### filter gmst predictions\n" if $verbose;

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
		system( "$bin/GeneMarkSTFiltering/filter.py --threads $cores --minUnsupportedLogOdds 0 --PROTHINT_PATH $bin/gmes/ProtHint/bin --outputFolder . $in_gmst_gff $in_tr_fasta $in_tr_gff $prot_db $in_genome 2> $workdir/filter_gmst.log" );

		system( "$bin/compare_intervals_exact.pl --f1 complete.gtf --f2 complete.gtf --trans --shared12 --out complete.id --original 1" );
		system( "$bin/select_by_transcript_id_from_gtf.pl complete.id complete.gtf complete_uniq.gtf" );
		system( "mv complete_uniq.gtf complete.gtf" );
	}
	else
	{
		print "# reusing HC gene filtering\n" if $verbose;
	}

	print "### filter gmst predictions ... done\n\n" if $verbose;
}
# -------------
sub GetFileNameFromPathName
{
	my $fpath = shift;

	die "error, file path is empty in GetFileNameFromPathName\n" if (!defined $fpath or ($fpath !~ /\S/));

	my $fname = '';

	if (( $fpath =~ /.+\/(.+?)\.gz$/ )or( $fpath =~ /.+\/(.+?)$/ ))
	{
		$fname = $1;
	}
	else
		{ die "error on file name parsing from file path: $fpath\n"; }

	die "error, file name is empty in GetFileNameFromPathName\n" if (!defined $fname or ($fname !~ /\S/));

	return $fname;
}
# -------------
sub PrepareProteinDBFile
{
	my $fpath = shift;

	print "### download protein db\n" if $verbose;

	die "error, protein file path is empty in PrepareProteinDBFile\n" if !$fpath;

	ChDir($workdir);
	MkDir("$workdir/data");
	MkDir("$workdir/arx");
	ChDir("$workdir/arx");

	my $fname = GetFileNameFromPathName($fpath);

	if ( CreateThis("$workdir/data/$fname") )
	{
		if ( CreateThis("$fname") )
		{
			FileFromPath($fpath);
		}
		StopIfNotFound($fname);
		system( "$bin/re_format_aa_fasta.pl $fname $workdir/data/$fname" ) and die "error in protein file parsing: $fname\n";
	}
	else
	{
		print "# reusing file 'data/$fname'\n" if $verbose;
	}

	StopIfNotFound("$workdir/data/$fname");

	if ( $clean and ( -e $fname ))
	{
		unlink $fname;
		print "# removed file $fname\n" if $verbose;
	}

	print "### download protein db ... done\n\n" if $verbose;

	return $fname;
}
# -------------
sub GeneMarkST
{
	print "### predict genes gmst\n" if $verbose;

	ChDir("$workdir/rnaseq");
	MkDir("$workdir/rnaseq/gmst");
	ChDir("$workdir/rnaseq/gmst");

	my $tseq = "$workdir/rnaseq/stringtie/transcripts_merged.fasta";
	StopIfNotFound($tseq);

	if( CreateThis("transcripts_merged.fasta.gff"))
	{
		system("$bin/gmst/gmst.pl --format GFF $tseq" );
		StopIfNotFound("transcripts_merged.fasta.gff");
	}
	else
	{
		print "# reusing gmst prediction: transcripts_merged.fasta.gff\n" if $verbose;
	}

	if ( CreateThis("genome_gmst.gtf"))
	{
		system( "$bin/gms2hints.pl --tseq $tseq --tgff transcripts_merged.fasta.gff --ggtf ../stringtie/transcripts_merged.gff --out genome_gmst.gtf --long --nodup --oneiso" );
		StopIfNotFound("genome_gmst.gtf");
	}

	if ( CreateThis("genome_gmst_for_HC.gtf"))
	{
		system( "$bin/gms2hints.pl --tseq $tseq --tgff transcripts_merged.fasta.gff --ggtf ../stringtie/transcripts_merged.gff --out genome_gmst_for_HC.gtf --nodup  --oneiso --min_cds 300");
		StopIfNotFound("genome_gmst_for_HC.gtf");
	}

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
	print "### create RNA-Seq hints\n" if $verbose;

	ChDir("$workdir/rnaseq");
	MkDir("$workdir/rnaseq/hints");
	ChDir("$workdir/rnaseq/hints");

	if ( CreateThis("hintsfile_merged.gff") )
	{
		my $manager = new Parallel::ForkManager( $cores );
		foreach my $set ( @rnaseq_sets )
		{
			print "# working on $set\n" if $verbose;

			my $in_bam = "$workdir/rnaseq/hisat2/mapping_". $set .".bam";
			my $out_hints = "bam2hints_". $set .".gff";

			if ( CreateThis($out_hints) )
			{
				StopIfNotFound($in_bam);
				$manager->start and next;
				system("$bin/bam2hints --intronsonly --in=$in_bam --out=$out_hints");
				$manager->finish;
			}
		}
		$manager->wait_all_children;

		if ( CreateThis("bam2hints_merged.gff") )
		{
			my @files = ();
			foreach my $set (@rnaseq_sets)
			{
				push @files, "bam2hints_". $set .".gff";
			}
		
			JoinIntronsGFF( "bam2hints_merged.gff", \@files );
		}

		if ( CreateThis("hintsfile_merged.gff"))
		{
			StopIfNotFound("$workdir/data/genome.fasta");
			StopIfNotFound("bam2hints_merged.gff");
			system( "$bin/filterIntronsFindStrand.pl $workdir/data/genome.fasta bam2hints_merged.gff --score > hintsfile_merged.gff" );
		}
	}
	else
	{
		print "# reusing hints file: hintsfile_merged.gff\n" if $verbose;
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
	print "### assemble transcripts\n" if $verbose;

	ChDir( "$workdir/rnaseq" );
	MkDir( "$workdir/rnaseq/stringtie" );
	ChDir( "$workdir/rnaseq/stringtie" );

	my $gff_list = '';

	if ( CreateThis("transcripts_merged.gff") )
	{
		foreach my $set ( @rnaseq_sets )
		{
			my $in_bam = "$workdir/rnaseq/hisat2/mapping_". $set .".bam";
			StopIfNotFound($in_bam);

			my $out_gff = "transcripts_". $set .".gff";
			$gff_list .= $out_gff ." ";

			print "# working on $set\n" if $verbose;
			if ( CreateThis($out_gff) )
			{
				system("stringtie -p $cores -o $out_gff $in_bam");
			}
			else
			{
				print "# reusing transcript assembly: $out_gff\n" if $verbose;
			}
		}

		if ( CreateThis("transcripts_merged.gff") )
		{
			system("stringtie --merge -o transcripts_merged.gff $gff_list");
		}
	}
	else
	{
		print "# reusing transcript assembly: transcripts_merged.gff\n" if $verbose;
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
	print "### map RNA-Seq\n" if $verbose;

	MkDir("$workdir/rnaseq");
	ChDir("$workdir/rnaseq");
	MkDir("$workdir/rnaseq/hisat2");
	ChDir("$workdir/rnaseq/hisat2");

	if ( !$bam and CreateThis("genome.1.ht2") )
	{
		print "# building genome index\n" if $verbose;
		system("hisat2-build --quiet -p $cores $workdir/data/genome.fasta genome");
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

		if ( CreateThis($out_bam) )
		{
			if ($bam)
			{
				system( "ln", "-sf", "$bam/$set.bam", "$out_bam" );
				print "# using BAM file $bam/$set.bam\n" if $verbose;
			}
			else
			{
				if ( CreateThis($out_sam) )
				{
					print "# mapping reads from $set\n" if $verbose;

					if ( -e $R )
					{
						system("hisat2 -x genome -U $R          --dta -p $cores -S $out_sam");
					}
					elsif ( -e $R1 and -e $R2 )
					{
						system("hisat2 -x genome -1 $R1 -2 $R2  --dta -p $cores -S $out_sam");
					}
					else
						{ die "error, reads are not found for mapping step: $set\n"; }

					print "# mapping reads from $set ... done\n" if $verbose;
				}
				else
				{
					print "# reusing SAM file $out_sam\n" if $verbose;
				}

				print "# from sam to sorted bam for $set\n" if $verbose;
				my $samtools_option = '';
				$samtools_option = " -m $SAMTOOLS_M " if $SAMTOOLS_M;
				system("samtools sort $samtools_option -o $out_bam -@ $cores $out_sam");
				print "# from sam to sorted bam for $set ... done\n" if $verbose;
			}
		}
		else
		{
			print "# reusing BAM file $out_bam\n" if $verbose;
		}

		if ( $clean and -e $out_sam and -e $out_bam) 
			{ unlink $out_sam or die "error on delete file: $out_sam $!\n"; }
	}

	print "### map RNA-Seq ... done\n\n" if $verbose;
}
# -------------
sub PrepareRNASeq
{
	print "### prepare RNA-Seq\n" if $verbose;

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

		if ((CreateThis($R1) or CreateThis($R2)) and CreateThis($R))
		{
			if ($local)
			{
				if ( -e "$local/$R1" and -e "$local/$R2" )
				{
					system( "ln", "-sf", "$local/$R1" );
					system( "ln", "-sf", "$local/$R2" );
				}
				elsif ( -e "$local/$R" )
				{
					system( "ln", "-sf", "$local/$R" );
				}
				else
					{ die "error, RNA-seq FASTQ file/s for library $set is not found at $local\n"; }
			}
			else
			{
				if ( CreateThis($sra) )
				{
					# 50G - is maximum allowed size of RNA-Seq SRA file to download
					# increase it for large SRA files
					system( "prefetch --max-size 50G $set" );
				}
				else
				{
					print "# reusing SRA file for $set\n" if $verbose;
				}

				StopIfNotFound("$workdir/rnaseq/reads/$set/$set.sra");

				system( "fastq-dump --split-3 $workdir/rnaseq/reads/$set/$set.sra" );

				if ((CreateThis($R1) or CreateThis($R2)) and CreateThis($R))
					{ die "error on unpacking SRA file $set\n"; }
			}
		}
		else
		{
			print "# reusing reads for $set\n" if $verbose;
		}

		if ( $clean and -e $sra ) { unlink $sra or die "error on delete file $sra: $!\n"; }
	}

	print "### prepare RNA-Seq ... done\n\n" if $verbose;
}
# -------------
sub MaskGenome
{
	my $fpath = shift;

	print "### mask genome\n" if $verbose;

	ChDir("$workdir/arx");

	if (!$fpath)
	{
		system( "touch ../data/repeats.gff" );
		print "# warning, no masking information was provided\n" if $verbose;
	}
	else
	{
		if ( CreateThis( "../data/repeats.gff" ) )
		{
			my $fname = FileFromPath($fpath);
			RepeatMaskerOutputToGFF( $fname, "../data/repeats.gff" );
			print "# repeat masking coordinates parsed from: $fname\n" if $verbose;

			if ($clean and ( -e $fname ))
			{
				unlink $fname;
				print "# removed file $fname\n" if $verbose;
			}
		}
		else
		{
			print "# reusing ../data/repeats.gff\n" if $verbose;
		}
	}

	ChDir("$workdir/data");

	StopIfNotFound( "genome.fasta" );

	if (!$fpath)
	{
		system( "ln", "-sf", "genome.fasta", "genome.softmasked.fasta" );
		print "# warning, no softmasking in file genome.softmasked.fasta\n" if $verbose;
	}
	else
	{
		if ( CreateThis( "genome.softmasked.fasta" ))
		{
			system("bedtools maskfasta -fi genome.fasta -bed repeats.gff -fo genome.softmasked.fasta -soft")
			and die "error on bedtools maskfasta\n";
		}
		else
		{
			print "# reusing genome.softmasked.fasta\n" if $verbose;
		}
	}

	StopIfNotFound( "genome.softmasked.fasta" );

	print "### mask genome ... done\n\n" if $verbose;
}
# -------------
sub CheckMasked
{
	my $MIN_RM = shift;

	my $value;

	print "### check genome masking\n" if $verbose;

	ChDir("$workdir/data");

	my $count_lower_case = 0;
	my $count_upper_case = 0;

	StopIfNotFound( "genome.softmasked.fasta" );

	my $txt = `$bin/probuild --stat --details --seq genome.softmasked.fasta`;

	if ( $txt =~ /SEQUENCE_atcg\s+(\d+)\s*/)
	{
		$count_lower_case = $1;
	}
	else { die "error in the code CheckRepeatsMasked\n"; }

	if ( $txt =~ /SEQUENCE_ACGT\s+(\d+)\s*/)
	{
		$count_upper_case = $1;
	}
	else { die "error in the code CheckRepeatsMasked\n"; }

	my $ratio = $count_lower_case/($count_upper_case + $count_lower_case);

	if ( $ratio < $MIN_RM )
	{
		$value = 0.03;
		print "# bases masked $ratio is below threshold $MIN_RM\n" if $verbose;
		print "# repeat penalty is set to: $value\n" if $verbose;
	}

	print "# bases atcg masked $count_lower_case, bases ATCG no masking $count_upper_case\n" if $verbose;
	print "# bases masked % ". sprintf("%.1f", 100.0*$ratio ) ."\n" if $verbose;

	print "### check genome masking ... done\n\n" if $verbose;

	return $value;
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

	print "### prepare genome\n" if $verbose;

	die "error, genome file path is empty in PrepareGenome\n" if !$fpath;

	ChDir($workdir);
	MkDir("$workdir/data");
	MkDir("$workdir/arx");
	ChDir("$workdir/arx");

	my $fname = '';

	if (CreateThis( "chr.names" ) or CreateThis( "../data/genome.fasta" ) or ( $softmask and CreateThis( "../data/genome.softmasked.fasta" )))
	{
		$fname = FileFromPath($fpath);
	}

	if ( CreateThis( "chr.names" ) )
	{
		 SetChrNames( $fname, "chr.names" );
	}
	else
	{
		print "# reusing file 'chr.names'\n" if $verbose;
	}

	if ( CreateThis( "../data/genome.fasta" ) )
	{
		system( "$bin/probuild --reformat_fasta --in $fname --out ../data/genome.fasta --uppercase 1 --letters_per_line 60 --include_sid chr.names --first_w --swap_sid" )
                and die "error on probuild in PrepareGenome\n";

		StopIfNotFound("../data/genome.fasta");
		print "# file 'data/genome.fasta' was created\n" if $verbose;
	}
	else
	{
		print "# reusing file 'data/genome.fasta'\n" if $verbose;
	}

	if ($softmask)
	{
		if( CreateThis( "../data/genome.softmasked.fasta" ) )
		{
			system( "$bin/probuild --reformat_fasta --in $fname --out ../data/genome.softmasked.fasta --uppercase 0 --letters_per_line 60 --include_sid chr.names --first_w --swap_sid" )
                        and die "error on probuild in PrepareGenome\n";;
	
			StopIfNotFound("../data/genome.softmasked.fasta");
			print "# file 'data/genome.softmasked.fasta' was created using softmask option\n" if $verbose;
		}
		else
		{
			print "# reusing file 'data/genome.softmasked.fasta'\n" if $verbose;
		}
	} 

	if ( $clean and ( -e $fname ))
	{
		unlink $fname;
		print "# removed file $fname\n" if $verbose;
	}

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

	my %h = (); # in case of the '--paper' mode
	            #   key is SeqID from FASTA defline
	            #   value is ChrID from FASTA defline
	            # in default case, doth key and value are SeqID

	my @arr = (); # SeqID's in this arrays are in the same order as in input FASTA file

	# collect SeqID's from FASTA deflines

	open( my $IN, $in ) or die "error on open file $in: $!\n";
	while(<$IN>)
	{
		if ( /^>/ )
		{
			my $seqid = '';
			my $chrid = '';

			if ( $paper ) 
			{
				if ( /mitochondrion|chloroplast|apicoplast|^>NW_/ )
				{
					print "warning, excluded FASTA record with defline: $_" if $warnings;
					next;
				}

				if ( !/Mus musculus/ and /^>(N[CT]_\S+) .* chromosome\s+(\S+)/ )
				{
					$seqid = $1;
					$chrid = $2;

					$chrid =~ s/,$//;
				}
				elsif ( /Mus musculus/ and /^>(NC_\S+) .* chromosome\s+(\S+)/ )
				{
					$seqid = $1;
					$chrid = $2;

					$chrid =~ s/,$//;
				}
				else
				{
					print "warning, excluded FASTA record with defline: $_" if $warnings;
					next;
				}
			}
			elsif ( /^>\s*(\S+)\s*/ )
			{
				$seqid = $1;
				$chrid = $seqid;
			}
			else
			{
				print "warning, excluded FASTA record with defline: $_" if $warnings;
				next;
			}

			# @arr is used to keep track the order of FASTA records in input file
			push @arr, $seqid;

			if ( ! exists $h{$seqid} )
			{
				$h{$seqid} = $chrid;
			}
			else
				{ die "error, FASTA SeqID duplication was detected: $seqid $_"; }

			print "# old-new seqid: ". $seqid ." ". $chrid ."\n" if $verbose;
		}
	}
	close $IN;

	if ( ! scalar @arr )
		{ die "error, FASTA selection was unsuccessful: $in\n"; }

	if ( my $result = ValueDuplicationInHash(\%h) )
		{ die "error, duplication was found in genome FASTA SeqID: $result\n"; }

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

	print "# loading file name from: $fpath\n" if $verbose;

	# extention '.gz' is excluded from file name
	my $fname = GetFileNameFromPathName($fpath);

	if ( $force or (! -e "$fname"))
	{
		if ( -e $fpath ) # if file is on local system - this should be absolute path
		{
			if ( $force or ( ! -e "$fname.gz" and ! -e $fname ))
			{
				system( "cp", $fpath, "./" ) and die "error on copy file from path: $fpath\n";
			}
		}
		elsif ( $fpath =~ /^https:\/\// or $fpath =~ /^http:\/\// or $fpath =~ /^ftp:\/\// )
		{
			if ( $force or ( ! -e "$fname.gz" and ! -e $fname ))
			{
				system( "wget", "-q", $fpath ) and die "error on wget file from path: $fpath\n";
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
		system( "gunzip", "-f", "$fname.gz" ) and die "error on gunzip\n";
	}

	die "error, file not found: $fname\n" if ( ! -e $fname );

	print "# loading file name from: ... done\n" if $verbose;

	return $fname;
}
# -------------
sub Start
{
	print "##### ETP running\n";
	my $date = localtime();
	print "# ". $date ."\n\n";
}
# -------------
sub Done
{
	my $date = localtime();
	print "\n# ". $date ."\n";
	print "##### ETP running ... done\n";
}
# -------------
sub StopIfNotFound
{
	my $fname = shift;

	die "error, file/forder name is empty\n" if !$fname;
	die "error, file/folder not found: $fname\n" if ( ! -e $fname );
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

	if ( !defined $fname or ( $fname !~ /\S/ ))
	{
		die "error, file name is empty in CreateThis\n";
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

	print "### extend hints\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	my $threshold = 4;

	if (!$gc)
	{
		if ( CreateThis("rna_conflicts.gff"))
		{
			system( "$bin/printRnaAlternatives.py pred_m/genemark.gtf r_hints_nonhc.gtf --otherIntrons prothint/prothint.gff --minIntronScore $threshold > rna_conflicts.gff" );
		}
	}
	else
	{
		if ( CreateThis("rna_conflicts_low.gff"))
		{
			system( "$bin/printRnaAlternatives.py pred_m_low/genemark.gtf r_hints_nonhc.gtf --otherIntrons prothint/prothint.gff --minIntronScore $threshold > rna_conflicts_low.gff" );
		}
		if ( CreateThis("rna_conflicts_medium.gff"))
		{
			system( "$bin/printRnaAlternatives.py pred_m_medium/genemark.gtf r_hints_nonhc.gtf --otherIntrons prothint/prothint.gff --minIntronScore $threshold > rna_conflicts_medium.gff" );
		}
		if ( CreateThis("rna_conflicts_high.gff"))
		{
			system( "$bin/printRnaAlternatives.py pred_m_high/genemark.gtf r_hints_nonhc.gtf --otherIntrons prothint/prothint.gff --minIntronScore $threshold > rna_conflicts_high.gff" );
		}
	}
	
	print "### extend hints ... done\n\n" if $verbose;
}
# -------------
sub SelectSupported
{
	my $wd = shift;

	print "### select supported genes\n" if $verbose;

	ChDir( "$workdir/$wd/nonhc" );

	# prothint.gff           from ProtHint and GMS-T seeds
	# hintsfile_merged.gff   from RNA-Seq mapping
	# prothint.gff           from ProtHint on nonhc regions with GM.hmm seeds
	# HC complete            from HC module
	# HC partial             from HC module

	system( "$bin/format_back.pl prothint/prothint.gff nonhc.trace > prothint_nonhc.gff" );
	system( "cat $hcc_genes $hcp_genes $prothint_hints $rnaseq_hints prothint_nonhc.gff > allHints.gff" );
	system( "$bin/selectSupportedSubsets.py ../genemark.gtf allHints.gff --fullSupport /dev/null --anySupport ../genemark_supported.gtf --noSupport /dev/null" );

	print "### select supported genes ... done\n\n" if $verbose;
}
# -------------
# function is similar to LoadPenaltyFromFile
# -------------
sub IsValidPenality
{
	my $fname = shift;

	my $status = 0;
	return $status if ( ! -e $fname );

	open( my $IN, $fname ) or die "error on open $fname: $!\n";
	my $line = <$IN>;
	if ( $line =~ /^\s*\d*\.*\d+\s*$/ )
	{
		$status = 1;
	}
	close $IN;

	return $status;
}
# -------------
sub EstinateMaskingPenalty
{
	my $wd = shift;

	print "### estimate masking penalty\n" if $verbose;

	ChDir("$workdir/$wd");
	MkDir("penalty");
	ChDir("penalty");

	my $estimated_penalty = 0;

	if (!$gc)
	{
		unlink "penalty.value" if !IsValidPenality("penalty.value");

		if ( CreateThis("penalty.value"))
		{	
			$estimated_penalty = `$bin/estimateMaskingPenalty.py  $hcc_genes  $workdir/data/genome.softmasked.fasta  ../model/output.mod  --GMES_PATH $bin/gmes  --threads $cores`;

			open( my $OUT, ">", "penalty.value" ) or die "error on open file penalty.value: $!\n";
			print $OUT "$estimated_penalty\n";
			close $OUT;
		}
	}
	else
	{
		if ( CreateThis("penalty.value.low"))
		{
			$estimated_penalty = `$bin/estimateMaskingPenalty.py  ../low/hc.gff  $workdir/data/genome.softmasked.fasta  ../low/output.mod  --GMES_PATH $bin/gmes  --threads $cores`;

			open( my $OUT, ">", "penalty.value.low" ) or die "error on open file penalty.value.low: $!\n";
			print $OUT "$estimated_penalty\n";
			close $OUT;
		}

		if ( CreateThis("penalty.value.medium"))
		{
			# medium GC
			$estimated_penalty = `$bin/estimateMaskingPenalty.py  ../medium/hc.gff  $workdir/data/genome.softmasked.fasta  ../medium/output.mod  --GMES_PATH $bin/gmes  --threads $cores`;

			open( my $OUT, ">", "penalty.value.medium" ) or die "error on open file penalty.value.medium: $!\n";
			print $OUT "$estimated_penalty\n";
			close $OUT;
		}

		if ( CreateThis("penalty.value.high"))
		{
			$estimated_penalty = `$bin/estimateMaskingPenalty.py  ../high/hc.gff  $workdir/data/genome.softmasked.fasta  ../high/output.mod  --GMES_PATH $bin/gmes  --threads $cores`;

			open( my $OUT, ">", "penalty.value.high" ) or die "error on open file penalty.value.high: $!\n";
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

	print "### combine hc and nonhc\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	if (!$gc)
	{
		StopIfNotFound("pred_m/genemark.gtf");
		system("$bin/format_back.pl  pred_m/genemark.gtf  nonhc.trace  > nonhc.gtf");
	}
	else
	{
		StopIfNotFound("pred_m_low/genemark.gtf");
		StopIfNotFound("pred_m_medium/genemark.gtf");
		StopIfNotFound("pred_m_high/genemark.gtf");

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

	print "### predict genes nonhc\n" if $verbose;

	my $mod        = $train_cfg{"model"};
	my $mod_low    = $train_cfg{"model_low"};
	my $mod_medium = $train_cfg{"model_medium"};
	my $mod_high   = $train_cfg{"model_high"};

	ChDir("$workdir/$wd/nonhc");

	# move to prothint section
	if ($gc)
	{
		# prothint folder structure differes between GC "0" and "1"
		# reusing folder from "0" in "1" 
		MkDir("prothint");
		system( "cat prothint_low/prothint.gff prothint_medium/prothint.gff prothint_high/prothint.gff > prothint/prothint.gff" );
		system( "cat prothint_low/evidence.gff prothint_medium/evidence.gff prothint_high/evidence.gff > prothint/evidence.gff" );
	}

	print "# creating hints for prediction step\n" if $verbose;

	# introns - overlap rnaseq and prothint
	system( "$bin/compare_intervals_exact.pl --f1 prothint/prothint.gff --f2 r_hints_nonhc.gtf --out evi_all.gff --shared12 --original 1 --intron --no" );

	# add start and stops from ProtHint-HMM-seeds and ProtHint-GMST-seeds to evidence
	system( "grep -P \"\tstart_codon\t\" prothint/evidence.gff >> evi_all.gff" );
	system( "grep -P \"\tstop_codon\t\"  prothint/evidence.gff >> evi_all.gff" );

	system( "grep -P \"\tstart_codon\t\" p_evi_nonhc.gtf >> evi_all.gff" );
	system( "grep -P \"\tstop_codon\t\"  p_evi_nonhc.gtf >> evi_all.gff" );

	# filter out hints overlapping HC-partial
	system( "$bin/filter_hints_partial.pl --hc hcp_hints_nonhc.gtf --hints evi_all.gff --out evi.gff" );

	# add HC-partial hints - excluding start codons
	system( "cat hcp_hints_nonhc.gtf | grep -v -P \"\tstart_codon\t\" >> evi.gff" );

	# from ExtendEvidence - how to handle in case of the update run?
	system( "cat rna_conflicts.gff >> evi.gff" )        if ( -e "rna_conflicts.gff" );
	system( "cat rna_conflicts_low.gff >> evi.gff" )    if ( -e "rna_conflicts_low.gff" );
	system( "cat rna_conflicts_medium.gff >> evi.gff" ) if ( -e "rna_conflicts_medium.gff" );
	system( "cat rna_conflicts_high.gff >> evi.gff" )   if ( -e "rna_conflicts_high.gff" );

	if (!$gc)
	{
		MkDir("pred_m");

		ChDir("pred_m");
		if ( CreateThis("genemark.gtf") or $train_cfg{"run"} )
		{
			print "# predicting using $mod\n" if $verbose;
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty --predict_with $mod --seq ../nonhc.fasta --cores $cores --verbose --evi ../evi.gff --max_gap 40000 --max_mask 40000 > loginfo");
		}
	}
	else
	{
		MkDir("pred_m_low");
		MkDir("pred_m_medium");
		MkDir("pred_m_high");

		ChDir("$workdir/$wd/nonhc/pred_m_low");
		if ( CreateThis("genemark.gtf") or $train_cfg{"run"} )
		{
			print "# predicting using $mod_low\n" if $verbose;
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_low --predict_with $mod_low --seq ../low.fasta --cores $cores --verbose --evi ../evi.gff  --max_gap 40000 --max_mask 40000 > loginfo");
		}

		ChDir("$workdir/$wd/nonhc/pred_m_medium");
		if ( CreateThis("genemark.gtf") or $train_cfg{"run"} )
		{
			print "# predicting using $mod_medium\n" if $verbose;
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_medium --predict_with $mod_medium --seq ../medium.fasta --cores $cores --verbose --evi ../evi.gff  --max_gap 40000 --max_mask 40000 > loginfo");
		}

		ChDir("$workdir/$wd/nonhc/pred_m_high");
		if ( CreateThis("genemark.gtf") or $train_cfg{"run"} )
		{
			print "# predicting using $mod_high\n" if $verbose;
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

	print "### ProtHint\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	if (!$gc)
	{
		MkDir("prothint");

		ChDir("prothint");
		if ( CreateThis("evidence.gff"))
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
		if ( CreateThis("evidence.gff"))
		{
			system("$bin/gmes/ProtHint/bin/prothint.py  --geneMarkGtf ../for_prothint_low/genemark.gtf ../low.fasta  $workdir/data/$prot 2> log");
		}

		ChDir("$workdir/$wd/nonhc/prothint_medium");
		if ( CreateThis("evidence.gff"))
		{
			system("$bin/gmes/ProtHint/bin/prothint.py  --geneMarkGtf ../for_prothint_medium/genemark.gtf ../medium.fasta  $workdir/data/$prot 2> log");
		}

		ChDir("$workdir/$wd/nonhc/prothint_high");
		if ( CreateThis("evidence.gff"))
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

	print "### predict genes for ProtHint\n" if $verbose;

	ChDir("$workdir/$wd/nonhc");

	my $mod        = $train_cfg{"model"};
	my $mod_low    = $train_cfg{"model_low"};
	my $mod_medium = $train_cfg{"model_medium"};
	my $mod_high   = $train_cfg{"model_high"};

	# introns - overlap rnaseq and prothint
	print "# creating hints for prediction: evi_ini.gff\n" if $verbose;

	system( "$bin/compare_intervals_exact.pl --f1 p_hints_nonhc.gtf --f2 r_hints_nonhc.gtf --out evi_ini.gff --shared12 --original 1 --intron --no" );
	system( "cat p_evi_nonhc.gtf >> evi_ini.gff" );
	system( "cat hcp_hints_nonhc.gtf >> evi_ini.gff" );

	if (!$gc)
	{
		MkDir("for_prothint");

		ChDir("for_prothint");
		if ( CreateThis("genemark.gtf"))
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty --predict_with $mod --seq ../nonhc.fasta --cores $cores --verbose --evi ../evi_ini.gff --max_gap 40000 --max_mask 40000 > loginfo");
		}
	}
	else
	{
		MkDir("for_prothint_low");
		MkDir("for_prothint_medium");
		MkDir("for_prothint_high");

		ChDir("$workdir/$wd/nonhc/for_prothint_low");
		if ( CreateThis("genemark.gtf"))
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_low --predict_with $mod_low --seq ../low.fasta --cores $cores --verbose --evi ../evi_ini.gff --max_gap 40000 --max_mask 40000 > loginfo");
		}

		ChDir("$workdir/$wd/nonhc/for_prothint_medium");
		if ( CreateThis("genemark.gtf"))
		{
			system("$bin/gmes/gmes_petap.pl --soft_mask 1000 --mask_penalty $penalty_medium --predict_with $mod_medium --seq ../medium.fasta --cores $cores --verbose --evi ../evi_ini.gff --max_gap 40000 --max_mask 40000 > loginfo");
		}

		ChDir("$workdir/$wd/nonhc/for_prothint_high");
		if ( CreateThis("genemark.gtf"))
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

	print "### training ab initio parameters\n" if $verbose;

	ChDir("$workdir/$wd");
	MkDir("$workdir/$wd/model");
	ChDir("$workdir/$wd/model");

	if ( CreateThis("output.mod"))
	{
		TranscriptIdFromRegions( "../hc_regions.gtf", "training.list" );
		system( "$bin/select_by_transcript_id_from_gtf.pl training.list $hcc_genes training.gtf" );

		if(!$bp)
		{
			system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
		}
		else
		{
			system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta --bp" );
		}
	}

	StopIfNotFound("$workdir/$wd/model/output.mod");

	$train_cfg{"model"} = "$workdir/$wd/model/output.mod";

	print "### training ab initio parameters ... done\n\n" if $verbose;
}
# -------------
sub TrainModelGC
{
	my $wd = shift;

	print "### training GC abinitio parametwers\n" if $verbose;

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
	if ( CreateThis("output.mod"))
	{
		system( "$bin/select_by_gene_id_from_gtf.pl ../hc/low.ids    $hcc_genes genes.gtf" );
		TranscriptIdFromRegions( "../hc_regions.gtf", "training.list" );
		system( "$bin/select_by_transcript_id_from_gtf.pl training.list genes.gtf training.gtf" );
		system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
	}

	ChDir("$workdir/$wd/medium");
	if ( CreateThis("output.mod"))
	{
		system( "$bin/select_by_gene_id_from_gtf.pl ../hc/medium.ids $hcc_genes genes.gtf" );
		TranscriptIdFromRegions( "../hc_regions.gtf", "training.list" );
		system( "$bin/select_by_transcript_id_from_gtf.pl training.list genes.gtf training.gtf" );
		system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
	}

	ChDir("$workdir/$wd/high");
	if ( CreateThis("output.mod"))
	{
		system( "$bin/select_by_gene_id_from_gtf.pl ../hc/high.ids   $hcc_genes genes.gtf" );
		TranscriptIdFromRegions( "../hc_regions.gtf", "training.list" );
		system( "$bin/select_by_transcript_id_from_gtf.pl training.list genes.gtf training.gtf" );
		system( "$bin/train_super.pl --hc training.gtf --dna $workdir/data/genome.fasta" );
	}

	StopIfNotFound("$workdir/$wd/low/output.mod");
	StopIfNotFound("$workdir/$wd/medium/output.mod");
	StopIfNotFound("$workdir/$wd/high/output.mod");

	$train_cfg{"model_low"}    = "$workdir/$wd/low/output.mod";
	$train_cfg{"model_medium"} = "$workdir/$wd/medium/output.mod";
	$train_cfg{"model_high"}   = "$workdir/$wd/high/output.mod";
	
	print "### training GC abinitio parametwers ... done\n\n" if $verbose;
}
# -------------
sub PrepareGenomeTraining
{
	my $wd = shift;

	# from complete and partial high confidence genes (GFF) and sequence (FASTA)

	my $margin = 50; # for non-hc regions

	print "### prepare genome sequence for training\n" if $verbose;

	ChDir("$workdir/$wd");

	if ( CreateThis("hc_regions.gtf"))
	{
		system( "$bin/create_regions.pl --hcc $hcc_genes --hcp $hcp_genes --out hc_regions.gtf --margin $margin" )
		and die "error on create_regions.pl";
	}

	MkDir("$workdir/$wd/nonhc");
	ChDir("$workdir/$wd/nonhc");
	MkDir("$workdir/$wd/nonhc/regions");
	ChDir("$workdir/$wd/nonhc/regions");

	if ( CreateThis("../nonhc.fasta"))
	{
		system( "$bin/probuild --cut nonhc --seq $workdir/data/genome.softmasked.fasta --regions ../../hc_regions.gtf --reverse --allow_x --trace ../nonhc.trace --soft --min_contig 120" )
		and die "error on probuild in PrepareGenomeTraining";

		JoinNonhs( ".", "../nonhc.fasta" );
	}

	ChDir("$workdir/$wd/nonhc");

	if ( CreateThis("hcp_hints_nonhc.gtf"))
	{
		# create hints from incomplete genes
		# remove start codons - if star exist
		system("$bin/gmes/rescale_gff.pl --out hcp_hints_nonhc.gtf --in $hcp_genes      --trace nonhc.trace");
		system("grep -v start_codon hcp_hints_nonhc.gtf > hcp.tmp ; mv hcp.tmp hcp_hints_nonhc.gtf");
	}

	if ( CreateThis("r_hints_nonhc.gtf"))
	{
		system("$bin/gmes/rescale_gff.pl --out r_hints_nonhc.gtf   --in $rnaseq_hints   --trace nonhc.trace");
	}

	if ( CreateThis("p_evi_nonhc.gtf"))
	{
		system("$bin/gmes/rescale_gff.pl --out p_evi_nonhc.gtf     --in $prothint_evi   --trace nonhc.trace");
	}

	if ( CreateThis("p_hints_nonhc.gtf"))
	{
		system("$bin/gmes/rescale_gff.pl --out p_hints_nonhc.gtf   --in $prothint_hints --trace nonhc.trace");
	}

	print "### prepare genome sequence for training ... done\n\n" if $verbose;
}
# -------------
sub PrepareHCregions
{
	my $wd = shift;

	print "### prepare HC regions training\n" if $verbose;

	ChDir("$workdir/$wd");
	MkDir("hc");
	ChDir("hc");
	MkDir("regions");
	ChDir("regions");

	if ( CreateThis("../hc.fasta"))
	{
		system( "$bin/probuild --cut hc --seq $workdir/data/genome.fasta --regions ../../hc_regions.gtf  --allow_x --trace ../hc.trace --soft --min_contig 120" )
		and die "error on parce hc regions";

		JoinFasta( ".", "../hc.fasta" );
	}

	print "### prepare HC regions training ... done\n\n" if $verbose;
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
# The coding of this "sub" was not finished.
# "Sub" should check all scripts and libs called from this package.
# -------------
sub CheckInstallation
{
	print "### Checking ETP installation\n";
	print "# Checking Perl setup\n";

	# to install modules on RedHat/CentOS
	# yum -y install perl-App-cpanminus.noarch
	# cpanm YAML::XS
	# ...

	# Here we are checking for Perl dependencies
	my @modules = ( "YAML::XS", "MCE::Mutex", "Thread::Queue", "Math::Utils", "Hash::Merge", "Parallel::ForkManager" );

	foreach my $module (@modules)
	{
		my $respond = `perl -M$module -e 1`;
		print "$respond\n" if ( $respond =~ /\S/ );
	}

	print "### Checking ETP installation ... done\n\n";
}
# -------------

