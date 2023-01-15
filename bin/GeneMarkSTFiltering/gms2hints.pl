#!/usr/bin/perl
# ---------------------
# Alex Lomsadze
# update 2022
# GaTech
#
# Move GMS-T predictions from transcripts to genome
# ---------------------

use strict;
use warnings;

use Getopt::Long qw( GetOptions );
use Data::Dumper;

# command line variables
my $in_genome_gtf = '';  # required
my $fix_single_with_strand = 0;
my $in_transcriptome_gff = ''; # required
my $in_transcriptome_sequence = ''; # required
my $move_start_to_lorf = 0;
my $move_partial_to_lorf = 0;
my $move_short_to_lorf = 0;
my $move_short_to_partial = 0;
my $nodup = 0;
my $oneiso = 0;
my $mask_th = 100;
my $out_hints_gtf = ''; # required
my $out_transcriptome_gff = '';

my $stringtie = 0; # optional - not used right now
my $gencode = 0; # optional - not used right now
my $min_5utr_for_complete = 0; # not used right now
my $in_genome_sequence = ''; # not used right now

my $min_cds = 0;

my $verbose = 0;
my $warnings = 0;
my $debug = 0;

Usage() if ( @ARGV < 1 );
ParseCMD();

my %data = ();

# most of transcript information is not changing - read ones and keep
# $data{$tid}{'t_strand'} and @{$data{$tid}{'exons'}}[3n+2] can change for some single exon genes from "+""-" to "."
#
# $data{$tid}{'seqid'}               - genome contig id
# $data{$tid}{'t_L'}                 - transcript location L on genome
# $data{$tid}{'t_R'}                 - transcript location R on genome
# $data{$tid}{'t_strand'}            - transcript strand on genome
# $data{$tid}{'gid'}                 - gene id
# $data{$tid}{'tid'}                 - transcript id
# $data{$tid}{'index'}               - order of transcript in gene by StringTie
# $data{$tid}{'num_exons'}           - number of exons in transcript
# $data{$tid}{'t_length'}            - cumulative exon length
# $data{$tid}{'t_cov'}               - transcript coverage from StringTie
# @{$data{$tid}{'exons'}}            - array of exons
# @{$data{$tid}{'exons'}}[3n+0]      - exon location L on genome  
# @{$data{$tid}{'exons'}}[3n+1]      - exon location R on genome
# @{$data{$tid}{'exons'}}[3n+2]      - exon strand on genome

# this status is set by script
# $data{$tid}{'ts_good'}             - TRUE FALSE
#                                       multi exon transcript is bad if no strand
#                                       single exon transcript is bad if strand and no alternative multi exon isoform

# transcript sequence is created from exon coordinates on genome and genome sequence
# transcript sequence id reverse complemented when transcript is on negative strand

# $data{$tid}{'has_cds'}             - TRUE if CDS was found, FALSE otherwice

# $data{$tid}{'cds_L'};              - CDS location L on transcript
# $data{$tid}{'cds_R'}               - CDS location R on transcript
# $data{$tid}{'cds_strand'}          - CDS strand on transcript
# $data{$tid}{'score'}               - logodd score of CDS
# $data{$tid}{'cds_length'}          - CDS length R-L+1
# $data{$tid}{'strand_match'}        - TRUE if CDS strand is not contradicting transcript strand

# $data{$tid}{'begin'}               - location of start - can be modified from prediction
# $data{$tid}{'end'}                 - location of end - can be modified from prediction
# $data{$tid}{'has_start'}           - has ATG start
# $data{$tid}{'has_stop'}            - has TAA/TAG/TGA stop
# $data{$tid}{'cds_complete'}        - both start and stop codon are found
# $data{$tid}{'5utr'}                - 5'utr length
# $data{$tid}{'3utr'}                - 3'utr length
# $data{$tid}{'LORF'}                - start codon location of LORF in CDS
# $data{$tid}{'nt_cds_masked'}       - CDS bases in lowcase letters
# $data{$tid}{'upstream_stop'}       - location of upstream (in 5UTR) stop triplet
# $data{$tid}{'is_LORF'}             - TRUE/FALSE if predicted start matches LORF 
# $data{$tid}{'is_shorter'}          - 
# $data{$tid}{'is_longer'}           - 

# @{$data{$tid}{'cds'}}, $L
# @{$data{$tid}{'cds'}}, $R
# @{$data{$tid}{'cds'}}, $strand
# @{$data{$tid}{'cds'}}, cds phase
# @{$data{$tid}{'cds'}}, cds type Initial, Internal, Terminal, Single, Partial

# $data{$tid}{'is_nodup'}            - TRUE if this is selected isoform from the duplicates
# $data{$tid}{'is_best_iso'}         - TRUE id this is selectd isoform

my %gdata = ();

# $gdata{$gid}{'isoforms'}          - number of isoforms in gene
# $gdata{$gid}{'isocds'}            - number of CDS isoforms in gene
# $gdata{$gid}{'strand'}            - strand; all transcripts from the same gene must have the same strand
# $gdata{$gid}{'has_strand'}        - TRUE if at least one isoform has strand
# @{$gdata{$gid}{'tid'}}, $tid      - array of tid's for this gene

my %tseq = ();

print "# running code in verbose mode\n" if $verbose;
print "# running code in warning mode\n" if $warnings;
print "# running code in debug mode\n" if $debug;

CheckGenomeGTF($in_genome_gtf, $stringtie, $gencode) if $debug; # just check - no other functionality
LoadGenomeGTF($in_genome_gtf);   # fill in %data and %gdata 
CheckExonsSortedInTranscripts() if $debug; # just check - no other functionality
CheckAndSetGeneStrand($fix_single_with_strand); # set in transcri 'ts_good' and in genes 'has_strand' 'strand'
ReportTranscriptIsoformHistogram() if $verbose;
LoadTranscriptomeGFF($in_transcriptome_gff);
SetUtrAndCompleteStatus();
UpdateCDSinfoFromSequence($in_transcriptome_sequence);
SetLORFs();
if($move_start_to_lorf or $move_partial_to_lorf or $move_short_to_lorf or $move_short_to_partial)
{
	MoveStarts();
	UpdateRecords();
	SetLORFs();
}
TransferCDStoDNA();
AddCDStypeLabels();
ReportStatus() if $debug;
ReduceCDSduplication() if $nodup;
ReduceAlternative() if $oneiso;
FilterOutRepeats();
ReportCDS($out_hints_gtf);
ReportTranscriptGFF($out_transcriptome_gff) if $out_transcriptome_gff;

# ============================
sub UpdateRecords
{
	print "### UpdateRecords update CDS related information after moving start codon\n" if $verbose;

	foreach my $tid ( keys %data )
	{
		my $ref = $data{$tid};

		next if ! $ref->{'has_cds'};

		if ( $ref->{'cds_strand'} eq '+' )
		{
			$ref->{'5utr'} = $ref->{'begin'} - 1;
			$ref->{'3utr'} = $ref->{'t_length'} - $ref->{'end'};
			$ref->{'cds_length'} = $ref->{'end'} - $ref->{'begin'} + 1;
			$ref->{'nt_cds_masked'} = CountLowcase( substr( $tseq{$tid}, $ref->{'begin'}, $ref->{'cds_length'} ));
		}
		elsif ( $ref->{'cds_strand'} eq '-' )
		{
			$ref->{'5utr'} = $ref->{'t_length'} - $ref->{'begin'};
			$ref->{'3utr'} = $ref->{'end'} - 1;
			$ref->{'cds_length'} = $ref->{'begin'} - $ref->{'end'} + 1;
			$ref->{'nt_cds_masked'} = CountLowcase( substr( $tseq{$tid}, $ref->{'end'}, $ref->{'cds_length'} ));
		}
		else
			{ die "strand is missing\n"; }

		if ( $ref->{'has_stop'} and $ref->{'has_start'} )
		{
			$ref->{'cds_complete'} = 1;
		}
		else
		{
			$ref->{'cds_complete'} = 0;
		}
	}

	print "### UpdateRecords update CDS related information after moving start codon ... done\n" if $verbose;
}
# ============================
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'ggtf=s'    => \$in_genome_gtf,
		'tgff=s'    => \$in_transcriptome_gff,
		'gseq=s'    => \$in_genome_sequence,
		'tseq=s'    => \$in_transcriptome_sequence,
		'out=s'     => \$out_hints_gtf,
		'tout=s'    => \$out_transcriptome_gff,
		'verbose'   => \$verbose,
		'debug'     => \$debug,
		'warnings'  => \$warnings,
		'min5utr'   => \$min_5utr_for_complete,
		'LORFStart' => \$move_start_to_lorf,
		'complete'  => \$move_partial_to_lorf,
		'long'      => \$move_short_to_lorf,
		'upmax'     => \$move_short_to_partial,
		'oneiso'    => \$oneiso,
		'nodup'     => \$nodup,
		'mask_th=f' => \$mask_th,
		'stringtie' => \$stringtie,
		'gencode'   => \$gencode,
		'fixss'     => \$fix_single_with_strand,
		'min_cds=i' => \$min_cds,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line\n" if( @ARGV > 0 );

	$verbose = 1 if $debug;

	die "error, --ggtf is not set\n" if !$in_genome_gtf;
	die "error, --tgff is not set\n" if !$in_transcriptome_gff;
	die "error, --tseq is not set\n" if !$in_transcriptome_sequence;
	die "error, --out is not set\n"  if !$out_hints_gtf;
}
# ----------------------------
sub Usage
{
	print "Usage:  $0  --ggtf []  --tgff []  --out [] \n";
	print "  --ggtf [] input GTF formatted file with transcript coordinates on genome: Genome GTF\n";
	print "  --tgff [] input GFF formatted file with CDS coordinates on transcript: Transcript GFF\n"; 
	print "  --tseq [] input FASTA formatted file with transcript sequences\n";
	print "  --out  [] output GTF formatted file with CDS coordinates on genome\n";
	print "\n";
	print "Optional:\n";
	print "  --LORFStart  move predicted start codon location to LORF\n";
	print "  --complete   move partial start predictions to LORF\n";
	print "  --long       move short start predictions to LORF\n";
	print "  --upmax      move start predictions to partial if possible\n";
	print "  --nodup      keep only one CDS from multiple duplicates\n";
	print "  --oneiso     keep only one isoform with longest CDS\n";
	print "  --tout []    output GFF formatted file with CDS coordinates on transcripts matching --out file\n";
	print "  --fixss      fix single transcript & single exon genes with strand\n";
	print "  --mask_th [$mask_th] keep only CDS with repeat coverage below threhold %\n";
	print "  --min_cds [$min_cds] report only CDS longer then\n";
	print "\n";
	print "Developer:\n";
#	print "  --stringtie\n";
#	print "  --gencode\n";
#	print "  --gseq []\n";
#	print "  --min5utr []\n";
	print "  --verbose\n";
	print "  --warnings\n";
	print "  --debug\n";

	exit 1;
}
# ----------------------------
sub FilterOutRepeats
{
	print "### FilterOutRepeats filter out CDS covered by repeats\n" if $verbose;

	my $count_masked = 0; 

	foreach my $tid ( keys %data )
	{
		next if ! $data{$tid}{'has_cds'};
		next if ( $nodup and !$data{$tid}{'is_nodup'} );
		next if ( $oneiso and !$data{$tid}{'is_best_iso'} );

		my $repeat_coverage = 100*$data{$tid}{'nt_cds_masked'}/$data{$tid}{'cds_length'};

		if ( $repeat_coverage >= $mask_th )
		{
			$data{$tid}{'has_cds'} = 0;
			$count_masked += 1;
		}
	}

	if($verbose)
	{
		print "# number of CDS removed by masking: $count_masked\n";
	}

	print "### FilterOutRepeats filter out CDS covered by repeats ... done\n" if $verbose;
}
# ----------------------------
sub ReduceAlternative
{
	print "### ReduceAlternative selecting one isoform representative\n" if $verbose;

	my %h = ();

	my $cds_alt_count = 0;

	my $length_cov_1_0 = 0;
	my $length_cov_0_1 = 0;

	foreach my $tid ( keys %data )
	{
		next if ! $data{$tid}{'has_cds'};
		next if ( $nodup and !$data{$tid}{'is_nodup'} );

		my $gid = $data{$tid}{'gid'};

		$h{$gid}{'count'} += 1;

		if ( $h{$gid}{'count'} == 1 )
		{
			$h{$gid}{'tid'} = $tid;
			$h{$gid}{'cds_length'} = $data{$tid}{'cds_length'};
			$h{$gid}{'t_cov'} = $data{$tid}{'t_cov'};

			$data{$tid}{'is_best_iso'} = 1;
		}
		else
		{
			if ( $data{$tid}{'cds_length'} > $h{$gid}{'cds_length'} )
			{
				$length_cov_1_0 += 1 if ( $data{$tid}{'t_cov'} < $h{$gid}{'t_cov'} );

				$data{ $h{$gid}{'tid'} }{'is_best_iso'} = 0;
				$data{$tid}{'is_best_iso'} = 1;

				$h{$gid}{'tid'} = $tid;
				$h{$gid}{'cds_length'} = $data{$tid}{'cds_length'};
				$h{$gid}{'t_cov'} = $data{$tid}{'t_cov'};
			}
			else
			{
				$length_cov_0_1 += 1 if ( $data{$tid}{'t_cov'} > $h{$gid}{'t_cov'} );

				$data{$tid}{'is_best_iso'} = 0;
			}

#			print "# warning, cds alt found: $tid\n" if $debug;
			$cds_alt_count += 1;
		}
	}

	if ($verbose)
        {
		print "# total cases with alt cds: $cds_alt_count\n";
		print "# length vs cov 1-0: $length_cov_1_0\n";
		print "# length vs cov 0-1: $length_cov_0_1\n";

		my %hist = ();
		foreach my $key (keys %h)
		{
			$hist{ $h{$key}{'count'} } += 1;
		}
		print "# histogram of cds alt\n";
		foreach my $key (sort{$a<=>$b} keys %hist)
		{
			print "# alt $key $hist{$key}\n";
		}
	}

	print "### ReduceAlternative selecting one isoform representative ... done\n" if $verbose; 
}
# ----------------------------
sub ReduceCDSduplication
{
	print "### ReduceCDSduplication selecting one CDS from duplicates\n" if $verbose;

	my %h = ();
	my $cds_dup_count = 0;

	foreach my $tid ( keys %data )
	{
		next if ( ! $data{$tid}{'has_cds'} );

		my $id = $data{$tid}{'seqid'} ."_";
		foreach my $key (@{$data{$tid}{'cds'}})
		{
			$id .= ($key ."_");
		}

		$data{$tid}{'CDS_id'} = $id;

		$h{$id}{'count'} += 1;

		if ( $h{$id}{'count'} == 1 )
		{
			$data{$tid}{'is_nodup'} = 1;
			$h{$id}{'best'} = $tid;
		}
		else
		{
			if ( $data{$tid}{'t_length'} > $data{ $h{$id}{'best'} }{'t_length'} )
			{
				$data{ $h{$id}{'best'} }{'is_nodup'} = 0;
				$data{$tid}{'is_nodup'} = 1;
				$h{$id}{'best'} = $tid;
			}
			else
			{
				$data{$tid}{'is_nodup'} = 0;
			}

#			print "# warning, cds dup found: $tid\n" if $debug;
			$cds_dup_count += 1;
		}
	}

	if ($verbose)
	{
		print "# total cases with cds duplication: $cds_dup_count\n";

		my %hist = ();
		foreach my $key (keys %h)
		{
			$hist{ $h{$key}{'count'} } += 1;
		}
		print "# histogram of cds duplications\n";
		foreach my $key (sort{$a<=>$b} keys %hist)
		{
			print "# dup $key $hist{$key}\n";
		}
	}

	print "### ReduceCDSduplication selecting one CDS from duplicates ... done\n" if $verbose;
}
# ----------------------------
sub AddComment
{
	my $tid = shift;

	my $str = '';

	my $ref = $data{$tid};

	if ( $ref->{'LORF'} == -1 )
	{
		$str = "noORF";
	}
	elsif ( $ref->{'is_LORF'} and ! $ref->{'upstream_stop'} )
	{
		$str = "LORF_NOUPSTOP";
	}
	elsif ( $ref->{'is_LORF'} and $ref->{'upstream_stop'} )
	{
		$str = "LORF_UPSTOP";
	}
	elsif ( $ref->{'is_shorter'} and ! $ref->{'upstream_stop'} )
	{
		$str = "sORF_NOUPSTOP";
	}
	elsif ( $ref->{'is_shorter'} and $ref->{'upstream_stop'} )
	{
		$str = "sORF_UPSTOP";
	}
	elsif ( $ref->{'is_longer'} )
	{
		$str = "upLORF";
	}
	else
	{
		die "$tid\n";
	}

	return $str;
}
# ----------------------------
sub ReportCDS
{
	my $fname = shift;

	print "### ReportCDS creating output GTF file\n" if $verbose;

	my $count_split_starts = 0;
	my $count_split_stops = 0;
	my $ignored_strand_mismatch = 0;

	open( my $OUT, ">", $fname ) or die "error on open file $fname; $!\n";
	foreach my $tid ( keys %data )
	{
		next if ! $data{$tid}{'has_cds'};

		# multi-CDS with strand mismatch 
		if ( (scalar @{$data{$tid}{'cds'}} > 5 ) and ! $data{$tid}{'strand_match'} )
		{
			$ignored_strand_mismatch += 1;
			next;
		}

		next if ( $nodup and !$data{$tid}{'is_nodup'} );
		next if ( $oneiso and !$data{$tid}{'is_best_iso'} );
		next if ( $data{$tid}{'cds_length'} < $min_cds );

		my $score = $data{$tid}{'score'};

		print $OUT "###\n";
                print $OUT ("# ". $tid ." ". $data{$tid}{'nt_cds_masked'} ." ". $data{$tid}{'cds_length'} ." ". $score ."\n");

		my $seqid = $data{$tid}{'seqid'};
		my $L = 0;
		my $R = 0;
		my $strand = $data{$tid}{'t_strand'};
		my $phase = 0;
		my $gid = $data{$tid}{'gid'};
		my $type = '';

		my $class = AddComment($tid);

		if ($strand eq '.')
		{
			$strand = $data{$tid}{'cds_strand'};
		}
		if ( ! $data{$tid}{'strand_match'} )
		{
			$strand = "-" if ($data{$tid}{'t_strand'} eq '+');
			$strand = "+" if ($data{$tid}{'t_strand'} eq '-');
		}

		my $status = '';
		if ( $data{$tid}{'cds_complete'} )
		{
			$status = 'complete';
		}
		else
		{
			$status = 'partial';
		}

		if ($strand eq '+')
		{
			if ( $data{$tid}{'has_start'} )
			{
				$L = $data{$tid}{'cds'}[0];
				$R = $data{$tid}{'cds'}[0] + 2;
				if ( $R > $data{$tid}{'cds'}[1] )
				{
					$R = $data{$tid}{'cds'}[1];
					print "# warning, split start codon: $tid ". ($R - $L + 1) ."\n" if $warnings;
					$count_split_starts += 1;
				}

				print $OUT  "$seqid\tgmst\tstart_codon\t$L\t$R\t$score\t$strand\t0\tgene_id \"$gid\"; transcript_id \"$tid\"; status \"$status\"; class \"$class\";\n";
			}
		}
		elsif ($strand eq '-')
		{
			if ( $data{$tid}{'has_stop'} )
			{
				$L = $data{$tid}{'cds'}[0];
				$R = $data{$tid}{'cds'}[0] + 2;
				if ( $R > $data{$tid}{'cds'}[1] )
				{
					$R = $data{$tid}{'cds'}[1];
					print "# warning, split stop codon: $tid ". ($R - $L + 1)."\n" if $warnings;
					$count_split_stops += 1;
				}
				print $OUT  "$seqid\tgmst\tstop_codon\t$L\t$R\t$score\t$strand\t0\tgene_id \"$gid\"; transcript_id \"$tid\"; status \"$status\"; class \"$class\";\n";
			}
		}

		if ( $data{$tid}{'num_cds'} == 1 )
		{
			$L = $data{$tid}{'cds'}[0];
			$R = $data{$tid}{'cds'}[1];
			$phase = $data{$tid}{'cds'}[3];
			$type = $data{$tid}{'cds'}[4];
			print $OUT  "$seqid\tgmst\tCDS\t$L\t$R\t$score\t$strand\t$phase\tgene_id \"$gid\"; transcript_id \"$tid\"; status \"$status\"; cds_type \"$type\"; class \"$class\";\n";
		}
		else
		{
			for( my $i = 0; $i < $data{$tid}{'num_cds'}; $i++ )
			{
				$L = $data{$tid}{'cds'}[$i *5 + 0];
				$R = $data{$tid}{'cds'}[$i *5 + 1];
				$phase = $data{$tid}{'cds'}[$i *5 + 3];
				$type = $data{$tid}{'cds'}[$i *5 + 4];
				print $OUT  "$seqid\tgmst\tCDS\t$L\t$R\t$score\t$strand\t$phase\tgene_id \"$gid\"; transcript_id \"$tid\"; status \"$status\"; cds_type \"$type\"; class \"$class\";\n";

				if ( $i <  $data{$tid}{'num_cds'} - 1 )
				{
					$L = $data{$tid}{'cds'}[$i *5 + 1] + 1;
					$R = $data{$tid}{'cds'}[($i+1) *5 + 0] - 1;
					print $OUT  "$seqid\tgmst\tintron\t$L\t$R\t$score\t$strand\t0\tgene_id \"$gid\"; transcript_id \"$tid\"; status \"$status\"; class \"$class\";\n";
				}
			}
		}

		if ($strand eq '+')
		{
                        if ( $data{$tid}{'has_stop'} )
                        {
				my $i = $data{$tid}{'num_cds'} - 1;
				$L = $data{$tid}{'cds'}[$i *5 + 1] - 2;
				$R = $data{$tid}{'cds'}[$i *5 + 1];
				if ( $L < $data{$tid}{'cds'}[$i *5 + 0] )
				{
					$L = $data{$tid}{'cds'}[$i *5 + 0];
					print "# warning, split stop codon: $tid ". ($R - $L + 1)."\n" if $warnings;
					$count_split_stops += 1;
				}
				print $OUT  "$seqid\tgmst\tstop_codon\t$L\t$R\t$score\t$strand\t0\tgene_id \"$gid\"; transcript_id \"$tid\"; status \"$status\"; class \"$class\";\n";
			}
		}
		elsif ($strand eq '-')
		{
			if ( $data{$tid}{'has_start'} )
			{
				my $i = $data{$tid}{'num_cds'} - 1;
				$L = $data{$tid}{'cds'}[$i *5 + 1] - 2;
				$R = $data{$tid}{'cds'}[$i *5 + 1];
				if ( $L < $data{$tid}{'cds'}[$i *5 + 0] )
				{
					$L = $data{$tid}{'cds'}[$i *5 + 0];
					print "# warning, split start codon: $tid ". ($R - $L + 1)."\n" if $warnings;
					$count_split_starts += 1;
				}
				print $OUT  "$seqid\tgmst\tstart_codon\t$L\t$R\t$score\t$strand\t0\tgene_id \"$gid\"; transcript_id \"$tid\"; status \"$status\"; class \"$class\";\n";
			}
		}
	}
	close $OUT;

	if ($verbose)
	{
		print "# split starts: $count_split_starts\n";
		print "# split stops: $count_split_stops\n";
		print "# removed mismatch: $ignored_strand_mismatch\n";
	}

	print "### ReportCDS creating output GTF file ... done\n" if $verbose;
}
# ----------------------------
# Print updated annotation on transcript level.
sub ReportTranscriptGFF
{
	my $fname = shift;

	print "# ReportTranscriptGFF creating output transcript level GFF file.\n" if $verbose;

	open( my $OUT, ">", $fname ) or die "error on open file $fname; $!\n";
	foreach my $tid ( keys %data )
	{
		my $ref = $data{$tid};
		
		next if ! $ref->{'has_cds'};
		next if ( $nodup and !$ref->{'is_nodup'} );
		next if ( $oneiso and !$ref->{'is_best_iso'} );

		print $OUT "$tid\tgmst\tCDS\t$ref->{'cds_L'}\t$ref->{'cds_R'}\t0\t$ref->{'cds_strand'}\t0\t.\n";
	}
	close $OUT;

	print "# ReportTranscriptGFF creating output transcript level GFF file ... done\n" if $verbose;
}
# ----------------------------
sub AddCDStypeLabels
{
	print "### AddCDStypeLabels add GeneMark.hmm stile CDS labels to CDS\n" if $verbose;

	my $count_single = 0;
	my $count_initial = 0;
	my $count_internal = 0;
	my $count_terminal = 0;
	my $count_partial = 0;
	my $count_all = 0;

	my $count_no_start = 0;
	my $count_no_stop = 0;

	foreach my $tid ( keys %data )
	{
		my $ref = $data{$tid};

		next if ! $ref->{'has_cds'};

		my $cds_length_processed = 0;

		if ( $ref->{'cds_strand'}  eq '-' )
		{
			$cds_length_processed = $ref->{'cds_length'};
		}

		for( my $i = 0; $i < $ref->{'num_cds'}; $i++ )
		{
			$count_all += 1;

			# set GFF phase
			if ( $ref->{'cds'}[$i * 5 + 2] eq '+' )
			{
				$ref->{'cds'}[$i * 5 + 3] = ( 3 - ($cds_length_processed % 3)) % 3;
				$cds_length_processed += ($ref->{'cds'}[$i * 5 + 1] - $ref->{'cds'}[$i * 5 + 0] + 1);
			}
			elsif ( $ref->{'cds'}[$i * 5 + 2] eq '-' )
			{
				$cds_length_processed -=( $ref->{'cds'}[$i * 5 + 1] - $ref->{'cds'}[$i * 5 + 0] + 1);
				$ref->{'cds'}[$i * 5 + 3] = ( 3 - ($cds_length_processed % 3)) % 3;
			}

			# set GeneMark.hmm type
			if ( $ref->{'num_cds'} == 1 )
			{
				if ( $ref->{'has_start'} and $ref->{'has_stop'} )
				{
					$ref->{'cds'}[$i * 5 + 4] = "Single";
					$count_single += 1;
				}
				else
				{
					$count_no_start += 1 if ! $ref->{'has_start'};
					$count_no_stop += 1 if ! $ref->{'has_stop'};
				}
			}
			else
			{
				if (( $i == 0 )and( $ref->{'has_start'} )and( $ref->{'cds'}[$i * 5 + 2] eq '+' ))
				{
					$ref->{'cds'}[$i * 5 + 4] = "Initial";
					$count_initial += 1;
				}
				elsif (( $i == 0 )and( $ref->{'has_stop'} )and( $ref->{'cds'}[$i * 5 + 2] eq '-' ))
				{
					$ref->{'cds'}[$i * 5 + 4] = "Terminal";
					$count_terminal += 1;
				}
				elsif (( $i == $ref->{'num_cds'} - 1 )and( $ref->{'has_stop'} )and( $ref->{'cds'}[$i * 5 + 2] eq '+' ))
				{
					$ref->{'cds'}[$i * 5 + 4] = "Terminal";
					$count_terminal += 1;
				}
				elsif (( $i == $ref->{'num_cds'} - 1 )and( $ref->{'has_start'} )and( $ref->{'cds'}[$i * 5 + 2] eq '-' ))
				{
					$ref->{'cds'}[$i * 5 + 4] = "Initial";
					$count_initial += 1;
				}
				elsif (  $i > 0 and ( $i < $ref->{'num_cds'} - 1 ))
				{
					$ref->{'cds'}[$i * 5 + 4] = "Internal";
					$count_internal += 1;
				}
				else
				{
					$count_no_start += 1 if ! $ref->{'has_start'};
					$count_no_stop += 1 if ! $ref->{'has_stop'};
				}
			}
		}
	}

	if ($debug)
	{
		print "# single: $count_single\n";
		print "# initial: $count_initial\n";
		print "# internal: $count_internal\n";
		print "# terminal: $count_terminal\n";
		print "# partial: ". ( $count_all - $count_single - $count_initial - $count_internal - $count_terminal ) ."\n";
		print "# all: $count_all\n";
		print "# no start: $count_no_start\n";
		print "# no stop: $count_no_stop\n";
	}

	print "### AddCDStypeLabels add GeneMark.hmm stile CDS labels to CDS ... done\n" if $verbose;
}
# ----------------------------
sub TransferCDStoDNA
{
	print "### TransferCDStoDNA transfering CDS to genome\n" if $verbose;

	my $count_strand_mismatch = 0;

	foreach my $tid ( keys %data )
	{
		my $ref = $data{$tid};

		next if ! $ref->{'has_cds'};

		# single exon
		if ( $ref->{'num_exons'} == 1 )
		{
			if (( $ref->{'begin'} < $ref->{'end'} )and( $ref->{'cds_strand'} eq '+' ))
			{
				push @{$ref->{'cds'}}, ($ref->{'exons'}[0] + $ref->{'begin'} - 1);
				push @{$ref->{'cds'}}, ($ref->{'exons'}[0] + $ref->{'end'} - 1);
			}
			elsif (( $ref->{'end'} < $ref->{'begin'} )and( $ref->{'cds_strand'} eq '-' ))
			{
				push @{$ref->{'cds'}}, ($ref->{'exons'}[0] + $ref->{'end'} - 1);
				push @{$ref->{'cds'}}, ($ref->{'exons'}[0] + $ref->{'begin'} - 1);
			}
			else
				{die "error, check the code\n";}

			push @{$ref->{'cds'}}, '.';
			push @{$ref->{'cds'}}, 0;  # CDS phase
			push @{$ref->{'cds'}}, "Partial"; # initialization - label is updated later
			$ref->{'num_cds'} = 1;

			if ( $ref->{'t_strand'} eq '.' )
			{
				$ref->{'cds'}[2] = $ref->{'cds_strand'};
				$ref->{'strand_match'} = 1;
			}
			elsif(( $ref->{'t_strand'} eq '+' )and( $ref->{'cds_strand'} eq '+' ))
			{
				$ref->{'cds'}[2] = '+';
				$ref->{'strand_match'} = 1;
			}
			elsif(( $ref->{'t_strand'} eq '+' )and( $ref->{'cds_strand'} eq '-' ))
			{
				$ref->{'cds'}[2] = '-';
				$ref->{'strand_match'} = 0;
				$count_strand_mismatch += 1;
#				print "# strand mismatch $tid\n" if $debug;
			}
			elsif(( $ref->{'t_strand'} eq '-' )and( $ref->{'cds_strand'} eq '+' ))
			{
				$ref->{'cds'}[0] = $ref->{'exons'}[0] + $ref->{'t_length'} - $ref->{'end'} + 1 - 1;
				$ref->{'cds'}[1] = $ref->{'exons'}[0] + $ref->{'t_length'} - $ref->{'begin'} + 1 - 1;
				$ref->{'cds'}[2] = '-';
				$ref->{'strand_match'} = 1;
			}
			elsif(( $ref->{'t_strand'} eq '-' )and( $ref->{'cds_strand'} eq '-' ))
			{
				$ref->{'cds'}[0] = $ref->{'exons'}[0] + $ref->{'t_length'} - $ref->{'begin'} + 1 - 1;
				$ref->{'cds'}[1] = $ref->{'exons'}[0] + $ref->{'t_length'} - $ref->{'end'} + 1 - 1;
				$ref->{'cds'}[2] = '+';
				$ref->{'strand_match'} = 0;
				$count_strand_mismatch += 1;
#				print "# strand mismatch $tid\n" if $debug;
			}
			else
				{ print "warning: $tid\n"; }
		}
		# if multi exon
		else
		{
			my $L = 0;
			my $R = 0;

			if (( $ref->{'t_strand'} eq '+' )and( $ref->{'cds_strand'} eq '+' ))
			{
				$L = $ref->{'begin'};
				$R = $ref->{'end'};
				$ref->{'strand_match'} = 1
			}
			elsif (( $ref->{'t_strand'} eq '-' )and( $ref->{'cds_strand'} eq '+' ))
			{
				$L = $ref->{'t_length'} - $ref->{'end'} + 1;
				$R = $ref->{'t_length'} - $ref->{'begin'} + 1;
				$ref->{'strand_match'} = 1
			}
			elsif (( $ref->{'t_strand'} eq '+' )and( $ref->{'cds_strand'} eq '-' ))
			{
				$ref->{'strand_match'} = 0;
				$count_strand_mismatch += 1;
				$L = $ref->{'end'};
				$R = $ref->{'begin'};
			}
			elsif (( $ref->{'t_strand'} eq '-' )and( $ref->{'cds_strand'} eq '-' ))
			{
				$ref->{'strand_match'} = 0;
				$count_strand_mismatch += 1;
				$L = $ref->{'t_length'} - $ref->{'begin'} + 1;
				$R = $ref->{'t_length'} - $ref->{'end'} + 1;
			}

			my $current_exon_length = 0;
                        my $L_found = 0;
                        my $R_found = 0;
                        my $L_out = 0;
                        my $R_out = 0;

			for( my $i = 0; $i < $ref->{'num_exons'}; $i++ )
			{
				$current_exon_length = ($ref->{'exons'}[$i * 3 + 1] - $ref->{'exons'}[$i * 3 + 0] + 1);

				if (( $current_exon_length < $L )and !$L_found )
				{
					$L -= $current_exon_length;
					$R -= $current_exon_length;
					next;
				}
						
				if (( $current_exon_length >= $L )and !$L_found )
				{
					$L_out = $ref->{'exons'}[$i * 3 + 0] + $L - 1;
					$L_found = 1;
				}
				elsif ( $L_found )
				{
					$L_out = $ref->{'exons'}[$i * 3 + 0];
				}

				if (( $current_exon_length >= $R )and !$R_found)
				{
					$R_out = $ref->{'exons'}[$i * 3 + 0] + $R - 1;
					$R_found = 1;
				}
				elsif (( $current_exon_length < $R )and !$R_found)
				{
					$R_out = $ref->{'exons'}[$i * 3 + 1];
				}
				elsif ( $R_found )
				{
					last;
				}

				push @{$ref->{'cds'}}, $L_out;
				push @{$ref->{'cds'}}, $R_out;
				if (( $ref->{'cds_strand'} eq '+' )and $ref->{'strand_match'} )
				{
					push @{$ref->{'cds'}}, $ref->{'t_strand'};
				}
				else
				{
					push @{$ref->{'cds'}}, "-" if ($ref->{'t_strand'} eq '+');
					push @{$ref->{'cds'}}, "+" if ($ref->{'t_strand'} eq '-');
				}
				push @{$ref->{'cds'}}, 0; # CDS phase
				push @{$ref->{'cds'}}, 'Partial';
				$ref->{'num_cds'} += 1;

				$L -= $current_exon_length;
				$R -= $current_exon_length;
			}
		}
	}

	print "# strand mismatches: $count_strand_mismatch\n" if $debug;

	print "### TransferCDStoDNA transfering CDS to genome ... done\n" if $verbose;
}
# ----------------------------
sub LoadSequnce
{
	my $fname = shift;
	my $ref = shift;

	my $id = '';

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		if ( /^>\s*(\S+)/ )
		{
			$id = $1;
			next;
		}

		my $seq = $_;
		$seq =~ s/\s//g;
		$seq =~ s/\d//g;

		$ref->{$id} .= $seq;
	}
}
# ----------------------------
sub LoadTranscriptomeGFF
{
	my $fname = shift;

	print "### LoadTranscriptomeGFF loading GeneMark.hmm GFF file: $fname\n" if $verbose;
	die "error, file name is empty\n" if !$fname;

	my $count_ignored_CDS = 0;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(my $line = <$IN>)
	{
		next if $line =~ /^#/;
		next if $line =~ /^\s*$/;

		if ( $line =~ /^(\S+)\t\S+\tCDS\t(\d+)\t(\d+)\t\S+\t([-+])\t/ )
		{
			my $tid = $1;
			my $L = $2;
			my $R = $3;
			my $strand = $4;

			my $logodd_score = 0;
			if ( $line =~ /logodd=(\S+)/ )
			{
				$logodd_score = $1;
			}

			die "error, seqid is not found in transcripts : $tid\n" if ( ! exists $data{$tid} );

			my $ref = $data{$tid};

			if ( exists $ref->{'has_cds'} and $ref->{'has_cds'} )
			{
				print "warning, more then one CDS found in transcript: $tid $ref->{'score'} $logodd_score\n" if $warnings;

				$count_ignored_CDS += 1;

				if ( $ref->{'score'} < $logodd_score )
				{
					$ref->{'cds_L'} = $L;
					$ref->{'cds_R'} = $R;
					$ref->{'cds_strand'} = $strand;
					$ref->{'score'} = $logodd_score;
					$ref->{'cds_length'} = $R - $L + 1;
					$ref->{'has_cds'} = 1;
					$ref->{'strand_match'} = 0; # initialization
				}
			}
			else
			{
				$ref->{'cds_L'} = $L;
				$ref->{'cds_R'} = $R;
				$ref->{'cds_strand'} = $strand;
				$ref->{'score'} = $logodd_score;
				$ref->{'cds_length'} = $R - $L + 1;
				$ref->{'has_cds'} = 1;
				$ref->{'strand_match'} = 0; # initialization
			}
		}
		else
			{die "error, unexpected line format found: $_";}
	}
	close $IN;

	if ($verbose)
	{
		print "# number of CDS ignored due to multi CDS per transcript: $count_ignored_CDS\n";
	}

	print "### LoadTranscriptomeGFF loading GeneMark.hmm GFF file: $fname ... done\n" if $verbose;
}
# ----------------------------
sub SetUtrAndCompleteStatus
{
	print "### SetUtrAndCompleteStatus setting UTR and complete status\n" if $verbose;

	# In old GeneMark.hmm GFF output the complete/partial gene status is missing
	# Complete status is estimated from location of predicted start and stop codons
	# This estimation is wrong in some cases
	# Some of these errors will be corrected later by parsing start/stops from sequence

	foreach my $tid (keys %data)
	{
		my $ref = $data{$tid};

		# transcripts without prediction
		$ref->{'has_cds'} = 0 if ! exists $ref->{'has_cds'};
		next if ! $ref->{'has_cds'};

		$ref->{'begin'} = -1;
		$ref->{'end'} = -1;
		$ref->{'has_start'} = 0;
		$ref->{'has_stop'} = 0;
		$ref->{'cds_complete'} = 0;
		$ref->{'5utr'} = -1;
		$ref->{'3utr'} = -1;

		my $L = $ref->{'cds_L'};
		my $R = $ref->{'cds_R'};
		
		if ( $ref->{'cds_strand'} eq '+' )
		{
			$ref->{'begin'} = $L;
			$ref->{'end'} = $R;
			$ref->{'5utr'} = $L - 1;
			$ref->{'3utr'} = $ref->{'t_length'} - $R;
		}
		elsif ( $ref->{'cds_strand'} eq '-' )
		{
			$ref->{'begin'} = $R;
			$ref->{'end'} = $L;
			$ref->{'5utr'} = $ref->{'t_length'} - $R;
			$ref->{'3utr'} = $L - 1;
		}
		else
			{ die "strand is missing\n"; }

		$ref->{'has_start'} = 1 if ( $ref->{'5utr'} > 2 );
		$ref->{'has_stop'}  = 1 if ( $ref->{'3utr'} > 2 );

		$ref->{'cds_complete'} = 1 if ( $ref->{'has_stop'} and $ref->{'has_start'} );
	}

	print "### SetUtrAndCompleteStatus setting UTR and complete status ... done\n" if $verbose;

	ReportStatus() if $debug;
}
# ----------------------------
sub ReportStatus
{
	my $num_good_ts = 0;
	my $num_bad_ts = 0;
	my $num_has_cds = 0;
	my $num_no_cds = 0;
	my $num_complete_cds = 0;
	my $num_partial_cds = 0;
	my $num_strand_match = 0;
	my $num_strand_mismatch = 0;
	my $num_no_start = 0;
	my $num_no_stop = 0;
	my $num_no_both = 0;

	foreach my $tid (keys %data)
	{
		my $ref = $data{$tid};

		if ( $ref->{'ts_good'} )
		{
			$num_good_ts += 1;
		}
		else
		{
			$num_bad_ts += 1;
		}

		if ( $ref->{'has_cds'} )
		{
			$num_has_cds += 1;
		}
		else
		{
			$num_no_cds += 1;
		}

		if ( $ref->{'has_cds'} )
		{
			if ( $ref->{'cds_complete'} )
			{
				$num_complete_cds += 1
			}
			else
			{
				$num_partial_cds += 1;
			}

			if ( $ref->{'strand_match'} )
			{
				$num_strand_match += 1;
			}
			else
			{
				$num_strand_mismatch += 1;
			}

			$num_no_start += 1 if !$ref->{'has_start'};
			$num_no_stop += 1 if !$ref->{'has_stop'};
			$num_no_both += 1 if ( !$ref->{'has_start'} and !$ref->{'has_stop'} );
		}
	}

	print "### Report\n";
	print "# good transcripts strand: $num_good_ts\n";
	print "# bad transcripts strand: $num_bad_ts\n";
	print "# has cds: $num_has_cds\n";
	print "# no cds $num_no_cds\n";
	print "# complete cds: $num_complete_cds\n";
	print "# partial cds $num_partial_cds\n";
	print "# strand match: $num_strand_match\n";
	print "# strand mismatch: $num_strand_mismatch\n";
	print "# no start codons: $num_no_start\n";
	print "# no stop codons: $num_no_stop\n";
	print "# no start and stop codons: $num_no_both\n";
	print "### Report ... done\n";
}
# ----------------------------
sub CountLowcase
{
	my $str = shift;

	my $count = () = $str =~ /[atcgn]/g;

	return $count;
}
# ----------------------------
sub UpdateCDSinfoFromSequence
{
	my $fname = shift;

	print "### UpdateCDSinfoFromSequence updating complete/partial CDS status from sequence: $fname\n" if $verbose;
	die "error, file name is empty\n" if !$fname;

	LoadSequnce( $fname, \%tseq );

	my $count_starts_added = 0;
	my $count_stops_added = 0;
	my $count_starts_removed = 0;
	my $count_stops_removed = 0;

	my $count_up_stop = 0;
	my $count_up_start = 0;

	foreach my $tid (keys %data)
	{
		my $ref = $data{$tid};

		next if ! $ref->{'has_cds'};

		$ref->{'nt_cds_masked'} = CountLowcase( substr( $tseq{$tid}, $ref->{'cds_L'}, $ref->{'cds_R'} - $ref->{'cds_L'} + 1) );

		$ref->{'upstream_stop'} = 0;
		$ref->{'LORF'} = -1;
		$ref->{'is_LORF'} = 0;

		my $start_codon = '';
		my $stop_codon = '';

		if ( $ref->{'cds_strand'} eq '+' )
		{
			$start_codon = uc( substr( $tseq{$tid}, $ref->{'cds_L'} - 1, 3 ));
			$stop_codon  = uc( substr( $tseq{$tid}, $ref->{'cds_R'} - 3, 3 ));
		}
		elsif ( $ref->{'cds_strand'} eq '-' )
		{
			$start_codon = RevComp( uc( substr( $tseq{$tid}, $ref->{'cds_R'} - 3, 3 )));
			$stop_codon  = RevComp( uc( substr( $tseq{$tid}, $ref->{'cds_L'} - 1, 3 )));
		}

		if ( $ref->{'has_stop'} == 1 )
		{
			if ( $stop_codon =~ /TAA|TAG|TGA/ )
			{
				; # as expected
			}
			else
			{
				# this can happen when sequence has "NNN"
				$ref->{'has_stop'} = 0;
				$ref->{'cds_complete'} = 0;
				$count_stops_removed += 1;
			}
		}
		else
		{
			if ( $stop_codon =~ /TAA|TAG|TGA/ )
			{
				$ref->{'has_stop'} = 1;
				$count_stops_added += 1;
#				print "# stop codon found with short UTR: $ref->{'3utr'}\n" if $debug;
			}
		}

		if ( $ref->{'has_start'} == 1 )
		{
			if ( $start_codon =~ /ATG/ )
			{
				; # as expected
			}
			else
			{
				# this can happen when sequence has "NNN"
				$ref->{'has_start'} = 0;
				$ref->{'cds_complete'} = 0;
				$count_starts_removed += 1;
			}
		}
		else
		{
			if ( $start_codon =~ /ATG/ )
			{
				$ref->{'has_start'} = 1;
				$count_starts_added += 1;
#				print "# start codon found with short UTR: $ref->{'5utr'}\n" if $debug;
			}
		}

		# find upstream stop and LORF
		if ( $ref->{'cds_strand'} eq '+' )
		{
			my $pos = $ref->{'cds_R'} - 6;
			my $found_up_ATG = 0;

			while( $pos >= 0 )
			{
				my $codon = uc( substr( $tseq{$tid}, $pos, 3 ));

				if ( $codon =~ /^(TAA|TAG|TGA)$/ )
				{
					$ref->{'upstream_stop'} = $pos + 1;
					$count_up_stop += 1;
					last;
				}

				if ( $codon eq "ATG" )
				{
					$ref->{'LORF'} = $pos + 1;
					$found_up_ATG += 1 if ( $ref->{'LORF'} < $ref->{'cds_L'} );
				}

				$pos -= 3;
			}

			$count_up_start += 1 if $found_up_ATG;
		}
		elsif ( $ref->{'cds_strand'} eq '-' )
		{
			my $pos = $ref->{'cds_L'} + 2;
			my $found_up_ATG = 0;

			while ( $pos <  $ref->{'t_length'} )
			{
				my $codon = RevComp( uc( substr( $tseq{$tid}, $pos, 3 )));

				if ( $codon =~ /^(TAA|TAG|TGA)$/ )
				{
					$ref->{'upstream_stop'} = $pos + 3;
					$count_up_stop += 1;
					last;
				}

				if ( $codon eq "ATG" )
				{
					$ref->{'LORF'} = $pos + 3;
					$found_up_ATG += 1 if ( $ref->{'cds_R'} < $ref->{'LORF'} );
				}

				$pos += 3;
			}

			$count_up_start += 1 if $found_up_ATG;

		}

		$ref->{'cds_complete'} = 1 if ( $ref->{'has_stop'} and $ref->{'has_start'} );
	}

	if ($debug)
	{
		print "# stops added: $count_stops_added\n";
		print "# stops removed: $count_stops_removed\n";
		print "# starts added: $count_starts_added\n";
		print "# starts removed: $count_starts_removed\n";
		print "# up stops found: $count_up_stop\n";
		print "# up starts found: $count_up_start\n";
	}

	print "### UpdateCDSinfoFromSequence updating complete/partial CDS status from sequence: $fname ... done\n" if $verbose;

	ReportStatus() if $debug;
}
# ----------------------------
sub SetLORFs
{
	print "### SetLORFs setting and collecting stat on LORFs\n" if $verbose;

	my $count_up_stop = 0;
	my $count_less_than_LORF = 0;
	my $count_more_than_LORF = 0;
	my $count_is_LORF = 0;
	my $count_is_LORF_dir = 0;
	my $count_is_LORF_rev = 0;
	my $count_no_LORF = 0;

	foreach my $tid ( keys %data )
	{
		my $ref = $data{$tid};

		next if ! $ref->{'has_cds'};

		if ( $ref->{'cds_strand'} eq "+" )
		{
			if ( $ref->{'LORF'} == -1 )
			{
				$ref->{'is_LORF'} = 0;
				$count_no_LORF += 1;
			}
			elsif ( $ref->{'LORF'} < $ref->{'begin'} )
			{
				$ref->{'is_LORF'} = 0;
				$ref->{'is_shorter'} = 1;
				$ref->{'is_longer'} = 0;
				$count_less_than_LORF += 1;
			}
			elsif ( $ref->{'LORF'} == $ref->{'begin'} )
			{
				$ref->{'is_LORF'} = 1;
				$ref->{'is_shorter'} = 0;
				$ref->{'is_longer'} = 0;
				$count_is_LORF += 1;
				$count_is_LORF_dir += 1;
			}
	       		else
			{
				$ref->{'is_LORF'} = 0;
				$ref->{'is_shorter'} = 0;
				$ref->{'is_longer'} = 1;
				$count_more_than_LORF += 1;
			}
		}
		elsif ( $ref->{'cds_strand'} eq "-" )
		{
			if ( $ref->{'LORF'} == -1 )
			{
				$ref->{'is_LORF'} = 0;
				$count_no_LORF += 1;
			}
			elsif ( $ref->{'LORF'} > $ref->{'begin'} )
			{
				$ref->{'is_LORF'} = 0;
				$ref->{'is_shorter'} = 1;
				$ref->{'is_longer'} = 0;
				$count_less_than_LORF += 1;
			}
			elsif ( $ref->{'LORF'} == $ref->{'begin'} )
			{
				$ref->{'is_LORF'} = 1;
				$ref->{'is_shorter'} = 0;
				$ref->{'is_longer'} = 0;
				$count_is_LORF += 1;
				$count_is_LORF_rev += 1;
			}
			else
			{
				$ref->{'is_LORF'} = 0;
				$ref->{'is_shorter'} = 0;
				$ref->{'is_longer'} = 1;
				$count_more_than_LORF += 1;
			}
		}

		if ( $ref->{'upstream_stop'} > 0 )
		{
			$count_up_stop += 1;
		}
	}

	if ($debug)
	{
		print "# ORFs with upstream stops: $count_up_stop\n";
		print "# CDS is shorter then LORF: $count_less_than_LORF\n";
		print "# CDS is longer then LORF: $count_more_than_LORF\n";
		print "# CDS is LORF: $count_is_LORF\n";
		print "# CDS is LORF on dir: $count_is_LORF_dir\n";
		print "# CDS is LORF on rev: $count_is_LORF_rev\n";
		print "# no start for ORF: $count_no_LORF\n";
	}

	print "### SetLORFs setting and collecting stat on LORFs ... done\n" if $verbose;
}
# ----------------------------
sub MoveStarts
{
	print "### MoveStarts change start codon location - cmd option based\n" if $verbose;

	if ( !$move_start_to_lorf and !$move_partial_to_lorf and !$move_short_to_lorf and !$move_short_to_partial )
	{
		print "# SetStarts change start codon location - was not turned on\n" if $verbose;
		return;
	}

	foreach my $tid ( keys %data )
	{
		my $ref = $data{$tid};

		next if ! $ref->{'has_cds'};

		if ( $move_start_to_lorf and !$move_partial_to_lorf and !$move_short_to_lorf and !$move_short_to_partial )
		{
			if ( $ref->{'LORF'} == -1 )
			{
				$ref->{'has_start'} = 0;
			}
			elsif ( $ref->{'cds_strand'} eq "+" )
			{
				$ref->{'begin'} = $ref->{'LORF'};
				$ref->{'has_start'} = 1;
			}
			elsif ( $ref->{'cds_strand'} eq "-" )
			{
				$ref->{'begin'} = $ref->{'LORF'};
				$ref->{'has_start'} = 1;
			}

			$ref->{'cds_complete'} = 1 if ( $ref->{'has_stop'} and $ref->{'has_start'} );
		}
		elsif ( !$move_start_to_lorf and $move_partial_to_lorf and !$move_short_to_lorf and !$move_short_to_partial )
		{
			if ( $ref->{'LORF'} == -1 )
			{
				$ref->{'has_start'} = 0;
			}
			elsif ( $ref->{'cds_strand'} eq "+" )
			{
				if ( $ref->{'begin'} < $ref->{'LORF'} )
				{
					$ref->{'begin'} = $ref->{'LORF'};
				}

				$ref->{'has_start'} = 1;
			}
			elsif ( $ref->{'cds_strand'} eq "-" )
			{
				if( $ref->{'LORF'} < $ref->{'begin'} )
				{
					$ref->{'begin'} = $ref->{'LORF'};
				}

				$ref->{'has_start'} = 1;
			}

			$ref->{'cds_complete'} = 1 if ( $ref->{'has_stop'} and $ref->{'has_start'} );
		}
		elsif ( !$move_start_to_lorf and !$move_partial_to_lorf and $move_short_to_lorf and !$move_short_to_partial )
		{
			if ( $ref->{'LORF'} == -1 )
			{
				$ref->{'has_start'} = 0;
			}
			elsif ( $ref->{'cds_strand'} eq "+" )
			{
				if ( $ref->{'LORF'} < $ref->{'begin'} )
				{
					$ref->{'begin'} = $ref->{'LORF'};
				}

				$ref->{'has_start'} = 1;
			}
			elsif ( $ref->{'cds_strand'} eq "-" )
			{
				if( $ref->{'begin'} < $ref->{'LORF'} )
				{
					$ref->{'begin'} = $ref->{'LORF'};
				}

				$ref->{'has_start'} = 1;
			}

			$ref->{'cds_complete'} = 1 if ( $ref->{'has_stop'} and $ref->{'has_start'} );
		}
		elsif ( !$move_start_to_lorf and !$move_partial_to_lorf and !$move_short_to_lorf and $move_short_to_partial )
		{
			 if ( $ref->{'LORF'} == -1 )
                        {
                                $ref->{'has_start'} = 0;
                        }
			elsif ( !$ref->{'upstream_stop'} )
			{
				$ref->{'has_start'} = 0;

				if ( $ref->{'cds_strand'} eq "+" )
				{
					$ref->{'begin'} = $ref->{'begin'} % 3;

					if ( ! ( $ref->{'end'} - $ref->{'begin'} + 1)%3 )
						{ die "error, CDS length is not 3n: $tid\n"; }
				}
				elsif ( $ref->{'cds_strand'} eq "-" )
				{
					$ref->{'begin'} = $ref->{'t_length'} - ($ref->{'t_length'} - $ref->{'begin'}) % 3;

					if ( ! ( $ref->{'begin'} - $ref->{'end'} + 1)%3 )
						{ die "error, CDS length is not 3n: $tid\n"; }
				}
			}

			$ref->{'cds_complete'} = 1 if ( $ref->{'has_stop'} and $ref->{'has_start'} );
		}
		else
			{ die "error, check cmd"; }
	}

	print "### MoveStarts change start codon location - cmd option based ... done\n" if $verbose;

	ReportStatus() if $debug;
}
# ----------------------------
sub RevComp
{
	my $str = shift;

	$str = reverse $str;
	$str =~ s/A/1/g;
	$str =~ s/T/A/g;
	$str =~ s/1/T/g;
	$str =~ s/C/2/g;
	$str =~ s/G/C/g;
	$str =~ s/2/G/g;
	
	return $str;
}
# ----------------------------
sub CheckExonsSortedInTranscripts
{
	print "### CheckExonsSortedInTranscripts checking exons in transcripts\n" if $verbose;

	foreach my $tid ( keys %data )
	{
		TranscriptLength($tid);

		my $ref = $data{$tid};

		die "error, no exons in transcript: $tid\n" if ! $ref->{'num_exons'};

		for( my $j = 0, my $i = 1 ; $i < $ref->{'num_exons'}; $i++, $j++)
		{
			if ( $ref->{'exons'}[$j*3 + 1] >= $ref->{'exons'}[$i*3 + 1] )
				{ die "error, exon coordiantes are not sorted: $tid\n"; }
		}
	}
	print "### CheckExonsSortedInTranscripts checking exons in transcripts ... done\n" if $verbose;
}
# ----------------------------
sub CheckAndSetGeneStrand
{
	my $fix = shift;

	print "### CheckAndSetGeneStrand checking strand in transcripts\n" if $verbose;

	# set ts_good - in transcrupts
	# set has_strand and strand - in genes

	my $count_bad_multi = 0;
	my $count_bad_single = 0;
	my $count_bad_single_fixed = 0;

	# loop over all multi-exon transcripts first
	# set strand based on the strand of multi exon transcripts
	foreach my $tid ( keys %data )
	{
		my $ref = $data{$tid};

                if ( $ref->{'num_exons'} > 1 )
		{
			# multi exon transcript without strand
			if ($ref->{'t_strand'} eq '.')
			{
				print "# warning, multi exon transcript without strand: $tid $ref->{'t_strand'}\n" if $warnings;
				$ref->{'ts_good'} = 0;
				$count_bad_multi += 1;
			}
			else
			{
				$ref->{'ts_good'} = 1;
				$gdata{$ref->{'gid'}}{'has_strand'} = 1;
			
				if ( ! $gdata{$ref->{'gid'}}{'strand'} )
				{
					$gdata{$ref->{'gid'}}{'strand'} = $ref->{'t_strand'};
				}
				elsif ( $gdata{$ref->{'gid'}}{'strand'} ne $ref->{'t_strand'} )
					{ die "error, gene has isoforms on different strands: $tid\n"; }
			}
		}
	}

	# check strand of single exon transcripts
	foreach my $tid ( keys %data )
	{
		my $ref = $data{$tid};

		if ( $ref->{'num_exons'} == 1 )
		{
			# single exon transcript with no alternative isoforms and strand
			if (( $ref->{'t_strand'} ne '.' )and( ! $gdata{$ref->{'gid'}}{'has_strand'} ))
			{
				print "# warning, single exon gene with strand: $tid $ref->{'t_strand'}\n" if $warnings;

				$count_bad_single += 1;

				if ( $fix )
				{
					$ref->{'t_strand'} = '.';
					$ref->{'exons'}[2] = '.';
					$gdata{$ref->{'gid'}}{'has_strand'} = 0;
					$ref->{'ts_good'} = 1;
					$count_bad_single_fixed += 1;
				}
				else
				{
					$ref->{'ts_good'} = 0;
				}
			}
			else
			{
				$ref->{'ts_good'} = 1;
			}
		}
	}

	if ($verbose)
	{
		print "# multi exon transcripts without strand: $count_bad_multi\n";
		print "# single exon transcript with strand: $count_bad_single\n";
		print "# single exon with fixed strand: $count_bad_single_fixed\n" if $fix;
	}

	print "### CheckAndSetGeneStrand checking strand in transcripts ... done\n" if $verbose;
}
# ----------------------------
sub ReportTranscriptIsoformHistogram
{
	print "### ReportTranscriptIsoformHistogram histogram isoforms count: isoforms-per-genes genes\n";

	my %isoform_count = ();
	my $count_genes = 0;
	my $count_transcripts = 0;

	foreach my $key (keys %gdata)
	{
		$isoform_count{ $gdata{$key}{'isoforms'} } += 1;
		$count_genes += 1;
		$count_transcripts += $gdata{$key}{'isoforms'};
	}

	print "# all transcripts (good+bad) are used in this histogram\n";
	print "# genes: $count_genes\n";
	print "# transcripts: $count_transcripts\n";
	print "# column_1 isoforms_per_gene ; column_2 number_of_genes_with_as_many_isoforms\n";
	foreach my $key ( sort{ $a <=> $b} keys %isoform_count )
	{
		print "# $key $isoform_count{$key}\n";
	}
	print "### ReportTranscriptIsoformHistogram histogram isoforms count: isoforms-per-genes genes ... done\n";
}
# ----------------------------
sub LoadGenomeGTF
{
	my $fname = shift;

	print "### LoadGenomeGTF loading exon/transcript from genome annotation GTF file: $fname\n" if $verbose;
	die "error, file name is empty\n" if !$fname;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		if ( /^(\S+)\t\S+\t(exon|transcript)\t(\d+)\t(\d+)\t\S+\t([-+.])\t\S+\tgene_id \"(\S+)\"; transcript_id \"(\S+)\";/ )
		{
			my $seqid = $1;
			my $type = $2;
			my $L = $3;
			my $R = $4;
			my $strand = $5;
			my $gid = $6;
			my $tid = $7;
			
			if ( $type eq "transcript" )
			{
				# data
				$data{$tid}{'seqid'} = $seqid;
				$data{$tid}{'t_L'} = $L;
				$data{$tid}{'t_R'} = $R;
				$data{$tid}{'t_strand'} = $strand;
				$data{$tid}{'gid'} = $gid;
				$data{$tid}{'tid'} = $tid;

				my $cov = 0;
				$cov = $1 if ( /cov \"(\S+)\";/ );
				$data{$tid}{'t_cov'} = $cov;

				$data{$tid}{'index'} = 0;
				if ( $tid =~ /^(.*)\.(\d+)$/ )
				{
					$data{$tid}{'index'} = $2 if ( $1 eq $gid );
				}

				$data{$tid}{'ts_good'} = 0; # initialization - value is set in CheckAndSetGeneStrand

				# gdata
				$gdata{$gid}{'isoforms'} += 1;
				push @{$gdata{$gid}{'tid'}}, $tid;
				$gdata{$gid}{'has_strand'} = 0; # initialization - value is set in CheckAndSetGeneStrand
				$gdata{$gid}{'strand'} = ''; # initialization - value is set in CheckAndSetGeneStrand
			}
			elsif ( $type eq "exon" )
			{
				$data{$tid}{'num_exons'} += 1;
				$data{$tid}{'t_length'} += ($R - $L + 1);
				push @{$data{$tid}{'exons'}}, $L;
				push @{$data{$tid}{'exons'}}, $R;
				push @{$data{$tid}{'exons'}}, $strand;
			}
		}
		else
			{ die "error, unexpected line format found: $_"; }
	}
	close $IN;

	print "### LoadGenomeGTF loading exon/transcript from genome annotation GTF file: $fname ... done\n" if $verbose;
}
# ----------------------------
sub TranscriptLength
{
	my $tid = shift;

	my $size = scalar @{$data{$tid}{'exons'}};
	die "error, check exon parsing for tid $tid\n" if ( $size %3 != 0 );

	print "warning, array size difference detected: $size/3 $data{$tid}{'num_exons'}\n" if ( $size/3 != $data{$tid}{'num_exons'} );

	my $length = 0;
	for( my $i = 0; $i < $size; $i += 3 )
	{
		$length += $data{$tid}{'exons'}[$i+1] - $data{$tid}{'exons'}[$i] + 1;
	}

	print "warning, transcript length difference detected: $data{$tid}{'t_length'} $length\n" if ( $data{$tid}{'t_length'} != $length );

	return $length;
}
# ----------------------------
sub CheckGenomeGTF
{
	# check for valid GTF file format
	# additional non-GTF rules are imposed

	my $fname = shift;
	my $stringtie = shift;
	my $gencode = shift;

	print "### CheckGenomeGTF checking GTF file format: $fname\n" if $verbose;
	die "error, file name is empty\n" if !$fname;

	# hash %h
	# tid - 
	#     - tid
	#     - gid
	#     - seqid
	#     - strand
	#     - L
	#     - R

	my %h = ();
	my $IN;

	# first pass through file
	open( $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		if ( /^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t\S+\tgene_id \"(\S+)\"; transcript_id \"(\S+)\";/ )
		{
			my $seqid = $1;
			my $label = $2;
			my $type = $3;
			my $L = $4;
			my $R = $5;
			my $strand = $6;
			my $gid = $7;
			my $tid = $8;

			if ( $stringtie )
			{
				if ( $label ne "StringTie" )
					{ die "error, StingTie label is not found: $_\n"; }
				if ( ! /cov \"(\S+)\";/ )
					{ die "error, cov is not found in transcript: $tid\n"; }
				if ( $type !~ /^(exon|transcript)$/ )
					{ die "error, unexpected type found in StingTie file: $_\n"; }
			}

			if ( $gencode )
			{
				if ( $label !~ /^(ENSEMBL|HAVANA)$/ )
					{ die "error, unexpected label in Gencode file: $_\n"; }
			}

			if ($L > $R)
				{ die "error, L is more than R: $L $R\n" };

			if ($L<1 or $R<1)
				{ die "error, L and R are not positive numnbers: $L $R\n" };

			if ( $type eq "transcript" )
			{
				if ( ! exists $h{$tid}{'tid'} )
				{
					$h{$tid}{'tid'} = $tid;
				}
				else
					{ die "error, transcript ID duplication found: $tid\n"; }

				$h{$tid}{'L'} = $L;
				$h{$tid}{'R'} = $R;
			}

			if ( ! exists $h{$tid}{'gid'} )
			{
				$h{$tid}{'gid'} = $gid;
			}
			elsif ( $h{$tid}{'gid'} ne $gid )
				{ die "error, one tid is connected to two gids: $tid $gid $h{$tid}{'gid'}\n"; }

			if ( ! exists $h{$tid}{'seqid'})
			{
				$h{$tid}{'seqid'} = $seqid;
			}
			elsif ( $h{$tid}{'seqid'} ne $seqid )
				{ die "error, one tid is connected to two seqids: $tid $seqid $h{$tid}{'seqid'}\n"; }

			if ( ! exists $h{$tid}{'strand'})
			{
				$h{$tid}{'strand'} = $strand;
			}
			elsif ( $h{$tid}{'strand'} ne $strand )
				{ die "error, one tid is connected to two strands: $tid $strand $h{$tid}{'strand'}\n"; }
		}
		else
			{ die "error, unexpected line format found: $_"; }
	}
	close $IN;

	# second pass through file
	open( $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		next if /^#/;
		next if /^\s*$/;

		if ( /^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t\S+\tgene_id \"(\S+)\"; transcript_id \"(\S+)\";/ )
		{
			my $seqid = $1;
			my $label = $2;
			my $type = $3;
			my $L = $4;
			my $R = $5;
			my $strand = $6;
			my $gid = $7;
			my $tid = $8;

			if ( $h{$tid}{'L'} > $L or $h{$tid}{'R'} < $R )
				{  die "error, feature is outside of transcript: $_\n"; }
		}
	}
	close $IN;

	print "### CheckGenomeGTF checking GTF file format: $fname ... done\n" if $verbose;
}
# ----------------------------

