#!/usr/bin/env perl
# --------------------------------------------
# Alex Lomsadze
# GaTech
# Last update 2021
#
# Create training file for GeneMark.hmm (ES style) from
#  * GTF
#  * genome
# --------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $VERSION = "v1_2022";

# --------------------------------------------
my $gtf = '';         # input GTF file
my $gseq = '';        # input file with genome sequence FASTA
my $out = '';         # output
my $margin = 1000;    # intergenic margin

my $verbose = '';
my $warnings = '';
my $debug = '';
# --------------------------------------------
Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();
# --------------------------------------------

my %genome = ();
LoadGenome( $gseq, \%genome );

my %gid2record = ();
LoadGTF( $gtf, \%gid2record );

# put one gene at a time into this array of arrays
my @record = ();

my $gcount = 0;

# output
open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";

foreach my $gid (sort{$a cmp $b } keys %gid2record)
{
	foreach my $line ( @{$gid2record{$gid}} )
	{
		AddLineToRecord( $line, \@record );
	}

	my $cds_length = CDSlength(\@record);

	if( ! $cds_length )
	{ 
		print "# warning, gene was ignored: $gid\n";
		@record = ();
		next;
	}

	$gcount += 1;


	my @cds = SelectCDSsorted(\@record);
	my @introns = CreateCdsIntrons(\@cds, $gid);
	CompareIntrons( \@record, \@introns );
	push (@record, @introns);

	print $OUT ("## seq_id: ". $record[0][0] ."\n");
	print $OUT ("## ". $gcount ."\t". $cds_length ."\t". (scalar @introns) ."\n");

	foreach my $feature (@record)
	{
		my $label = '';
		my $len = $feature->[4] - $feature->[3] + 1;
		my $Lph = -1;
		my $Rph = -1;
		my $seq = '';

		my $signal = '';
		my $stop_label = '';

		if ( $feature->[2] eq "CDS" )
		{
#			print "# processing CDC\n" if $verbose;

			my $phase = $feature->[7];

			if ( $feature->[6] eq '+' )
			{
				if    ( $phase == 0 ) { $Lph = 1; }
				elsif ( $phase == 1 ) { $Lph = 3; }
				elsif ( $phase == 2 ) { $Lph = 2; }

				$phase = (3 - ( 3 - $phase + $len )%3)%3;
				if    ( $phase == 0 ) { $Rph = 3; } 
				elsif ( $phase == 1 ) { $Rph = 2; }
				elsif ( $phase == 2 ) { $Rph = 1; }
			}
			else
			{
				if    ( $phase == 0 ) { $Rph = 1; }
				elsif ( $phase == 1 ) { $Rph = 3; }
				elsif ( $phase == 2 ) { $Rph = 2; }

				$phase = (3 - ( 3 - $phase + $len )%3)%3;
				if    ( $phase == 0 ) { $Lph = 3; }
				elsif ( $phase == 1 ) { $Lph = 2; }
				elsif ( $phase == 2 ) { $Lph = 1; }
			}

			$label = LabelForCDC($feature->[8]);

			$seq = substr( $genome{$feature->[0]}, $feature->[3] -1, $len );
			$seq = RevComp($seq) if ( $feature->[6] eq '-' );

			if ( $feature->[6] eq '+' )
			{
				if    ( $Lph == 1 ) {;}
				elsif ( $Lph == 2 ) { $seq = "n". $seq; }
				elsif ( $Lph == 3 ) { $seq = "nn". $seq; }

				if    ( $Rph == 1 ) { $seq .= "nn"; }
				elsif ( $Rph == 2 ) { $seq .= "n"; }
				elsif ( $Rph == 3 ) {;}
			}
			else
			{
				if    ( $Rph == 1 ) {;}
				elsif ( $Rph == 2 ) { $seq = "n". $seq; }
				elsif ( $Rph == 3 ) { $seq = "nn". $seq; }

				if    ( $Lph == 1 ) { $seq .= "nn"; }
				elsif ( $Lph == 2 ) { $seq .= "n"; }
				elsif ( $Lph == 3 ) {;}
			}

			PrintLineOut( "# $label", $gcount, $cds_length, $feature->[3], $feature->[4], $len, $feature->[6], $Lph, $Rph, $seq );

			# intergenic
			if (( $label eq "Initial" )or( $label eq "Single" )or( $label eq "Terminal" ))
			{
				$seq = '';

				if (( $label eq "Initial" )or( $label eq "Single" ))
				{
					$seq  = substr( $genome{$feature->[0]}, $feature->[3] -1 -$margin, $margin );
				}
				if (( $label eq "Single" )or( $label eq "Terminal" ))
				{
					$seq .= substr( $genome{$feature->[0]}, $feature->[4], $margin );
				}
				PrintLineOut( "# Intergenic", $gcount, $cds_length, $feature->[3], $feature->[4], $margin, $feature->[6], $Lph, $Rph, $seq );
			}

			# start
			if (( $label eq "Initial" )or( $label eq "Single" ))
			{
				if ( $feature->[6] eq '+' )
				{
					$signal = substr( $genome{$feature->[0]}, $feature->[3] -1,    3 );
					$seq    = substr( $genome{$feature->[0]}, $feature->[3] -1 -6, 3 +6 +3 );
				}
				else
				{
					$signal = substr( $genome{$feature->[0]}, $feature->[4] -1 -2,    3 );
					$signal = RevComp($signal);
					$seq    = substr( $genome{$feature->[0]}, $feature->[4] -1 -2 -3, 3 +6 +3 );
					$seq    = RevComp($seq);
				}

				if ($signal eq "ATG")
				{
					PrintLineOut( "# START", $gcount, $cds_length, $feature->[3], $feature->[4], $len, $feature->[6], $Lph, $Rph, $seq );
				}
				else { print "# $signal start $gid\n"; }
			}

			# stop
			if (( $label eq "Terminal" )or( $label eq "Single" ))
			{
				if ( $feature->[6] eq '+' )
				{
					$stop_label = LabelForStop( substr( $genome{$feature->[0]}, $feature->[4] -1 -2,     3 ));
					$seq                      = substr( $genome{$feature->[0]}, $feature->[4] -1 -2 - 3, 3 +6 +3 );
				}
				else
				{
					$stop_label = LabelForStop( RevComp(substr( $genome{$feature->[0]}, $feature->[3] -1, 3 )));
					$seq                              = substr( $genome{$feature->[0]}, $feature->[3] -1 -6 , 3 +6 +3 );
					$seq = RevComp($seq);
				}

				if ( $stop_label )
				{	
					PrintLineOut( "# $stop_label", $gcount, $cds_length, $feature->[3], $feature->[4], $len, $feature->[6], $Lph, $Rph, $seq );
				}
				else { print "# $seq stop $gid\n"; }
			}
		}

		if ( $feature->[2] eq "Intron" )
		{
#			print "# processing intron\n" if $verbose;
			
			my $phase = $feature->[7];

			if    ( $phase == 0 ) { $Lph=0; $Rph=0; }
			elsif ( $phase == 1 ) { $Lph=1; $Rph=1; }
			elsif ( $phase == 2 ) { $Lph=2; $Rph=2; }

			$seq = substr( $genome{$feature->[0]}, $feature->[3] -1, $len );
			$seq = RevComp($seq) if ( $feature->[6] eq '-' );
			PrintLineOut( "# Intron", $gcount, $cds_length, $feature->[3], $feature->[4], $len, $feature->[6], $Lph, $Rph, $seq );

			# Donor
			if ( $feature->[6] eq '+' )
			{
				$signal = substr( $genome{$feature->[0]}, $feature->[3] -1,    2 );
				$seq    = substr( $genome{$feature->[0]}, $feature->[3] -1 -3, 2 +3 +4 );
			}
			else
			{
				$signal = substr( $genome{$feature->[0]}, $feature->[4] -1 -1,    2 );
				$signal = RevComp($signal);
				$seq    = substr( $genome{$feature->[0]}, $feature->[4] -1 -1 -4, 2 +3 +4 );
				$seq    = RevComp($seq);
			}

			if ($signal eq "GT")
			{
				PrintLineOut( "# DON", $gcount, $cds_length, $feature->[3], $feature->[4], $len, $feature->[6], $Lph, $Rph, $seq );
			}
			else { print "# $signal don $gid\n"; }

			# Acceptor

			if ( $feature->[6] eq '+' )
			{
				$signal = substr( $genome{$feature->[0]}, $feature->[4] -1 -1,     2 );
				$seq    = substr( $genome{$feature->[0]}, $feature->[4] -1 -1 -18, 2 +18 +1 );
			}
			else
			{
				$signal = substr( $genome{$feature->[0]}, $feature->[3] -1,    2 );
				$signal = RevComp($signal);
				$seq    = substr( $genome{$feature->[0]}, $feature->[3] -1 -1, 2 +18 +1 ); 
				$seq    = RevComp($seq);
			}

			if ( $signal eq "AG" )
			{
				PrintLineOut( "# ACC", $gcount, $cds_length, $feature->[3], $feature->[4], $len, $feature->[6], $Lph, $Rph, $seq );
			}
			else { print "# $signal acc $gid\n"; }
		}		
	}

	@record = ();
}

close $OUT;

exit 0;
# ====================
my $in_gff3;
my $out_gff3;
my $tseq_out;
my $gseq_in;
my $no_gc;
my $get_sites;
my $tgff_out;
my $min_intron = 1;
my $adjust_cds;
my %info = ();
my $line;
my %meta;


# --------------------------------------------
sub PrintLineOut
{
	my @arr = @_;

	my $out_str = join( "\t", @arr );
	$out_str .= "\n";

	print $OUT $out_str;
}
# --------------------------------------------
sub LabelForCDC
{
	my $str = shift;
	my $label = '';

	if    ( $str =~ /\"Initial\"/ )  { $label = "Initial"; }
	elsif ( $str =~ /\"Internal\"/ ) { $label = "Internal"; }
	elsif ( $str =~ /\"Terminal\"/ ) { $label = "Terminal"; }
	elsif ( $str =~ /\"Single\"/ )   { $label = "Single"; }
	else { die "error, unexpected label found: $str\n"; }

	return $label;
}
# --------------------------------------------
sub LabelForStop
{
	my $str = shift;
	my $label = '';

	if    ($str eq "TAA") { $label = "STOP_TAA"; }
	elsif ($str eq "TAG") { $label = "STOP_TAG"; }
	elsif ($str eq "TGA") { $label = "STOP_TGA"; }

	return $label;
}
# --------------------------------------------
sub LoadGTF   #inuse
{
	my $fname = shift;
	my $ref = shift;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while( my $line = <$IN> )
	{
		# skip empty lines
		next if( $line =~ /^\s*$/ );
		# skip comments
		next if( $line =~ /^#/ );

		# get gene id	
		my $gid = GetGid($line);

		push @{$ref->{$gid}}, $line;
	}
	close $IN;

	print Dumper($ref) if $debug;
}
# --------------------------------------------
sub GetGid   #inuse
{
	my $str = shift;
	my $gene_id = '';

	if ( $str =~ /gene_id \"(\S+?)\"/ )
	{
		$gene_id = $1;
	}

	print "$gene_id\n" if $debug;

	return $gene_id;
}
# --------------------------------------------
sub GC
{
	my $s = shift;

	my $count_gc = 0;
	my $count_at = 0;

	$count_gc = () = $s =~ m/[gcGC]/g;
	$count_at = () = $s =~ m/[atAT]/g;

	my $gc = 0;

	if ( $count_at + $count_at == 0 )
	{
		print "# warning, no ATCG letters in sequence\n";
	}
	else
	{
		$gc = int( 100 * $count_gc / ($count_gc + $count_at));
	}

	return $gc;
}
# --------------------------------------------
sub LoadGenome
{
	my $fname = shift; # source
	my $ref = shift; # destination

	my $seqid = '';

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(my $line = <$IN>)
	{
		next if ( $line =~ /^#/ );
		next if ( $line =~ /^\s*$/ );

		if ( $line =~ /^\s*>/ )
		{
			if ( $line =~ /^>\s*(\S+)\s+/ )
			{
				$seqid = $1;
				print "# seqid: $seqid\n" if $verbose;
			}
			else
				{ die "error, unexpected defline format found: $line\n"; }
		}
		else
		{
			if ( ! $seqid )
				{ die "error, seqid is missing\n"; }

			$line =~ s/[\s0-9]//g;

			$line = uc $line;

			$ref->{$seqid} .= $line;
		}
	}
	close $IN;

	if ($verbose)
	{
		print "# genome sequence from file: $fname\n";
		print "# number of sequences: ". (scalar keys %{$ref}) ."\n";
	}
}
# --------------------------------------------
sub GetLength
{
	my $ref = shift;

	my $len = 0;

	foreach my $entry (@{$ref})
	{
		$len += ($entry->[4] - $entry->[3] + 1)
	}

	return $len;
}
# --------------------------------------------
sub CheckCanonical
{
	my $tid = shift;

	my $is_canonical = 1;

	$is_canonical = 0 if !$meta{$tid}{'is_cod'};

	$is_canonical = 0 if ! exists $meta{$tid}{'has_supported_start'};
	$is_canonical = 0 if ! exists $meta{$tid}{'has_supported_stop'};
	$is_canonical = 0 if ! exists $meta{$tid}{'has_supported_introns'};
	$is_canonical = 0 if ! exists $meta{$tid}{'has_gap'};

	if ( $is_canonical )
	{
		$is_canonical = 0 if !$meta{$tid}{'has_supported_start'};
		$is_canonical = 0 if !$meta{$tid}{'has_supported_stop'};
		$is_canonical = 0 if !$meta{$tid}{'has_supported_introns'};
		$is_canonical = 0 if !$meta{$tid}{'has_gap'};
	}

	$meta{$tid}{'is_canonical'} = $is_canonical;
}
# --------------------------------------------
sub CheckComplete
{
	my $tid = shift;

	my $is_complete = 1;

	$is_complete = 0 if !$meta{$tid}{'is_cod'};

	$is_complete = 0 if ! exists $meta{$tid}{'has_supported_start'};
	$is_complete = 0 if ! exists $meta{$tid}{'has_supported_stop'};
	$is_complete = 0 if ! exists $meta{$tid}{'change_in_phase'};

	if ( $is_complete )
	{
		$is_complete = 0 if !$meta{$tid}{'has_supported_start'};
		$is_complete = 0 if !$meta{$tid}{'has_supported_stop'};
		$is_complete = 0 if $meta{$tid}{'change_in_phase'};
	}

	$meta{$tid}{'is_complete'} = $is_complete;
}
# --------------------------------------------
sub CalcUTRs
{
	my $cds_ref = shift;
	my $exon_ref = shift;
	my $tid = shift;

	my %cds_to_exon_index = ();
	my %exon_to_cds_index = ();

	my $total_exons = scalar @{$exon_ref};
	my $total_cds = scalar @{$cds_ref};

	my $i = 1; # index+1 on exon arr
	my $j = 1; # index+1 on CDS arr

	while( $i < $total_exons + 1 )
	{
		if ( $exon_ref->[$i - 1][4] <  $cds_ref->[0][4] )
		{
			$exon_to_cds_index{$i} = 0;
			$i += 1;
		}
		elsif ( $exon_ref->[$i - 1][4] >= $cds_ref->[0][4] )
		{
			$exon_to_cds_index{$i} = 1;
			$cds_to_exon_index{1} = $i;
			$j += 1;
			$i += 1;
			last
		}
	}

	while( $j < $total_cds + 1)
	{
		if (  $exon_ref->[$i -1][3] == $cds_ref->[$j -1][3] )
		{
			$exon_to_cds_index{$i} = $j;
			$cds_to_exon_index{$j} = $i;
			$i+=1;
			$j+=1;
		}
		else
		{
			$i+=1;
			last;
		}
	}

	while( $i < $total_exons + 1 )
	{
		$exon_to_cds_index{$i} = 0;
		$i += 1;
	}

	my @arr = (sort{$a<=>$b} keys %cds_to_exon_index);

	# ----
	my $utr_L = 0;
	my $utr_R = 0;

	my $first_cds_idx = $arr[0];
	my $last_cds_idx = $arr[-1];
	my $first_exon_match = $cds_to_exon_index{$first_cds_idx};
	my $last_exon_match = $cds_to_exon_index{$last_cds_idx};

	for( my $idx = 1; $idx < $first_exon_match; $idx += 1 )
	{
		$utr_L += ( $exon_ref->[$idx -1][4] - $exon_ref->[$idx -1][3] + 1 );
	}
	$utr_L += ( $cds_ref->[$first_cds_idx -1][3] - $exon_ref->[$first_exon_match -1][3] );

	for( my $idx = $total_exons; $idx > $last_exon_match; $idx -= 1 )
	{
		$utr_R += ( $exon_ref->[$idx -1][4] - $exon_ref->[$idx -1][3] + 1 );
	}
	$utr_R += ( $exon_ref->[$last_exon_match -1][4] - $cds_ref->[$last_cds_idx -1][4] );

	if ( $cds_ref->[0][6] eq '+' )
	{
		$meta{$tid}{'utr5_length'} = $utr_L;
		$meta{$tid}{'utr3_length'} = $utr_R;
	}
	elsif ( $cds_ref->[0][6] eq '-' )
	{
		$meta{$tid}{'utr5_length'} = $utr_R;
		$meta{$tid}{'utr3_length'} = $utr_L;
	}
}
# --------------------------------------------
sub CheckInFrameStops
{
	my $tid = shift;

	my $utr5  = $meta{$tid}{'utr5_length'};
	my $utr3  = $meta{$tid}{'utr3_length'};
	my $total = $meta{$tid}{'t_length'};
	my $cds   = $meta{$tid}{'cod_length'};

	if ( $utr5 + $utr3 + $cds != $total )
	{
		print Dumper($meta{$tid});
		die "error, in parsing CDS on transcript level\n";
	}

	for( my $pos = $utr5; $pos < $utr5 + $cds - 3; $pos += 3 )
	{
		my $codon = substr( $meta{$tid}{'seq'}, $pos, 3 );
		if ($codon =~ /TAA|TAG|TGA/)
		{
			print "# in frame stop $meta{$tid}{'tid'}\n";
		}
	}
}
# --------------------------------------------
sub AdjustCDS
{
	my $ref = shift;
	my $tid = shift;

	$meta{$tid}{'change_in_phase'} = 0;

	if ( $ref->[0][6] eq '+' )
	{
		if ( $ref->[0][7] != 0 )
		{
			if ( $ref->[0][3] + $ref->[0][7] <= $ref->[0][4] )
			{
				$ref->[0][3] += $ref->[0][7];
				$ref->[0][7] = 0;
			}
			else
			{
				print "# warning, short cds". $meta{$tid}{'tid'} ."\n";
			}
			$meta{$tid}{'change_in_phase'} = 1;
		}
	}
	elsif ( $ref->[-1][6] eq '-' )
	{
		if ( $ref->[-1][7] != 0 )
		{
			if ( $ref->[-1][4] - $ref->[-1][7] >= $ref->[-1][3] )
			{
				$ref->[-1][4] -= $ref->[-1][7];
				$ref->[-1][7] = 0;
			}
			else
			{
				print "# warning, short cds". $meta{$tid}{'tid'} ."\n";
			}
			$meta{$tid}{'change_in_phase'} = 1;
		}
	}

	$meta{$tid}{'cod_length'} = GetLength($ref);

	my $mod = $meta{$tid}{'cod_length'} % 3;

	if ( $mod != 0 )
	{
		if ( $ref->[-1][6] eq '+' )
		{
			if ( $ref->[-1][4] - $mod >= $ref->[-1][3] )
			{
				$ref->[-1][4] -= $mod;
			}
			else
			{
				print "# warning, short cds". $meta{$tid}{'tid'} ."\n";
			}
			$meta{$tid}{'change_in_phase'} = 1;
		}
		elsif ( $ref->[0][6] eq '-' )
		{
			if ( $ref->[0][3] + $mod <= $ref->[0][4] )
			{
				$ref->[0][3] += $mod;
			}
			else
			{
				print "# warning, short cds". $meta{$tid}{'tid'} ."\n";
			}
			$meta{$tid}{'change_in_phase'} = 1;
		}
	}

	$meta{$tid}{'cod_length'} = GetLength($ref);
}
# -------------------------------------------
sub EnrichRecord
{
	# ref on array of arrays of GFF values from one gene
	my $ref = shift;

	my $use_seq = shift;

	# return enriched record
	my @new_rec = ();

	# with features
	my @gene = ();
	my @transcripts = ();
	my @exons = ();
	my @cds = ();
	my @introns = ();
	my @start_codon = ();
	my @stop_codon = ();

	# Find "gene" line and put it into @gene
	@gene = GetGene($ref);
	if ( (scalar @gene) == 0 ) { print Dumper($ref); die "error, no gene in record\n"; }
	my $gene_id = GetGeneID(\@gene);

	# prepare mrna lines
	@transcripts = SelectTranscripts($ref);
	if ( (scalar @transcripts) == 0 ) { print Dumper($ref); die "error, no transcript in record\n"; }
	AddCountToAttrSimple(\@transcripts);

	my %mrna = GetMrnaIDs(\@transcripts);
	SplitByMRNA( $ref, \%mrna );

	# strt output array
	push @new_rec, @gene;

	foreach my $key ( keys %mrna )
	{
		@exons = SelectExonsSorted($mrna{$key});
		@cds = SelectCDSsorted($mrna{$key});

		$meta{$key}{'gid'} = $gene_id;
		$meta{$key}{'tid'}= $key;
		$meta{$key}{'t_length'} = GetLength(\@exons);
		$meta{$key}{'cod_length'} = GetLength(\@cds);

		ParseAnnotationQuality($mrna{$key}, $key);
		ParsePseudo($mrna{$key}, $key);
		ParseGap($mrna{$key}, $key);

		foreach my $entry (@exons) { $entry->[8] = "Parent=". $key .";"; }
		foreach my $entry (@cds)   { $entry->[8] = "Parent=". $key .";"; }

		AddCountToAttr(\@exons);
		CreateTseq(\@exons, \%genome, $key) if $tseq_out;

		if ( @cds > 0 )
		{
			$meta{$key}{'is_cod'} = 1;
			
			AdjustCDS(\@cds, $key) if $adjust_cds;
			CalcUTRs( \@cds, \@exons, $key );

			@introns = CreateCdsIntrons(\@cds, $key);
			CompareIntrons( $ref, \@introns ) if $warnings;
			@start_codon = CreateStartCodon(\@cds, $key);
			@stop_codon = CreateStopCodon(\@cds, $key);

			foreach my $entry (@introns)      { $entry->[8] = "Parent=". $key .";"; }
			foreach my $entry (@start_codon)  { $entry->[8] = "Parent=". $key .";"; }
			foreach my $entry (@stop_codon)   { $entry->[8] = "Parent=". $key .";"; }

			if ( $use_seq )
			{
				my $are_splice_sites_supported =  AddSpliceSites( \@introns, \%genome );
				if ( $are_splice_sites_supported )
					{ $meta{$key}{'has_supported_introns'} = 1; }
				else
					{ $meta{$key}{'has_supported_introns'} = 0; }
				
				my $codon = '';

				$codon = AddStartStopSites( \@start_codon, \%genome );
				if ($codon eq 'ATG')
					{ $meta{$key}{'has_supported_start'} = 1; }
				else
					{ $meta{$key}{'has_supported_start'} = 0; }

				$codon = AddStartStopSites( \@stop_codon, \%genome );
				if ($codon =~ /TAA|TAG|TGA/)
					{ $meta{$key}{'has_supported_stop'} = 1; }
				else
					{ $meta{$key}{'has_supported_stop'} = 0; }

				CheckInFrameStops($key) if $tseq_out;
				CheckCanonical($key);
				CheckComplete($key);
			}

			AddLabelsToCDS(\@cds, $key);
			AddCountToAttr(\@cds);
			AddCountToAttr(\@introns);
			AddCountToAttrSemiReverse(\@start_codon, 1);
			AddCountToAttrSemiReverse(\@stop_codon, 0);

			push @new_rec, [ @{$mrna{$key}[0]} ];
			push @new_rec, @exons;
			push @new_rec, @cds;
			push @new_rec, @introns;
			push @new_rec, @start_codon;
			push @new_rec, @stop_codon;
		}
		else
		{
			$meta{$key}{'is_cod'} = 0;
			push @new_rec, [ @{$mrna{$key}[0]} ];
			push @new_rec, @exons;
		}
	}

	return @new_rec;
}
# --------------------------------------------
sub CompareIntrons
{
	my $annot = shift;
	my $calc = shift;

	my %h = ();

	# collect all annotated introns
	foreach my $entry ( @{$annot} )
	{
		if ( $entry->[2] =~ /[Ii]ntron/ )
		{
			my $key = $entry->[0] ."_". $entry->[3] ."_". $entry->[4] ."_". $entry->[6];
			$h{$key} = 0;
		}
	}

	$info{"annot_introns_in_input"} += scalar (keys %h);
	$info{"calc_introns_match_annot"} = 0;
	$info{"calc_introns_mismatch_annot"} = 0;

	foreach my $entry ( @{$calc} )
	{
		my $key = $entry->[0] ."_". $entry->[3] ."_". $entry->[4] ."_". $entry->[6];

		if ( exists $h{$key} )
		{
			$info{"calc_introns_match_annot"} += 1;
		}
		else
		{
			$info{"calc_introns_mismatch_annot"} += 1;
		}
	}

	$info{"annot_introns_mismatch_calc"} = $info{"annot_introns_in_input"} - $info{"calc_introns_match_annot"};

#	print Dumper(\%info);
}
# --------------------------------------------
sub AddStartStopSites
{
	my $ref = shift;
	my $h_genome = shift;

	my $scodon = '';

	foreach my $entry ( @{$ref} )
	{
		if ($warnings)
		{
			if ( $entry->[2] !~ /[Ss]tart_codon/ and $entry->[2] !~ /[Ss]top_codon/ )
				{ die "error, start/stop type is expected: $entry->[2]\n"; }
		}

		my $CODON = '';

		if ( $entry->[4] - $entry->[3] + 1 == 3 )
		{
			$CODON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 3 );
			$scodon = $CODON;
		}
		else
		{
			if ( $entry->[4] - $entry->[3] + 1 == 1 )
			{
				$CODON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 1 );
				$scodon .= $CODON;
			}
			elsif ( $entry->[4] - $entry->[3] + 1 == 2 )
			{
				$CODON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 2 );
				$scodon .= $CODON;
			}
			else
				{die;}
		}

		if ( $entry->[6] eq '+' )
		{
			;
		}
		elsif ( $entry->[6] eq '-' )
		{
			$CODON = RevComp($CODON);
		}
		else
			{ die "error, strand is required for sequence extraction:"; }

		$entry->[8] .= ("site_seq=". $CODON .";");
	}

	if ( $ref->[0][6] eq '-' )
	{
		$scodon = RevComp($scodon);
	}

	return uc($scodon);
}
# --------------------------------------------
sub AddSpliceSites
{
	my $ref = shift;
	my $h_genome = shift;

	my $is_supported = 1;

	foreach my $entry ( @{$ref} )
	{
		if ($warnings)
		{
			if ( $entry->[2] !~ /[Ii]ntron/ and $entry->[2] !~ /gap/ )
				{ die "error, intron type is expected: $entry->[2]\n"; }
		}

		my $DON = '';
		my $ACC = '';

		my $DONsite = '';
		my $ACCsite = '';

		if ( $entry->[6] eq '+' )
		{
			$DON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 2 );
			$ACC = substr( $h_genome->{$entry->[0]}, $entry->[4] -2, 2 );

			if ($get_sites)
			{
				$DON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1 -3,  2 +3  +4 );
				$ACC = substr( $h_genome->{$entry->[0]}, $entry->[4] -2 -18, 2 +18 +1 );
			}
		}
		elsif ( $entry->[6] eq '-' )
		{
			$ACC = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 2 );
			$DON = substr( $h_genome->{$entry->[0]}, $entry->[4] -2, 2 );

			$DON = RevComp($DON);
			$ACC = RevComp($ACC);

			if ($get_sites)
			{
				$ACCsite = substr( $h_genome->{$entry->[0]}, $entry->[3] -1 -1, 2 +18 +1);
				$DONsite = substr( $h_genome->{$entry->[0]}, $entry->[4] -2 -4, 2 +3  +4);

				$DONsite = RevComp($DONsite);
				$ACCsite = RevComp($ACCsite);
			}
		}
		else
			{ die "error, strand is required for sequence extraction:"; }

		$entry->[8] .= ("site_seq=". $DON ."_". $ACC .";");

		if ($get_sites)
		{
			$entry->[8] .= ("site=". $DONsite ."_". $ACCsite .";");
		}

		if ((uc($DON) ne 'GT')and(uc($DON) ne 'GC'))
		{
			$is_supported = 0;
		}

		if (uc($ACC) ne 'AG')
		{
			$is_supported = 0;
		}
	}

	return $is_supported;
}
# --------------------------------------------
sub RevComp
{
	my $s = shift;

	$s = reverse($s);

	$s =~ s/A/1/g;
	$s =~ s/T/A/g;
	$s =~ s/1/T/g;

	$s =~ s/C/2/g;
	$s =~ s/G/C/g;
	$s =~ s/2/G/g;

	$s =~ s/a/3/g;
	$s =~ s/t/a/g;
	$s =~ s/3/t/g;

	$s =~ s/c/4/g;
	$s =~ s/g/c/g;
	$s =~ s/4/g/g;

	return $s;
}
# --------------------------------------------
sub AddCountToAttrSimple
{
	my $ref = shift;

	my $size = scalar @{$ref};

	for( my $i = 0; $i < $size; $i += 1 )
	{
		$ref->[$i][8] .= ";" if ( $ref->[$i][8] !~ m/;\s*$/ );
		$ref->[$i][8] .= "count=". ($i+1) ."_". $size .";";
	}
}
# --------------------------------------------
sub AddCountToAttr
{
	my $ref = shift;

	my $size = scalar @{$ref};

	return if ( $size == 0 );

	if ( $ref->[0][6] eq "+" )
	{
		for( my $i = 0; $i < $size; $i += 1 )
		{
			$ref->[$i][8] .= ";" if ( $ref->[$i][8] !~ m/;$/ );
			$ref->[$i][8] .= "count=". ($i+1) ."_". $size .";";
		}
	}
	elsif ( $ref->[0][6] eq "-" )
	{
		for( my $i = $size -1; $i >= 0; $i -= 1 )
		{
			$ref->[$i][8] .= ";" if ( $ref->[$i][8] !~ m/;$/ );
			$ref->[$i][8] .= "count=". ($i+1) ."_". $size .";";
		}
	}
}
# --------------------------------------------
sub AddCountToAttrSemiReverse
{
	my $ref = shift;
	my $is_start = shift;

	my $size = scalar @{$ref};

	return if ( $size == 0 );

	if ( $size == 1 )
	{
		$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
		$ref->[0][8] .= "count=1_1;";
	}
	elsif ( $size == 2 )
	{
		$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
		$ref->[1][8] .= ";" if ( $ref->[1][8] !~ /;$/ );

		if (( $is_start and ($ref->[0][6] eq "+" )) or ( !$is_start and ($ref->[0][6] eq "-" )))
		{
			$ref->[0][8] .= "count=1_2;";
			$ref->[1][8] .= "count=2_2;";
		}
		elsif (( !$is_start and ($ref->[0][6] eq "+" )) or ( $is_start and ($ref->[0][6] eq "-" )))
		{
			$ref->[0][8] .= "count=2_2;";
			$ref->[1][8] .= "count=1_2;";
		}
	}
}
# --------------------------------------------
sub AddCountToAttrInHashMRNA
{
	my $ref = shift;

	my $size = scalar ( keys %{$ref} );
	my $i = 0;

	foreach my $key ( keys %{$ref} )
	{
		$ref->{$key}[0][8] .= ";" if ( $ref->{$key}[0][8] !~ m/;$/ );
		$ref->{$key}[0][8] .= "count=". ($i+1) ."_". $size .";";

		$i += 1;
	}
}
# --------------------------------------------
sub AddLabelsToCDS
{
	my $ref = shift;
	my $tid = shift;

	my $size = scalar @{$ref};

	if ( $size == 1 )
	{
		$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
		if ( !$meta{$tid}{'has_supported_start'} or !$meta{$tid}{'has_supported_stop'} )
		{
			$ref->[0][8] .= "cds_type=Partial;";
		}
		else
		{
			$ref->[0][8] .= "cds_type=Single;";
		}
	}
	elsif ( $size > 1 )
	{
		if ( $ref->[0][6] eq "+" )
		{
			$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
			if ( ! $meta{$tid}{'has_supported_start'} )
			{
				$ref->[0][8] .= "cds_type=Partial;";
			}
			else
			{
				$ref->[0][8] .= "cds_type=Initial;";
			}

			$ref->[$size -1][8] .= ";" if ( $ref->[$size -1][8] !~ /;$/ );
			if ( ! $meta{$tid}{'has_supported_stop'} )
			{
				$ref->[$size -1][8] .= "cds_type=Partial;";
			}
			else
			{
				$ref->[$size -1][8] .= "cds_type=Terminal;";
			}
		}
		elsif ( $ref->[0][6] eq "-" )
		{
			$ref->[$size -1][8] .= ";" if ( $ref->[$size -1][8] !~ /;$/ );
			if ( ! $meta{$tid}{'has_supported_start'} )
			{
				$ref->[$size -1][8] .= "cds_type=Partial;";
			}
			else
			{
				$ref->[$size -1][8] .= "cds_type=Initial;";
			}

			$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
			if ( ! $meta{$tid}{'has_supported_stop'} )
			{
				$ref->[0][8] .= "cds_type=Partial;";
			}
			else
			{
				$ref->[0][8] .= "cds_type=Terminal;";
			}
		}

		for( my $i = 1; $i < $size -1; $i += 1 )
		{
			$ref->[$i][8] .= ";" if ( $ref->[$i][8] !~ /;$/ );
			$ref->[$i][8] .= "cds_type=Internal;";
		}
	}
}
# --------------------------------------------
sub CreateStopCodon
{
	# ref on array of arrays
	my $ref = shift;
	my $tid = shift;

	# output
	my @arr = ();

	# first CDS from the left
	my @current = @{$ref->[0]};

	if ( $current[6] eq '+' )
	{
		my $idx = (scalar @{$ref}) - 1;

		# last CDS from the left
		@current = @{$ref->[$idx]};

		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "stop_codon";
			$current[3] = $current[4] - 3 + 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			# last CDS from the left

			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "stop_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 2;

				push @arr, [ @current ];

				# before last CDS from left
				@current = @{$ref->[$idx - 1]};

				$current[2] = "stop_codon";
				$current[3] = $current[4];
				$current[7] = 0;

				push @arr, [ @current ];
                        }
			elsif ( $current[4] - $current[3] + 1 == 1 )
			{
				$current[2] = "stop_codon";
				$current[3] = $current[4];
				$current[7] = 1;

				push @arr, [ @current ];

				@current = @{$ref->[$idx - 1]};

				$current[2] = "stop_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 0;

				push @arr, [ @current ];
			}

			if ($warnings)
			{
				print "warning, split stop codon detected: $meta{$tid}{'tid'}\n";
				print Dumper(\@arr) if $debug;
			}
		}
	}
	elsif ( $current[6] eq '-' )
	{
		# first CDS from the left

		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "stop_codon";
			$current[4] = $current[3] + 3 - 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "stop_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 2;

				push @arr, [ @current ];

				@current = @{$ref->[1]};

				$current[2] = "stop_codon";
				$current[4] = $current[3];
				$current[7] = 0;

				push @arr, [ @current ];
			}
			elsif ( $current[4] - $current[3] + 1 == 1 )
			{
				$current[2] = "stop_codon";
				$current[4] = $current[3];
				$current[7] = 1;

				push @arr, [ @current ];

				@current = @{$ref->[1]};

				$current[2] = "stop_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 0;

				push @arr, [ @current ];
			}

			if ($warnings)
			{
				print "warning, split stop codon detected: $meta{$tid}{'tid'}\n";
				print Dumper(\@arr) if $debug;
			}
		}
	}

	@arr = sort{ $a->[3] <=> $b->[3] } @arr;

	return @arr;
}
# --------------------------------------------
sub CreateStartCodon
{
	# ref on array of array - values of CDS from one mrna - sorted
	my $ref = shift;
	my $tid = shift;

	# output array of arrays
	my @arr = ();

	# this is first CDS on the left side
	my @current = @{$ref->[0]};

	if ( $current[6] eq '+' )
	{
		# complete start codon
		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "start_codon";
			$current[4] = $current[3] + 3 - 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			if ( @{$ref} < 2 )
				{ die "error, not enough data to position start codon\n"; }

			# to do : check strand

			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "start_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 0;

				push @arr, [ @current ];

				# mode to second CDS from left
				@current = @{$ref->[1]};

				$current[2] = "start_codon";
				$current[4] = $current[3];
				$current[7] = 1;

				push @arr, [ @current ];
			}
			elsif ( $current[4] - $current[3] + 1 == 1 )
			{
				$current[2] = "start_codon";
				$current[4] = $current[3];
				$current[7] = 0;

				push @arr, [ @current ];

				# mode to second CDS from left
				@current = @{$ref->[1]};

				$current[2] = "start_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 2;

				push @arr, [ @current ];
			}

			if ($warnings)
			{
				print "warning, split start codon detected: $meta{$tid}{'tid'}\n";
				print Dumper(\@arr) if $debug;
			}
		}
	}
	elsif ( $current[6] eq '-' )
	{
		my $idx = (scalar @{$ref}) - 1;

		# this is last CDS from the left
		@current = @{$ref->[$idx]};

		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "start_codon";
			$current[3] = $current[4] - 3 + 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			if ( @{$ref} < 2 )
				{ die "error, not enough data to position start codon\n"; }

			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "start_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 0;

				push @arr, [ @current ];

				# before last CDS from left
				@current = @{$ref->[$idx - 1]};

				$current[2] = "start_codon";
				$current[3] = $current[4];
				$current[7] = 1;

				push @arr, [ @current ];
                        }
                        elsif ( $current[4] - $current[3] + 1 == 1 )
                        {
                                $current[2] = "start_codon";
                                $current[3] = $current[4];
                                $current[7] = 0;

                                push @arr, [ @current ];

				# before last CDS from left
                                @current = @{$ref->[$idx - 1]};

                                $current[2] = "start_codon";
                                $current[3] = $current[4] - 1;
                                $current[7] = 2;

                                push @arr, [ @current ];
                        }

			if ($warnings)
			{
				print "warning, split start codon detected: $meta{$tid}{'tid'}\n";
				print Dumper(\@arr) if $debug;
			}
		}
	} 

	@arr = sort{ $a->[3] <=> $b->[3] } @arr;

	return @arr;
}
# --------------------------------------------
sub PrintRecord
{
	my $ref = shift;

	foreach my $entry ( @{$ref} )
	{
		print $OUT   $entry->[0] ."\t". $entry->[1] ."\t". $entry->[2] ."\t". $entry->[3] ."\t". $entry->[4] ."\t". $entry->[5] ."\t". $entry->[6] ."\t". $entry->[7];
		if ( defined $entry->[8] )
		{
			print $OUT "\t". $entry->[8];
		}
		print $OUT  "\n";
	}
}
# --------------------------------------------
sub CreateCdsIntrons
{
	# ref on array of arrys - CDS values from GFF file - one mrna - sorted
	my $ref = shift;
	my $tid = shift;

	my $size = scalar @{$ref};

	# ref on arry of arrays with introns - output
	my @arr = ();

	# two CDS minimum for intron deriviation
	return @arr if ($size < 2);

	my $i = 0;
	my $j = 1;

	while( $j < $size)
	{
		# two exons must be on the same strand

		if ( $ref->[$i][6] eq $ref->[$j][6] )
		{
			my @current = @{$ref->[$i]};

			$current[2] = "Intron";
			$current[3] = $ref->[$i][4] + 1;
			$current[4] = $ref->[$j][3] - 1;

			if ( $ref->[$i][6] eq "+" )
			{
				$current[7] = (3 - $ref->[$j][7])%3;
			}
			elsif ( $ref->[$i][6] eq "-" )
			{
				$current[7] = (3 - $ref->[$i][7])%3;
			}

			push @arr, [@current];
		}
		else
		{
			print "warning, oposite strand CDS were detected: intron is not assigned in such cases\n" if $warnings;
		}

		$i += 1;
		$j += 1;
	}

	return @arr;
}
# --------------------------------------------
sub SelectTranscripts
{
	my $ref = shift;
	my @arr = ();

	my %labels = ( "mRNA" => 1, "transcript" => 1);
	$labels{"nc_primary_transcript"} = 1;
       	$labels{"lnc_RNA"} = 1;
	$labels{"unconfirmed_transcript"} = 1;
	$labels{"pseudogene"} = 1;
	$labels{"tRNA"} = 1;
	$labels{"rRNA"} = 1;
	$labels{"ncRNA"} = 1;
	$labels{"snRNA"} = 1;
	$labels{"snoRNA"} = 1;
	$labels{"pre_miRNA"} = 1;
	$labels{"antisense_lncRNA"} = 1;
	$labels{"miRNA_primary_transcript"} = 1;
	$labels{"miRNA"} = 1;
	$labels{"pseudogenic_transcript"} = 1;
	$labels{"ncRNA_gene"} = 1;
	$labels{"scRNA"} = 1;
	$labels{"transcript_region"} = 1;
	$labels{"antisense_RNA"} = 1;

	foreach my $entry (@{$ref})
	{
		if (( exists $labels{ $entry->[2] } )or( $entry->[2] =~ /^\S_gene_segment$/ ))
		{
			push @arr, [ @{$entry} ];
		}
	}

	if ( scalar @arr == 0 )
	{
		print Dumper($ref);
		die "error, lines with transcript info not found\n";
	}

	return @arr;
}
# --------------------------------------------
sub SelectExonsSorted
{
	my $ref = shift;
	my @arr = ();

	foreach my $entry (@{$ref})
	{
		if (( $entry->[2] eq "exon" )or( $entry->[2] eq "pseudogenic_exon" )or( $entry->[2] eq "miRNA_" ))
		{
			push @arr, [ @{$entry} ];
		}
	}

	@arr = sort{ $a->[3] <=> $b->[3] } @arr;

	return @arr;
}
# --------------------------------------------
sub SelectCDSsorted
{
	my $ref = shift;
	my @arr = ();

	foreach my $entry (@{$ref})
	{
		if ( $entry->[2] eq "CDS" )
		{
			if ( $entry->[6] !~ /^[+-]$/ )
				{ die "error, CDS strand value is missing: $entry->[6]\n"; } 

			push @arr, [ @{$entry} ];
		}
	}

	# is it possible to have trans-splicing from different chromosomes?
	
	@arr = sort{ $a->[0] cmp $b->[0] || $a->[3] <=> $b->[3] } @arr;

	if ( $warnings and (@arr > 0))
	{
		my $i = 0;
		my $j = 1;
		my $size = scalar @arr;

		while( $j < $size )
		{
			if ( $arr[$i][4] >= $arr[$j][3] )
			{
				print "warning, two CDS from the same mRNA overlap: $arr[$i][3] .. $arr[$i][4]  $arr[$j][3] .. $arr[$j][4]\n"; 
			}

			$i += 1;
			$j += 1;
		}

		my $strand = $arr[0][6];
		my $seqid = $arr[0][0];

		foreach my $current (@arr)
		{
			if ( $strand ne $current->[6] )
			{
				print "warning, two strands detected in one mRNA: $strand  $current->[6]\n";
				last;
			}	
		}

		foreach my $current (@arr)
		{
			if ( $seqid ne $current->[0] )
			{
				print "warning, two seqid detected in one mRNA: $seqid  $current->[0]\n";
				last;
			}
		}
	}

	return @arr;
}
# --------------------------------------------
sub SplitByMRNA
{
	# ref array of arrays with GGF values of one gene
	my $ref = shift;
	# ref on hash of arrays - with GFF values separated by mrna ID
	my $h_mrna = shift;

	my %labels = ( "mRNA" => 1, "transcript" => 1 );
	$labels{"nc_primary_transcript"} = 1;
	$labels{"lnc_RNA"} = 1;
	$labels{"unconfirmed_transcript"} = 1;
	$labels{"gene"} = 1;
	$labels{"stop_codon_redefined_as_selenocysteine"} = 1;
	$labels{"pseudogene"} = 1;
	$labels{"tRNA"} = 1;
	$labels{"rRNA"} = 1;
	$labels{"ncRNA"} = 1;
	$labels{"snRNA"} = 1;
	$labels{"snoRNA"} = 1;
	$labels{"pre_miRNA"} = 1;
	$labels{"antisense_lncRNA"} = 1;
	$labels{"miRNA_primary_transcript"} = 1;
	$labels{"ncRNA_gene"} = 1;
	$labels{"miRNA"} = 1;
	$labels{"pseudogenic_transcript"} = 1;
	$labels{"scRNA"} = 1;
	$labels{"transcript_region"} = 1;
	$labels{"antisense_RNA"} = 1;
	$labels{"transposable_element_gene"} = 1;

	foreach my $entry (@{$ref})
	{
		next if ( exists $labels{ $entry->[2] } );
		next if ( $entry->[2] =~ /^\S_gene_segment$/ );

		if ( $entry->[8] =~ /Parent=(\S+?);/ or $entry->[8] =~ /Parent=(\S+)/)
		{
			my @parent_set = split( ',', $1 );

			foreach my $value (@parent_set)
			{
				if ( exists $h_mrna->{$value} )
				{
					push @{$h_mrna->{$value}}, [ @{$entry} ];
				}
				else
				{
					print Dumper($ref) if $debug;
					die "error, feature not in mRNA:\n$value\n";
				}
			}
		}
		else
			{ die "error, line without the Parent in:\n$entry->[8]\n"; }
	}
}
# --------------------------------------------
sub GetMrnaIDs
{
	# ref on array of arrays - GFF values one gene
	my $ref = shift;

	# one or many mRNA per record 
	my %mRNA = ();

	my %labels = ( "mRNA" => 1, "transcript" => 1 );
	$labels{"nc_primary_transcript"} = 1;
	$labels{"lnc_RNA"} = 1;
	$labels{"unconfirmed_transcript"} = 1;
       	$labels{"pseudogene"} = 1;
	$labels{"tRNA"} = 1;
	$labels{"rRNA"} = 1;
	$labels{"ncRNA"} = 1;
	$labels{"snRNA"} = 1;
	$labels{"snoRNA"} = 1;
	$labels{"pre_miRNA"} = 1;
	$labels{"antisense_lncRNA"} = 1;
	$labels{"miRNA_primary_transcript"} = 1;
	$labels{"lnc_RNA"} = 1;
	$labels{"miRNA"} = 1;
	$labels{"pseudogenic_transcript"} = 1;
	$labels{"scRNA"} = 1;
	$labels{"transcript_region"} = 1;
	$labels{"antisense_RNA"} = 1;

	foreach my $entry (@{$ref})
	{
		if (( exists $labels{ $entry->[2] } ) or ( $entry->[2] =~ /^\S_gene_segment$/))
		{
			if ( $entry->[8] =~ /ID=(\S+?);/ )
			{
				my $ID = $1;

				if ( ! exists $mRNA{$ID} )
				{
					push @{$mRNA{$ID}}, [ @{$entry} ];
				}
				else
					{ die "error, mRNA ID duplication found: $ID\n"; }
			}
			else
				{ die "error, mRNA ID field not found in: $entry->[8]\n"; }
		}
	}

	if ( scalar (keys %mRNA ) == 0 )
	{
		print Dumper($ref); # if $debug;
		die "error, no mRNA in record\n";
	}

	return %mRNA;
}
# --------------------------------------------
sub GetGene
{
	my $ref = shift;
	
	my @arr = ();

	foreach my $entry (@{$ref})
	{
		if (( $entry->[2] eq "gene" )or( $entry->[2] eq "ncRNA_gene" )or( $entry->[2] eq "pseudogene_" )or( $entry->[2] eq "transposable_element_gene" ))
		{
			 push @arr, [@{$entry}];
		}
	}

	return @arr;
}
# --------------------------------------------
sub GetGeneID
{	
	my $ref = shift;

	my $gene_id = '';

	foreach my $entry (@{$ref})
	{
		if ( !$gene_id )
		{
			if ( $entry->[8] =~ /ID=(\S+?);/ ) 
			{
				$gene_id = $1;
			}
			else
				{ die "error, gene ID field not found in: $entry->[8]\n"; }
		}
		else
			{ die "error, gene entry duplication was detected in record: $gene_id\n"; }
	}

	if ( ! $gene_id )
		{ die "error, gene id is missing\n"; }

	return $gene_id;
}
# --------------------------------------------
sub CDSlength   #inuse
{
	# ref on array of arrays of GFF values from one gene
	my $ref = shift;

	my $cds_length = 0;

	my $partial = 0;

	foreach my $entry (@{$ref})
	{
		if ( $entry->[2] eq "CDS" )
		{
			$cds_length += ($entry->[4] - $entry->[3] + 1);

			if ( $entry->[8] =~ /\"Partial\"/ )
			{
				$partial = 1;
			}
		}
	}

	$cds_length = 0 if $partial;

	return $cds_length;
}
# --------------------------------------------
sub CheckForValidGFF   #inuse
{
	my $str = shift;
	my $ref = shift;

	# seqid 1
	if ( $ref->[0] !~ /^\w+$/ )
		{ die "error, unexpected seqid format found:\n$str\n"; }

	# start 4
	if ( $ref->[3] !~ /^\d+$/ )
		{ die "error, unexpected start format found:\n$str\n"; }

	# end 5
	if ( $ref->[4] !~ /^\d+$/ )
		{ die "error, unexpected end format found:\n$str\n"; }

	# start <= end
	if ( $ref->[3] > $ref->[4] )
		{ die "error, start is more than end:\n$str\n"; }

	# strand 6
	if ( $ref->[6] !~ /^[+-]$/ )
		{ die "error, wrong strand value:\n$str\n"; }

	# phase 7
	if ( $ref->[7] !~ /^[012]$/ )
		 { die "error, wrong phase value:\n$str\n"; }
}
# --------------------------------------------
sub AddLineToRecord   #inuse
{
	# str - line from GFF file
	# ref on array of arrays - put all split lines from one gene here
	my $str = shift;
	my $ref = shift;

	chomp $str;

	my @arr = split( '\t', $str );

	my $size = @arr;

	if ( $size != 8 and $size != 9 )
		{ die "error, unexpected number of TABs found:\n$str\n"; }

	CheckForValidGFF( $str, \@arr ) if $debug;

	push @{$ref}, [@arr];
}
# ------------------------------------------------
sub CheckBeforeRun   #inuse
{
	die "error, file not found: option --in $gtf\n" if( ! -e $gtf );
	die "error, file not found: option --gseq $gseq\n" if( ! -e $gseq );
	die "error, output file name matches input file: $gtf $gseq $out\n" if (( $out eq $gtf )or( $out eq $gseq));
}
# ------------------------------------------------
sub ParseCMD   #inuse
{
	my $opt_results = GetOptions
	(
		'gtf=s'    => \$gtf,
		'gseq=s'   => \$gseq,
		'out=s'    => \$out,

		'verbose'       => \$verbose,
		'warnings'      => \$warnings,
		'debug'         => \$debug,
        );

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line\n" if( @ARGV > 0 );

	$verbose = 1 if $debug;
	$warnings = 1 if $debug;
}
# ------------------------------------------------
sub Usage   #inuse
{
        print qq(
Usage: $0 --gtf [name] --gseq [name] --out [name]

  --verbose
  --warnings
  --debug

Version: $VERSION

);
	exit 1;
}
# ------------------------------------------------

