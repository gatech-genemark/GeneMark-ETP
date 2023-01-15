#!/usr/bin/env perl
# --------------------------------------------
# Alex Lomsadze
# GaTech
# January 2022
#
# This is ETP variant of old script
# --------------------------------------------

use strict;
use warnings;

use Getopt::Long;

# --------------------------------------------
if ( $#ARGV == -1 ) { print PrintUsage(); exit 1; }

my $in = '';
my $out = '';
my $list = '';
my $v = '';
my $col = 1;
my $swap = 0;

my $result = GetOptions
(
  'in=s'   => \$in,
  'out=s'  => \$out,
  'list=s' => \$list,
  'verbose'=> \$v,
  'col=i'  => \$col,
  'swap'   => \$swap,
);

if ( !$result ) { print "error on cmd"; print PrintUsage(); exit 1; }
if ( !$in or !$out or !$list ) { die "error on cmd: --out or --in or --list are missing\n"; }
if (( $in eq $out ) or ( $list eq $out)) { die "error on cmd: --out equals --in or --list\n"; }

# --------------------------------------------
my %h;
LoadList( $list, \%h );

# --------------------------------------------
my $count_in = 0;
my $count_out = 0;

my %allowed =
(
	"gene" => 1,
	"mRNA" => 1,
	"CDS" => 1,
	"exon" => 1,
	"start_codon" => 1,
	"stop_codon" => 1,

	"Repeat" => 1,
	"repeat" => 1,
	"similarity" => 1,
	"match" => 1,

	"five_prime_UTR" => 1,
	"three_prime_UTR" => 1,
	"intron" => 1,
	"Intron" => 1,
	"transposable_element" => 1,

	"pseudogenic_transcript" => 1,
	"piRNA" => 1,
	"lincRNA" => 1,
	"miRNA_primary_transcript" => 1,
	"nc_primary_transcript" => 1,
	"primary_transcript" => 1,
	"pseudogenic_rRNA" => 1,
	"pseudogenic_tRNA" => 1,
	"scRNA" => 1,
	"guide_RNA" => 1,
	"SRP_RNA" => 1,
	"telomerase_RNA" => 1,

	"lnc_RNA" => 1,
	"antisense_lncRNA" => 1,
	"transcript_region" => 1,
	"antisense_RNA" => 1,
	"transposable_element_gene" => 1,

	"ncRNA" => 1,
	"circular_ncRNA" => 1,
	"snoRNA" => 1,
	"pseudogene" => 1,
	"snRNA" => 1,
	"tRNA" => 1,
	"rRNA" => 1,
	"pre_miRNA" => 1,
	"miRNA" => 1,

	"ncRNA_gene" => 1,
	"unconfirmed_transcript" => 1,
	"C_gene_segment" => 1,
	"D_gene_segment" => 1,
	"J_gene_segment" => 1,
	"V_gene_segment" => 1,

	"transcript" => 1,
	"tRNA_gene" => 1,
	"lincRNA_gene" => 1,
	"lincRNA" => 1,
	"miRNA_gene" => 1,
	"miRNA" => 1,
	"RNase_P_RNA" => 1,
	"RNase_MRP_RNA" => 1,
	"Y_RNA" => 1,

	"stop_codon_redefined_as_selenocysteine" => 1,
	"gap" => 1,

	"golden_path" => 1,
	"repeat_region" => 1,
	"TSS" => 1,
	"exon_junction" => 1,
	"origin_of_replication" => 1,
	"TF_binding_site" => 1,
	"RNAi_reagent" => 1,
	"tandem_repeat" => 1,
	"BAC_cloned_genomic_insert" => 1,
	"orthologous_to" => 1,
	"paralogous_to" => 1,
	"polypeptide" => 1,
	"polypeptide_region" => 1,
	"chromosome_breakpoint" => 1,
	"region" => 1,
	"transposable_element_insertion_site" => 1,
	"PCR_product" => 1,
	"polyA_site" => 1,
	"oligo" => 1,
	"rescue_region" => 1,
	"delins" => 1,
	"syntenic_region" => 1,
	"orthologous_region" => 1,
	"regulatory_region" => 1,
	"modified_RNA_base_feature" => 1,
	"sequence_alteration" => 1,
	"point_mutation" => 1,
	"sgRNA" => 1,
	"insulator" => 1,
	"chromosome_band" => 1,
	"deletion" => 1,
	"insertion" => 1,
	"insertion_site" => 1,
	"MNV" => 1,
	"complex_substitution" => 1,
	"mature_protein_region" => 1,
	"sequence_variant" => 1,
	"protein_binding_site" => 1,
	"DNA_motif" => 1,

	"protein" => 1,
	"transposon_fragment" => 1,
	"pseudogenic_exon" => 1,
	"uORF" => 1,

	"chromosome" => 1,
	"biological_region" => 1,

	"sequence_feature" => 1,
);

my %exclude = (
	"chromosome_breakpoint" => 1,
);

my %out_only = (

	"mRNA" => 1,
	"gene" => 1,
	"exon" => 1,
	"CDS" => 1,
	"start_codon" => 1,
	"stop_codon" => 1,
	"transcript" => 1,
	"gap" => 1,
	"pseudogene" => 1,
	"ncRNA" => 1,
	"snoRNA" => 1,
	"snRNA" => 1,
	"tRNA" => 1,
	"rRNA" => 1,
	"pre_miRNA" => 1,
	"lnc_RNA" => 1,
	"antisense_lncRNA" => 1,
	"transcript_region" => 1,
	"antisense_RNA" => 1,
	"transposable_element_gene" => 1,
	"repeat" => 1,
	"ncRNA_gene" => 1,
	"miRNA" => 1,
	"pseudogenic_transcript" => 1,
	"scRNA" => 1,
	"V_gene_segment" => 1,
	"miRNA_primary_transcript" => 1,
	"piRNA" => 1,
	"pseudogenic_tRNA" => 1,
	"primary_transcript" => 1,
	"pseudogenic_rRNA" => 1,
);

my %found;

open( my $IN, $in ) or die "Can't open $in: $!\n";
open( my $OUT, ">", $out ) or die "Can't open $out: $!\n";

while ( my $record = <$IN> )
{
	next if ( $record =~ /^\s*$/ );

	++$count_in;

	if ( $record =~ /^##gff-version\s+3/ )
	{
		print $OUT $record;
		++$count_out;
		next;
	}

	if ( $record =~ /^##sequence-region\s+(\S+)\s*/ )
	{
		if ( exists $h{$1} )
		{
			print $OUT $record;
			 ++$count_out;
			next;
		}
	}

	if ( $record =~/^(\S+)\t[^\t]+\t(\S+)\t/ )
	{
		my $id = $1;
		my $type = $2;

		die "error, unexpected type was found: $type\n" if ( ! exists $allowed{$type} );
		next if ( ! exists $out_only{$type} );
		next if ( exists $exclude{$type} );

		if ( exists $h{$id} )
		{
			$record =~ s/^$id/$h{$id}/;

			if ( $record =~ /ath-miR/ )
			{
				print "# warning, excludeed record $record\n";
				next;
			}

			if ( $record =~ /FBtr00840|FBtr03077/ )
			{
				if ( $record =~ /FBtr0084066|FBtr0084067|FBtr0084068|FBtr0084069|FBtr0084070|FBtr0084071|FBtr0084072|FBtr0084073|FBtr0084074|FBtr0084075|FBtr0084077|FBtr0084078|FBtr0084079|FBtr0084080|FBtr0084081|FBtr0084082|FBtr0084083|FBtr0084084|FBtr0084085|FBtr0084060|FBtr0084061|FBtr0084062|FBtr0084063|FBtr0084064|FBtr0084065|FBtr0114359|FBtr0114360|FBtr0114361|FBtr0307759|FBtr0307760|FBtr0343765/ )
				{
					if ( $record =~ /\t\+\t/ )
					{
						$record =~ s/\t\+\t/\t\-\t/;

						print "# warning, modified record with $record\n";
					}
				}
			}

			print $OUT $record;
			++$count_out;
			$found{$id} +=1 ;
		}
		else
		{
			# skip this line - seqid is not in the list
			# print "$id $type\n";
		}
	}
	else
	{
		if ( $record =~ /\t/ )
		{
			print $record;
		}
	}
}

close $IN;
close $OUT;

if ($list)
{
	CheckAllFound( \%found );
}

print "$count_in in $in file\n";
print "$count_out in $out file\n";

# --------------------------------------------
sub CheckAllFound
{
	my $ref = shift;

	for my $key (keys %$ref )
	{
		if ( ! exists $ref->{$key} )
		{
			print "warning, record was not found: $ref->{$key}\n";
		}
	}	
}
# --------------------------------------------
sub LoadList
{
	my $fname = shift;
	my $ref = shift;
	
	die "error, file name is empty\n" if (!$fname);
	
	open( my $F, $fname ) or die "error on  open $fname: $!\n";
	while( my $line = <$F> )
	{
		next if ( $line =~ /^\s*$/ );
		next if ( $line =~ /^#/ );
		
		if ( $line =~ /^\s*(\S+)\s+(\S+)\s*/ )
		{
			my $col_a = $1;
			my $col_b = $2;

			if ( $col == 1 )
			{
				;
			}
			elsif ( $col == 2 )
			{
				$col_a = $2;
				$col_b = $1;
			}
			else
				{ die "error, unsupported value for column was provided: $col\n"; }

			if ( exists $ref->{$col_a} )
			{
				die "error, duplicated entry found in the list: $fname\n";
			}

			if ( $swap )
			{
				$ref->{$col_a} = $col_b;
			}
			else
			{
				$ref->{$col_a} = $col_a;
			}
		}
		else
		{
			print "warning, unexpected format was found in $fname: $line";
		}
	}
	close $F;

	if ( scalar keys %$ref < 1 )
	{
		die "error, list is empty: $fname\n";
	}

	if ( $v )
	{
		for my $key (keys %$ref)
		{
			print $key ." ". $ref->{$key} ."\n";
		}
	}
}
# --------------------------------------------
sub PrintUsage
{
	my $txt = "Usage: $0  --in <file name>  --out <file name>  --list <file name>

  This program takes as input GFF formatted anotation and
  outputs GFF records only for sequence ID's listed on the --list file.

  List file may have two columns
  --list <file name> : name_1 name_2
  This specifies which column to use as
  --col  <column>    :   1  or  2
  --swap             : swap one name to another
";
	return $txt;
}
# --------------------------------------------

