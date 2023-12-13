#!/usr/bin/env perl
# -----------------
# Alex Lomsadze
# 2023
#
# check GeneMark-ETP instalation
# ----------------

use strict;
use warnings;

print "### Checking ETP installation\n";

CheckPerlModules();
CheckGeneMark_hmm();
Check_p_GeneMark_hmm();
CheckExe();
CheckPython3();
CheckThirdParty();
RunHmm();
RunWithProtHint();

print "### Checking ETP installation ... done\n\n";

# ----------------
sub RunWithProtHint
{
	print "# Checking prediction with ProtHint\n";

	my $path = "bin/gmes/GeneMark-E-tests/EP";
	if ( -d $path )
	{
		mkdir "$path/test";
		chdir "$path/test";

		if ( -e "genemark.gtf" )
		{
			unlink "genemark.gtf";
		}

		system("../../../gmes_petap.pl --seq ../input/genome.fasta --EP --dbep ../input/proteins.fasta --cores=8 --max_intergenic 10000 --mask_penalty 0 > loginfo 2>&1 ");

		if ( ! -e "genemark.gtf" )
			{ die "error, EP run test failed\n"; }

		system("../../../compare_intervals_exact.pl  --f1 genemark.gtf  --f2 ../output/genemark.gtf  --v" );	
	}
	else
		{ print "warning, EP test folder is missing: $path\n"; }

	print "# Checking prediction with ProtHint ... done\n";
}
# ----------------
sub RunHmm
{
	print "# Checking prediction by hmm\n";

	my $path = "bin/gmes/GeneMark-E-tests/GeneMark.hmm";
	if ( -d $path )
	{
		mkdir "$path/test";

		if ( -e "$path/test/genemark.gtf" )
		{
			unlink "$path/test/genemark.gtf";
		}

		system( "bin/gmes/gmhmme3 -o $path/test/genemark.gtf -m $path/input/athaliana.mod -f gtf $path/input/sequence.fasta" );

		if ( ! -e "$path/test/genemark.gtf" )
			{ die "error, GeneMark.hmm run test failed\n"; }

		system( "bin/compare_intervals_exact.pl --f1 $path/output/genemark.gff3 --f2 $path/test/genemark.gtf --v" );
	}
	else
		{ print "warning, GeneMark.hmm test folder is missing: $path\n"; }

	print "# Checking prediction by hmm ... done\n";
}
# ----------------
sub CheckThirdParty
{
	print "# Checking third party tools\n";

	if ( ! -d "tools" )
		{ die "error, folder 'tools' not found\n"; }

	chdir "tools";
	system("./check_ETP_tools.pl");
	chdir "../";

	print "# Checking third party tools ... done\n";
}
# ----------------
sub CheckPython3
{
	print "# Checking Python3\n";

	my $message = `/usr/bin/env python3 --version`;

	if ( $message !~ /^Python 3/ )
		{ die "error, check for python3\n $message"; }

	print "# Checking Python3 ... done\n";
}
# ----------------
sub CheckExe
{
	print "# Checking exec\n";

	my $message = `"bin/probuild"`;
	if ( $message !~ /Version: 2/ )
		{ die "error, check bin/probuild\n $message"; }

	$message = `bin/bam2hints -h`;
	if ( $message !~ /Usage:/ )
		{ die "error, check bin/bam2hints\n $message"; }

	$message = `bin/gmes/ProtHint/dependencies/diamond version`;
	if ( $message !~ /version/ )
		{ die "error, check bin/gmes/ProtHint/dependencies/diamond\n $message"; }

	$message = `bin/gmes/ProtHint/dependencies/spaln 2>&1 | grep version`;
	if ( $message !~ /version/ )
		{ die "error, check bin/gmes/ProtHint/dependencies/spaln\n $message"; }

	$message = `bin/gmes/ProtHint/dependencies/spaln_boundary_scorer 2>&1`;
	if ( $message !~ /Usage:/ )
		{ die "error, check bin/gmes/ProtHint/dependencies/spaln_boundary_scorer\n $message"; }

	$message = `bin/gmes/Gibbs3`;
	if ( $message !~ /Gibbs 3/ )
		{ die "error, check bin/gmes/Gibbs3\n $message"; }

	print "# Checking exec ...  done\n";
}
# ----------------
sub Check_p_GeneMark_hmm
{
	print "# Checking P GeneMark.hmm\n";

	if (( ! -e "bin/gmst/gmhmmp" )or( ! -x "bin/gmst/gmhmmp"))
	{
		die "error, file gmhmmp is not an executable or does not exist\n";
	}

	my $message = `bin/gmst/gmhmmp`;

	if ( $message =~ /License key/ )
	{
		die "error, installation key is missing or expired\n";
	}

	if ( $message =~ /cannot execute/ )
	{
		die "error, mismatch between the binary and OS was detected\n";
	}

	print "# Checking P GeneMark.hmm ... done\n";
}
# ----------------
sub CheckGeneMark_hmm
{
	print "# Checking GeneMark.hmm\n";

	if (( ! -e "bin/gmes/gmhmme3" )or( ! -x "bin/gmes/gmhmme3"))
	{
		die "error, file gmhmme3 is not an executable or does not exist\n";
	}

	my $message = `"bin/gmes/gmhmme3"`;

	if ( $message =~ /License key/ )
	{
		die "error, installation key is missing or expired\n";
	}

	if ( $message =~ /cannot execute/ )
	{
		die "error, mismatch between the binary and OS was detected\n";
	}
	
	print "# Checking GeneMark.hmm ... done\n";
}
# ----------------
sub CheckPerlModules
{
        print "# Checking for Perl modules\n";

        my @modules = (
		"Cwd",
		"Data::Dumper",
		"File::Path",
		"File::Spec",
		"File::Temp",
		"FindBin",
		"Getopt::Long",
		"Hash::Merge",
		"List::Util",
		"MCE::Mutex",
		"Math::Utils",
		"Parallel::ForkManager",
		"Statistics::LineFit",
		"Storable",
		"Thread::Queue",
		"YAML",
		"YAML::XS",
		"threads",
	);

        foreach my $module (@modules)
        {
                my $respond = `perl -M$module -e 1`;
                print "$respond\n" if ( $respond =~ /\S/ );
        }

        print "# Checking for Perl modules ... done\n";
}
# ----------------

