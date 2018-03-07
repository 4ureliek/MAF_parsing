#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta - https://github.com/4ureliek
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  parsing a .maf (multispecies alignment file) to get total alignment length for the blocks containing the species of interest
##########################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::AlignIO;

my $version = "1.0";
my $changelog = "
#	- v1.0 = 21 Jul 2015
\n";

my $usage = "\nUsage [v$version]: 
	perl MAF_microdel--1--get-gaps.pl -in <dir> -sp <X> [-v] [-chlog] [-h]
	
	This script parse .maf file(s) (multispecies alignment file) and get alignment length for each species.
	It can filter for blocks containing info for X species (see -sp).
	Not needed for the pipeline MAFmicrodel but can be useful.
	
    MANDATORY ARGUMENTS:	
     -in     => (STRING) directory containing input .maf file(s). 
     
    OPTIONAL ARGUMENTS
     -sp     => (INT)    Number of species retained to filter blocks that have no info for any species
                         Files needs to be previously grepped for species of interest beforehand if this option is chosen
                         For example: grep '^\$\|maf\|^a\|hg19\|panTro4\|ponAbe2\|rheMac3\' chr_all.maf > chr_all.primates.maf
     -v      => (BOOL)   verbose mode, make the script talks to you
     -v      => (BOOL)   print version if only option
     -chlog  => (BOOL)   print change log (updates)
     -h|help => (BOOL)   print this usage\n\n";  


my ($in,$nbsp,$help,$v,$chlog);
GetOptions ('in=s' => \$in, 'sp=s' => \$nbsp, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script MAF_microdel--1--get-gaps.pl version $version\n\n" if ((! $in) && (! $nbsp) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die "\nERROR: Please provide an input file (see usage below)\n$usage" if (! $in);
die "\n Directory with input files ($in) does not exist\n\n" if (! -e $in);
die "\n Value for -sp ($nbsp) needs to be an integer\n\n" if (($nbsp) && ($nbsp !~ /[0-9]+/));

$in =~ s/\/$//; #remove / at the end of directory
if ($v) {
	print STDERR "\n --- MAF_microdel--1--get-aln-len.pl version $version started, with:\n";
	print STDERR "      - Input files in $in\n";
	print STDERR "      - Number of species per block: $nbsp\n" if ($nbsp);
	print STDERR "\n --- Looping through file(s) to obtain alignment length:\n";
}

my @in = `ls $in/*.maf | grep -v $in/*.maf.*.maf`;
FILE: foreach my $file (@in) {
	chomp ($file);
	next FILE unless ($file);
	#####################################################
	# Get MAF files containing only blocks with these X species if not done yet
	#####################################################
	print STDERR " --- $file\n" if ($v);
	my %totlen = ();
	my @spIDs = ();
	my $outfile = "$file.$nbsp.maf";
	unless (-e $outfile) {
		print STDERR "      - $outfile does not exists => get new MAF file containing only blocks with these $nbsp species\n" if ($v);
		open (my $infh, "<$file") or confess  "\t    ERROR - Failed to open to read $file $!\n";
		open (my $outfh, ">$outfile") or confess  "\t    ERROR - Failed to open to write $outfile $!\n";
		my $c = 0;
		my $score;
		my $prevline;
		my @currblock = ();
		my %checksp = ();
		my $processed = 0;
		my $ID;
		while (<$infh>) {
			chomp(my $currline = $_);
			print $outfh "$currline\n" if (substr($currline,0,1) eq "#");
			if ($currline =~ /^a/) { #beginning of a block => check number of species
				$c = 0; #counter => initiate number of species per block to 0 since it is beginning of the block
				%checksp = (); #reinitialize this as well
				@currblock = (); #reiniate list since it is beginning of block
				push(@currblock,$currline);
			}
			if ($currline =~ /^s/) { #s lines = alignement lines ie sequence infos
				push(@currblock,$currline);
				my @name = split(/\./,$currline);
				$ID = $name[0];
				$c++ unless ($checksp{$ID});
				if ($c == $nbsp) { #ie reach number of species asked for
					print $outfh "\n";
					foreach my $seq (@currblock) {
						print $outfh "$seq\n"; #print the s lines of the the whole block
					}
					$processed++;
					print STDERR "          ..$processed blocks with $nbsp species\n" if (($v) && ($processed =~ /^[1-9]00+$/)); #just to have an idea of progression
				}
				$checksp{$ID}=1;
			}
		}
		$file = $outfile;
		close ($infh);
		close ($outfh);
	} else {
		print STDERR "      - $file.$nbsp.maf exists, = MAF file containing only blocks with $nbsp species\n" if ($v);
		$file = $outfile;
	}


	#####################################################
	# put sequences in table + get lengths
	#####################################################
	print STDERR "      - Getting length of alignment for each species...\n" if ($v);
	my $alignio = Bio::AlignIO->new(-file => $file, -format => 'maf');
	my $nb = 0;
	while(my $aln = $alignio->next_aln()){
		print STDERR "          ..$nb blocks done\n" if (($v) && ($nb =~ /^[1-9]00+$/));
		my %seqs = ();
		my $length;
		foreach my $seq ($aln->each_seq() ) {
			my @name = split(/\./,$seq->display_id);
			my $ID = $name[0];	
			push (@spIDs,$ID) if ($nb == 0); #do that only for first block to get species
			my @sequence = split(//,$seq->seq);
			$seqs{$ID} = \@sequence;
			$length = $seq->length;
			$totlen{$ID}+=$length;
		}
	}
	
	#AMOUNT FILE
	my $amount = "$file"."_amount.txt";
	open(my $amoutfh, ">$amount") or confess "\t    ERROR - can not create text file $amount $!\n";
	foreach my $spec (sort keys %totlen) {
		print $amoutfh "$spec\t$totlen{$spec}\n";
	}
	close ($amoutfh);
}

#####################################################
print STDERR " --- Script is done\n\n" if ($v);
exit;

