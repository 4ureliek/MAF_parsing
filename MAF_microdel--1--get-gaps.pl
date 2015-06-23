#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta - https://github.com/4ureliek
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  parsing a .maf (multispecies alignment file) to list gaps to get small deletions.
#            This script is step 1/2 and it will list the gaps.
#            The .maf files needs to be grepped for species of interest before running this script.
##########################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::AlignIO;

my $version = "2.4";
my $changelog = "
#	- v1.0 = 22 Dec 2011. Other script named MAF_to_fasta_Gapfreq_SpeBranches_9spec.pl.
#				Was doing the same thing, but problem was situation like that:
#					
#					specie 1 ATGCATGCATGCATGCATGCATGCAT----------------------ATGCATGCATGCATGCATGCATGC
#					specie 2 ATGCATGCATGC----------------------ATGCATGCATGCATGCATGCATGCATGCATGCATGC
#					                    | lineage spe | shared | lineage spe |
#
# 					But these should be counted as 2 indep gaps probably => I probably underestimated lineage specific gaps	
#	
#	- v2.0 = Jul 2013
#				No fasta step
#				Gap analysis made with bedtools
#				Added part of script to filter blocks that don't have the 11 species
#	- v2.1 = 12 Jul 2013
#				Now this script just lists gaps => allow to run that on // on all files, plus it is a general one
#	- v2.2 = 07 May 2015
#				Nb species in command line + couple of other little things
#	- v2.3 = 16-17 Jun 2015
#				Changelog, usage, options
#				Bug fix in the sub get_amounts (there were no outputs)
#	- v2.4 = 18 Jun 2015
#				Bug fix in printing gaps
\n";

my $usage = "\nUsage [v$version]: 
	perl MAF_microdel--1--get-gaps.pl -in <dir> -sp <X> [-v] [-chlog] [-h]
	
	This script parse .maf file(s) (multispecies alignment file) and list gaps for all species.
	It is step 1/2 to obtain small deletions (step 2 = MAF_microdel--2--analyze-gaps-XXX.pl).
	
    MANDATORY ARGUMENTS:	
     -in     => (STRING) directory containing input .maf file(s). 
                         They needs to be previously grepped for species of interest before running this script.
                         For example: grep '^\$\|maf\|^a\|hg19\|panTro4\|ponAbe2\|rheMac3\' chr_all.maf > chr_all.primates.maf
     -sp     => (INT)    number of species retained, to filter blocks that have no info for any species

    OPTIONAL ARGUMENTS     
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
die "\nERROR: -in or -sp not set (see usage below)\n$usage" if ((! $in) || (! $nbsp));
die "\n Directory with input files ($in) does not exist\n\n" if (! -e $in);
die "\n Value for -sp ($nbsp) needs to be an integer\n\n" if ($nbsp !~ /[0-9]+/);

$in =~ s/\/$//; #remove / at the end of directory
if ($v) {
	print STDERR "\n --- MAF_microdel--1--get-gaps.pl version $version started, with:\n";
	print STDERR "      - Input files in $in\n";
	print STDERR "      - Number of species per block: $nbsp\n\n";
}

print STDERR " --- Gaps in alignement are being listed:\n" if ($v);
my @in = `ls $in/*.maf`;
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
	print STDERR "      - Putting sequences into table + get gap coordinates in a file per species\n" if ($v);
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

		#####################################################
		# Get gaps and print their coordinates
		#####################################################
		my %start = ();
		my %gap_nb = ();
		for (my $n=1;$n<$length;$n++) {
			for (my $i=0;$i<=$#spIDs;$i++) {
				#nt for each species at each position is ($seqs{$spIDs[$i]}->[$n]
				#print in a file gap start-end couples, one per line, one file per species [to keep track]
				my $gaps_small = "$file"."_$spIDs[$i].gaps.1-30.bed";
				my $gaps_big = "$file"."_$spIDs[$i].gaps.31-x.bed";
				open(my $gapsmallfh,">>$gaps_small") or confess "\t    ERROR - can not open file $gaps_small $!";
				open(my $gapbigfh,">>$gaps_big") or confess "\t    ERROR - can not open file $gaps_big $!";
				my $chr = $file; #bug fix v2.4
				$chr =~ s/^([A-Za-z1-9]+)\..*$/$1/;
				my $region = $chr."_block.".$nb;
				if (($seqs{$spIDs[$i]}->[$n] eq "-") && ($seqs{$spIDs[$i]}->[$n-1] ne "-")) { #this is a gap opening
					$start{$spIDs[$i]} = $n+1;
				}
				if (($seqs{$spIDs[$i]}->[$n-1] eq "-") && ($seqs{$spIDs[$i]}->[$n] ne "-") && ($start{$spIDs[$i]})) { 
				#this is a gap ending; unless gap is on begining of alignement - I don't want to count these since I don't know their length
					my $end = $n+1;
					my $len = $end - $start{$spIDs[$i]}; #gap length
					$gap_nb{$spIDs[$i]}++;
					print $gapsmallfh "$region\t$start{$spIDs[$i]}\t$end\t$gap_nb{$spIDs[$i]}\t.\t+\t$len\t$length\n" if ($len <= 30);
					print $gapbigfh "$region\t$start{$spIDs[$i]}\t$end\t$gap_nb{$spIDs[$i]}\t.\t+\t$len\t$length\n" if ($len > 30);
				}
				close $gapsmallfh;
				close $gapbigfh;
			}
		}
		$nb++;
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
