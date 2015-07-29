#!/usr/bin/perl -w
##########################################################
# SUBROUTINES
# FOR SUBSET OF SCRIPTS = MAF_microdel--2--analyze-gaps-XXX.pl
##########################################################
# Author  :  Aurelie Kapusta - https://github.com/4ureliek
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  parsing a .maf (multispecies alignment file) to list gaps to get deletions 1 to 30, based on gaps in alignments.
#            This script is step 2/2 and it will analyze gaps previously listed
#			    -> separate micro deletions in two categories: all and 1-30nt.
#            It is very stringent - a gap is considered specific to a species or group of species if not shared by any others.
#
# BEDTOOLS NOTES:
#	consider same gap if overlap is >90% (ie flexibility of 1nt per 10nt of gaps) => option -f 0.80
#	this min overlap is for BOTH sense, ie if BIG gap in one specie, and a small one in the other, not same gap => option -r [only for Intersect]
#	options -wa -wb allow to keep the whole thing (ie => join; keep original entry of a or b or both)
#	to do a subtraction that keeps only shared stuff, ie exclude => -v means keep only stuff that DON'T overlap
##########################################################
package MAFmicrodel;
use strict;
use warnings;
use Carp;

my $version = "3.1";

my $changelog = "
# UPDATES for MAF_microdel--2--analyze-gaps-XXX.pl
#	- v1.0 = 22 Dec 2011
#				New script, to correct the problem coming from a situation like that:
#					
#					specie 1 ATGCATGCATGCATGCATGCATGCAT----------------------ATGCATGCATGCATGCATGCATGC
#					specie 2 ATGCATGCATGC----------------------ATGCATGCATGCATGCATGCATGCATGCATGCATGC
#					                    | lineage spe | shared | lineage spe |
#
# 					But these should be counted as 2 indep gaps probably => I probably underestimated lineage specific gaps	
#	
#	- v2.0 = Jul 2013
#				No fasta step
#				Gap analysis made with bedtools [so much FASTER]
#				Added part of script to filter blocks that don't have the X species
#	- v2.1 = 12 Jul 2013
#				Now this script = just analyze gaps (would be specific to an analysis, other script is not)
#	- v3.0 = 07 May 2015 + 18-24 Jun 2015
#				Write the ones for birds + added cow in mammals
#               Subroutines + MAFmicrodel.pm
#				Bug fix
#               Use of -f 0.8 instead of 0.9 => was making stuff specific when visually they were not
#               Changelog, usage, options (+ added -concat, -out and -bed options)
#               Use of intersectBed and not subtractBed, to be able to use -r with the -f 0.80
#	- v3.1 = 08 Jul 2015
#				Branches of chiroptera/rodents/cow are not really resolved
\n";

my $usage = "\nUsage [v$version]: 
    perl MAF_microdel--2--analyze-gaps-XXX.pl -in <_*.gaps.1-30.bed> -sp <X> [-aln <file.maf>] [-out <path>] [-bed <path>] [-v] [-chlog] [-h]
	
    This script is step 2/2 and it will analyze gaps previously listed by MAF_microdel--1--get-gaps.pl
	
    REQUIREMENTS
    BEDtools is required for this pipeline to run (Quinlan AR and Hall IM, 2010. Bioinformatics)
    Available at: https://github.com/arq5x/bedtools2
       Notes:
        - Two gaps between two species are considered to be specific if they overlap less than 80%
             => option -f 0.80 -v of intersectBed
        - This min overlap is for BOTH sense
             (if BIG gap in one specie, and a small one in the other, not same gap)
             => option -r of intersectBed
        - To be able to use -r, intersectBed -v instead of subtractBed is used
	
    MANDATORY ARGUMENTS:	
     -in     => (STRING) Directory containing the files = _*.gaps.1-30.bed
                         Should be one file per species (outputs of MAF_microdel--1--get-gaps.pl)
                         Should contain *_amount.txt file(s) (output of MAF_microdel--1--get-gaps.pl). 
                         If not, you need to also specify -aln
     -sp     => (INT)    number of species retained
                         
    OPTIONAL ARGUMENTS
     -concat => (BOOL)   Choose this if there are several files for each species, e.g. chr1, chr2, etc.
                         They will be concatenated prior to running everything
     -aln    => (STRING) Original alignment file that was the input of the first script, 
                         in case you don't have the *_amount.txt (outputs of MAF_microdel--1--get-gaps.pl)
                         Typically: ALL.birds.maf.8.aln
     -out    => (STRING) path of a directory where to write output files (default = directory with input files)
     -bed    => (STRING) if BEDTools are not in your path, use this option to provide path of BEDtools bin directory
     -v      => (BOOL)   verbose mode, make the script talks to you
     -v      => (BOOL)   print version if only option
     -chlog  => (BOOL)   print change log (updates)
     -h|help => (BOOL)   print this usage\n\n"; 


#----------------------------------------------------------------------------
# Extract info from the amount file(s), or from aln files
# my $path = MAFmicrodel::check_options($bedtools,$in,$out,$nbsp,$aln,$concat,$help,$chlog,$v);
#----------------------------------------------------------------------------
sub check_options {
	my ($bedtools,$in,$out,$nbsp,$aln,$concat,$help,$chlog,$v) = @_;
	my $path;
	($out eq "na")?($path=$in):($path = $out);

	#check step to see if mandatory argument is provided + if help/changelog
	die "\n MAF_microdel--2--analyze-gaps-XXX.pl version $version\n\n" if (($in eq "na") && ($nbsp eq "na") && ($help eq "na") && ($chlog eq "na") && ($v));
	die $changelog if ($chlog ne "na");
	die $usage if ($help ne "na");
	die "\nERROR: -in and -sp need to be provided (see usage below)\n$usage" if (($in eq "na") || ($nbsp eq "na"));
	die "\n Directory with input files ($in) does not exist\n\n" if (! -e $in);
	die "\n Alignment file ($aln) does not exist\n\n" if (($aln ne "na") && (! -e $aln));
	die "\n Value for -sp ($nbsp) needs to be an integer\n\n" if ($nbsp !~ /[0-9]+/);
	$bedtools =~ s/\/$//; #remove / at the end of directory
	die "\n $bedtools does not contain BEDtools binaries?\n\n" if (($bedtools ne "na") && (! -e "$bedtools/intersectBed"));

	$in =~ s/\/$//; #remove / at the end of directory
	$out =~ s/\/$//; #remove / at the end of directory
	`mkdir $out` unless (-e $out);
	
	if ($v) {
		print STDERR "\n --- Script MAF_microdel--2--analyze-gaps-XXX.pl version $version started, with:\n";
		print STDERR "      - Input files in $in\n";
		print STDERR "      - Several input files per species -> will be concatenated\n" if ($concat ne "na");
		print STDERR "      - Output directory: $path\n";
		print STDERR "      - Alignment file(s): $aln\n" if ($aln ne "na");
		print STDERR "      - Number of species per block: $nbsp\n\n";
		print STDERR "      - BEDtools path: $bedtools\n" if ($bedtools ne "na");
	}
	return ($path,$in,$out,$bedtools);
}

#----------------------------------------------------------------------------
# Extract info from the amount file(s), or from aln files
# my ($totlen,$spIDs) = get_amounts($in,$aln,$nbsp,$path,$v); #totlen = hash and spIDs = list
#----------------------------------------------------------------------------
sub get_amounts {
	my ($in,$aln,$nbsp,$path,$v) = @_;
	print STDERR " --- Getting total length of alignement + list of species...\n" if ($v);
	my $totlen = ();
	my $spIDs = ();
	my @amounts = `ls $in/*_amount.txt`;
	if (exists $amounts[0]) {
		print STDERR "      - amount file(s) exit(s), extracting amounts...\n" if ($v);
		my $c = 0;
		foreach my $amtfile (@amounts) {
			chomp ($amtfile);
			print STDERR "      - from $amtfile\n" if ($v);
			open(my $amountfh, "<$amtfile") or confess "\t    ERROR - can not open to read $amtfile $!\n";
			while(<$amountfh>) {
				chomp (my $line = $_);
				my($ID,$len) = split(/\t/,$line);
				($totlen->{$ID})?($totlen->{$ID}+=$len):($totlen->{$ID}=$len);
				push (@{$spIDs},$ID) if ($c < $nbsp); #get species, only one time though
				$c++;
			}
			close $amountfh;
		}
	} elsif ($aln ne "na") {	
		my $alncheck = 0;
		print STDERR "      - amount file(s) not found, getting amount from $aln\n" if ($v);
		my $alignio = Bio::AlignIO->new(-file => $aln, -format => 'maf');
		while (my $aln = $alignio->next_aln()){
			foreach my $seq ($aln->each_seq() ) {
				my @name = split(/\./,$seq->display_id);
				my $len = $seq->length;
				my $ID = $name[0];
				($totlen->{$ID})?($totlen->{$ID}+=$len):($totlen->{$ID}=$len);
				push (@{$spIDs},$ID) if ($alncheck == 0); #get species, only one time though	
			}
			$alncheck = 1;
		}	
		my $amtfile = "$path/_amount.txt";
		open(my $amountfh, ">$amtfile") or confess "\t    ERROR - can not open to write $amtfile $!\n";
		foreach my $spec (sort keys %{$totlen}) {
			print $amountfh "$spec\t$totlen->{$spec}\n";
		}
		close $amountfh;
	} else {
		print STDERR "      - option -aln not set but amount.txt files not found\n" if ($v);
		confess "\t    ERROR - option -aln not set but amount.txt files not found $!\n";
	}
	print STDERR "      - done\n" if ($v);
	print STDERR "      - Species studied are @{$spIDs}\n\n" if ($v);
	return ($totlen,$spIDs);
}

#----------------------------------------------------------------------------
# Concat files if set, or test if it should be
# concat_gaps($spIDs,$in,$concat,$v);
#----------------------------------------------------------------------------
sub concat_gaps {
	my ($spIDs,$in,$concat,$v) = @_;	
	print STDERR " --- Concatenating .bed files for each species into 1...\n" if (($concat ne "na") && ($v));
	SPECIES: foreach my $sp (@{$spIDs}) {
		my @list = `ls $in/*$sp.gaps.1-30.bed`;
		die "\t    ERROR - Several input files *.gaps.1-30.bed found for $sp, did you forget to set -concat?\n" if (($concat eq "na") && ($list[1]));
		 if (-e "$in/_$sp.gaps.1-30.bed") {
		 	print STDERR "     - skipping $sp ($in/_$sp.gaps.1-30.bed exists)\n" if ($v);
		 	next SPECIES;
		 }
		`cat $in/*$sp.gaps.1-30.bed > $in/_$sp.gaps.1-30.bed`;
	}
	print STDERR "      - done\n\n" if ($v);
	return;
}

#----------------------------------------------------------------------------
# use 2 lists, species that have to share the gap, and species that have to not shared the gap [this is very stringent and probably not biological]
# $files = split_gaps(\@spIDs,"no","gaps.1-30","all",$in,$path,$files,$bedtools,$v); 
#----------------------------------------------------------------------------
sub split_gaps {
	my ($to_loop,$to_sub,$ID,$type,$in,$path,$files,$bedtools,$v) = @_;
	#ID=>#define previous outputs as input here
	my ($intersectBed,$subtractBed) = ("intersectBed","subtractBed");
	($intersectBed,$subtractBed) = ("$bedtools/intersectBed","$bedtools/subtractBed") if ($bedtools ne "na");
	my $ok = 1;
	my %noshared = ();
	foreach my $specie (@{$to_loop}) {
		print STDERR "        -> species = $specie; ID = $ID, type = $type\n" if ($v);
		#DEFINE FILE NAMES
		my $shared     = "$path/_$specie.gaps.1-30.shared.$type.bed"; #final output will shared gaps
		my $not_shared = "$path/_$specie.gaps.1-30.not-shared.$type.bed"; #new output with these gaps above subtracted
		#LOOP
		my $sp = 0;
		my $indels_in;
		($to_sub eq "no")?($indels_in = "$in/_$specie.$ID.bed"):($indels_in = "$path/_$specie.$ID.bed"); #files previously generated; all gaps that are not shared (besides first time this is used)
		my $indel_in_temp = $indels_in; #to avoid rw of $indel_in during the looping stuff
		print STDERR "             -> loop in list of species of the same group to get shared gaps\n" if ($v);
		foreach my $otherspec (@{$to_loop}) { #loop on all the species, to intersect one by one the current $specie with the $otherspecie
			my $indels_out_temp;
			if ($otherspec ne $specie) { #useless to intersect if same species
				$indels_out_temp = "$shared.temp.$sp";
				my $indels_to_check;
				($to_sub eq "no")?($indels_to_check = "$in/_$otherspec.$ID.bed"):($indels_to_check = "$path/_$otherspec.$ID.bed"); #files previously generated
				print STDERR "              intersectBed -a $indel_in_temp -b $indels_to_check -f 0.80 -r -wa > $indels_out_temp\n" if ($v);
				`$intersectBed -a $indel_in_temp -b $indels_to_check -f 0.80 -r -wa > $indels_out_temp`;
				# now this output is the input for next round
				$indel_in_temp = $indels_out_temp;
			}
			$sp++;
			#check previous output; if empty, it means there won't be shared gaps between the species
			#But still need the files to exist so do not exit the loop
			unless ($otherspec eq $specie) {
				my $lines = 0;
				open(my $prev,"<", $indels_out_temp) or warn "           ERROR - can not open to read file $indels_out_temp $!";
				CHECKLOOP: while(<$prev>) {
					$lines++;
					last CHECKLOOP if ($lines > 0);
				}
				close $prev;
				$ok = 0 if ($lines < 1);
				if (($ok == 0) && (! $noshared{$type})) {
					print STDERR "           WARN: No shared gaps between $type species\n" if ($v);
					$noshared{$type}=1;
				}
			}	
		}
		# If there are shared ones:
		# - keep file with shared gaps between all species
		`mv $indel_in_temp $shared`; # rename to keep the last file
		`rm -f $shared.temp*`; # remove temp intermediate files
		# - GET FILE FOR NEXT STEPS => SUBTRACT SHARED STUFF TO GET NON SHARED STUFF (will be input for next time this subroutine is used)
		print STDERR "             -> getting input for the next round\n" if ($v);
		print STDERR "              subtractBed -a $indels_in -b $shared > $not_shared\n" if ($v);
		`$subtractBed -a $indels_in -b $shared > $not_shared`; #no need flexibility here, same spec => same coords	

		# - NOW FILTER TO GET SHARED THAT ARE SPECIFIC (super stringent since here can't be convergence) unless no need
		#   here I could sub the original files, but would be longer. Only do that if there are no "not-shared" for that species (out of groups)
		unless (($to_sub eq "no") || (($noshared{$type}) && ($noshared{$type}==1))) {
			my $i = 0;
			my $indel_in_temp2 = $shared; #to avoid rw of $indel_in
			print STDERR "             -> loop in complementary list of species to remove gaps shared with these\n" if ($v);
			foreach my $spec (@{$to_sub}) {
				my $indel_to_sub = "$path/_$spec.$ID.bed"; #stuff are going to be subtracted
				$indel_to_sub = "$in/_$spec.gaps.1-30.bed" if (! -e "$path/_$spec.$ID.bed"); #sub the original gaps if needed
				my $indels_out_temp2 = "$shared.temp.$i"; #files previously generated = shared but not filtered
				print STDERR "                 subtractBed -a $indel_in_temp2 -b $indel_to_sub -f 0.80 -r > $indels_out_temp2\n" if ($v);
				`$intersectBed -a $indel_in_temp2 -b $indel_to_sub -f 0.80 -r -v -wa > $indels_out_temp2`; #flexibility here, b/c not same species. This syntax means subtracting b from a, return a.
				$indel_in_temp2 = $indels_out_temp2; #now this output is the input for next round
				$i++;
			}
			#      keep file with shared gaps between all species
			`mv $indel_in_temp2 $shared`; # rename to keep the last file => really shared ones now
			#`rm -f $shared.temp*`; # remove temp intermediate files	
		} 
		push(@{$files},$shared);
	}
	print STDERR "       ..done\n\n" if ($v);
	return($files);
}

#----------------------------------------------------------------------------
# Get specific gaps, not shared with anything else
# $files = MAFmicrodel::spe_gaps($spIDs,$in,$path,$files,$bedtools,$v);
#----------------------------------------------------------------------------
sub spe_gaps {
	my ($spIDs,$in,$path,$files,$bedtools,$v) = @_;
	print STDERR "     - Get gaps species specific..\n" if ($v);
	my $intersectBed = "intersectBed";
	$intersectBed = "$bedtools/intersectBed" if ($bedtools ne "na");
	foreach my $specie (@{$spIDs}) {
		print STDERR "        -> species = $specie\n" if ($v);
		my $indels_out = "$path/_$specie.gaps.1-30.specific.bed";
		my $indels_in = "$in/_$specie.gaps.1-30.bed"; #files previously generated; all gaps for this species
	
		#LOOP to subtract all other gaps
		my $i = 0;
		my $indel_in_temp = $indels_in; #to avoid rw of $indel_in
		foreach my $otherspec (@{$spIDs}) {
			unless ($otherspec eq $specie) {
				my $indel_to_sub = "$in/_$otherspec.gaps.1-30.bed"; #here I can use all gaps, doesn't change anything
				my $indels_out_temp = "$indels_out.temp.$i";
				print STDERR "              intersectBed -a $indel_in_temp -b $indel_to_sub -f 0.80 -r -wa -v > $indels_out_temp\n" if ($v);
				`$intersectBed -a $indel_in_temp -b $indel_to_sub -f 0.80 -r -wa -v > $indels_out_temp`; #flexibility here, b/c not same species. This syntax means subtracting b from a, return a.
				$indel_in_temp = $indels_out_temp; #now this output is the input for next round
				$i++;
			}	
		}
		#      keep file with only gaps that are specific to current species
		system "mv $indel_in_temp $indels_out"; # rename to keep the last file
		system "rm -f $indels_out.temp*"; # remove temp intermediate files
		push(@{$files},$indels_out);
	}
	print STDERR "       ..done\n\n" if ($v);
	return($files);
}

#----------------------------------------------------------------------------
# print the amounts of gaps
# print_amounts(($path,$totlen,$files,$v);
#----------------------------------------------------------------------------
sub print_amounts {
	my ($path,$totlen,$files,$v) = @_;
	print STDERR " --- Analyzing shared/specific gap files (located in $path)\n" if ($v);

	my $out = "$path/__microdeletions_RESULTS.tab";
	open(my $outfh,">$out") or confess "\t    ERROR - an not open to write $out $!";
	print $outfh "Total length (nt) of alignement per species (should be the same!):\n";

	#tot length is stored in hash = my %totlen (key = spec so values should same for all);
	my $tot = 0;
	foreach my $spec (sort keys %{$totlen}) {
		print $outfh "$spec\t$totlen->{$spec}\n";
		$tot = $totlen->{$spec};
	}
	print $outfh "File_name\tmicrodel(number)\tmicrodel(nt)\tmicrodel(%_of_aln)\n\n";
	foreach my $f (@{$files}) {
		#0				1		2	3	4	5	6			7
		#chr1_block.4	start	end	nb	.	+	gap_len		aln_len(block)
		open(my $infh, "<$f") or confess "\t    ERROR - can not open to read file $f $!";
		my $totgaplen;
		my %gapcount;
		while (<$infh>) {
			chomp (my $line = $_);
			my @line = split(/\t/,$line);
			my $gaplen = $line[2]-$line[1];
			my $ID = $line[0]."_".$line[3];
			$totgaplen += $gaplen;
			$gapcount{$ID}++;
		}
		close $infh;
		my $gaps = 0;
		foreach my $gap (sort keys %gapcount)  {
			$gaps++;
		}
		my $gap_per = $totgaplen / $tot * 100;	
		print $outfh "$f\t$gaps\t$totgaplen\t$gap_per\t\n";	
	}
	close $outfh;
	print STDERR "     -> file $out\n" if ($v);
	print STDERR "       ..done\n\n" if ($v);
	return;
}

