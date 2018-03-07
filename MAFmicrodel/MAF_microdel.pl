#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta - https://github.com/4ureliek
# email   :  4urelie.k@gmail.com
# PURPOSE :  parsing a .maf (multispecies alignment file) to list gaps to get deletions 1 to 30, based on gaps in alignments.
#            It is quite stringent - a gap is considered specific to a species or group of species if not shared by ANY others.
##########################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::AlignIO;
use Acme::Tools qw(minus);

#-------------------------------------------------------------------------------
#--- LOAD AND CHECK ------------------------------------------------------------
#-------------------------------------------------------------------------------
my $VERSION = "4.0";
my $SCRIPTNAME = "MAF_microdel.pl";

my $CHANGELOG;
set_chlog();
sub set_chlog {
	$CHANGELOG = "
\nUPDATES for $SCRIPTNAME
 - v4.0 = XX Mar 2018
          Big update = a newick tree can be loaded, to make the gap processing dynamic
          which means that there is no need to edit the code anymore
          So it is easier to have only one script for the whole process:
          Thus, this version is a merge of the 3 steps:
             - grep of the species of interest
             	Typically: grep '^\$\|maf\|^a\|hg19\|panTro4\|ponAbe2\|rheMac3\' chr_all.maf > chr_all.primates.maf
             - MAF_microdel--1--get-gaps.pl v2.8
             - MAF_microdel--2--analyze-gaps-XXX.pl & MAFmicrodel.pm v3.1
\n";
	return 1;
}

my $USAGE;
set_usage();
sub set_usage {
	$USAGE = "\nUsage [v$VERSION]: 
    perl $SCRIPTNAME -m <dir_with_maf> -o <output_dir> -t <tree.nwk>
                        [-b <path>] [-v] [-l] [-h]
	
    PLEASE CITE: 
    Kapusta, Suh & Feschotte (2017) PNAS (doi: 10.1073/pnas.1616702114)  
    And link to $SCRIPTNAME, vX.X, https://github.com/4ureliek/MAF_parsing 
     	
    REQUIREMENTS
    BEDtools is required for this pipeline to run (Quinlan AR and Hall IM, 2010. Bioinformatics)
    Available at: https://github.com/arq5x/bedtools2
    
    STEPS OF THE SCRIPT:
    1. Rewrite the maf files, with only the lines of the species present in the tree (set with -t)
       Like all files, these will be in the directory set with -o
    
    2. Extract gap coordinates from these maf files; one file per species per original maf file
       Gap files will be as follow:
          [0]        [1]          [2]        [3]             [4]  [5]       [6]    [7]
	      file_block start_of_gap end_of_gap gap_nb_in_block  .  \"strand\" length length_of_block
	   
       If a gap file is missing for even just one species they will be all deleted and then 
       re obtained; but if they exist for all species, the total length of gaps will be loaded
       from the *_amount.txt file
	
    3. Cat the gap files if several input maf files; ln -s otherwise
    
    4. Gap analysis, using BEDtools:
       Notes on BEDtools command lines:
        - Two gaps between two species are considered to be specific if they overlap less than 80%
          (ie flexibility of 1nt per 10nt of gaps)
             => option -f 0.80 -v of intersectBed
                These 2 gaps would be shared (17nt of overlap is > 80% of the gap)
                   specie 1 ATGCATGCATGCATGCA--------------------ATGCATGCATGCATGCATGCATGC
                   specie 2 ATGCATGCATGCAT--------------------GCAATGCATGCATGCATGCATGCATGC      
                But these 2 would be specific (15nt of overlap is < 80% of the gap)
                   specie 1 ATGCATGCATGCATGCA--------------------ATGCATGCATGCATGCATGCATGC
                   specie 2 ATGCATGCATGCAT------------------ATGCAATGCATGCATGCATGCATGCATGC            
        - This min overlap is reciprocal
             (if BIG gap in one specie, and a small one in the other, likely not from the same deletion event)
             => option -r of intersectBed
        - To be able to use -r, intersectBed -v instead of subtractBed is used
          (-v means keep only stuff that don't overlap)
	
    MANDATORY ARGUMENTS:	
     -m,--maf  => (STRING) Directory containting the .maf files to parse
     -o,--out  => (STRING) Directory where to write the files = _*.gaps.1-30.bed
                           Will also contain *_amount.txt file(s)
                           And the analyzed gaps
     -t,--tree => (STRING) Species tree, in newick format. This will guide the gap analysis.
                           Thus, the species IDs NEED to match between the .maf files and the tree,
                           and ONLY the species of interest can be in the tree
                         
    OPTIONAL ARGUMENTS
     -f,--f    => (INT)    Set -f to X to filter out blocks < X nt
                           Default: 50 (30 nt indels won't be in small blocks...)
     -b,--bed  => (STRING) If BEDTool binaries are not in your path, 
                           use this option to provide path of BEDtools bin directory
     -v        => (BOOL)   Verbose mode, make the script talks to you
     -v        => (BOOL)   Print version if only option
     -l,--log  => (BOOL)   Print change log (updates)
     -h,--help => (BOOL)   Print this usage\n\n"; 
	return 1;	
}

my ($MAF,$OUT,$TREE,$BEDTOOLS,$V,$HELP,$CHLOG);
my $LEN = 50;
GetOptions ('m=s'  => \$MAF, 
            'o=s'  => \$OUT, 
            't=s'  => \$TREE, 
            'b=s'  => \$BEDTOOLS,
            'f=s'  => \$LEN, 
            'l'    => \$CHLOG, 
            'h'    => \$HELP, 
            'v'    => \$V);

my ($OUTTBED,$SUBBED);
check_options();
sub check_options {
	#CHECK OPTIONS:
	die "\n $SCRIPTNAME version $VERSION\n\n" if (! $MAF && ! $OUT && ! $TREE && ! $HELP && ! $CHLOG && $V);
	die $CHANGELOG if ($CHLOG);
	die $USAGE if ($HELP);
	die "\nERROR: -m, -o and -t need to be provided (see usage below)\n$USAGE" if (! $MAF || ! $OUT || ! $TREE);
	die "\n Tree file ($TREE) does not exist\n\n" if (! -e $TREE);
	die "\n Directory with input files ($MAF) does not exist?\n\n" if (! -e $MAF);
	$MAF =~ s/\/$//;
	$OUT =~ s/\/$//;
	unless (-e $OUT) {
		mkdir $OUT or confess "     ERROR - can not mkdir $OUT $!";
	}
	if ($BEDTOOLS) {	
		$BEDTOOLS =~ s/\/$//;
		($OUTTBED,$SUBBED) = ("$BEDTOOLS/intersectBed","$BEDTOOLS/subtractBed");
		confess "\n $BEDTOOLS does not contain BEDtools binaries?\n\n" if (! -e $OUTTBED || ! -e $SUBBED);
	} else {
		($OUTTBED,$SUBBED) = ("intersectBed","subtractBed");
		my $OUTt = `which $OUTTBED`;
		chomp($OUTt);
		confess "\n -bed not set, but intersectBed not in your path?\n\n" if (! -e $OUTt);
	}

	#LOG:
	if ($V) {
		print STDERR "\n --- Script MAF_microdel--2--analyze-gaps-XXX.pl version $VERSION started, with:\n";
		print STDERR "      - Input maf files: $MAF\n";
		print STDERR "      - Output directory: $OUT\n";
		print STDERR "      - Corresponding tree = $TREE\n";
		print STDERR "      - BEDtools path: $BEDTOOLS\n" if ($BEDTOOLS);
	}
	return 1;
}


#-------------------------------------------------------------------------------
#--- LOAD & PARSE TREE ---------------------------------------------------------
#-------------------------------------------------------------------------------
#Load Newick
my $NWCK;
load_newick_tree();

#Parse Newick & load species IDs in hash
my %NWCK = ();
my $ITER;
my $LEAF = 0;
my %SPIDS = ();
my $MAX;
until ($ITER) {
	parse_newick_tree();
}
my %GROUP = ();
get_groups_from_ids();

my %PAR = ();
get_group_parents();

#Species to consider:
my @SPIDS = @{$SPIDS{$MAX}};
my $NBSP = scalar(@SPIDS);


#-------------------------------------------------------------------------------
#--- REWRITE THE MAF FILES WITH SPECIES FROM THE TREE --------------------------
#-------------------------------------------------------------------------------
my $CONCAT;
my @RWMAF = ();
grep_species_from_MAF();


#-------------------------------------------------------------------------------
#--- DEAL WITH THE GAPS --------------------------------------------------------
#-------------------------------------------------------------------------------
my %TOTLEN = ();
print STDERR " --- Processing maf files:\n" if ($V);
print STDERR "     1. Rewrite maf files to keep only blocks\n" if ($V);
print STDERR "       - with all $NBSP species\n" if ($V); 
print STDERR "       - and of at least $LEN nt\n" if ($V);
print STDERR "     2. Listing gaps for each species\n" if ($V);
foreach my $maf (@RWMAF) {
	print STDERR "       ..from $maf\n" if ($V);
	#Get MAF files containing only blocks with these X species if not done yet
	my $rwmaf = "$maf.$NBSP.$LEN.maf";
	filter_blocks($maf,$rwmaf) unless (-e $rwmaf);
	#Now extract gaps and write the files if needed, extract amounts if not
	my $exist = check_if_gaps($rwmaf);
	if ($exist eq "y") {
		print STDERR "       - gap files exist for all species => load amounts instead\n" if ($V);		
		load_amounts($rwmaf);
	} else {
		print STDERR "       - at least one species was missing the gap file => reobtain them\n" if ($V);		
		extract_and_write_gaps($rwmaf);
		#Write the amount files to keep track - even if TOTLEN is saved; unless files exist
		write_amounts($rwmaf);
	}	
}

print STDERR "     3. Several input maf files => concatenate the gap files\n" if ($CONCAT && $V);
print STDERR "     3. Only one maf file per species => rename them\n" if (! $CONCAT && $V);
#Concatenate the gap files if relevant; ln -s if not
concat_gaps();


#-------------------------------------------------------------------------------
#--- ANALYZE THE GAPS ----------------------------------------------------------
#-------------------------------------------------------------------------------
print STDERR " --- Analyzing gaps:\n" if ($V);
print STDERR " --- Obtaining shared gaps\n" if ($V);
my @FILES = (); #list that will contain the files to process: 
my ($ANC,$CURR); #for file names

#Now loop through the levels
foreach my $lvl (sort { $b <=> $a } keys %NWCK) { 
	my @sp = @{$SPIDS{$lvl}};
	my $par = $PAR{$lvl};
	if ($lvl == $MAX) {
		print STDERR "     GROUP $lvl => @sp (root)\n";	
		#top level is the "all species" => all shared gaps => not orientable		
		$ANC = "gaps.1-30";
		$CURR = $lvl;
		split_gaps($SPIDS{$lvl},"no"); 		
	} else {
		my @parsp = @{$SPIDS{$par}};
		print STDERR "     GROUP $lvl => @sp (parent = $par => @parsp)\n";	
		#Each time, set $ANC and $CURR, with
		#   $ANC = suffix file name from the first round of "split_gaps" of the last common ancestor group processed
		#   $CURR = label of this group in output file names => label the species in the list of interest
		#Use the parent for the "not shared"
		$ANC = "gaps.1-30.not-shared.".$par;
		$CURR = $lvl;
		#the use of 'minus' avoids having to retype species names all the time; second list is subtracted from the first one
		my @not_sp = minus($SPIDS{$MAX},$SPIDS{$lvl});		
		split_gaps($SPIDS{$lvl},\@not_sp);
	}
}

#Finally, species specific gaps
print STDERR " --- Obtaining species specific gaps\n" if ($V);
spe_gaps();

#Then print the data:
print_final();

print STDERR " --- Script is done\n\n" if ($V);
exit;


#-------------------------------------------------------------------------------
#--- SUBROUTINES ---------------------------------------------------------------
#-------------------------------------------------------------------------------
sub load_newick_tree {
	$NWCK = `cat $TREE`;
	chomp($NWCK);
	$NWCK =~ s/\n//;
	print STDERR " --- Tree loaded:\n" if ($V);
	print STDERR "       $NWCK\n" if ($V);
	return 1;
}

#----------------------------------------------------------------------------
sub parse_newick_tree {	
	if ($NWCK =~ /\(([\w\-_,]+?)\)/) {
		$NWCK{$LEAF} = $NWCK;
		$NWCK{$LEAF} =~ s/^.*\(([\w\-_,]+?)\).*$/$1/;		
		get_spids(); #load list of species IDs of this leaf
		$NWCK =~ s/(^.*)\(([\w\-_,]+?)\)(.*$)/$1$LEAF$3/;
		$MAX = $LEAF;
		$LEAF++;		
	} else {
		$ITER = 1;
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub get_spids {
	if ($NWCK{$LEAF} =~ /,/) {
		my @ids = ();
		@ids = split(",",$NWCK{$LEAF});
		foreach my $id (@ids) {
			if ($SPIDS{$id}) {
				my @id = @{$SPIDS{$id}};
				push (@{$SPIDS{$LEAF}},@id);
			} else {
				push (@{$SPIDS{$LEAF}},$id);
			}
		}
	} else {
		my $id = $NWCK{$LEAF};
		$id = $SPIDS{$id} if ($SPIDS{$id});
		$SPIDS{$LEAF}=$id;
	}
	return 1;
}

#----------------------------------------------------------------------------
sub get_groups_from_ids {
	while (my($key, $value) = each %SPIDS) {
		$GROUP{@{$value}}=$key;
	}
	return 1;
}

#----------------------------------------------------------------------------
sub get_group_parents {	
	$PAR{$MAX}="root"; #max lvl = all species
	foreach my $lvl (sort { $a <=> $b } keys %NWCK) { 
		my @sp = @{$SPIDS{$lvl}};
		#parent is the last group that contained the species ids in its list
		for (my $i = $lvl+1; $i < $MAX+1; $i++) {
			my @thsp = @{$SPIDS{$i}};
			my $found;
			$found = check_array(\@sp,\@thsp);
			$PAR{$lvl}=$GROUP{@thsp} if ($found);
			last if ($found);
		}
	}
	return 1;
}

#----------------------------------------------------------------------------
sub check_array {
    my ( $test, $source ) = @_;
    my %exists = map { $_ => 1 } @$source;
    foreach my $ts ( @{$test} ) {
        return if ! $exists{$ts};
    }
    return 1;
}

#----------------------------------------------------------------------------
sub grep_species_from_MAF {
	print STDERR " --- Rewriting maf files with only the species of interest\n" if ($V);
	my $grepsp = "^\$\\|maf\\|^a\\|".join("\\|",@SPIDS);	
	print STDERR "       Grep regex = $grepsp\n";
	my @mafs = `ls $MAF | grep .maf`;
	foreach my $maf (@mafs) {
		chomp($maf);
		my $rwmaf = "$OUT/$maf";
		push (@RWMAF,$rwmaf);
		next if (-e $rwmaf); #avoid rewriting if exists
		`grep "$grepsp" $MAF/$maf > $rwmaf`;
	}
	$CONCAT = 1 if (scalar (@mafs) > 1);
	return 1;
}

#----------------------------------------------------------------------------
sub filter_blocks {
	my $maf = shift;
	my $rwmaf = shift;
	open (my $fhi,"<",$maf) or confess  "    ERROR - Failed to open to read $maf $!\n";
	open (my $fho,">",$rwmaf) or confess  "    ERROR - Failed to open to write $rwmaf $!\n";
	my $c = 0;
	my $score;
	my $prevline;
	my @currblock = ();
	my %checksp = ();
	my $processed = 0;
	my $id;
	my $blocklen;
	while(defined(my $currline = <$fhi>)) {
		chomp($currline);
		print $fho "$currline\n" if (substr($currline,0,1) eq "#");
		if ($currline =~ /^a/) { #beginning of a block => check number of species
			$c = 0; #counter => initiate number of species per block to 0 since it is beginning of the block
			%checksp = (); #reinitialize this as well
			@currblock = (); #reiniate list since it is beginning of block
			push(@currblock,$currline);
			$blocklen = 0; #reinitialize
		}
		#s lines = alignement lines ie sequence infos; also need to either have $c = 0 (first species), or the block len > the len set in options
		if (($currline =~ /^s/) && (($c == 0) || ($blocklen > $LEN))) { 		
			push(@currblock,$currline);
			my @name = split(/\./,$currline);
			$id = $name[0];
			my @sequence = split('\s+',$currline);
			$blocklen = length($sequence[6]); #this will be done on the first line only
			$c++ unless ($checksp{$id});
			if ($c == $NBSP) { #ie reach number of species asked for
				print $fho "\n";
				foreach my $seq (@currblock) {
					print $fho "$seq\n"; #print the s lines of the the whole block
				}
				$processed++;
#				print STDERR "          ..$processed blocks with $NBSP species\n" if (($V) && ($processed =~ /^[1-9]00+$/)); #just to have an idea of progression
			}
			$checksp{$id}=1;
		}
	}
	close ($fhi);
	close ($fho);	
	return 1;
}

#----------------------------------------------------------------------------
sub check_if_gaps {
	my $rwmaf = shift;
	foreach my $sp (@SPIDS) {
		if (! -e $rwmaf."_$sp.gaps.1-30.bed") {	
			#at least one of the species is missing:
			# => delete all the gap files - otherwise gaps will get appended if any exists
			`rm -Rf $rwmaf*gaps.1-30.bed`;
			return "n";
		}
	}
	return "y";
}

#----------------------------------------------------------------------------
sub load_amounts {
	my $rwmaf = shift;
	my @amounts = `ls $rwmaf*_amount.txt`;
	foreach my $f (@amounts) {
		chomp($f);
		open(my $fh, "<",$f) or confess "     ERROR - can not open to read $f $!\n";
		while(defined(my $l = <$fh>)) {
			chomp $l;
			my($id,$len) = split(/\t/,$l);
			$TOTLEN{$id}+=$len;
		}
		close $fh;
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub extract_and_write_gaps {
	my $rwmaf = shift;	
	my $alignio = Bio::AlignIO->new(-file => $rwmaf, -format => 'maf');
	my $nb = 0;
	while(my $aln = $alignio->next_aln()){
#		print STDERR "          ..$nb blocks done\n" if (($v) && ($nb =~ /^[1-9]00+$/));
		my %seqs = ();
		my $length;
		foreach my $seq ($aln->each_seq() ) {
			my @name = split(/\./,$seq->display_id);
			my $id = $name[0];	
			my @sequence = split(//,$seq->seq);
			$seqs{$id} = \@sequence;
			$length = $seq->length;
			$TOTLEN{$id}+=$length;
		}
		# Get gaps and print their coordinates - unless the files exist
		my %start = ();
		my %gap_nb = ();
		for (my $n=1;$n<$length;$n++) {
			foreach my $sp (@SPIDS) {
				#nt for each species at each position is ($seqs{$sp}->[$n]
				#print in a file gap start-end couples, one per line, one file per species
				my $gaps = $rwmaf."_$sp.gaps.1-30.bed";
				open(my $fh,">>",$gaps) or confess "    ERROR - can not open file $gaps $!";
				my $chr = $rwmaf;
				$chr =~ s/^([A-Za-z1-9]+)\..*$/$1/;
				my $region = $chr."_block.".$nb;
				if ($seqs{$sp}->[$n] eq "-" && $seqs{$sp}->[$n-1] ne "-") { #this is a gap opening
					$start{$sp} = $n+1;
				}
				if ($seqs{$sp}->[$n-1] eq "-" && $seqs{$sp}->[$n] ne "-" && $start{$sp}) { 
				#this is a gap ending; unless gap is on begining of alignement - I don't want to count these since I don't know their length
					my $end = $n+1;
					my $len = $end - $start{$sp}; #gap length
					 if ($len <= 30) {
						$gap_nb{$sp}++;
						print $fh "$region\t$start{$sp}\t$end\t$gap_nb{$sp}\t.\t+\t$len\t$length\n";
					}	
				}
				close $fh;
			}
		}
		$nb++;
	}
	return 1;
}

#----------------------------------------------------------------------------
sub write_amounts {
	my $rwmaf = shift;
	my $amt = "$rwmaf"."_amount.txt";
	open(my $fh, ">$amt") or confess "    ERROR - can not create text file $amt $!\n";
	foreach my $sp (sort keys %TOTLEN) {
		print $fh "$sp\t$TOTLEN{$sp}\n";
	}
	close ($fh);
	return 1;
}

#-------------------------------------------------------------------------------
sub concat_gaps {
	foreach my $sp (@SPIDS) {
		my $out = "$OUT/_$sp.gaps.1-30.bed";
		 if ($CONCAT) {
			`cat $OUT/*$sp.gaps.1-30.bed > $out`;
		} else {
			my $file = `ls $OUT | col | grep $sp.gaps.1-30.bed | grep -v "^_"`	;
			chomp($file);
			`ln -s $file $out` unless (-e $out);
		}
	}
	return 1;
}

#----------------------------------------------------------------------------
sub split_gaps {
	my ($to_loop,$to_sub) = @_;
	my @to_loop = @{$to_loop};
	my $ok = 1;
	my %noshared = (); #check hash for shared gaps, if none
	foreach my $specie (@to_loop) {
		print STDERR "        $specie\n" if ($V);
		#DEFINE FILE NAMES
		my $shared     = "$OUT/_$specie.gaps.1-30.shared.$CURR.bed"; #final output will shared gaps
		my $not_shared = "$OUT/_$specie.gaps.1-30.not-shared.$CURR.bed"; #new output with these gaps above subtracted
		#LOOP
		my $sp = 0;
		my $i_in;
		if ($to_sub eq "no") {
			$i_in = "$OUT/_$specie.$ANC.bed"; #all gaps in this case
		} else {
			$i_in = "$OUT/_$specie.$ANC.bed"; #will be longer but now hard to automatize what group is the "parent"
		}
		my $i_in_temp = $i_in; #to avoid rw of $i_in during the looping stuff
		print STDERR "              (shared gaps in $shared; not shared in $not_shared)\n" if ($V && ! $noshared{$CURR});
		SP: foreach my $othersp (@to_loop) { 
			my $i_out_temp;
			if ($othersp ne $specie) { #useless to intersect if same species
				$i_out_temp = "$shared.temp.$sp";
				my $i_to_check = "$OUT/_$othersp.$ANC.bed";
				#Avoid doing the intersections if no gap shared between all - but need the files
				if ($noshared{$CURR}) {
					$i_in_temp = $i_out_temp;
					`touch $i_in_temp`;
					last SP;
				} else {
					print STDERR "              intersectBed -a $i_in_temp -b $i_to_check -f 0.80 -r -wa > $i_out_temp\n" if ($V);
					`$OUTTBED -a $i_in_temp -b $i_to_check -f 0.80 -r -wa > $i_out_temp`;
					# now this output is the input for next round
					$i_in_temp = $i_out_temp;
					#check this i_in_temp file; if empty, it means there won't be shared gaps between the species
					#So no need to intersectBed with the others
					my $lines = 0;
					open(my $prev,"<", $i_in_temp) or warn "           ERROR - can not open to read file $i_in_temp $!";
					while(<$prev>) {
						$lines++;
						last if ($lines > 0);
					}
					close $prev;
					$ok = 0 if ($lines < 1);
					if ($ok == 0) {
						print STDERR "           WARN: No shared gaps between $specie and all species of the group $CURR\n" if ($V);
						print STDERR "           => skipping the rest of the intersections\n" if ($V);
						$noshared{$CURR}=1;
						last SP;
					}
				}	
			}
			$sp++;	
		}
		# If there are shared ones:
		# - keep file with shared gaps between all species
		`mv $i_in_temp $shared`; # rename to keep the last file
		# - GET FILE FOR NEXT STEPS => SUBTRACT SHARED STUFF TO GET NON SHARED STUFF (will be input for next time this subroutine is used)
		print STDERR "           -> getting input for the next round\n" if ($V);
		print STDERR "              subtractBed -a $i_in -b $shared > $not_shared\n" if ($V);
		`$SUBBED -a $i_in -b $shared > $not_shared`; #no need flexibility here, same spec => same coords	

		# - NOW FILTER TO GET SHARED THAT ARE SPECIFIC (super stringent since here can't be convergence) unless no need
		#   here I could sub the original files, but would be longer. Only do that if there are no "not-shared" for that species (out of groups)
		unless ($to_sub eq "no" || $noshared{$CURR}) {
			my $i = 0;
			my $i_in_temp2 = $shared; #to avoid rw of $i_in
			print STDERR "           -> loop in complementary list of species to remove gaps shared with these\n" if ($V);
			foreach my $spec (@{$to_sub}) {
				my $i_to_sub = "$OUT/_$spec.$ANC.bed"; #stuff are going to be subtracted
				$i_to_sub = "$OUT/_$spec.gaps.1-30.bed" if (! -e "$OUT/_$spec.$ANC.bed"); #sub the original gaps if needed
				my $i_out_temp2 = "$shared.temp.$i"; #files previously generated = shared but not filtered
				print STDERR "              intersectBed -a $i_in_temp2 -b $i_to_sub -f 0.80 -r -v -wa > $i_out_temp2\n" if ($V);
				`$OUTTBED -a $i_in_temp2 -b $i_to_sub -f 0.80 -r -v -wa > $i_out_temp2`; 
				     #flexibility here, b/c not same species. This syntax means subtracting b from a, return a.
				$i_in_temp2 = $i_out_temp2; #now this output is the input for next round
				$i++;
			}
			#      keep file with shared gaps between all species
			`mv $i_in_temp2 $shared` unless ($i_in_temp2 eq $shared); # rename to keep the last file if needed (= if several files)
			#`rm -f $shared.temp*`; # remove temp intermediate files	
		} 
		push(@FILES,$shared);
	}
#	`rm -f $OUT/*.temp*`; # remove temp intermediate files
	return 1;
}

#----------------------------------------------------------------------------
sub spe_gaps {
	foreach my $sp (@SPIDS) {
		print STDERR "        $sp\n" if ($V);
		my $i_out = "$OUT/_$sp.gaps.1-30.specific.bed";
		my $i_in = "$OUT/_$sp.gaps.1-30.bed"; #files previously generated; all gaps for this species
	
		#LOOP to subtract all other gaps
		my $i = 0;
		my $i_in_temp = $i_in; #to avoid rw of $i_in
		foreach my $othersp (@SPIDS) {
			unless ($othersp eq $sp) {
				my $i_to_sub = "$OUT/_$othersp.gaps.1-30.bed"; #here I can use all gaps, doesn't change anything
				my $i_out_temp = "$i_out.temp.$i";
				print STDERR "              intersectBed -a $i_in_temp -b $i_to_sub -f 0.80 -r -wa -v > $i_out_temp\n" if ($V);
				`$OUTTBED -a $i_in_temp -b $i_to_sub -f 0.80 -r -wa -v > $i_out_temp`; 
			     #=> flexibility here, b/c not same species. This syntax means subtracting b from a, return a.
				$i_in_temp = $i_out_temp; #now this output is the input for next round
				$i++;
			}	
		}
		#      keep file with only gaps that are specific to current species
		system "mv $i_in_temp $i_out"; # rename to keep the last file
		system "rm -f $i_out.temp*"; # remove temp intermediate files
		push(@FILES,$i_out);
	}
}

#----------------------------------------------------------------------------
sub print_final {
	print STDERR " --- Analyzing shared/specific gap files (located in $OUT)\n" if ($V);
	my $out = "$OUT/__microdeletions_RESULTS.tab";
	open(my $fh,">$out") or confess "     ERROR - can not open to write $out $!";
	print $fh "#Total length of alignement per species:\n";
    print $fh "#(Should be the same - if not, the maf files may have several lines per species)\n";
    print $fh "#species\tlen(nt)\n";
	#tot length is stored in hash = my %totlen (key = spec so values should same for all);
	my $tot = 0;
	foreach my $spec (sort keys %TOTLEN) {
		print $fh "$spec\t$TOTLEN{$spec}\n";
		$tot = $TOTLEN{$spec};
	}
	print $fh "\n";
	print $fh "#Species groups are as follow:\n";
	print $fh "#group_id\tspecies_list\n";
	foreach my $lvl (sort { $b <=> $a } keys %NWCK) { 
		my @sp = @{$SPIDS{$lvl}};
		print $fh "$lvl\t@sp\n";	
	}
	print $fh "\n";
	print $fh "File_name\tmicrodel(number)\tmicrodel(nt)\tmicrodel(%_of_aln)\n\n";
	foreach my $f (@FILES) {
		#0				1		2	3	4	5	6			7
		#chr1_block.4	start	end	nb	.	+	gap_len		aln_len(block)
		open(my $ifh, "<$f") or confess "     ERROR - can not open to read file $f $!";
		my $totgaplen = 0;
		my %gapcount;
		while (<$ifh>) {
			chomp (my $line = $_);
			my @line = split(/\t/,$line);
			my $gaplen = $line[2]-$line[1];
			my $id = $line[0]."_".$line[3];
			$totgaplen += $gaplen;
			$gapcount{$id}++;
		}
		close $ifh;
		my $gaps = 0;
		foreach my $gap (sort keys %gapcount)  {
			$gaps++;
		}
		my $gap_per = 0;
		$gap_per = $totgaplen / $tot * 100 if ($totgaplen > 0);	
		print $fh "$f\t$gaps\t$totgaplen\t$gap_per\t\n";	
	}
	close $fh;
	print STDERR " --- results in file $out\n" if ($V);
	return;
}



