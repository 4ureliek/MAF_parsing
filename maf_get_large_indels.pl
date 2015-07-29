#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below
# email   :  4urelie.k@gmail.com
# PURPOSE :  extract big deletions/insertions from MAF file from UCSC (multiz)
##########################################################
# UPDATES
#	- v1.0 = 2 Jun 2014
#	- v1.1 = 3 Jun 2014
#			 Store info of the deletions in human genome so I can go back and check later
#   - v1.2 = 4 Jun 2014
#			 Added context, suggestion of Carson: threads->create({ 'context' => 'scalar' }, ...
#	- v1.3 = 10 Jun 2014
#			 Correction deference the $finished flag in while test ($finished was passed in as a reference).
#			 Correction of opening AND closing FH in subs to avoid Semaphore to stall
#			 Output only a certain size of deletions, default 10kbs after removing lc hg19
#			 Filter out Ns in hg19   
#	- v2.0 = 03 Feb 2015
#            first species in aln can be set
#	- v2.1 = 02 Jun 2015
#			 Little debugging 
#              - option -len usage was not matching the script
#              - trying to have threads returning (it was printing the outputs and was done but threads were still started, loop was wrong...)
#	- v2.2 = 22 Jul 2015
#            Forgot to collect -sp value
#	- v3.0 = 23-28 Jul 2015
#            Merge with the maf_get_empty_data.pl script to extract with hg19 coordinates, both the C lines and empty blocks.
#            Because thechnically, they are both potential large deletions.
#            Also, output both length with and without lower cases in hg19
#            Filter on no lc length

#  TO DO: add the R plots subroutine  
#######################################################
#always load forks before anything else
use forks;
use forks::shared;
#load the rest
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Cwd qw(cwd);

#keep STDOUT and STDERR from buffering
select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately

my $version = "3.0";
my $usage = "\nUsage [v$version]: 
    perl <scriptname.pl> -dir <dir_with_alignments> [-sp <SP>] [-lenc <X>] [-lene <X>] [-cpu <X>] [-v] [-h|help]
    	
    Outputs:
      I.  List of sizes of all consecutive empty data for a given species (no data on the browser)
	      Only when the empty data is interrupting a scaffold (note that it could correspond to a misassembly)
	      Empty data = no \"s\" line. However, there could be insertion (non aligning bases): in that case,
	      the large indel will be printed only if 80% of its lower case length is > insertion length
	      
      II. List of sizes of all consecutive empty blocks for a given species (C lines in the maf = continuous:
	      \"C -- the sequence before and after is contiguous implying that this region was either deleted in the 
	      source or inserted in the reference sequence. The browser draws a single line or a '-' in base mode in these blocks.\"

      In both cases, outputs will be:
	      [0]     [1]       [2]     [3]                 [4]                  [5]            [6]               [7]
	      chr(SP) start(SP) end(SP) chr/scaffold(Query) DeletionCoord(Query) length_with_lc length_without_lc length_insertion

	      lc = lower cases. In column 6, lower cases in reference are removed 
	           Length in col 6 would be the minimum deletion len (underestimated when ancient TEs, 
	           but this will remove all TE insertions in the reference that would not be a deletion).
	
    MANDATORY ARGUMENT:	
    -dir  (STRING) => directory with all maf files. There can be sub directories in it too.
	 
    OPTIONAL ARGUMENTS
    -sp   (STRING) => species ID that is the first in the alignment = reference species
                      default = hg19 (multiz46 and 100ways)
    -lenc (STRING) => minimum length for a gap to be reported (with lowercases not included)
                      for continuous data. Default = 1000.
    -lene (STRING) => minimum length for a gap to be reported (with lowercases not included)
                      for empty data. Default = 10000.
    -out (STRING)  => output directory. Default = <dir>_large_indels
                      /!\\ This directory will be deleted at each run if it was previously generated
    -cpu  (INT)    => max number of parallel jobs running (CPUs used). 
                      Useless if higher than number of files.
                      Default = 1 (e.g. no threading)
                      Note that files are put up in memory and a lot is stored, 
                      so don't start too many threads if don't have a lot of memory.
    -v    (BOOL)   => verbose mode, make the script talks to you
    -v    (BOOL)   => if only option, returns version
    -h    (BOOL)   => print this help
    -help (BOOL)   => print this help\n\n";
	
##########################################################################################################
# MAIN
##########################################################################################################
#Deal with options
my $sp = "hg19";
my $cpu = 1;
my $min_lenc = 1000;
my $min_lene = 10000;
my ($dir,$data,$help,$v);
GetOptions ('dir=s' => \$dir, 'sp=s' => \$sp, 'lenc=s' => \$min_lenc, 'lene=s' => \$min_lene, 'out=s' => \$data, 'cpu=s' => \$cpu, 'h' => \$help, 'help' => \$help, 'v' => \$v);

#check step to see if mandatory argument is provided + if help
die "\n Version = $version\n\n" if ((! $dir) && ($v));	
die $usage if ((! $dir) || ($help));
die "\n ERROR: $dir does not exit?\n\n" if (! -e $dir);

#avoid / at the end of paths + check blast location provided
$dir = $1 if ($dir =~ /^(.*)\/$/);

#print "log"
my $path = cwd($dir);
$data = $dir."_large_indels" unless ($data);
if ($v) {
	print STDERR "\n --- Script maf_get_empty.pl started (v$version), with following parameters:\n";
	print STDERR "      - Directory containing input maf files = $path/$dir\n";
	print STDERR "      - First species in aln blocks set to: $sp\n";
	print STDERR "      - Only large indel of length without lowercases > $min_lenc will be in outputs for continuous lines\n";
	print STDERR "      - Only large indel of length without lowercases > $min_lene will be in outputs for empty data\n";
	print STDERR "        (also, 80% of the length without lower case will need to be > insertion in the query if any)\n";
	print STDERR "      - Output files will be located in $path/$data\n";
	print STDERR "      - Max number of CPUs used = $cpu\n";
	print STDERR "\n";
	print STDERR " --- Now running to obtain potential large deletions\n";
}


#Initialize and open Thread stuff
# For deletion part
my @maf_list :shared; #shared list of files on which actions will be done

#get maf files list
print STDERR "     - Getting list of .maf files\n" if ($v);
@maf_list = `find $dir -type f -name \"*.maf\"` or confess "\nERROR (main): can't list files in $dir $!\n"; 
die "\nERROR (main): no .maf files in $dir?\n$usage\n" unless (@maf_list);

#directory with outputs; clean if exists
print STDERR "     - Creating folder for outputs:\n" if ($v);
if (-e $data) {
	print STDERR "       $data previously created, deleting...\n" if ($v);
	system "rm -Rf $data"; 
}	
print STDERR "       mkdir $data\n" if ($v);
system "mkdir $data";
$data = "$data/raw";
print STDERR "       mkdir $data\n" if ($v);
system "mkdir $data"; 

#start threads and feed subroutine
print STDERR "     - Starting $cpu threads\n" if ($v);
print STDERR "     - Looping through maf files...\n\n" if ($v);
for(my $i = 1; $i < $cpu; $i++){
    threads->create({ 'context' => 'scalar' }, \&thread_deletions, \@maf_list, \$data, \$sp, \$min_lene, \$min_lenc, \$v);   
}
thread_deletions(\@maf_list, \$data, \$sp, \$min_lene, \$min_lenc, \$v);

#clean open threads
print STDERR "\n --- Cleaning all threads\n" if ($v);
my $total = 0;
foreach my $thr (threads->list){
    my ($count) = $thr->join();
    $total+=$count;
}
$total=$total+1;
#Just check that files were processed and same number of position files generated
print STDERR "     => $total maf files processed\n" if ($v);

#concat output files
print STDERR "\n --- Concatenating output files\n" if ($v);
concat_species($data,"cont_del");
concat_species($data,"empty_del");

print STDERR "\n --- Script DONE\n\n" if ($v);
exit;


##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# sub thread_run = what will be done on each file
# thread_deletions(\@maf_list, \$data, \$sp, \$min_lene, \$min_lenc, \$v);
#----------------------------------------------------------------------------
sub thread_deletions {
	my ($maf_list,$data,$spe,$min_lene,$min_lenc,$v) = @_;
    my $c = 0;
    my $maf_list_nb;
    FILE: while(my $maf = shift @{$maf_list}){
		next FILE unless ($maf);
		chomp ($maf);
		$maf_list_nb = @{$maf_list};
		print STDERR "       STARTING: $maf (thr ".threads->tid().") [$maf_list_nb files to process]...\n" if ($$v);
		
		#Screen the MultiZ aln for continuous blocks without the species in the list, keep only when same scaffold before and after though.       	
		my %ref = ();
		my %spinfo = ();
		my %empty_data = ();
		my %insert = ();
		my %refinfo = ();
		my %del_data = ();		
		my $blk=0; #initialize block counts
		my %delid = (); #initialize deletion counts
		my $id = "none";
		my $refchr;	
		open(my $fh, "<", $maf) or warn "\nERROR (sub thread_deletions): could not open $maf $!\n";
		LINE: while (<$fh>) {
			chomp (my $l = $_);
			next LINE unless ($l =~ /\w/);
			my @ll = split(/\s+/,$l);
			if ($ll[0] eq "a") {
				$blk++; #count blocks
				next LINE;
			}
			#next lines are going through the species in that block
			my ($src,$start,$size,$srcsize) = ($ll[1],$ll[2],$ll[3],$ll[5]);
			my ($sp,$chr) = ($src,$src); #retrieve assembly name; this way just to be safe for un-initialized stuff
			($sp,$chr) = ($1,$2) if $src =~ /^(.+?)\.(.+?)$/;			
			if ($ll[0] eq "s") { #s line = alignement lines ie sequence infos.
				#First one is going to be reference species => get length and memorize stuff
				if ($sp eq "$$spe") {
					my $seq = $ll[6];
					($ref{'chr'}{$blk},$ref{'st'}{$blk},$ref{'len'}{$blk}) = ($chr,$start,$size);
					$seq =~ s/-|N|n//g; #remove gaps, Ns
					$ref{'lenlc'}{$blk} = length($seq);
					$seq =~ s/a|t|g|c//g; #remove lower case => no repeats
					$ref{'lennolc'}{$blk} = length($seq);
				} #Other species: 
				else {								
					my @infos = ($chr,$blk,$start); #infos specific to this line that will be stored so that when species is not seen for 1 or more block it can be detected
					#Check if species previously seem, if empty at least for one block and if same chr => continuous
					unless (($spinfo{$sp}->[0]) && ($spinfo{$sp}->[1] < $blk-1) && ($spinfo{$sp}->[0] eq $chr)) {
						$spinfo{$sp} = \@infos;
						next LINE;
					}
					#Carry one only if there was contiguous empty data
					my $firstempty = $spinfo{$sp}->[1]+1; #=> the first block where it was not seen
					my ($len,$lenlc,$lennolc,$ins) = (0,0,0,0);
					for (my $i = $firstempty; $i < $blk; $i++) { #increment length of reference
						$lenlc += $ref{'lenlc'}{$i};
						$lennolc += $ref{'lennolc'}{$i};
						$len += $ref{'len'}{$i};
					}
					#Now store info about this deletion, unless < minlen stuff or insertion > 80% of deletion
					unless (($lennolc < $$min_lene) && (($insert{$sp}{$firstempty-1}) && ($insert{$sp}{$firstempty-1} > ($lennolc*80/100)))) {
						$delid{$maf}{$src}{'empty_del'}++;
						$id = $src."#".$start."#indelempty_".$delid{$maf}{$src}{'empty_del'}; #Will start at 1
						$empty_data{$sp}{$id}{'chr'}=$ref{'chr'}{$firstempty};
						$empty_data{$sp}{$id}{'st'}=$ref{'st'}{$firstempty};
						$empty_data{$sp}{$id}{'en'}=$ref{'st'}{$firstempty}+$len;#length of the reference
						$empty_data{$sp}{$id}{'lenlc'}=$lenlc;
						$empty_data{$sp}{$id}{'lennolc'}=$lennolc;
						$empty_data{$sp}{$id}{'Qst'}=$spinfo{$sp}->[2];									
						$empty_data{$sp}{$id}{'ins'}=$insert{$sp}{$firstempty-1} if ($insert{$sp}{$firstempty-1}); #If inserted sequence, block just before the first empty block will have the info
					}
					undef $insert{$sp}; #reinitialize anyway, no need to keep in memory
					$spinfo{$sp} = \@infos; #reinitialize
				}
			} 
			#Make sure there are no "I" lines when empty data => load them with block ID to be able to recover it
			elsif (($ll[0] eq "i") && ($ll[-2] eq "I")) {
				$insert{$sp}{$blk}=$ll[-1]; #remember that there was inserted sequence	
			}
			#now if line is not species info and reference => check if e, look for C. Same as deletions, need to keep the start in memory.
			elsif (($ll[0] eq "e") && ($ll[-1] eq "C")) {
				#Now simply increment reflen as deletion
				#Not first line => increment stuff
				if ((exists $del_data{$sp}) && (exists $del_data{$sp}{$id})) {
					$del_data{$sp}{$id}{'lenlc'}+=$ref{'lenlc'}{$blk};
					$del_data{$sp}{$id}{'lennolc'}+=$ref{'lennolc'}{$blk};
					$del_data{$sp}{$id}{'en'}+=$ref{'len'}{$blk}; #end in hg19 (will be incremented until the end of the missing data)
				} else { #If first line, then also memorize the start in ref sp and set the rest
					$delid{$maf}{$src}{'cont_del'}++;
					$id = $src."#".$start."#indelcont_".$delid{$maf}{$src}{'cont_del'}; #Will start at 1, increment at each new del for each src
					$del_data{$sp}{$id}{'lenlc'}=$ref{'lenlc'}{$blk};
					$del_data{$sp}{$id}{'lennolc'}=$ref{'lennolc'}{$blk};
					$del_data{$sp}{$id}{'chr'}=$ref{'chr'}{$blk};
					$del_data{$sp}{$id}{'st'}=$ref{'st'}{$blk}; #start in hg19 of same block
					$del_data{$sp}{$id}{'en'}=$ref{'st'}{$blk}+$ref{'len'}{$blk};				
				}
			}	
		}			
		#now print all deletion data for this file for all species
		print_del($data,$maf,\%empty_data,$min_lene,"empty_del"); #NB: filtering done before storage for these, to save memory
		print_del($data,$maf,\%del_data,$min_lenc,"cont_del");		
		print STDERR "        ...DONE: $maf (thr ".threads->tid().")\n" if ($$v);
		$c++;
	}
	print STDERR "\n       ==> thread ".threads->tid()." returning\n"  if ($$v);
    print STDERR "           => size of list of files still to process = $maf_list_nb [should be 0]\n" if ($$v);
    print STDERR "           => Number of files preocessed by this thread = $c\n" if ($$v);
	return ($c);
}
	
#----------------------------------------------------------------------------
# print_del (called by sub thread_deletions)
# print_del($data,$maf,\%del_data,$min_len,"cont_del");
# print_del($data,$maf,\%empty_data,$min_len,"empty_del");
#----------------------------------------------------------------------------
sub print_del {
	my ($data,$maf,$del,$minlen,$type) = @_;
	$maf =~ s/.*\/(.*)$/$1/;
	foreach my $sp (sort keys %{$del}) {
		my $out = "$$data/$maf.$sp.$type.tab";		
		open (my $fh,">",$out) or warn "\n      ERROR (sub thread_deletion; sub print_del): can't open to write $out $!\n"; 	
		ID: foreach my $id (keys %{$del->{$sp}}) {
			my ($src,$start,$delid) = split(/#/,$id);
			my ($srcsp,$srcchr) = ($1,$2) if $src =~ /^(.+?)\.(.+?)$/;			
 			$start = $del->{$sp}{$id}{'Qst'} if ($type eq "empty_del");
 			my $lenlc = $del->{$sp}{$id}{'lenlc'};
			my $lennolc = $del->{$sp}{$id}{'lennolc'};
			$del->{$sp}{$id}{'ins'} = "na" unless ($del->{$sp}{$id}{'ins'});
 			print $fh "$del->{$sp}{$id}{'chr'}\t$del->{$sp}{$id}{'st'}\t$del->{$sp}{$id}{'en'}\t$srcsp\t$srcchr\t$start\t$delid\t$lenlc\t$lennolc\t$del->{$sp}{$id}{'ins'}\n" if ($lennolc > $$minlen);	
		}
		close $fh;
		
		#delete empty files
		my $lines = 0;
		open(my $fhout,"<", $out) or warn "\n      ERROR (sub thread_deletion; sub print_del): can't open to read $out $!";
		CHK: while(<$fhout>) {
			$lines++;
			last CHK if ($lines > 0);
		}
		close $fhout;			
		unlink $out if ($lines == 0);
	}
}

#----------------------------------------------------------------------------
# concat_species 
# concat_species($data,"cont_del");
# concat_species($data,"empty_del");
#----------------------------------------------------------------------------
sub concat_species {
	my ($data,$type) = @_;	
	my $path = cwd($data);
	my $fullpath = "$path/$1" if ($data =~ /^(.*)\/raw$/);
	
	my @check = glob("$data/*.$type.tab");
	if (@check) {	
		my @dels = `ls $data/*.$type.tab`;
		foreach my $file (@dels) {
			chomp ($file);
			$file =~ s/.*\/(.*)$/$1/;
			my $sp = $1 if $file =~ /^.*\.maf\.(.+)\.$type\.tab$/;
			system "cat $data/*$sp.$type.tab > $fullpath/$sp.$type.tab";
		}
	} else {
		 warn "     WARN: no data for $type for any species\n"
	}	
}
