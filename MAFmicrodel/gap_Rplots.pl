#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  Density graphs of gaps
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Statistics::R; #required to make plots through R; comment this if not needed

my $version = "1.1";

my $changelog = "
#	- v1.0 = 13 Jul 2015
#	- v1.1 = 16 Jul 2015
#            Change at the ls, to match names of the MAFmicrodel scripts as well
#            Error, there shouldn't have had a +1 for the length. Uh-Oh.
\n";

my $usage = "\nUsage [v$version]: 
	perl <gap_density.pl> -in <file> [-dir] [-pdf] [-type <plottype>] [-v] [-chlog] [-h]
	
	This script will process gap coordinates and return command lines to plot in R (density or boxplot)
	Command lines will be for 8 graps on one image
	.pdf will be output with one file per input file if -pdf is set
	
    MANDATORY ARGUMENTS:	
     -in     => (STRING) RepeatMasker output .out OR RepeatMasker output .align
                         Typically, process first the .align to correct %div for higher mutation rate at CpG sites,
                         and with Kimura 2-Parameter divergence metric
                         To do so, run: RepeatMasker/util/calcDivergenceFromAlign.pl genome.fa.align -a genome.fa.align.k2p.noCpG
                         IMPORTANT: if 2 repeats were merged in the final .out, here they will be treated separately.
                         Therefore, amount of DNA masked by several repeats will be much higher than if the .out is parsed.
                                               
     OPTIONAL ARGUMENTS
     -dir    => (BOOL)   chose this if -in is in fact a directory containing a bunch of .bed files                       
                         Note that the directory should contain no other file with .bed than the ones to parse
     -pdf    => (BOOL)   To get pdfs (images) of the plots
     -type   => (STRING) Type of plotting:
                         -type den for density plots
                         -type box for boxplot with noth [default]
     -v      => (BOOL)   verbose mode, make the script talks to you
     -v      => (BOOL)   print version if only option
     -chlog  => (BOOL)   print change log (updates)
     -h|help => (BOOL)   print this usage\n\n";      


################################################################################
# Get arguments/options / checks
################################################################################
my $type = "box";
my ($in,$dir,$help,$pdf,$chlog,$v);
GetOptions ('in=s' => \$in, 'dir' => \$dir, 'pdf' => \$pdf, 'type=s' => \$type, 'h' => \$help, 'help' => \$help, 'chlog' => \$chlog, 'v' => \$v);

#check step to see if mandatory argument is provided + if help
die "\n Script gap_density.pl version $version\n\n" if ((! $in) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die "\n please provide input file(s)\n $usage" if (! $in);
die "\n please set -type as \"box\" or \"den\" (see usage)" if (($type ne "box") && ($type ne "den"));

#avoid / at the end of path in case it's a directory
$in = $1 if (($dir) && ($in =~ /^(.*)\/$/));

#"log"
print STDERR "\n --- Script gap_density.pl started (v$version)\n" if ($v);
print STDERR "      - Directory containing input files = $in (-dir chosen)\n" if (($dir) && ($v));
print STDERR "      - Input file = $in\n" if ((! $dir) && ($v));
die "\n $in does not exist?\n\n" if (($in) && (! -e $in));
print STDERR "      - R command lines will be written in a file per input file\n" if (($pdf) && ($v));
print STDERR "      - Rplots will be output in pdf format (one file per input file)\n" if (($pdf) && ($v));
print STDERR "         -> boxplots with notch\n" if (($type eq "box") && ($v));
print STDERR "         -> density lines\n" if (($type eq "den") && ($v));

################################################################################
# MAIN
################################################################################
print STDERR "\n --- Now extracting lengths from file(s)\n" if ($v);
print STDERR "     - Getting list of files from $in...\n" if (($dir) && ($v));
#Get list of input files if -dir, or just load the one from -in
my @files = ();
($dir)?(@files = `ls $in/*.specific.bed* | grep -v pdf | grep -v len`):(push(@files,$in));
my $path = path($in);
print STDERR "       ..done\n" if (($dir) && ($v));

($dir)?($dir = "y"):($dir = "n");
my @glist = ();
FILE: foreach my $f (@files) {
	chomp($f);
	next FILE unless ($f); #ls is weird sometimes, blank values
	$f = $in."/".filename($f) if ($dir eq "y");
	$f = $path."/".filename($f) if ($dir eq "n");
	print STDERR "     -> $f..\n" if ($v);	
	print STDERR "        ..skipped, does not exist?\n" if ((! -e $f) && ($v));
	next FILE unless (-e $f);
	my $out = $f.".len.txt";
	open my $ifh, "<", $f or confess "\nERROR (Main): can't open to read $f $!\n";	
	open my $ofh, ">", $out or confess "\nERROR (Main): can't open to write $out $!\n";	
	LINE: while(<$ifh>) {
		chomp (my $line = $_);	
		# Typically, gap files _*.gaps.specific.bed.big.tab from DelGet pipeline
		# region                                                start{$spIDs[$i]}	$end	$gap_nb{$spIDs[$i]}	.	+	$len	$length
		# ./Del/Deletions.0/_ExtractAlign/reg10-142.fa.align.fa	4923				5095	29					.	+	173		11094
		# ./Del/Deletions.0/_ExtractAlign/reg10-142.fa.align.fa	8033				8067	41					.	+	35		11094

		# Typically, gap files _*.gaps.1-30.specific.bed from MAFmicrodel pipeline
		# Block													start{$spIDs[$i]}	$end	$gap_nb{$spIDs[$i]}	.	+	$len	$length
		# maf_12mammals/chr10.12mammals.maf.12.maf_block.75		76					78		3					.	+	2		251

		#here, just need the start and end to make density plots (need to recalculate length)
		my @l = split('\t',$line);
		my ($st,$en) = ($l[1],$l[2]);
		my $len = $en - $st; #No +1 here because for the 1nt gaps I didn't put same value for start and end.
		print $ofh "$len\n";
	}
	close $ifh;
	close $ofh;
	push (@glist,$out);
}
print STDERR "     ..done\n\n" if ($v);

#print R command lines 
print STDERR " --- Printing R command lines\n" if ($v); 
get_Rcmdlines($in,\@glist,$dir,$type);

#get plots through R
print STDERR " --- Getting R plots\n" if (($v) && ($pdf)); 
get_Rplots(\@glist,$path,$type) if ($pdf);

print STDERR " --- Script done -> see files in $path\n\n" if ($v); 
exit;


##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# get a filename from a full path
# my $name = filename($filename);
#----------------------------------------------------------------------------
sub filename {
	my($file) = shift;
	$file =~ s/.*\/(.*)$/$1/;
	return $file;
}

#----------------------------------------------------------------------------
# from a filename keep only its path - independently of the current dir
# my $path = path($filename);
#----------------------------------------------------------------------------
sub path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
# print R command lines to plot coverages
# get_Rcmdlines($in,\@glist,$dir,$type);
#----------------------------------------------------------------------------
sub get_Rcmdlines {
	my ($in,$list,$dir,$type) = @_;	
	my $out = $in.".R.$type.txt";
	if ($dir eq "y") {
		open my $rfh, ">$out" or confess "ERROR (sub get_Rcmdlines): Failed to open to write $out $!\n";
		print $rfh "setwd(\"\") #fill with location of files on computer\n\n";
		close $rfh;
	}
	#now loop
	my $nb = 1;
	my $i = 1;
	foreach my $f (@{ $list }) {
		$out = $f.".R.txt" if ($dir eq "n");
		open my $rfh, ">>$out" or confess "ERROR (sub get_Rcmdlines): Failed to open to write $out $!\n";
		my $fname = filename($f);		
		my $sp = filename($fname);
		$sp = $1 if ($fname =~ /.+_(.+?)\.gaps\./);
		print $rfh "\npar(mfrow=c(4,2))\n" if (($nb == 1) && ($dir eq "y"));
		print $rfh "setwd(\"\") #fill with location of files on computer\n\n" if ($dir eq "n");
		print $rfh "dat<-read.table(\"$f\", sep=\"\\t\", header=FALSE)\n";
		if ($type eq "den") {			
			print $rfh "d<-density(dat\$V1)\n";
			print $rfh "plot(d, col=\"blue\", pch = 18, main = \"$sp\", ylab = \"Density of x=deletion size (nt)\")\n";
		} elsif ($type eq "box") {
			print $rfh "boxplot(dat\$V1,notch = T)\n";
			print $rfh "title(\"$sp\")\n";
		}
		$i++;	
		$nb = 0 if ($nb == 8);
		$nb++;
		close $rfh;
	}	
#Table[!is.na(Table$Col),])
}


#----------------------------------------------------------------------------
# get R plots
# get_Rplots(\@glist,$path,$type) if ($pdf);
#----------------------------------------------------------------------------
sub get_Rplots {
	my ($del_list,$dir,$type) = @_;
	#Start R bridge
	my $R = Statistics::R->new();
	$R->startR;
	#$R->send(q`setwd(".")`);
	#Now loop
	TAB: foreach my $tab (@$del_list){
		chomp $tab;		
		#Check if file is empty
		if (-z $tab) {
			print "     $tab empty, skipped\n" if ($v);
			next TAB;
		}	
		#plot
		$R->send(qq`dat <- read.table("$tab",header=FALSE)`);
		my $out = $tab.".".$type.".pdf"; 
		my $sp = filename($tab);
		$sp = $1 if ($tab =~ /.+_(.+?)\.gaps\./);
		$R->run(qq`pdf("$out")`);
		if ($type eq "den") {
			$R->send(qq`d<-density(dat\$V1)`);
			$R->send(qq`plot(d, col="blue", pch = 18, main = "$sp", ylab = "Density of x=deletion size (nt)")`);
		} elsif ($type eq "box") {
			$R->send(qq`boxplot(dat\$V1,notch = T)`);
			$R->send(qq`title("$sp")`);
		}
		$R->run(q`dev.off()`);
	}	
	#End R bridge
	$R->stopR() ;
}
	
# 		# PDF Format  
# 		# The 'pdf' device allows for PDF output.  
# 		#  
# 		# Height and Width: Dimensions of the image in inches  
# 		# Onefile: If true, it plots multiple figures on a page. If false,  
# 		#      each figure gets its own numbered page.  
# 		# Family: The font family used to render text. You can use most fonts  
# 		#     available to your computer.  
# 		# Paper: The type of paper determines the pdf size.  
# 		# Pointsize: Font size.  
# 		pdf(file='PDF Output.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12) 






