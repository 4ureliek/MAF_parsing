#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  Plot graphs of gaps, for outputs of script maf_get_large_indels.pl
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Statistics::R; #required to make plots through R; comment this if not needed and creates errors

my $version = "2.1";

my $changelog = "
#	- v1.0 = 13 Jul 2015
#            Script named gap_Rplots.pl
#	- v2.0 = 29-30 Jul 2015
#			 Adapt to outputs of maf_get_large_indels.pl
#            ecd plots
#	- v2.1 = 31 Jul 2015
#			 summary file + avoid negative values in the .len files (shouldn't be any though)
\n";

my $usage = "\nUsage [v$version]: 
	perl <maf_get_large_indels_Rplot.pl> -in <file> [-dir] [-type <plottype>] [-xlim <X,Y>] [-stats] [-pdf] [-v] [-chlog] [-h]
	
	This script will process lengths of gaps from output files of the script maf_get_large_indels.pl
	(without lower cases in reference and minus insertion in Query if relevant)
	and return command lines to plot in R (ecd, boxplot or density)
	Command lines will be for 8 graps on one image
	Plots in pdf format will be output if -pdf is set
	
    MANDATORY ARGUMENTS:	
     -in     => (STRING) RepeatMasker output .out OR RepeatMasker output .align
                         Typically, process first the .align to correct %div for higher mutation rate at CpG sites,
                         and with Kimura 2-Parameter divergence metric
                         To do so, run: RepeatMasker/util/calcDivergenceFromAlign.pl genome.fa.align -a genome.fa.align.k2p.noCpG
                         IMPORTANT: if 2 repeats were merged in the final .out, here they will be treated separately.
                         Therefore, amount of DNA masked by several repeats will be much higher than if the .out is parsed.
                                               
     OPTIONAL ARGUMENTS
     -dir    => (BOOL)   Chose this if -in is in fact a directory containing a bunch of .bed files                       
                         Note that the directory should contain no other file with .bed than the ones to parse
     -type   => (STRING) Type of plot:
                         -type ecd for empirical cumulative distributions (log scale for x axis) [default]
                         -type box for boxplot with notches, no outliers on plots
                         -type den for density plots
     -xlim   => (STRING) To set xlim=c(xmin,xmax) in R, for ecd and density plots only. 
                         Typically: -xlim 1000,50000 will plot from size 1 kb to 50 kb
                         [Default = smallest and biggest values for each file]
     -stats  => (BOOL)   To get an additional output with numbers, max, median, averages for each species
     -pdf    => (BOOL)   To get pdfs (images) of the plots (one per input file). 
                         Requires Statistics::R (can be commented if creates issues)
                         Requires to install calibrate. As sudo:
                         	R
                         	install.packages(\"calibrate\")
                         	quit(\"default\")
     -v      => (BOOL)   verbose mode, make the script talks to you
     -v      => (BOOL)   print version if only option
     -chlog  => (BOOL)   print change log (updates)
     -h|help => (BOOL)   print this usage\n\n";      


################################################################################
# Get arguments/options / checks
################################################################################
my $type = "ecd";
my ($in,$dir,$xlim,$help,$stats,$pdf,$chlog,$v);
GetOptions ('in=s' => \$in, 'dir' => \$dir, 'pdf' => \$pdf, 'stats' => \$stats, 'type=s' => \$type, 'xlim=s' => \$xlim, 'h' => \$help, 'help' => \$help, 'chlog' => \$chlog, 'v' => \$v);

#check step to see if mandatory argument is provided + if help
die "\n Script gap_density.pl version $version\n\n" if ((! $in) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die "\n please provide input file(s)\n $usage" if (! $in);
die "\n please set -type as  \"ecd\", \"box\" or \"den\" (see usage)" if (($type ne "ecd") && ($type ne "box") && ($type ne "den"));

#avoid / at the end of path in case it's a directory
$in = $1 if (($dir) && ($in =~ /^(.*)\/$/));

#"log"
print STDERR "\n --- Script maf_get_large_indels_Rplot.pl started (v$version)\n" if ($v);
print STDERR "      - Directory containing input files = $in (-dir chosen)\n" if (($dir) && ($v));
print STDERR "      - Input file = $in\n" if ((! $dir) && ($v));
die "\n $in does not exist?\n\n" if (($in) && (! -e $in));
print STDERR "      - R command lines will be printed in $in.R.$type.txt, for:\n" if ($v);
print STDERR "         -> empirical cumulative distributions (log scale for x axis)\n" if (($type eq "ecd") && ($v));
print STDERR "         -> boxplots with notch and without outliers\n" if (($type eq "box") && ($v));
print STDERR "         -> density lines\n" if (($type eq "den") && ($v));
print STDERR "      - with x axis minimum anx maximum values = $xlim\n" if (($xlim) && ($v));
print STDERR "      - with x axis minimum anx maximum values = min and max of each dataset\n" if ((! $xlim) && ($v));
print STDERR "      - plots will be output in pdf format (one per input file)\n" if (($pdf) && ($v));
print STDERR "      - There will be an additional output with some stats = $in.summary.tab with\n" if (($stats) && ($v));
$xlim = "min,max" unless ($xlim);

################################################################################
# MAIN
################################################################################
print STDERR "\n --- Getting list of files from $in...\n" if (($dir) && ($v));
#Get list of input files if -dir, or just load the one from -in
my @files = ();
($dir)?(@files = `ls $in/*_del.tab | grep -v _del.tab.pdf | grep -v _del.tab.len.txt`):(push(@files,$in)); #grep -v in case previous outputs of this script
my $path = path($in);
print STDERR "       ..done\n" if (($dir) && ($v));

print STDERR "\n --- Now extracting lengths from file(s)\n" if ($v);
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
		# ref_chr   ref_st      ref_end     Query      Q_chr    Q_location   Q_indel_ID     len_lc   len_no_lc   insertion_len_in Q
        # chr11     48901944    48909865    gorGor3    chr11    49685897     indel#22008    7921     7195       na

		#here, just need the length - use no lc, remove insertion if relevant
		my @l = split('\t',$line);
		my ($len,$ins) = ($l[8],$l[9]);
		$len = $len - $ins unless ($ins eq "na");
		print $ofh "$len\n" unless ($len < 1000); #fake bug fix here, because cases of insertions > deletion length in ref were printed
	}
	close $ifh;
	close $ofh;
	push (@glist,$out);
}
print STDERR "     ..done\n\n" if ($v);

#Get stats if relevant
if ($stats) {
	print STDERR " --- Getting summary values\n" if ($v); 
	get_stats($in,\@glist);
}

#print R command lines 
print STDERR " --- Printing R command lines\n" if ($v); 
get_Rcmdlines($in,\@glist,$dir,$type,$xlim);

#get plots through R
print STDERR " --- Getting R plots\n" if (($v) && ($pdf)); 
get_Rplots(\@glist,$path,$type,$xlim) if ($pdf);

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
# Get Stats / summary file
# get_stats($in,\@glist);
#----------------------------------------------------------------------------
sub get_stats {
	my ($in,$list) = @_;	
	my $out = $in.".summary.tab";
	`rm $out` if (-e $out);
	open my $ofh, ">>$out" or confess "ERROR (sub get_stats): Failed to open to write $out $!\n";
	print $ofh "#species\tcounts\tsum\taverage\tmedian\tminimum\tmaximum\n\n";
	close $ofh;
	foreach my $f (@{ $list }) {
		my $sp = filename($f);
		$sp =~ s/(.+?)\..+?_del\.tab\.len\.txt/$1/;
		my $sum = 0;
		my @values = ();
		open my $ifh, "<", $f or confess "\nERROR (sub get_stats): can't open to read $f $!\n";	
		while(<$ifh>) {
			chomp(my $val = $_);
			push(@values,$val);
			$sum+=$val;
		}
		close $ifh;
		#get values
		@values = sort {$a cmp $b} @values;
		my $min = $values[0];
		my $max = $values[-1];
		my $count = @values;
		my $avg = $sum/$count;
		my $median = median(\@values);
		#now print
		open my $ofh, ">>$out" or confess "ERROR (sub get_stats): Failed to open to write $out $!\n";
		print $ofh "$sp\t$count\t$sum\t$avg\t$median\t$min\t$max\n";
		close $ofh;
	}	
	return;
}

sub median {
	my ($array_ref) = @_; 
	my $count = scalar @$array_ref;
	my @array = sort { $a <=> $b } @$array_ref; 
	if ($count % 2) { 
		return $array[int($count/2)]; 
	} else { 
		return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
	} 
} 

#----------------------------------------------------------------------------
# print R command lines to plot coverages
# get_Rcmdlines($in,\@glist,$dir,$type,$xlim);
#----------------------------------------------------------------------------
sub get_Rcmdlines {
	my ($in,$list,$dir,$type,$xlim) = @_;	
	my $out = $in.".R.$type.txt";
	if ($dir eq "y") {
		open my $rfh, ">$out" or confess "ERROR (sub get_Rcmdlines): Failed to open to write $out $!\n";
		print $rfh "\nlibrary(calibrate)\n\n";
		print $rfh "setwd(\"\") #fill with location of files on computer\n\n";		
		close $rfh;
	}
	#now loop
	my $nb = 1;
	my $i = 1;
	foreach my $f (@{ $list }) {
		$out = $f.".R.txt" if ($dir eq "n");
		open my $rfh, ">>$out" or confess "ERROR (sub get_Rcmdlines): Failed to open to write $out $!\n";
		my $sp = filename($f);		
		$sp = $1 if ($sp =~ /(.+?)\..+?_del\.tab/);		
		print $rfh "\npar(mfrow=c(4,2))\n" if (($nb == 1) && ($dir eq "y"));
		print $rfh "\nlibrary(calibrate)\n\n" if ($dir eq "n");	
		print $rfh "setwd(\"\") #fill with location of files on computer\n\n" if ($dir eq "n");		
		print $rfh "$sp<-read.table(\"$f\", header=FALSE)\n";
		if ($type eq "ecd") {
			print $rfh "nb <- length($sp\$V1)\n";
			print $rfh "if (nb > 1) max<-sort($sp\$V1,partial=nb-1)[nb] else max<-$sp\$V1\n"; 
			print $rfh "if (nb > 1) min<-sort($sp\$V1)[1] else min<-$sp\$V1\n";
			print $rfh "plot(ecdf($sp\$V1), log=\"x\", xlim=c($xlim), main=paste(\"$sp (n=\",nb,\")\"), xlab=\"indel size = length in reference without lower cases (nt)\", ylab=\"Empirical cumulative distribution\",lwd=1,cex.lab=1.2,cex.axis=1.2,cex.main=1.2,cex.sub=0.8)\n";						
			print $rfh "plot(ecdf(max),add=T)\n";
			print $rfh "textxy(max,1,max,cex=1,pos=3)\n";
			print $rfh "textxy(min,1,min,cex=1,pos=3)\n";
		} elsif ($type eq "den") {
			print $rfh "nb <- length($sp\$V1)\n";
			print $rfh "d<-density(dat\$V1)\n";
			print $rfh "plot(d, col=\"blue\", pch = 18, main=paste(\"$sp (n=\",nb,\")\"), ylab = \"Density of x=deletion size (nt)\"),xlim=c($xlim)\n";
		} elsif ($type eq "box") {
			print $rfh "boxplot(dat\$V1, notch=T, outline=F)\n";
			print $rfh "title(paste(\"$sp (n=\",nb,\")\"))\n";
		}
		$i++;	
		$nb = 0 if ($nb == 8);
		$nb++;
		close $rfh;
	}	
	#Table[!is.na(Table$Col),])
	return;
}

#----------------------------------------------------------------------------
# get R plots
# get_Rplots(\@glist,$path,$type,$xlim) if ($pdf);
#----------------------------------------------------------------------------
sub get_Rplots {
	my ($list,$dir,$type,$xlim) = @_;
	#Start R bridge
	my $R = Statistics::R->new();
	$R->startR;
	$R->send(qq`library(calibrate)`);
	#$R->send(q`setwd(".")`);
	#Now loop
	FILE: foreach my $f (@$list){
		chomp $f;		
		#plot
		my $sp = filename($f);
		$sp = $1 if ($sp =~ /(.+?)\..+?_del\.tab/);	
		$R->send(qq`$sp <- read.table("$f",header=FALSE)`);
		my $out = $f.".".$type.".pdf"; 
		$R->run(qq`pdf("$out")`);	
		if ($type eq "ecd") {
			$R->send(qq`nb <- length($sp\$V1)`);
			$R->send(qq`if (nb > 1) max<-sort($sp\$V1,partial=nb-1)[nb] else max<-$sp\$V1`); 
			$R->send(qq`if (nb > 1) min<-sort($sp\$V1)[1] else min<-$sp\$V1`);			
			$R->send(qq`plot(ecdf($sp\$V1), log="x", xlim=c($xlim), main=paste("$sp (n=",nb,")"), xlab="indel size = length in reference without lower cases (nt)", ylab="Empirical cumulative distribution",lwd=1,cex.lab=1.2,cex.axis=1.2,cex.main=1.2,cex.sub=0.8)`);
			$R->send(qq`plot(ecdf(max),add=T)`);
			$R->send(qq`textxy(max,1,max,cex=1,pos=3)`);
			$R->send(qq`textxy(min,1,min,cex=1,pos=3)`);	
		} elsif ($type eq "den") {
			$R->send(qq`nb<-length(dat\$V1)`);
			$R->send(qq`d<-density(dat\$V1)`);
			$R->send(qq`plot(d, col=\"blue\", pch = 18, main=paste("$sp (n=",nb,")"), ylab = "Density of x=deletion size (nt)"),xlim=c($xlim)`);
		} elsif ($type eq "box") {
			$R->send(qq`boxplot(dat\$V1,notch = T,outline=F)`);
			$R->send(qq`title(paste("$sp (n=",nb,")")))`);
		}
		$R->run(q`dev.off()`);
	}	
	#End R bridge
	$R->stopR() ;
	return;
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






