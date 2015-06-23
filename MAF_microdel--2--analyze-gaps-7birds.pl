#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta - https://github.com/4ureliek
# email   :  4urelie.k@gmail.com
# version :  see MAFmicrodel.pm
# purpose :  see MAFmicrodel.pm
# usage   :  see MAFmicrodel.pm
##########################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::AlignIO;
use vars qw($BIN);
use Cwd 'abs_path';
BEGIN { 	
	$BIN = abs_path($0);
	$BIN =~ s/(.*)\/.*$/$1/;
	unshift(@INC, "$BIN");
}
use MAFmicrodel;

##########################################################################################################
# I. Get options and prep stuff
##########################################################################################################
my $v;
my ($in,$out,$nbsp,$aln,$concat,$bedtools,$help,$chlog) = ("na","na","na","na","na","na","na","na");
GetOptions ('in=s' => \$in, 'sp=s' => \$nbsp, 'concat' => \$concat, 'aln=s' => \$aln, 'bed=s' => \$bedtools, 'out=s' => \$out, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

my $path = MAFmicrodel::check_options($bedtools,$in,$out,$nbsp,$aln,$concat,$help,$chlog,$v);

# Get total length of alignement and list of species, using $amount files
my ($totlen,$spIDs) = MAFmicrodel::get_amounts($in,$aln,$nbsp,$path,$v); #totlen = hash and spIDs = list

# Concatenate files, or test if it should be
MAFmicrodel::concat_gaps($spIDs,$in,$concat,$v);


##########################################################################################################
# II. Now deal with the gaps. THIS IS THE VARIABLE PART.
##########################################################################################################
print STDERR "\n --- Analyzing gaps [from UCSC MultiZ 100 ways, 7 birds]
        According to the following phylogeny:

                                  |-------- falPer1     Falco peregrinus = Peregrine falcon 
                       AUSTRALAVES|  FINCH
                               |--|     |-- taeGut2     Taeniopygia guttata = Zebra finch
                               |  |  |--|
                               |  |  |  |-- geoFor1     Geospiza fortis = Medium ground finch 
                        PASSERA|  |--|        
                            |--|     |----- melUnd1     Melopsittacus undulatus = Budgerigar (common pet parakeet)
                            |  |  PSITTACOPASSERAE 
                            |  |
                     NEOAVES|  |----------- colLiv1     Columbia livia = Rock pigeon              
                         |--|        
                         |  |           |-- anaPla1     Anas platyrynchos = Mallard duck              
                       --|  |-----------|
                         |  GALLO       |-- galGal4     Gallus gallus = Chicken   
                         |
                         |----------------- anoCar2     Anolis carolinensis = Lizard (outgroup)\n\n\n" if ($v);

			  
# DEFINE LOOPS => possible to extend to more species more easily
#NOTE : improve here with subtraction of lists
# use Acme::Tools qw(minus);
# my @diff = minus(\@a, \@bl);

my @aves = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1","anaPla1","galGal4");
my @not_aves = ("anoCar2");

my @pass = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1");
my @not_pass = ("anaPla1","galGal4","anoCar2");
my @gallo = ("anaPla1","galGal4");
my @not_gallo = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1","anoCar2");

my @austr = ("falPer1","taeGut2","geoFor1","melUnd1");
my @not_austr = ("colLiv1","anaPla1","galGal4","anoCar2");

my @psitta = ("taeGut2","geoFor1","melUnd1");
my @not_psitta = ("falPer1","colLiv1","anaPla1","galGal4","anoCar2");

my @finch = ("taeGut2","geoFor1");
my @not_finch = ("melUnd1","falPer1","colLiv1","anaPla1","galGal4","anoCar2");


my $files = ();
print STDERR "     - subtract all gaps shared by all species, including outgroup...\n" if ($v);
print STDERR "        (basically it means clean the previous files that have a lot of noise, since it contains all other gaps from species in aln that are removed now)\n" if ($v);
$files = MAFmicrodel::split_gaps($spIDs,"no","gaps.1-30","all",$in,$path,$files,$bedtools); 
print STDERR "       ..done\n" if ($v);

print STDERR "     - split gaps between ones shared by all species (NEOAVES) but not with outgroup [not orientable]..\n" if ($v);
$files = MAFmicrodel::split_gaps(\@aves,\@not_aves,"gaps.1-30.not-shared.all","neoaves",$in,$path,$files,$bedtools); 
print STDERR "       ..done\n" if ($v);

print STDERR "     - split gaps between ones shared by all PASSERA (\"pass\") but not with any others\n" if ($v);
$files = MAFmicrodel::split_gaps(\@pass,\@not_pass,"gaps.1-30.not-shared.neoaves","pass",$in,$path,$files,$bedtools);
print STDERR "       ..done\n" if ($v);
print STDERR "     - split gaps between ones shared by all GALLO (\"gallo\") but not with any others\n" if ($v);
$files = MAFmicrodel::split_gaps(\@gallo,\@not_gallo,"gaps.1-30.not-shared.neoaves","gallo",$in,$path,$files,$bedtools);
print STDERR "       ..done\n" if ($v);	

print STDERR "     - split gaps between ones shared by all AUSTRALAVES (\"austr\") but not with any others\n" if ($v);
$files = MAFmicrodel::split_gaps(\@austr,\@not_austr,"gaps.1-30.not-shared.pass","austr",$in,$path,$files,$bedtools);
print STDERR "       ..done\n" if ($v);

print STDERR "     - split gaps between ones shared by all PSITTACOPASSERAE but not with any others\n" if ($v);
$files = MAFmicrodel::split_gaps(\@psitta,\@not_psitta,"gaps.1-30.not-shared.austr","psitta",$in,$path,$files,$bedtools);
print STDERR "       ..done\n" if ($v);

print STDERR "     - split gaps between ones shared by all FINCH but not with any others\n" if ($v);
$files = MAFmicrodel::split_gaps(\@finch,\@not_finch,"gaps.1-30.not-shared.psitta","finch",$in,$path,$files,$bedtools);
print STDERR "       ..done\n" if ($v);



##########################################################################################################
# III. Now get species spe gaps and print
##########################################################################################################
$files = MAFmicrodel::spe_gaps($spIDs,$in,$path,$files,$bedtools,$v);
MAFmicrodel::print_amounts($path,$totlen,$files,$v);
print STDERR " --- Script is done\n\n";
exit;

