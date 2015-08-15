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

my $path;
($path,$in,$out,$bedtools) = MAFmicrodel::check_options($bedtools,$in,$out,$nbsp,$aln,$concat,$help,$chlog,$v);

# Get total length of alignement and list of species, using $amount files
my ($totlen,$spIDs) = MAFmicrodel::get_amounts($in,$aln,$nbsp,$path,$v); #totlen = hash and spIDs = list

# Concatenate files, or test if it should be
MAFmicrodel::concat_gaps($spIDs,$in,$concat,$v);


##########################################################################################################
# II. Now deal with the gaps. THIS IS THE VARIABLE PART.
##########################################################################################################
print STDERR "\n --- Analyzing gaps [from file birds_reestimated_archosaurRef.13sp.maf]
        According to the following phylogeny:

                                  |-------- falPer1     Falco peregrinus = Peregrine falcon 
                       AUSTRALAVES|  FINCH
                               |--|     |-- taeGut2     Taeniopygia guttata = Zebra finch
                               |  |  |--|
                               |  |  |  |-- geoFor1     Geospiza fortis = Medium ground finch 
                        NEOAVES|  |--|        
                            |--|     |----- melUnd1     Melopsittacus undulatus = Budgerigar (common pet parakeet)
                            |  |  PASSERIMORPHAE 
                            |  |
                  NEOGNATHAE|  |----------- colLiv1     Columbia livia = Rock pigeon              
                         |--|        
                         |  |           |-- anaPla1     Anas platyrynchos = Mallard duck              
                         |  |-----------|
                     AVES|  GALLO       |-- galGal4     Gallus gallus = Chicken  
                      |--|  
      ARCHOSAURIFORMES|  |----------------- strCam0     Struthio camelus = Ostrich                  
                   |--|   	     
                   |  |              |----- allMis2     Alligator mississippiensis = American alligator
                   |  |--------------|  
   ARCHOSAUROMORPHA|  CROCS          |  |-- ghaGan1     Gavialis gangeticus = Gharial
                |--|                 |--|
                |  |                    |-- croPor2     Crocodylus porosus = Crocodile
                |  |     
                |  |----------------------- chrPic1     Chrysemys picta bellii = Painted turtle                  
                |
                |-------------------------- anoCar2     Anolis carolinensis = Lizard (outgroup)\n\n\n";

			  
# DEFINE LOOPS => possible to extend to more species more easily
#NOTE : improve here with subtraction of lists
# use Acme::Tools qw(minus);
# my @diff = minus(\@a, \@bl);

my @archo1 = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1","anaPla1","galGal4","strCam0","allMis2","ghaGan1","croPor2","chrPic1");
my @not_archo1 = ("anoCar2");

my @archo2 = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1","anaPla1","galGal4","strCam0","allMis2","ghaGan1","croPor2");
my @not_archo2 = ("chrPic1","anoCar2");

my @crocs = ("allMis2","ghaGan1","croPor2");
my @not_crocs = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1","anaPla1","galGal4","strCam0","chrPic1","anoCar2");
my @aves = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1","anaPla1","galGal4","strCam0");
my @not_aves = ("allMis2","ghaGan1","croPor2","chrPic1","anoCar2");

my @neognathae = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1","anaPla1","galGal4");
my @not_neognathae = ("strCam0","allMis2","ghaGan1","croPor2","chrPic1","anoCar2");

my @neoaves = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1");
my @not_neoaves = ("anaPla1","galGal4","strCam0","allMis2","ghaGan1","croPor2","chrPic1","anoCar2");
my @gallo = ("anaPla1","galGal4");
my @not_gallo = ("falPer1","taeGut2","geoFor1","melUnd1","colLiv1","strCam0","allMis2","ghaGan1","croPor2","chrPic1","anoCar2");

my @austr = ("falPer1","taeGut2","geoFor1","melUnd1");
my @not_austr = ("colLiv1","anaPla1","galGal4","strCam0","allMis2","ghaGan1","croPor2","chrPic1","anoCar2");

my @passeri = ("taeGut2","geoFor1","melUnd1");
my @not_passeri = ("falPer1","colLiv1","anaPla1","galGal4","strCam0","allMis2","ghaGan1","croPor2","chrPic1","anoCar2");

my @finch = ("taeGut2","geoFor1");
my @not_finch = ("melUnd1","falPer1","colLiv1","anaPla1","galGal4","strCam0","allMis2","ghaGan1","croPor2","chrPic1","anoCar2");



my $files = ();
print STDERR "     - subtract all gaps shared by all species, including outgroup...\n";
print STDERR "        (basically it means clean the previous files that have a lot of noise, since it contains all other gaps from species in aln that are removed now)\n";
$files = MAFmicrodel::split_gaps($spIDs,"no","gaps.1-30","all",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by all species (\"archo1\") but not with outgroup [not orientable]..\n";	
$files = MAFmicrodel::split_gaps(\@archo1,\@not_archo1,"gaps.1-30.not-shared.all","archo1",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by all ARCHOSAURIFORMES (\"archo2\" = aves+crocs) but not with any others\n";	
$files = MAFmicrodel::split_gaps(\@archo2,\@not_archo2,"gaps.1-30.not-shared.archo1","archo2",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by all AVES but not with any others\n";
$files = MAFmicrodel::split_gaps(\@aves,\@not_aves,"gaps.1-30.not-shared.archo2","aves",$in,$path,$files,$bedtools);
print STDERR "     - split gaps between ones shared by all CROCS but not with any others\n";
$files = MAFmicrodel::split_gaps(\@crocs,\@not_crocs,"gaps.1-30.not-shared.archo2","crocs",$in,$path,$files,$bedtools);

print STDERR "     - split gaps between ones shared by all NEOGNATHAE but not with any others\n";
$files = MAFmicrodel::split_gaps(\@neognathae,\@not_neognathae,"gaps.1-30.not-shared.aves","neognathae",$in,$path,$files,$bedtools);

print STDERR "     - split gaps between ones shared by all NEOAVES (\"pass\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@neoaves,\@not_neoaves,"gaps.1-30.not-shared.neognathae","neoaves",$in,$path,$files,$bedtools);
print STDERR "     - split gaps between ones shared by all GALLO (\"gallo\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@gallo,\@not_gallo,"gaps.1-30.not-shared.neognathae","gallo",$in,$path,$files,$bedtools);

print STDERR "     - split gaps between ones shared by all AUSTRALAVES (\"austr\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@austr,\@not_austr,"gaps.1-30.not-shared.neoaves","austr",$in,$path,$files,$bedtools);
print STDERR "     - split gaps between ones shared by all PASSERIMORPHAE but not with any others\n";
$files = MAFmicrodel::split_gaps(\@passeri,\@not_passeri,"gaps.1-30.not-shared.neoaves","passeri",$in,$path,$files,$bedtools);

print STDERR "     - split gaps between ones shared by all FINCH but not with any others\n";
$files = MAFmicrodel::split_gaps(\@finch,\@not_finch,"gaps.1-30.not-shared.passeri","finch",$in,$path,$files,$bedtools);



##########################################################################################################
# III. Now get species spe gaps and print
##########################################################################################################
$files = MAFmicrodel::spe_gaps($spIDs,$in,$path,$files,$bedtools,$v);
MAFmicrodel::print_amounts($path,$totlen,$files,$v);
print STDERR " --- Script is done\n\n" if ($v);
exit;




