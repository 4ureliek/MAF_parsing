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
print STDERR "\n --- Analyzing gaps [from UCSC MultiZ 100 ways, 18 reptiles]
        According to the following phylogeny: 
				
							FALCONS			  |-- falChe1 	Falco cherrug = Saker falcon
								  |-----------|               
								  |           |-- falPer1	Falco peregrinus = Peregrine falcon
								  |
								  |     |-------- pseHum1	Pseudopodoces humilis = Ground tit = Tibetan ground jay = Hume's ground-tit
		  	 PASSERA/AUSTRALAVES  |     |
							   |--|     |  |----- ficAlb2	Ficedula albicollis = Collared flycatcher
							   |  PSITTA|  |
							   |  |  |--|  | |--- taeGut2	Taeniopygia guttata = Zebra finch
							   |  |  |  |--| |
							   |  |--|     |-| |- geoFor1	Geospiza fortis = Medium ground finch
				       PASSERIMORPHA |       |-|
				      PASSERA  |     |         |- zonAlb1	Zonotrichia albicollis = White throated sparrow
							|--|     |
							|  |     |     |----- melUnd1	Melopsittacus undulatus = Budgerigar
							|  |     |-----|
			    AVES/NEOAVES|  |  PARROTS  |  |-- amaVit1	Amazona vittata = Parrot   
					  |-----|  |           |--|
					  |     |  |              |-- araMac1	Ara macao = Scarlet macaw
				      |     |  |
					  |     |  |----------------- colLiv1	Columbia livia = Rock pigeon
	 ARCHOSAURIFORMES |     |            
				   |--|     |              |----- anaPla1	Anas platyrynchos = Mallard duck
				   |  |     |--------------|
				   |  |   GALLO            |----- galGal4	Gallus gallus = Chicken  
     ARCHOSOMORPHA |  |                   
				|--|  |-------------------------- allMis1	Alligator mississippiensis = American alligator
				|  |                   TURTLES                
				|  |                       |----- pelSin1	Pelodiscus sinensis = Soft-shell turtle 
				|  |-----------------------|                   
				|                          |  |-- cheMyd1	Chelonia mydas = Green seaturtle
				|                          |--|
				|                             |-- chrPic1	Chrysemys picta bellii = Painted turtle
				|
				|-------------------------------- anoCar2	Anolis carolinensis = Lizard
\n\n";

		  
# DEFINE LOOPS => possible to extend to more species more easily
#NOTE : improve here with subtraction of lists
# use Acme::Tools qw(minus);
# my @diff = minus(\@a, \@bl);

my @archo1 = ("falPer1","falChe1","ficAlb2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1","colLiv1","anaPla1","galGal4","taeGut2","allMis1","chrPic1","cheMyd1","pelSin1");
my @not_archo1 = ("anoCar2");

my @turtles = ("chrPic1","cheMyd1","pelSin1");
my @not_turtles = ("falPer1","falChe1","ficAlb2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1","colLiv1","anaPla1","galGal4","taeGut2","allMis1","anoCar2");
my @turtles2 = ("chrPic1","cheMyd1");
my @not_turtles2 = ("falPer1","falChe1","ficAlb2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1","colLiv1","anaPla1","galGal4","taeGut2","allMis1","pelSin1","anoCar2");

my @archo2 = ("falPer1","falChe1","ficAlb2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1","colLiv1","anaPla1","galGal4","taeGut2","allMis1");
my @not_archo2 = ("chrPic1","cheMyd1","pelSin1","anoCar2");

my @aves = ("falPer1","falChe1","ficAlb2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1","colLiv1","anaPla1","galGal4","taeGut2");
my @not_aves = ("allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");

my @gallo = ("anaPla1","galGal4");
my @not_gallo = ("falPer1","falChe1","ficAlb2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1","colLiv1","taeGut2","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");
my @passera = ("falPer1","falChe1","ficAlb2","taeGut2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1","colLiv1");
my @not_passera = ("anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");

my @austr = ("falPer1","falChe1","ficAlb2","taeGut2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1");
my @not_austr = ("colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");

my @passeri = ("ficAlb2","taeGut2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1");
my @not_passeri = ("falPer1","falChe1","colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");
my @falcons = ("falPer1","falChe1");
my @not_falcons = ("ficAlb2","taeGut2","zonAlb1","geoFor1","pseHum1","melUnd1","amaVit1","araMac1","colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");

my @psitta = ("ficAlb2","taeGut2","zonAlb1","geoFor1","pseHum1");
my @not_psitta = ("melUnd1","amaVit1","araMac1","falPer1","falChe1","colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");
my @parrots = ("melUnd1","amaVit1","araMac1");
my @not_parrots = ("ficAlb2","taeGut2","zonAlb1","geoFor1","pseHum1","falPer1","falChe1","colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");

my @parrots2 = ("amaVit1","araMac1");
my @not_parrots2 = ("melUnd1","ficAlb2","taeGut2","zonAlb1","geoFor1","pseHum1","falPer1","falChe1","colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");

my @finch1 = ("ficAlb2","taeGut2","zonAlb1","geoFor1");
my @not_finch1 = ("pseHum1","melUnd1","amaVit1","araMac1","falPer1","falChe1","colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");
my @finch2 = ("taeGut2","zonAlb1","geoFor1");
my @not_finch2 = ("ficAlb2","pseHum1","melUnd1","amaVit1","araMac1","falPer1","falChe1","colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");
my @finch3 = ("zonAlb1","geoFor1");
my @not_finch3 = ("taeGut2","ficAlb2","pseHum1","melUnd1","amaVit1","araMac1","falPer1","falChe1","colLiv1","anaPla1","galGal4","allMis1","chrPic1","cheMyd1","pelSin1","anoCar2");


my $files = ();
print STDERR "     - subtract all gaps shared by all species, including outgroup...\n";
print STDERR "        (basically it means clean the previous files that have a lot of noise, since it contains all other gaps from species in aln that are removed now)\n";
$files = MAFmicrodel::split_gaps($spIDs,"no","gaps.1-30","all",$in,$path,$files,$bedtools); 
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all species (\"archo1\") but not with outgroup [not orientable]..\n";	
$files = MAFmicrodel::split_gaps(\@archo1,\@not_archo1,"gaps.1-30.not-shared.all","archo1",$in,$path,$files,$bedtools); 
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all TURTLES but not with any others\n";	
$files = MAFmicrodel::split_gaps(\@turtles,\@not_turtles,"gaps.1-30.not-shared.archo1","turtles",$in,$path,$files,$bedtools); 
print STDERR "       ..done\n";	
print STDERR "     - split gaps between ones shared by 2 TURTLES (\"turtles2\")but not with any others\n";	
$files = MAFmicrodel::split_gaps(\@turtles2,\@not_turtles2,"gaps.1-30.not-shared.turtles","turtles2",$in,$path,$files,$bedtools); 
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all ARCHOSAURIFORMES (\"archo2\" = aves+crocs) but not with any others\n";	
$files = MAFmicrodel::split_gaps(\@archo2,\@not_archo2,"gaps.1-30.not-shared.archo1","archo2",$in,$path,$files,$bedtools); 
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all (NEO)AVES but not with any others\n";
$files = MAFmicrodel::split_gaps(\@aves,\@not_aves,"gaps.1-30.not-shared.archo2","aves",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all PASSERA but not with any others\n";
$files = MAFmicrodel::split_gaps(\@passera,\@not_passera,"gaps.1-30.not-shared.aves","passera",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	
print STDERR "     - split gaps between ones shared by all GALLO (\"gallo\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@gallo,\@not_gallo,"gaps.1-30.not-shared.aves","gallo",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all AUSTRALAVES (\"austr\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@austr,\@not_austr,"gaps.1-30.not-shared.pass","austr",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all PASSERIMORPHA (\"passeri\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@passeri,\@not_passeri,"gaps.1-30.not-shared.austr","passeri",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	
print STDERR "     - split gaps between ones shared by all FALCONS but not with any others\n";
$files = MAFmicrodel::split_gaps(\@falcons,\@not_falcons,"gaps.1-30.not-shared.austr","falcons",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all PSITTACOPASSERAE (\"psitta\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@psitta,\@not_psitta,"gaps.1-30.not-shared.passeri","psitta",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	
print STDERR "     - split gaps between ones shared by all PARROTS but not with any others\n";
$files = MAFmicrodel::split_gaps(\@parrots,\@not_parrots,"gaps.1-30.not-shared.passeri","parrots",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	
print STDERR "     - split gaps between ones shared by 2 PARROTS (\"parrots2\") (but not with any others\n";
$files = MAFmicrodel::split_gaps(\@parrots2,\@not_parrots2,"gaps.1-30.not-shared.parrots","parrots2",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	

print STDERR "     - split gaps between ones shared by all 3 FINCH + FLYCATCHER (\"finch1\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@finch1,\@not_finch1,"gaps.1-30.not.psitta.psitta","finch1",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	
print STDERR "     - split gaps between ones shared by the 3 FINCH (\"finch2\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@finch2,\@not_finch2,"gaps.1-30.not.psitta.finch1","finch2",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	
print STDERR "     - split gaps between ones shared by 2 FINCH (\"finch3\") but not with any others\n";
$files = MAFmicrodel::split_gaps(\@finch3,\@not_finch3,"gaps.1-30.not.psitta.finch2","finch3",$in,$path,$files,$bedtools);
print STDERR "       ..done\n";	



##########################################################################################################
# III. Now get species spe gaps and print
##########################################################################################################
$files = MAFmicrodel::spe_gaps($spIDs,$in,$path,$files,$bedtools,$v);
MAFmicrodel::print_amounts($path,$totlen,$files,$v);
print STDERR " --- Script is done\n\n";
exit;



