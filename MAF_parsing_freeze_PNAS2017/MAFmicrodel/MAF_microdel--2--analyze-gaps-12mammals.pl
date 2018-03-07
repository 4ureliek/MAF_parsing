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
print "\n --- Analyzing gaps (MultiZ 100 Way - 12 species extracted),
        According to the following phylogeny:

                             |-- hg19
                          |--|             2 Primates 
         SUPRAPRIMATES    |  |-- rheMac3
                    |-----|
                    |     |  |-- mm10
                    |     |--|             2 Rodents
        BOREOTHERIA |        |-- rn5
                 |--|
                 |  |        |-- canFam3
                 |  |  |-----|             2 Carnivora
                 |  |  |     |-- felCat5
        EUTHERIA |  |  |                   LAURASIATHERIA
              |--|  |  |     |-- myoLuc2					
              |  |  |--|-----|             2 Bats
              |  |     |     |-- pteVam1
              |  |     |
              |  |     |-------- bosTau7    1 Bovidae
              |  |
            --|  |  |----------- loxAfr3										
              |  |--|                      AFROTHERIA			
              |     |----------- echTel2
              |  
              |----------------- monDom5 = Opossum (outgroup)\n\n\n";

			  
# DEFINE LOOPS => possible to extend to more species more easily
#NOTE : improve here with subtraction of lists
# use Acme::Tools qw(minus);
# my @diff = minus(\@a, \@bl);

my @eutherian = ("hg19","rheMac3","mm10","rn5","canFam3","felCat5","myoLuc2","pteVam1","bosTau7","loxAfr3","echTel2");
my @not_eutherian = ("monDom5");

my @boreo = ("hg19","rheMac3","mm10","rn5","canFam3","felCat5","myoLuc2","pteVam1","bosTau7");
my @not_boreo = ("loxAfr3","echTel2","monDom5");

my @afro = ("loxAfr3","echTel2");
my @not_afro = ("hg19","rheMac3","mm10","rn5","canFam3","felCat5","myoLuc2","pteVam1","bosTau7","monDom5");

my @supra = ("hg19","rheMac3","mm10","rn5");
my @not_supra = ("canFam3","felCat5","myoLuc2","pteVam1","bosTau7","loxAfr3","echTel2","monDom5");

my @laura = ("canFam3","felCat5","myoLuc2","pteVam1","bosTau7");
my @not_laura = ("hg19","rheMac3","mm10","rn5","loxAfr3","echTel2","monDom5");

my @prim = ("hg19","rheMac3");
my @not_prin = ("mm10","rn5","canFam3","felCat5","myoLuc2","pteVam1","bosTau7","loxAfr3","echTel2","monDom5");

my @rod = ("mm10","rn5");
my @not_rod = ("hg19","rheMac3","canFam3","felCat5","myoLuc2","pteVam1","bosTau7","loxAfr3","echTel2","monDom5");

my @carn = ("canFam3","felCat5");
my @not_carn = ("hg19","rheMac3","mm10","rn5","myoLuc2","pteVam1","bosTau7","loxAfr3","echTel2","monDom5");

my @bats = ("myoLuc2","pteVam1");
my @not_bats = ("hg19","rheMac3","mm10","rn5","canFam3","felCat5","bosTau7","loxAfr3","echTel2","monDom5");


my $files = ();
print STDERR "     - subtract all gaps shared by all species, including outgroup...\n";
print STDERR "        (basically it means clean the previous files that have a lot of noise, since it contains all other gaps from species in aln that are removed now)\n";
$files = MAFmicrodel::split_gaps($spIDs,"no","gaps.1-30","all",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by all species (EUTHERIAN) but not with outgroup [not orientable]..\n";	
$files = MAFmicrodel::split_gaps(\@eutherian,\@not_eutherian,"gaps.1-30.not-shared.all","eutherian",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by BOREOTHERIAN but not with any others\n";
$files = MAFmicrodel::split_gaps(\@boreo,\@not_boreo,"gaps.1-30.not-shared.eutherian","boreo",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by AFROTHERIAN but not with any others\n";
$files = MAFmicrodel::split_gaps(\@afro,\@not_afro,"gaps.1-30.not-shared.eutherian","afro",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by SUPRAPRIMATES but not with any others\n";
$files = MAFmicrodel::split_gaps(\@supra,\@not_supra,"gaps.1-30.not-shared.boreo","supra",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by LAURASIATHERIAN but not with any others\n";
$files = MAFmicrodel::split_gaps(\@laura,\@not_laura,"gaps.1-30.not-shared.boreo","laura",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by PRIMATES but not with any others\n";
$files = MAFmicrodel::split_gaps(\@prim,\@not_prin,"gaps.1-30.not-shared.supra","prim",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by RODENTS but not with any others\n";
$files = MAFmicrodel::split_gaps(\@rod,\@not_rod,"gaps.1-30.not-shared.supra","rod",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by CARNIVORA but not with any others\n";
$files = MAFmicrodel::split_gaps(\@carn,\@not_carn,"gaps.1-30.not-shared.laura","carn",$in,$path,$files,$bedtools); 

print STDERR "     - split gaps between ones shared by CHIROPTERA (bats) but not with any others\n";
$files = MAFmicrodel::split_gaps(\@bats,\@not_bats,"gaps.1-30.not-shared.laura","bats",$in,$path,$files,$bedtools); 



##########################################################################################################
# III. Now get species spe gaps and print
##########################################################################################################
$files = MAFmicrodel::spe_gaps($spIDs,$in,$path,$files,$bedtools,$v);
MAFmicrodel::print_amounts($path,$totlen,$files,$v);
print STDERR " --- Script is done\n\n" if ($v);
exit;
