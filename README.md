MAFmicrodel
=====
Scripts to obtain micro deletions from a multi species alignment in maf format

A gap will be considered specific to a branch only if it is specific to the species considered.

	For examples:
       gap i is shared between all species
       gap ii and iii are specific to species 2 and 1 respectively (not enough overlap)
       gap iv is shared between species 3 and 4

                        i        ii            iii          iv     
        |-- specie 1 ATG-ATGCATGCATGCATGCATGCAT----------ATGCATGCATGCATGCATGCATGC
     |--|
     |  |-- specie 2 ATG-ATGCATGC----------------ATGCATGCATGCATGCATGCATGCATGCATGC
    -|
     |  |-- specie 3 ATG-ATGCATGCATGCATGCATGCATGCATGCATGCATG-----------GCATGCATGC
     |--|
        |-- specie 4 ATG-ATGCATGCATGCATGCATGCATGCATGCATGCATGG-----------CATGCATGC

IMPORTANT: the .maf file need to be grepped for species of interest beforehand

	Typically: grep '^\$\|maf\|^a\|hg19\|panTro4\|ponAbe2\|rheMac3\' chr_all.maf > chr_all.primates.maf

You will also need to make your own version of the script MAF_microdel--2--analyze-gaps-XXX.pl, by modifying the part II (open the file in a text editor). Use the existing ones as templates to edit the lists (the tree is not necessary, it is just fyi).


perl MAF_microdel--1--get-gaps.pl
=====
    USAGE
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
     -h|help => (BOOL)   print this usage


perl MAF_microdel--2--analyze-gaps-XXX.pl
=====
    USAGE
    perl MAF_microdel--2--analyze-gaps-XXX.pl -in <_*.gaps.1-30.bed> -sp <X> [-aln <file.maf>] [-out <path>] [-bed <path>] [-v] [-chlog] [-h]
	
    This script is step 2/2 and it will analyze gaps previously listed by MAF_microdel--1--get-gaps.pl
	
    REQUIREMENTS
    BEDtools is required for this pipeline to run (Quinlan AR and Hall IM, 2010. Bioinformatics)
    Available at: https://github.com/arq5x/bedtools2
       Notes:
        - Two gaps in two species are considered to be shared if they overlap more than 90% 
             (flexibility of 1nt per 10nt of gaps) 
             => option -f 0.90 of intersectBed
        - This min overlap is for BOTH sense
             (if BIG gap in one specie, and a small one in the other, not same gap)
             => option -r of intersectBed
	
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
     -h|help => (BOOL)   print this usage