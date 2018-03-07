## MAF_parsing_freeze_PNAS2017
This folder contains the code used in [this publication:](http://www.pnas.org/content/114/8/E1460)
`Kapusta, Suh & Feschotte (2017) PNAS doi: 10.1073/pnas.1616702114`

MAF_microdel.pl v4.0 gives identical results, but is WAY easier to run.

## MAFmicrodel
From a multi species alignment in maf format and a newick tree, this script outputs microdeletions (1 to 30nt) of the various branches.
(See its own README for more details)[https://github.com/4ureliek/MAF_parsing/blob/master/MAFmicrodel/README.md]

## maf_get_large_indels.pl
	USAGE:
	perl <scriptname.pl> -dir <dir_with_alignments> [-sp <SP>] [-lenc <X>] [-lene <X>] [-cpu <X>] [-v] [-h|help]
    	
	Outputs:
	I.  List of sizes of all consecutive empty data for a given species (no data on the browser)
	    Only when the empty data is interrupting a scaffold (note that it could correspond to a misassembly)
	    Empty data = no \"s\" line. However, there could be insertion (non aligning bases): in that case,
	    the large indel will be printed only if lower case length - insertion length > lene
	      
	II. List of sizes of all consecutive empty blocks for a given species (C lines in the maf = continuous:
	    "C -- the sequence before and after is contiguous implying that this region was either deleted in the 
	    source or inserted in the reference sequence. The browser draws a single line or a '-' in base mode in these blocks."

	In both cases, outputs will be:
	      [0]     [1]       [2]     [3]                 [4]                  [5]            [6]               [7]
	      chr(SP) start(SP) end(SP) chr/scaffold(Query) DeletionCoord(Query) length_with_lc length_without_lc length_insertion

	      lc = lower cases. In column 6, lower cases in reference are removed 
	           Length in col 6 would be the minimum deletion len (underestimated when ancient TEs, 
	           but this will remove all TE insertions in the reference that would not be a deletion).
