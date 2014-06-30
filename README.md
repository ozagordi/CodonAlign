## Codon Alignment
Standalone version of [Codon Alignment](http://www.hiv.lanl.gov/content/sequence/CodonAlign/codonalign.html)
from Los Alamos National Laboratory.

### Description
The staff at [LANL](http://www.lanl.gov/) kindly provided me with the
standalone version of their tool for translation aware alignment. The original
tool is written in perl and is able to read sequences in (almost?) any format.
I packaged a minimal standalone version that works only with fasta files by
writing a wrapper in python for the core program `codonalign`.

### Usage
It needs python (2.7) and perl. Only Biopython is needed as external dependency
(though it would be easy to do without). Running `ca_wrapper -h` gives the
usual help
	[user@host CodonAlign]$ ca_wrapper.py -h
	usage: ca_wrapper.py [-h] [-f FILENAME] [-r FRAME] [-s FS_COMP]

	Run codon align as in
	http://www.hiv.lanl.gov/content/sequence/CodonAlign/codonalign.html

	optional arguments:
	  -h, --help            show this help message and exit
	  -f FILENAME, --file FILENAME
	                        input file in fasta format <sample.fasta>
	  -r FRAME, --frame FRAME
	                        frame [1, 2 or 3], otherwise best is determined
	  -s FS_COMP, --fs_comp FS_COMP
	                        max n of compensating frameshifts <5>