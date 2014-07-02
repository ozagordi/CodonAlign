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

******************************************************************************

### Licensing
Code is released under the GNU GPL version 2.0. For the script `codonalign`,
obtained from LANL, please consider the following statement as well.

	The script has been approved by LANL/LANS for open source release.
	Include the following copyright and disclaimer. [JM 2012]

	Copyright 2012.  Los Alamos National Security, LLC. This material was produced under
	U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL),
	which is operated by Los Alamos National Security, LLC for the U.S. Department of
	Energy. The U.S. Government has rights to use, reproduce, and distribute this software.
	NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS
	OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
	modified to produce derivative works, such modified software should be clearly marked,
	so as not to confuse it with the version available from LANL.

	Additionally, this program is free software; you can redistribute it and/or modify it
	under the terms of the GNU General Public License as published by the Free Software
	Foundation; version 2.0 of the License. Accordingly, this program is distributed in
	the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
	warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
	Public License for more details.
