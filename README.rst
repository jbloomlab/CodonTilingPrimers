==========================
CodonTilingPrimers
==========================

The *create_primers.py* script can be used to create *NNN* primers that tile the codons of a gene in both the forward and reverse direction. You can use this primers to make codon mutagenesis libraries. You might want to order these primers in the form of an `IDT 96-well plate`_.

This script was written by `Jesse Bloom`_, and can be downloaded `on GitHub`_.


Running the script
-------------------------
The script can be run with Python 2. It has been tested with version 2.6, but will probably work with other versions.

To run the script, create an input file with the format described below. Make sure this input file and the files that it references into the same directory as the script, and then run the script. This script is distributed with an example input file named *EN72-HA_infile.txt*. To run that example, use the following command::

    python create_primers EN72-HA_infile.txt

This will create the output file *EN72-HA_primers.txt* which lists all of the primers. You can create a similar input file for your analysis.


Format of the input file
--------------------------

The input file is in text format. Any blank lines or lines beginning with *#* are ignored (the latter allows for comment lines). All other lines should in the format *key value* where *key* is one of the strings listed below and *value* is the value assigned to that key. The required *key* values are:

* *primerlength* the length of the primers to create. An *NNN* will be placed in the middle of the primer. The primer length must be odd so that there will be the same number of flanking nucleotides on each side of the *NNN*.

* *sequencefile* is the name of a file giving the sequence for which we are designing the primers. This file should only contain the sequence, and should not have any headers or other content. For the sequence, make the 5' and 3' ends that you do not want to mutate in lower case. Make the portion of the coding sequence that you want to tile with *NNN* in upper case. Typically, for example, you would not want to mutate the initial start codon, so the first *atg* would be lower case. You must have at least *(primerlength - 3) / 2* nucleotides in lower case at each end of the upper case sequence that you are mutating. This is because at least this much flanking sequence is needed to design primers of the indicated length.

* *outfile* is the name of the output file that we create which lists all of the primers. This file is overwritten if it already exists.

* *primerprefix* is the prefix attached to the primer names.

* *firstcodon* is the number of the first codon being mutated. This is used for naming the primers.


Example input file
-------------------

Here is an example input file::

    # input file to the create_primers.py script of CodonTilingPrimers
    primerlength 37
    sequencefile EN72-HA.txt
    outfile EN72-HA_primers.txt
    primerprefix EN72
    firstcodon 2


Output of the script
---------------------

The result of running this script is the file specified by *outfile*. It lists the primers to be ordered. All of the forward primers are have names which are the prefix listed by *primerprefix*, then *-for-mut*, then the codon number starting with *firstcodon*. The reverse primers are named similarly, but with the *for* replaced by *rev*. The forward primers are grouped in sets of 96 (for ordering in 96-well plates), as are the reverse primers Here are the first few lines of the output of an example output file::

    Plate 1
    EN72-for-mut2, taattctattaatcatgNNNACTATCATTGCTTTGAG
    EN72-for-mut3, ttctattaatcatgAAGNNNATCATTGCTTTGAGCTA
    EN72-for-mut4, tattaatcatgAAGACTNNNATTGCTTTGAGCTACAT
    EN72-for-mut5, taatcatgAAGACTATCNNNGCTTTGAGCTACATTTT
    EN72-for-mut6, tcatgAAGACTATCATTNNNTTGAGCTACATTTTCTG
    EN72-for-mut7, tgAAGACTATCATTGCTNNNAGCTACATTTTCTGTCT
    EN72-for-mut8, AGACTATCATTGCTTTGNNNTACATTTTCTGTCTGGT
    EN72-for-mut9, CTATCATTGCTTTGAGCNNNATTTTCTGTCTGGTTCT
    EN72-for-mut10, TCATTGCTTTGAGCTACNNNTTCTGTCTGGTTCTCGG

Here are the last few lines of the same file::

    EN72-rev-mut556, CACCTAATGTTGCCTTTNNNGCAGGCCCACATGATGA
    EN72-rev-mut557, TTGCACCTAATGTTGCCNNNTTGGCAGGCCCACATGA
    EN72-rev-mut558, ATGTTGCACCTAATGTTNNNTTTTTGGCAGGCCCACA
    EN72-rev-mut559, CAAATGTTGCACCTAATNNNGCCTTTTTGGCAGGCCC
    EN72-rev-mut560, ATGCAAATGTTGCACCTNNNGTTGCCTTTTTGGCAGG
    EN72-rev-mut561, caAATGCAAATGTTGCANNNAATGTTGCCTTTTTGGC
    EN72-rev-mut562, actcaAATGCAAATGTTNNNCCTAATGTTGCCTTTTT
    EN72-rev-mut563, tacactcaAATGCAAATNNNGCACCTAATGTTGCCTT
    EN72-rev-mut564, taatacactcaAATGCANNNGTTGCACCTAATGTTGC
    EN72-rev-mut565, tactaatacactcaAATNNNAATGTTGCACCTAATGT
    EN72-rev-mut566, aattactaatacactcaNNNGCAAATGTTGCACCTAA



.. _`on GitHub`: https://github.com/jbloom/CodonTilingPrimers
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`IDT 96-well plate`: http://www.idtdna.com/pages/products/dna-rna/96-and-384-well-plates
