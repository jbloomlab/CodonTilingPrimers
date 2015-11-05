==========================
CodonTilingPrimers
==========================

The *create_primers.py* script can be used to create *NNN* primers that tile the codons of a gene in both the forward and reverse direction. You can use this primers to make codon mutagenesis libraries. You might want to order these primers in the form of an `IDT 96-well plate`_.

This script was written by `Jesse Bloom`_, and can be downloaded `on GitHub`_.

Edited by Adam Dingens Nov 2015 to generate primers of differing lengths to all have a similar melting temperature (Tm). Notes on edits below. 
This script first makes an ORIGINAL primer of specified length  according to the infile (usually 37 bps). 
If the ORIGINAL primer has a Tm of greater than MaxPrimerTm, then nucleotides are trimmed off one by one (first 5', then 3', then 5' etc) until the Tm is less than MaxPrimerTm. Note that this could be over a degree less than the MaxPrimerTm. 
If the ORIGINAL primer has a Tm less than MinPrimerTm, then nucelotides are added one by one (first 3', then 5', then 3' etc) until the Tm is over MinPrimerTm. Note that this could be over a degree more than the MinPrimerTm
If the ORIGINAL primer has a Tm of less than MaxPrimerTm but greater than MinPrimerTm, it is not altered. 
The primers are constrained to be between MinPrimerlength and MaxPrimerLength bps long. The Tm of some MaxPrimerLength primers may not be > MinPrimerTemp, and the Tm of some MinPrimerLength primers may not be < MaxPrimerTm.

An infile (see below), MinPrimerTm, MaxPrimerTm, MinPrimerLength, and MaxPrimerLength are given as command line arguments, in that order.  

The  Tm_NN command of the MeltingTemp Module of Biopython (http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html) is used to calculate Tm of primers. 
This calculation is based on nearest neighbor thermodynamics. nucelotides labeled N are given average values in the Tm calculation. 
It is possible to vary salt concentration and other addatives if needed.


Running the script
-------------------------
The script can be run with Python 2. It has been tested with version 2.6, but will probably work with other versions.

To run the script, create an input file with the format described below, and give the infile (see below), MinPrimerTm, MaxPrimerTm, MinPrimerLength, and MaxPrimerLength as arguments in that order. Make sure this input file and the files that it references into the same directory as the script, and then run the script. This script is distributed with an example input file named *EN72-HA_infile.txt*. To run that example, use the following command::

    python create_primers.py EN72-HA_infile.txt 60 61 25 51 

This will create the output file *EN72-HA_primers.txt* which lists all of the primers, each with Tms around 60-61C, with a minimum length of 25 bp and a maximum length of 51 bp. 
Note that there are primers with Tms less that MinPrimerTm and more than MaxPrimerTm due to either length limitations or instances when trimming a primer resulted in a large increase or decrease in Tm past the MinPrimerTm or MaxPrimerTm respectively. 
You can create a similar input file for your analysis.


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
    startprimerlength 37
    sequencefile EN72-HA.txt
    outfile EN72-HA_primers.txt
    primerprefix EN72
    firstcodon 2


Output of the script
---------------------

The result of running this script is the file specified by *outfile*. It lists the primers to be ordered. All of the forward primers are have names which are the prefix listed by *primerprefix*, then *-for-mut*, then the codon number starting with *firstcodon*. The reverse primers are named similarly, but with the *for* replaced by *rev*. The forward primers are grouped in sets of 96 (for ordering in 96-well plates), as are the reverse primers Here are the first few lines of the output of an example output file::

    
	Plate 1
	EN72-for-mut2, ggggataattctattaatcatgNNNACTATCATTGCTTTGAGCTACA
	EN72-for-mut3, gggataattctattaatcatgAAGNNNATCATTGCTTTGAGCTACATTTTC
	EN72-for-mut4, ataattctattaatcatgAAGACTNNNATTGCTTTGAGCTACATTTTCTGT
	EN72-for-mut5, ctattaatcatgAAGACTATCNNNGCTTTGAGCTACATTTTCTGT
	EN72-for-mut6, ttaatcatgAAGACTATCATTNNNTTGAGCTACATTTTCTGTCTGG
	EN72-for-mut7, atgAAGACTATCATTGCTNNNAGCTACATTTTCTGTCTGG
	EN72-for-mut8, gAAGACTATCATTGCTTTGNNNTACATTTTCTGTCTGGTTCT
	EN72-for-mut9, CTATCATTGCTTTGAGCNNNATTTTCTGTCTGGTTCTC
	EN72-for-mut10, ATTGCTTTGAGCTACNNNTTCTGTCTGGTTCTCG

Here are the last few lines of the same file::

    
	EN72-rev-mut556, CTAATGTTGCCTTTNNNGCAGGCCCACATG
	EN72-rev-mut557, CCTAATGTTGCCNNNTTGGCAGGCCC
	EN72-rev-mut558, TTGCACCTAATGTTNNNTTTTTGGCAGGCCC
	EN72-rev-mut559, AATGTTGCACCTAATNNNGCCTTTTTGGCAGG
	EN72-rev-mut560, CAAATGTTGCACCTNNNGTTGCCTTTTTGGC
	EN72-rev-mut561, caAATGCAAATGTTGCANNNAATGTTGCCTTTTTGG
	EN72-rev-mut562, cactcaAATGCAAATGTTNNNCCTAATGTTGCCTTTTTG
	EN72-rev-mut563, acactcaAATGCAAATNNNGCACCTAATGTTGCC
	EN72-rev-mut564, taatacactcaAATGCANNNGTTGCACCTAATGTTGC
	EN72-rev-mut565, ttaattactaatacactcaAATNNNAATGTTGCACCTAATGTTGCCT
	EN72-rev-mut566, tttttaattactaatacactcaNNNGCAAATGTTGCACCTAATGTTG




.. _`on GitHub`: https://github.com/jbloom/CodonTilingPrimers
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`IDT 96-well plate`: http://www.idtdna.com/pages/products/dna-rna/96-and-384-well-plates
