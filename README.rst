==========================
CodonTilingPrimers
==========================

The *create_primers.py* script can be used to create *NNN* primers that tile the codons of a gene in both the forward and reverse direction. You can use this primers to make codon mutagenesis libraries. You might want to order these primers in the form of an `IDT 96-well plate`_.

This script was written by Adam Dingens and `Jesse Bloom`_, and can be downloaded `on GitHub`_.

The script takes command line arguments; for a listing of how to provide the arguments, type the following to get the help message::

    python create_primers.py -h

For instance, the file ``EN72_HA_primers.txt`` included in the repository can be generated with::

    python create_primers.py EN72-HA.txt EN72 2 En72-HA_primers.txt

There are a variety of optional parameters specifing primer length and melting temperature constraints; the default values for these optional parameters are displayed when you run the program with the ``-h`` option to get the help message.

You can also adjust the optional parameters described in the help message, such as::
	
    python create_primers.py EN72-HA.txt EN72 2 EN72-HA_primers.txt --minprimertm 65 --maxprimertm 66

The script works as follows:

    1) For each codon, it first makes an ORIGINAL primer of the length specified by ``--startprimerlength``

    2) If the original primer has a melting temperature (Tm) greater than the value specified by ``--maxprimertm``, then nucleotides are trimmed off one by one (first from the 5' end, then the 3' end, then the 5' end again, etc) until the melting temperature is less than ``--maxprimertm`` or the length is reduced to ``--minlength``.

    3) If the original primer has a Tm greater than ``--minprimertm``, then nucleotides are added one-by-one (first to the 3' end, then the 5' end, then the 3' end again, etc) until the melting temperature is greater than ``--minprimertm`` or the length reaches ``--maxlength``.

    4) Note that because the primers are constrained to be between ``--minprimerlength`` and ``--maxprimerlength``, the Tm may not always fall between ``--minprimertm`` and ``--maxprimertm``. This can also happen if a primer initially exceeds ``--maxprimertm`` but the first trimming that drops it below this value also drops it below ``--minprimertm``, or vice-versa if the primer is being extended to increase its melting temperature.

The  *Tm_NN* command of the *MeltingTemp* Module of *Biopython* (http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html) is used to calculate Tm of primers. 
This calculation is based on nearest neighbor thermodynamics; nucleotides labeled ``N`` are given average values in the Tm calculation. 

The result of running this script is the file specified by ``outfile``. It lists the primers. All of the forward primers are have names which are the prefix specified by ``primerprefix``, then ``-for-mut``, then the codon number starting with ``firstcodon``. The reverse primers are named similarly, but with the ``for`` replaced by ``rev``. The forward primers are grouped in sets of 96 (for ordering in 96-well plates), as are the reverse primers Here are the first few lines of the output of an example output file::

    
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
