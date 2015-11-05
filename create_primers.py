"""Script for creating gene assembly primers.
    
Primers containing an NNN codon at their center are tiled along a gene.

Jesse Bloom, 2013. 

Edited by Adam Dingens Nov 2015 to generate primers of differing lengths to all have a Tm of ~60C. Notes on edits below. 
This script first makes an ORIGINAL primer of specified length  according to the infile (usually 37 bps). 
If the ORIGINAL primer has a Tm of greater than 61C, then nucleotides are trimmed off one by one (first 5', then 3', then 5' etc) until the Tm is less than 61C (which could result in less than 60C primer).
If the ORIGINAL primer has a Tm less than 60C, then nucelotides are added one by one (first 3', then 5', then 3' etc) until the Tm is over 60C.
If the ORIGINAL primer has a Tm of less than 61C but greater than 60C, it is not altered. 
The primers are constrained to be between 25 and 51 bps long. Some 51 bp primers may not be > 60C, and some 25 bp primers may not be < 61C.


The  Tm_NN command of the MeltingTemp Module of Biopython (http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html) is used to calculate Tm of primers. 
This calculation is based on nearest neighbor thermodynamics. nucelotides labeled N are given average values in the Tm calculation. 
It is possible to vary salt concentration and other addatives if needed."""


import os
import sys
import math
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import re



def ReverseComplement(seq):
    """Returns reverse complement of sequence. Preserves upper/lower case."""
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', 'n':'n'}
    rc = [d[nt] for nt in seq]
    rc.reverse()
    return ''.join(rc)


def CreateMutForOligos(seq, primerlength, prefix, firstcodon):

    """Creates oligos to tile a gene and introduce NNN at each codon.

    *seq* : sequence of the gene. The gene itself should be upper case. The
    flanking regions and start / stop codons should be lower case. 
    All upper case codons are randomized. The length of the lower
    case sequences at each end must be >= (primerlength - 3) / 2.0

    *primerlength* : length of primers. Must be an odd number, so that equal length
    flanking on each side.

    *prefix* : string prefix attached to primer names.

    *firstcodon* : number assigned to first codon in primer name.

    Tiles primers across the gene in the forward direction. The primers
    are all of length primerlength with NNN at the middle codon.
    Note that only upper case letters are randomized.
    Primers are named as follows:

    "%s-for-mut%d" % (prefix, i) -> 5' tiling primers, where i = 2, 3, ...
    In other words, the indices cover codons 2 and up.

    Returns a list of all these primers as *(name, sequence)* 2-tuples.
    """
    n = len(seq)
    assert primerlength % 2 == 1, "primer length not odd"
    flanklength = (primerlength - 3) // 2
    upperseq = ''.join([nt for nt in seq if nt.istitle()])
    assert upperseq in seq, "upper case nucleotides not substring"
    assert len(upperseq) % 3 == 0, "length of upper case not multiple of 3"
    startupper = seq.index(upperseq)
    if startupper < flanklength:
        raise ValueError("not enough 5' lower case flanking nucleotides")
    if n - len(upperseq) - startupper < flanklength:
        raise ValueError("not enough 3' lower case flanking nucleotides")
    ncodons = len(upperseq) // 3
    primers = []
    for icodon in range(ncodons):
        i = startupper + icodon * 3
        primer = "%sNNN%s" % (seq[i - flanklength : i], seq[i + 3 : i + 3 + flanklength])
        name = "%s-for-mut%d" % (prefix, firstcodon + icodon)
        primers.append((name, primer))
    return primers


def CreateMutForOligosVarLength(seq, primerlength, prefix, firstcodon):

    """Creates oligos to tile a gene and introduce NNN at each codon.

    *seq* : sequence of the gene. The gene itself should be upper case. The
    flanking regions and start / stop codons should be lower case. 
    All upper case codons are randomized. The length of the lower
    case sequences at each end must be >= (primerlength - 3) / 2.0

    *primerlength* : length of primers. Must be an odd number, so that equal length
    flanking on each side.

    *prefix* : string prefix attached to primer names.

    *firstcodon* : number assigned to first codon in primer name.

    Tiles primers across the gene in the forward direction. The primers
    are all of length primerlength with NNN at the middle codon.
    Note that only upper case letters are randomized.
    Primers are named as follows:

    "%s-for-mut%d" % (prefix, i) -> 5' tiling primers, where i = 2, 3, ...
    In other words, the indices cover codons 2 and up.

    Returns a list of all these primers as *(name, sequence)* 2-tuples.
    """
    uppertemplimit = float(sys.argv[3])
    lowertemplimit = float(sys.argv[2])
    maxlength = float(sys.argv[5])
    minlength = float(sys.argv[4])
    n = len(seq)
    assert primerlength % 2 == 1, "primer length not odd"
    initial_flanklength = (primerlength - 3) // 2
    upperseq = ''.join([nt for nt in seq if nt.istitle()])
    assert upperseq in seq, "upper case nucleotides not substring"
    assert len(upperseq) % 3 == 0, "length of upper case not multiple of 3"
    startupper = seq.index(upperseq)
    if startupper < initial_flanklength:
        raise ValueError("not enough 5' lower case flanking nucleotides")
    if n - len(upperseq) - startupper < initial_flanklength:
        raise ValueError("not enough 3' lower case flanking nucleotides")
    ncodons = len(upperseq) // 3
    primers = []
    for icodon in range(ncodons):
        i = startupper + icodon * 3
        primer = "%sNNN%s" % (seq[i - initial_flanklength : i], seq[i + 3 : i + 3 + initial_flanklength]) 
        name = "%s-for-mut%d" % (prefix, firstcodon + icodon)
        primerseq = Seq(primer)
        
        TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
        print name
        print primer
        print "Original primer has Tm of:"
        print TmNN
        add_3 = True
        minus_5 = True
        flank5 = flank3 = initial_flanklength
        if float(TmNN) > float(uppertemplimit):
            while float(TmNN) > float(uppertemplimit) and len(primer) > minlength:
                if minus_5:
                    flank5 -= 1
                    primer = "%sNNN%s" % (seq[i - (flank5) : i], seq[i + 3 : i + 3 + flank3])
                    minus_5 = False
                else:
                    flank3 -= 1
                    primer =  "%sNNN%s" % (seq[i - (flank5) : i], seq[i + 3 : i + 3 + flank3])

                    minus_5 = True
                if startupper < flank5:
                    raise ValueError("not enough 5' lower case flanking nucleotides")
                if n - len(upperseq) - startupper < flank3:
                    raise ValueError("not enough 3' lower case flanking nucleotides") #not sure if this is correct!

                
                primerseq = Seq(primer)
                TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
                primerlength= len(primer)
            print "Redesigned %s has a length of %s and a Tm of %s and the sequence is:" % (name, primerlength, TmNN)
            print primer
            print "\n"
        else: 
            if float(TmNN) < float(lowertemplimit):
                while float(TmNN) < float(lowertemplimit) and len(primer) < maxlength:
                    if add_3:
                        flank3 += 1
                        primer = "%sNNN%s" % (seq[i - (flank5) : i], seq[i + 3 : i + 3 + flank3])
                        add_3 = False
                    else:
                        flank5 +=1
                        primer = "%sNNN%s" % (seq[i - (flank5) : i], seq[i + 3 : i + 3 + flank3])
                        add_3 = True    
                    primerseq = Seq(primer)
                    TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
                    primerlength= len(primer)
                    if startupper < flank5:
                        raise ValueError("not enough 5' lower case flanking nucleotides")
                    if n - len(upperseq) - startupper < flank3:
                        raise ValueError("not enough 3' lower case flanking nucleotides") 
                print "Redesigned %s has a length of %s and a Tm of %s and the sequence is:" % (name, primerlength, TmNN)
                print primer
                print "\n"

            else:
                print "Original %s bp primer %s has a Tm of %s. This is between 60 and 61C and doesn't require any trimming. The sequence is:" % (primerlength, name, TmNN)
                print primer
                print "\n"








 




        primers.append((name, primer))
    #print primers 
    return primers



def main():
    args = sys.argv[1 : ]
    uppertemplimit = float(sys.argv[3])
    lowertemplimit = float(sys.argv[2])
    maxlength = float(sys.argv[5])
    minlength = float(sys.argv[4])
    if len(args) != 5:
        raise IOError("Script must be run with five arguments specifying the file name, minimum primer Tm, maximum primer Tm, minimum primer length, and maximum prmer length.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Cannot find infile of %s" % infilename)
    inlines = [line for line in open(infilename).readlines() if not line.isspace() and line[0] != '#']
    d = {}
    for line in inlines:
        entries = line.split(None, 1)
        if len(entries) != 2:
            raise ValueError("Line does not specify a key and a value:\n%s" % line)
        (key, value) = (entries[0].strip(), entries[1].strip())
        if key in d:
            raise ValueError("Duplicate key of %s" % key)
        d[key] = value
    primerlength = int(d['startprimerlength'])

    if (primerlength <=3 ) or (primerlength % 2 == 0):
        raise ValueError("Does not appear to be valid primer length: %d" % primerlength)
    print "Designing primers first of  length %d, then trimming or adding to have %s < Tm <%s but primer length between %s and %s bps" % (primerlength, lowertemplimit, uppertemplimit, minlength, maxlength)
    sequencefile = d['sequencefile']
    if not os.path.isfile(sequencefile):
        raise IOError("Cannot find sequencefile %s" % sequencefile)
    sequence = open(sequencefile).read()
    sequence = sequence.replace(' ', '')
    sequence = sequence.replace('\n', '')
    print "Read a sequence of length %d from %s:\n%s" % (len(sequence), sequencefile, sequence)
    outfile = d['outfile']
    primerprefix = d['primerprefix']
    firstcodon = int(d['firstcodon'])
    print "The primers will be named with the prefix %s, and the first codon numbered as %d." % (primerprefix, firstcodon)

    # Design forward mutation primers
    mutforprimers = CreateMutForOligosVarLength(sequence, primerlength, primerprefix, firstcodon)
    print "Designed %d mutation forward primers." % len(mutforprimers)


    
    # Design reverse mutation primers
    mutrevprimers = [(name.replace('for', 'rev'), ReverseComplement(seq)) for (name, seq) in mutforprimers]
    print "Designed %d mutation reverse primers." % len(mutrevprimers)
   
    # Print out all of the primers
    primers = mutforprimers + mutrevprimers
    print "This gives a total of %d primers." % len(primers)


    print "\nNow writing these primers to %s" % outfile
    iplate = 1
    f = open(outfile, 'w')
    for primers in [mutforprimers, mutrevprimers]:
        f.write("\r\nPlate %d\r\n" % iplate)
        n_in_plate = 0
        for (name, primer) in primers:
            f.write("%s, %s\r\n" % (name, primer))
            n_in_plate += 1
            if n_in_plate == 96:
                n_in_plate = 0
                iplate += 1
                f.write("\r\nPlate %d\r\n" % iplate)
        if n_in_plate:
            iplate += 1



main() # run the main program
