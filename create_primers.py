"""Script for creating gene assembly primers.
    
Primers containing an NNN codon at their center are tiled along a gene.

Jesse Bloom, 2013."""


import os
import sys
import math



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


def main():
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be run with one argument specifying the file name.")
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
    primerlength = int(d['primerlength'])
    if (primerlength <=3 ) or (primerlength % 2 == 0):
        raise ValueError("Does not appear to be valid primer length: %d" % primerlength)
    print "Designing primers of length %d" % primerlength
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
    mutforprimers = CreateMutForOligos(sequence, primerlength, primerprefix, firstcodon)
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
