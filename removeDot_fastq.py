import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
for title, seq, qual in FastqGeneralIterator(sys.stdin):
    print("@%s\n%s\n+\n%s" % (title, seq.replace(".", "N"), qual))