from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO
import sys
fasta1=sys.argv[1]
fasta2=sys.argv[2]
try:
    name1=sys.argv[3]
except:
    name1=None

try:
    name2=sys.argv[4]
except:
    name2=None

handle = open(fasta1, "rU")
fasta1_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

handle = open(fasta2, "rU")
fasta2_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()


matrix = matlist.blosum62



if name1 is None and name2 is None:
    for name1 in fasta1_dict:
        seq_fasta1=fasta1_dict[name1]
        for name2 in fasta2_dict:
            seq_fasta2=fasta2_dict[name2]
            for a in pairwise2.align.localms(seq_fasta1, seq_fasta2, 2, -1, -.5, -.1):
                print name1,name2
                print format_alignment(*a)

elif name1 is not None and name2 is not None:
    seq_fasta1=fasta1_dict[name1]
    seq_fasta2=fasta2_dict[name2]
    for a in pairwise2.align.localms(seq_fasta1, seq_fasta2, 2, -1, -.5, -.1):
        print name1,name2
        print format_alignment(*a)
        
else:
    "Prrr"
