import os
import io
import re
import subprocess

import Bio.SeqIO
from Bio.Data.IUPACData import protein_letters
import pandas as pd
from createmutationfile import createMutFile

def createmutfiles(inputfastas, outputfilenames, parameter):
    for inputfasta, outputfilename in zip(inputfastas, outputfilenames):
        createMutFile(inputfasta, outputfilename, parameter)

if __name__ == '__main__':
    inputfastas = ["data/spikenuc0523_1.fasta", "data/spikenuc0523_2.fasta"]
    outputfilenames = ['result/multtest_1.csv', 'result/multtest_2.csv']
    wildtype = "data/wildtype_sequence.fasta"
    ref_seq = Bio.SeqIO.read(wildtype, 'fasta')
    parameter = {'min_length': 3800,
                 'max_length': 3900,
                 'max_ambig': 100,
                 'ref_seq': ref_seq,
                 'refprotname': wildtype,
                 'mafft': "C:\Program Files\mafft-win\mafft.bat",
                 'max_muts': 150,
                 'site_offset': 330
                 }
    createmutfiles(inputfastas, outputfilenames, parameter)