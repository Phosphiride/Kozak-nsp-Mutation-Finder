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
    parameter = {'min_length': 29500,
                 'max_length': 30000,
                 'max_ambig': 100,
                 'ref_seq': ref_seq,
                 'refprotname': wildtype,
                 'mafft': "C:/Program Files/mafft-win/mafft.bat",
                 'max_muts': 100000,
                 'site_offset': 10055,  # nsp5: 10055; nsp12: 13442
                 'exclude_ambig': True,
                 'align_size': 1000
                 }
    createmutfiles(inputfastas, outputfilenames, parameter)