import os
import io
import re
import subprocess

import Bio.SeqIO
from createmutationfile import createMutFile

def createmutfiles(inputfastas, outputfilenames, unaggoutputs, countryaggs, parameter):
    for inputfasta, outputfilename, unaggoutputs, countryaggs in zip(inputfastas, outputfilenames, unaggoutputs, countryaggs):
        createMutFile(inputfasta, outputfilename, unaggoutputs, countryaggs, parameter)

if __name__ == '__main__':
    inputfastas = ["data/multsingle_1.fasta", "data/multsingle_2.fasta"]
    outputfilenames = ['result/multtest_1.csv', 'result/multtest_2.csv']
    unaggoutputs = ['result/unagg_1.csv', 'result/unagg_2.csv']
    countryaggs = ['result/country_1.csv', 'result/country_2.csv']
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
    createmutfiles(inputfastas, outputfilenames, unaggoutputs, countryaggs, parameter)