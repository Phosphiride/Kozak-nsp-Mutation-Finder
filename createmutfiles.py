import os
import io
import re
import subprocess
from glob import glob

import Bio.SeqIO
from createmutationfile import createMutFile

def createmutfiles(inputfastas, outputfilenames, unaggoutputs, countryaggs, parameter):
    for inputfasta, outputfilename, unaggoutput, countryagg in zip(inputfastas, outputfilenames, unaggoutputs, countryaggs):
        createMutFile(inputfasta, outputfilename, unaggoutput, countryagg, parameter)

def outputname(inputfastas):
    name_list = []
    count = 1
    for file in inputfastas:
        name_list.append(f'result/multtest_{count}.csv')
        count += 1
    return name_list

def unagg(inputfastas):
    name_list = []
    count = 1
    for file in inputfastas:
        name_list.append(f'result/unagg_{count}.csv')
        count += 1
    return name_list

def country(inputfastas):
    name_list = []
    count = 1
    for file in inputfastas:
        name_list.append(f'result/country_{count}.csv')
        count += 1
    return name_list

if __name__ == '__main__':
    #inputfastas = ["data/20211104_gisaid_genomes.fasta", "data/20210613_gisaid_genomes.fasta"]
    inputfastas = glob("./data/multfolder/*.fasta")
    #outputfilenames = ['result/multtest_1.csv', 'result/multtest_2.csv']
    outputfilenames = outputname(inputfastas)
    #unaggoutputs = ['result/unagg_1.csv', 'result/unagg_2.csv']
    unaggoutputs = unagg(inputfastas)
    #countryaggs = ['result/country_1.csv', 'result/country_2.csv']
    countryaggs = country(inputfastas)
    wildtype = "data/GISAID_nsp5.fasta"
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