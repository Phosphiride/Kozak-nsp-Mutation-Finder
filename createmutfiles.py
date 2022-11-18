import os
import io
import re
import subprocess
from glob import glob
from datetime import datetime

import Bio.SeqIO
from createmutationfile import createMutFile
#from test import createMutFile

def createmutfiles(inputfastas, outputfilenames, unaggoutputs, countryaggs, dateaggs, parameter):
    for inputfasta, outputfilename, unaggoutput, countryagg, dateagg in zip(inputfastas, outputfilenames, unaggoutputs, countryaggs, dateaggs):
        createMutFile(inputfasta, outputfilename, unaggoutput, countryagg, dateagg, parameter)

def createfolders():
    filedate = datetime.today().strftime('%Y%m%d%H%M%S')
    #outputpath = f'result/{filedate}/output'                   #windows
    #countrypath = f'result/{filedate}/country'
    #unaggpath = f'result/{filedate}/unagg'
    #datepath = f'result/{filedate}/date'
    outputpath = f'/work/data/kozak_data/{filedate}/output'                    #server
    countrypath = f'/work/data/kozak_data/{filedate}/country'
    unaggpath = f'/work/data/kozak_data/{filedate}/unagg'
    datepath = f'/work/data/kozak_data/{filedate}/date'
    for path in [outputpath, countrypath, unaggpath, datepath]:
        os.makedirs(path)
    return outputpath, countrypath, unaggpath, datepath

def outputname(inputfastas, outputpath, countrypath, unaggpath, datepath):
    output_list = []
    country_list = []
    unagg_list = []
    date_list = []
    count = 1
    for file in inputfastas:
        output_list.append(f'{outputpath}/output_{count}.csv')
        unagg_list.append(f'{countrypath}/country_{count}.csv')
        country_list.append(f'{unaggpath}/unagg_{count}.csv')
        date_list.append(f'{datepath}/date_{count}.csv')
        count += 1
    return output_list, country_list, unagg_list, date_list

if __name__ == '__main__':
    #inputfastas = ["data/20211104_gisaid_genomes.fasta", "data/20210613_gisaid_genomes.fasta"]
    #inputfastas = glob("./data/work_1/*.fasta")
    inputfastas = glob("/work/data/kozak_data/allnuc0521_split/*.fasta")
    outputpath, countrypath, unaggpath, datepath = createfolders()

    #outputfilenames = ['result/multtest_1.csv', 'result/multtest_2.csv']
    outputfilenames, countryaggs, unaggoutputs, dateaggs = outputname(inputfastas, outputpath, countrypath, unaggpath, datepath)
    #unaggoutputs = ['result/unagg_1.csv', 'result/unagg_2.csv']
    #countryaggs = ['result/country_1.csv', 'result/country_2.csv']
    wildtype = "data/GISAID_nsp5.fasta"
    ref_seq = Bio.SeqIO.read(wildtype, 'fasta')
    parameter = {'min_length': 900,
                 'max_length': 1100,
                 'max_ambig': 100000,
                 'ref_seq': ref_seq,
                 'refprotname': wildtype,
                 #'mafft': "C:/Program Files/mafft-win/mafft.bat",              #home pc
                 'mafft': "/home/tools/mafft-7.505-with-extensions/core/mafft",  # server
                 'max_muts': 10000,
                 'site_offset': 10055,  # nsp5: 10055; nsp12: 13442
                 'exclude_ambig': True,
                 'align_size': 100000
                 }
    createmutfiles(inputfastas, outputfilenames, unaggoutputs, countryaggs, dateaggs, parameter)