import os
import pandas

def fastasplit(inputfasta, outputpath, outputfastas, file_name, parameter):
    splitfile = None
    linesperfile = parameter['linesperfile']

    isExist = os.path.exists(outputpath)

    count = 0

    if not isExist:
        os.makedirs(outputpath)

    with open(inputfasta) as input:
        for lineno, line in enumerate(input):
            if line.startswith('>hCoV-19/'):
                print(line)

                if count % parameter['linesperfile'] == 0:
                    if splitfile:
                        splitfile.close()
                    split_filename = outputfastas.format(file_name, file_name, int(count / linesperfile + 1))
                    splitfile = open(split_filename, 'w')
                    print(f'creating {split_filename}')

                count += 1
            splitfile.write(line)

    if not splitfile.closed:
        splitfile.close()

    print(f'Finished creating {int(count / linesperfile)} files.')

def metasplit(inputmeta):
    meta_db = p

def getmetaname(inputname):
    name_split = inputname.split('.')
    meta_name = f"{name_split[0]}.metadata.tsv"
    file_name = name_split[0].split('/')
    return meta_name, file_name[2]


if __name__ == '__main__':
    inputfasta = "data/test2/1660403438836.sequences.fasta"
    inputmeta, file_name = getmetaname(inputfasta)
    outputpath = f"data/{file_name}"
    outputfastas = "data/{}/{}_{}.fasta"
    parameter = {'linesperfile': 3}

    fastasplit(inputfasta, outputpath, outputfastas, file_name, parameter)
