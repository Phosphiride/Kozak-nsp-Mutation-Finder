import os
import re

def dnatest(test):
    out = re.search("^[ATCGMRWSYKVHDBN]+$", test)

    if out:
        return True

    else:
        return False

def fastasplit(inputfasta, outputpath, outputfastas, file_name, parameter):
    splitfile = None
    linesperfile = parameter['linesperfile']

    isExist = os.path.exists(outputpath)

    count = 0

    if not isExist:
        os.makedirs(outputpath)

    with open(inputfasta) as inputfasta_fp:
        line = inputfasta_fp.readline()
        while line:
            if dnatest(line) == True and header_type == 1:
                splitfile.write(line)

            elif line.startswith('>NSP5|hCoV-19'):
                header_type = 1
                if count % parameter['linesperfile'] == 0:
                    if splitfile:
                        splitfile.close()
                    split_filename = outputfastas.format(file_name, file_name, int(count / linesperfile + 1))
                    splitfile = open(split_filename, 'w')
                    print(f'creating {split_filename}')

                count += 1
                splitfile.write(line)

            else:
                header_type = 0

            line = inputfasta_fp.readline()

    if not splitfile.closed:
        splitfile.close()


        '''for lineno, line in enumerate(input):
            if line.startswith('>NSP5|hCoV-19'):

                if count % parameter['linesperfile'] == 0:
                    if splitfile:
                        splitfile.close()
                    split_filename = outputfastas.format(file_name, file_name, int(count / linesperfile + 1))
                    splitfile = open(split_filename, 'w')
                    print(f'creating {split_filename}')

                count += 1
                splitfile.write(line)
                
                if dnatest(line) == True:
                    splitfile.write(line)



    if not splitfile.closed:
        splitfile.close()

    print(f'Finished creating {int(count / linesperfile)} files.')'''

def getfilename(inputname):
    name_split = inputname.split('.')
    file_name = name_split[0].split('/')
    return file_name[1]


if __name__ == '__main__':
    inputfasta = "data/spike_nsp5.fasta"
    file_name = getfilename(inputfasta)
    outputpath = f"data/{file_name}"
    outputfastas = "data/{}/{}_{}.fasta"
    parameter = {'linesperfile': 4}

    fastasplit(inputfasta, outputpath, outputfastas, file_name, parameter)
