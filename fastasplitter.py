import os

def splitfasta(inputfasta, outputpath, outputfastas, parameter):
    splitfile = None

    isExist = os.path.exists(outputpath)

    count = 0

    if not isExist:
        os.makedirs(outputpath)

    with open(inputfasta) as input:
        for lineno, line in enumerate(input):
            if splitfile:
                splitfile.close()
            split_filename = outputfastas.format(int(lineno/parameter['linesperfile'] + 1))
            splitfile = open(split_filename, 'w')
            print(f'creating {split_filename}')
        splitfile.write(line)
        count += 1
    if splitfile:
        splitfile.close()
        count += 1

    print(f'Finished creating {count} files.')


if __name__ == '__main__':
    inputfasta = "data/spikenuc0523_1.fasta"
    outputpath = "data/spikenuc0523_1"
    outputfastas = "data/spikenuc0523_1/spikenuc0523_1_{}.fasta"
    parameter = {'linesperfile': 2500}

    splitfasta(inputfasta, outputpath, outputfastas, parameter)
