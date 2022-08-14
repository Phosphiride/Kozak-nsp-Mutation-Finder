import getopt
import sys


def extractMutation( input_file_name , output_file_name , start_line_no , number_of_lines):

    line_no = 0
    line_count = 0
    with open (input_file_name ) as input_file_fp:


        line = input_file_fp.readline()
        while line :
            line_no = line_no +1
            if line_no == start_line_no:
                output_file_fp = open(output_file_name , 'w')
                line_count  = line_count+1
                while line and line_count <= number_of_lines:
                    output_file_fp.write(line)
                    line = input_file_fp.readline()
                    line_count = line_count+1
                output_file_fp.close()
                input_file_fp.close()
                return
            else:
                continue

if __name__ == "__main__":

    options, remainder = getopt.getopt(sys.argv[1:], '', ['startlineno=',
                                                             'numberoflines=',
                                                          'inputfile=',
                                                          "outputfile="

                                                             ])
    for opt , arg in options:
        if opt in ['--inputfile']:
            input_file_name = arg
        elif opt in ['--outputfile']:
            output_file_name = arg
        elif opt in ['--startlineno']:
            start_line_no = int(arg)
        elif opt in ['--numberoflines']:
            number_of_lines = int(arg)
        else:
            print(' usage : python3 ExtractMutation.py --inputfile [inputfile_name] --outputfile [ output file name] --startlineno [start line number] --numberoflines [ number of lines to extract]')
            exit(1)

    extractMutation(input_file_name , output_file_name , start_line_no , number_of_lines)