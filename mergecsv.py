import os
import pandas as pd
from glob import glob

def mergecsv(inputfolder, outputpath, parameter):
    if parameter['datatype'] == 'output':
        mergeframe = pd.DataFrame(columns=['gene site', 'genome site', 'wt nt', 'mutant nt', 'aa site', 'wt aa',
                                                  'mutant aa', 'total seq.', 'count', 'n_countries', 'frequency'])
        createoutput(inputfolder, outputpath, mergeframe, parameter)
    if parameter['datatype'] == 'country':
        mergeframe = pd.DataFrame(columns=['gene site', 'genome site', 'wt nt', 'mutant nt', 'aa site', 'wt aa',
                                                  'mutant aa', 'country', 'count', 'n_countries', 'frequency'])
        createoutput(inputfolder, outputpath, mergeframe, parameter)
    if parameter['datatype'] == 'date':
        mergeframe = pd.DataFrame(columns=['gene site', 'genome site', 'wt nt', 'mutant nt', 'aa site', 'wt aa',
                                                  'mutant aa', 'total seq.', 'date', 'count', 'n_countries', 'frequency'])
        createoutput(inputfolder, outputpath, mergeframe, parameter)
    if parameter['datatype'] == 'unagg':
        mergeframe = pd.DataFrame(columns=['gene site', 'genome site', 'wt nt', 'mutant nt', 'aa site', 'wt aa',
                                                  'mutant aa', 'total seq.', 'country', 'date'])
        createoutput(inputfolder, outputpath, mergeframe, parameter)

def createoutput(inputfolder, outputpath, mergeframe, parameter):
    for file in inputfolder:
        wk_df = pd.read_csv(file)
        filtered_wk = wk_df[wk_df['aa site'].isin(parameter['aasite'])]
        exception_wk = wk_df[(wk_df['count'] >= parameter['except']) & ~wk_df['aa site'].isin(parameter['aasite'])]
        mergeframe = pd.concat([mergeframe, filtered_wk, exception_wk])
    mergeframe.to_csv(outputpath, index=False)
    print(f'Finished writing outputs to {outputpath}')



if __name__ == '__main__':
    #inputpath = 'result/20220819001526/country'
    #inputpath = '/work/data/kozak_data/result_1-30/output'
    inputpath = '/work/data/kozak_data/result_1-30/country'
    #inputpath = '/work/data/kozak_data/result_1-30/date'
    #inputpath = '/work/data/kozak_data/result_1-30/unagg'

    inputfolder = glob(f'{inputpath}/*.csv')

    #inputfolder = glob('result/20220819001526/country/*.csv')
    #inputfolder = glob('result/20220819001526/date/*.csv')
    #inputfolder = glob('result/20220819001526/unagg/*.csv')

    outputpath = f'{inputpath}/country_1-30.csv'
    parameter = {#'aasite': [15, 17, 21],
                 'aasite': [15, 17, 21, 41, 45, 49, 50, 54, 55, 66, 70, 75, 83, 88, 89, 90, 96, 105, 108, 129, 132, 135,
                            140, 141, 142, 143, 144, 145, 163, 164, 165, 166, 167, 168, 172, 184, 186, 187, 188, 189,
                            190, 191, 192, 205, 212, 213, 220, 221, 234, 246, 248, 252, 253, 255, 260, 266, 274, 303],
                 'except': 500,
                 #'datatype': 'output'
                 'datatype': 'country'
                 #'datatype': 'date'
                 #'datatype': 'unagg'
                 }
    mergecsv(inputfolder, outputpath, parameter)

'''To do:
    1. Make empty pandas dataframe for merging files
    2. For loop iterating through all files in inputfolder
    3. Make pd dataframe housing csv
    4. For loop iterating through all rows in csv dataframe > pd.query instead???
    5. If aasite in dataframe = one of the aasite values in parameter list, add row to merge files dataframe
    6. If aasite not found in parameter list but count exceeds 500, also add to merge files dataframe
    7. Export merge files dataframe as another csv, then process by hand I guess.'''
