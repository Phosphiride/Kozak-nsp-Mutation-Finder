import pandas as pd
from glob import glob

def mergedf(inputcsv, parameter):
    unagg_df = pd.DataFrame(columns=parameter['column_names'])
    total_seq = 0

    for file in inputcsv:
        wk_df = pd.read_csv(file)
        df_sum = wk_df.iloc[-1, wk_df.columns.get_loc('total seq.')]
        total_seq += df_sum
        wk_df = wk_df[wk_df['aa site'].isin(parameter['aasite'])]
        unagg_df = pd.concat([unagg_df, wk_df])

    return unagg_df, total_seq

def dfagg(unagg_df, total_seq, outputpath):
    final_df = pd.DataFrame()

    for amino in parameter['aasite']:
        temp_df = unagg_df[(unagg_df['aa site'] == amino)]
        final_df = pd.concat([final_df, temp_df])

    #final_df.to_csv('result/test1.csv')

    #agg_func = {'count': ['sum']}
    #grouplist = ['gene site', 'genome site', 'wt nt', 'mutant nt', 'aa site', 'wt aa', 'mutant aa']

    #final_df = final_df.groupby(grouplist).aggregate(count=pd.NamedAgg('country', 'count'), n_countries=pd.NamedAgg('country', 'nunique'))\
    #    .reset_index() .sort_values('count', ascending=False).assign(frequency=lambda x: x['count'] / total_seq)

    #group_by_col_names = ['gene site', 'genome site', 'wt nt', 'mutant nt', 'aa site', 'wt aa', 'mutant aa']
    '''agg_func = {
        'count': ['sum']
    }

    df4 = final_df.groupby(grouplist).agg(agg_func)\
        .assign(total = lambda x: total_seq).assign(frequency = lambda x: x['count'] / total_seq)\
        .sort_values('gene site', ascending=True)'''

    #df4.to_csv(outputpath)
    final_df.to_csv(outputpath, index=False)
    print(f'Finished writing output to {outputpath}.')


if __name__ == '__main__':
    #inputpath = 'result/output_1-95'
    inputpath = '/work/data/kozak-data/raw_1-95_11202022'
    inputcsv = glob(f'{inputpath}/*.csv')
    outputpath = f'{inputpath}/raw_1-95_11202022.csv'

    parameter = {  #'aasite': [15, 17, 21],
        #'aasite': [15, 17, 21, 41, 45, 49, 50, 54, 55, 66, 70, 75, 83, 88, 89, 90, 96, 105, 108, 129, 132, 135,
        #           140, 141, 142, 143, 144, 145, 163, 164, 165, 166, 167, 168, 172, 184, 186, 187, 188, 189,
        #           190, 191, 192, 205, 212, 213, 220, 221, 234, 246, 248, 252, 253, 255, 260, 266, 274, 303],
        'aasite': [75],
        'column_names': ['description', 'row', 'gene site', 'genome site', 'wt nt', 'mutant nt', 'aa site', 'wt aa',
                         'mutant aa', 'country', 'date', 'total seq.']
        #'except': 500
        }
    unagg_df, total_seq = mergedf(inputcsv, parameter)
    dfagg(unagg_df, total_seq, outputpath)