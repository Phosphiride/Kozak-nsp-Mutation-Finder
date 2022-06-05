import os
import io
import re
import subprocess

import Bio.SeqIO
from Bio.Data.IUPACData import protein_letters
import pandas as pd

def createMutFile(inputfasta, outputfilename, parameter):
    sequences = list(Bio.SeqIO.parse(inputfasta, 'fasta'))      #parse sequences into python
    print(f'Read {len(sequences)} sequences.')

    seq_df = (
        pd.DataFrame({'seqrecord': sequences})
        .assign(description = lambda x: x['seqrecord'].map(lambda rec: rec.description),
                country = lambda x: x['description'].str.split('|').str[-1],
                host = lambda x: x['description'].str.split('|').str[6].str.strip(),
                length = lambda x: x['seqrecord'].map(len),
                n_ambiguous=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('X') + rec.seq.count('x')))
    )

    max_length = parameter['max_length']
    min_length = parameter['min_length']

    print(f'Only keeping sequences with lengths between {min_length} and {max_length}.')

    seq_df = (
        seq_df.assign(valid_length = lambda x: (min_length <= x['length']) & (x['length'] <= max_length))
    )

    seq_df = seq_df.query('valid_length')

    max_ambig = parameter["max_ambig"]

    print(f'Filtering sequences with greater than {max_ambig} ambiguous bases.')

    seq_df = (
        seq_df.assign(excess_ambiguous = lambda x: x['n_ambiguous'] > max_ambig)
    )

    temp_file = "temp/human_full-length-seq.fasta"

    Bio.SeqIO.write(seq_df['seqrecord'].tolist(), temp_file, 'fasta')
    ref_seq = parameter['ref_seq']
    rbd_df = alignment(50000, seq_df, parameter["refprotname"], ref_seq, parameter['mafft'])

    rbd_df = rbd_df.query('n_ambiguous == 0').query('n_gaps == 0')
    assert rbd_df['all_valid_nts'].all()
    print(f'Retained {len(rbd_df)} sequences')

    refseq_str = str(ref_seq.seq)
    ref_df = max_muts(rbd_df, refseq_str, parameter['max_muts'])

    site_offset = parameter['site_offset']

    write_output(ref_df, outputfilename, site_offset, refseq_str)


def alignment(chunksize, spikes_df, refprotfile, refseq, mafft):
    aligned_rbds = []

    for i in range(0, len(spikes_df), chunksize):
        spikes_file = os.path.join("temp/",
                                   f"human_full-length_spikes_{i + 1}-to-{i + chunksize}.fasta")
        print(f"Writing spikes {i + 1} to {i + chunksize} to {spikes_file}")
        _ = Bio.SeqIO.write(spikes_df['seqrecord'].tolist()[i: i + chunksize], spikes_file, 'fasta')
        print('Now aligning these sequences...')
        # cmds = ['mafft', '--auto', '--thread', str(config['max_cpus']),
        #        '--keeplength', '--addfragments', spikes_file, refprotfile]

        cmds = [mafft, '--auto', '--thread', '-1',
                '--keeplength', '--addfragments', spikes_file, refprotfile]

        res = subprocess.run(cmds, capture_output=True)
        if res.returncode:
            raise RuntimeError(f"Error in alignment:\n{res.stderr}")
        else:
            print('Alignment complete.\n')
            with io.StringIO(res.stdout.decode('utf-8')) as f:
                iseqs = list(Bio.SeqIO.parse(f, 'fasta'))

                # What mafft created are all lowercase sequence
                for iseq in iseqs:
                    iseq.seq = iseq.seq.upper()

                # remove reference sequence, which should be first in file
                print(iseqs[0].seq)
                print(iseqs[0].description)
                assert iseqs[0].seq == refseq.seq and iseqs[0].description == refseq.description
                iseqs = iseqs[1:]
                assert len(iseqs) == min(chunksize, len(spikes_df) - i)
                aligned_rbds += iseqs

    assert len(aligned_rbds) == len(spikes_df)

    rbd_df = (
        pd.DataFrame({'seqrecord': aligned_rbds})
        .assign(description=lambda x: x['seqrecord'].map(lambda rec: rec.description),
                country=lambda x: x['description'].str.split('|').str[-1],
                host=lambda x: x['description'].str.split('|').str[6].str.strip(),
                length=lambda x: x['seqrecord'].map(len),
                n_ambiguous=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('X') + rec.seq.count('x')),
                n_gaps=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('-')),
                all_valid_nts=lambda x: x['seqrecord'].map(lambda rec: re.fullmatch(f"[{protein_letters}]+",
                                                                                    str(rec.seq)) is not None),
                )
    )

    assert all(rbd_df['length'] == len(refseq))

    return rbd_df

def max_muts(rbd_df, refseq_str, max_muts):
    rbd_df = (
        rbd_df.assign(seq = lambda x: x['seqrecord'].map(lambda rec: str(rec.seq)),
                      n_mutations = lambda x: x['seq'].map(lambda s: sum(x != y for x, y in zip(s, refseq_str))))
    )

    rbd_df = rbd_df.query('n_mutations <= @max_muts')

    return rbd_df

def write_output(rbd_df, outputfile, site_offset, refseq_str):

    records = []
    for tup in rbd_df[['seq', 'country']].itertuples():
        for isite, (mut, wt) in enumerate(zip(tup.seq, refseq_str), start=1):
            if mut != wt:
                records.append((isite, isite + site_offset, wt, mut, tup.country))

    muts_df = (pd.DataFrame.from_records(records,
                                         columns=['isite', 'site', 'wildtype', 'mutant', 'country'])
               .groupby(['isite', 'site', 'wildtype', 'mutant'])
               .aggregate(count=pd.NamedAgg('country', 'count'),
                          n_countries=pd.NamedAgg('country', 'nunique'))
               .reset_index()
               .sort_values('count', ascending=False)
               .assign(frequency=lambda x: x['count'] / len(rbd_df))
               )

    print(f'Writing mutation counts to {outputfile}')
    muts_df.to_csv(outputfile, index = False)


if __name__ == '__main__':
    inputfasta = "data/spikenuc0523_1.fasta"
    outputfilename = 'result/test_2.csv'
    wildtype = "data/wildtype_sequence.fasta"
    ref_seq = Bio.SeqIO.read(wildtype, 'fasta')
    parameter = {'min_length': 3800,
                 'max_length': 3900,
                 'max_ambig': 100,
                 'ref_seq': ref_seq,
                 'refprotname': wildtype,
                 'mafft': "C:\Program Files\mafft-win\mafft.bat",
                 'max_muts': 150,
                 'site_offset': 330
                 }

    createMutFile(inputfasta, outputfilename, parameter)

