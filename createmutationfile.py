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

    alignment(50000, seq_df, parameter["refprotname"], parameter["ref_seq"])

def alignment(chunksize, spikes_df, refprotfile, refseq):
    aligned_rbds = []

    for i in range(0, len(spikes_df), chunksize):
        spikes_file = os.path.join("temp/",
                                   f"human_full-length_spikes_{i + 1}-to-{i + chunksize}.fasta")
        print(f"Writing spikes {i + 1} to {i + chunksize} to {spikes_file}")
        _ = Bio.SeqIO.write(spikes_df['seqrecord'].tolist()[i: i + chunksize], spikes_file, 'fasta')
        print('Now aligning these sequences...')
        # cmds = ['mafft', '--auto', '--thread', str(config['max_cpus']),
        #        '--keeplength', '--addfragments', spikes_file, refprotfile]

        cmds = [parameter["mafft"], '--auto', '--thread', '-1',
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

    return aligned_rbds


if __name__ == '__main__':
    inputfasta = "data/spikenuc0523_1.fasta"
    outputfilename = 'result/test_1.csv'
    wildtype = "data/wildtype_sequence.fasta"
    ref_seq = Bio.SeqIO.read(wildtype, 'fasta')
    parameter = {'min_length': 3800,
                 'max_length': 3900,
                 'max_ambig': 100,
                 'ref_seq': ref_seq,
                 'refprotname': wildtype,
                 'mafft': "C:\Program Files\mafft-win\mafft.bat"
                 }

    createMutFile(inputfasta, outputfilename, parameter)

