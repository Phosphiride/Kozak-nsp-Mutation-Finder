import os
import io
import re
import subprocess

import Bio.SeqIO
import Bio.Seq
from Bio.Data.IUPACData import unambiguous_dna_letters, protein_letters
import pandas as pd


def createMutFile(inputfasta, outputfilename, parameter):
    sequences = list(Bio.SeqIO.parse(inputfasta, 'fasta'))  # parse sequences into python
    print(f'Read {len(sequences)} sequences.')

    seq_df = (
        pd.DataFrame({'seqrecord': sequences})
        .assign(description=lambda x: x['seqrecord'].map(lambda rec: rec.description),
                country=lambda x: x['description'].str.split('/').str[1],
                #host=lambda x: x['description'].str.split('|').str[6].str.strip(),
                length=lambda x: x['seqrecord'].map(len),
                n_ambiguous=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('X') + rec.seq.count('x')))
    )

    max_length = parameter['max_length']
    min_length = parameter['min_length']

    print(f'Only keeping sequences with lengths between {min_length} and {max_length}.')

    seq_df = (
        seq_df.assign(valid_length=lambda x: (min_length <= x['length']) & (x['length'] <= max_length))
    )

    seq_df = seq_df.query('valid_length')

    max_ambig = parameter["max_ambig"]

    print(f'Filtering sequences with greater than {max_ambig} ambiguous bases.')

    seq_df = (
        seq_df.assign(excess_ambiguous=lambda x: x['n_ambiguous'] > max_ambig)
    )

    temp_file = "temp/human_full-length-seq.fasta"

    Bio.SeqIO.write(seq_df['seqrecord'].tolist(), temp_file, 'fasta')
    ref_seq = parameter['ref_seq']
    gen_df = alignment(50000, seq_df, parameter["refprotname"], ref_seq, parameter['mafft'])

    gen_df = gen_df.query('n_ambiguous == 0').query('n_gaps == 0')
    assert gen_df['all_valid_nts'].all()
    print(f'Retained {len(gen_df)} sequences')

    refseq_str = str(ref_seq.seq)
    refseq_aa_str = str(ref_seq.translate)
    ref_df = max_muts(gen_df, refseq_str, parameter['max_muts'])

    site_offset = parameter['site_offset']

    write_output(ref_df, outputfilename, site_offset, refseq_str, refseq_aa_str)


def alignment(chunksize, genes_df, refprotfile, refseq, mafft):
    aligned_genes = []

    for i in range(0, len(genes_df), chunksize):
        gene_file = os.path.join("temp/",
                                   f"human_full-length_spikes_{i + 1}-to-{i + chunksize}.fasta")
        print(f"Writing genes {i + 1} to {i + chunksize} to {gene_file}")
        _ = Bio.SeqIO.write(genes_df['seqrecord'].tolist()[i: i + chunksize], gene_file, 'fasta')
        print('Now aligning these sequences...')
        # cmds = ['mafft', '--auto', '--thread', str(config['max_cpus']),
        #        '--keeplength', '--addfragments', gene_file, refprotfile]

        cmds = [mafft, '--auto', '--thread', '-1',
                '--keeplength', '--addfragments', gene_file, refprotfile]

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
                assert len(iseqs) == min(chunksize, len(genes_df) - i)
                aligned_genes += iseqs

    assert len(aligned_genes) == len(genes_df)

    gen_df = (
        pd.DataFrame({'seqrecord': aligned_genes})
        .assign(description=lambda x: x['seqrecord'].map(lambda rec: rec.description),
                country=lambda x: x['description'].str.split('/').str[1],
                #host=lambda x: x['description'].str.split('/').str[6].str.strip(),
                length=lambda x: x['seqrecord'].map(len),
                n_ambiguous=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('X') + rec.seq.count('x')),
                n_gaps=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('-')),
                all_valid_nts=lambda x: x['seqrecord'].map(lambda rec: re.fullmatch(f"[{unambiguous_dna_letters}]+",
                                                                                    str(rec.seq)) is not None),
                )
    )

    if parameter['exclude_ambig'] == True:
        before_exclude = len(gen_df)
        gen_df = gen_df[gen_df['all_valid_nts'] == True]
        after_exclude = len(gen_df)
        print(f'Retained {after_exclude} rows out of {before_exclude}.')

    assert all(gen_df['length'] == len(refseq))

    gen_df = (
        gen_df.assign(all_valid_prot = lambda x: x['seqrecord'].map(lambda rec: rec.seq.translate()))
    )

    return gen_df


def max_muts(gen_df, refseq_str, max_muts):
    gen_df = (
        gen_df.assign(seq=lambda x: x['seqrecord'].map(lambda rec: str(rec.seq)),
                      n_mutations=lambda x: x['seq'].map(lambda s: sum(x != y for x, y in zip(s, refseq_str))))
    )

    gen_df = gen_df.query('n_mutations <= @max_muts')

    return gen_df


def write_output(gen_df, outputfile, site_offset, refseq_str, refseq_aa_str):
    records = []
    for tup in gen_df[['seq', 'country']].itertuples():
        for isite, (mut_nt, wt_nt) in enumerate(zip(tup.seq, refseq_str), start=1):
            if mut_nt != wt_nt:
                records.append((isite, isite + site_offset, wt_nt, mut_nt))

    for tup in gen_df[['all_valid_prot', 'country']].itertuples():
        for aasite, (mut_aa, wt_aa) in enumerate(zip(tup.all_valid_prot, refseq_aa_str), start=1):
            if mut_aa != wt_aa:
                records.append((aasite, wt_aa, mut_aa, tup.country))

    muts_df = (pd.DataFrame.from_records(records,
                                         columns=['gene site', 'genome site', 'wt nt', 'mutant nt','aa site', 'wt aa', 'mutant aa', 'country', 'count'])
               .groupby(['gene site', 'genome site', 'wt nt', 'mutant nt', 'aa site', 'wt aa', 'mutant aa'])
               .aggregate(count=pd.NamedAgg('country', 'count'),
                          n_countries=pd.NamedAgg('country', 'nunique'))
               .reset_index()
               .sort_values('count', ascending=False)
               .assign(frequency=lambda x: x['count'] / len(gen_df))
               )

    print(f'Writing mutation counts to {outputfile}')
    muts_df.to_csv(outputfile, index=False)


if __name__ == '__main__':
    inputfasta = "data/20210613_gisaid_genomes.fasta"
    outputfilename = 'result/test_13_aa.csv'
    wildtype = "data/GISAID_nsp5.fasta"
    ref_seq = Bio.SeqIO.read(wildtype, 'fasta')
    parameter = {'min_length': 29500,
                 'max_length': 30000,
                 'max_ambig': 100,
                 'ref_seq': ref_seq,
                 'refprotname': wildtype,
                 'mafft': "C:/Program Files/mafft-win/mafft.bat",
                 'max_muts': 100000,
                 'site_offset': 10055,     #   nsp5: 10055; nsp12: 13442
                 'exclude_ambig': True
                 }

    createMutFile(inputfasta, outputfilename, parameter)
