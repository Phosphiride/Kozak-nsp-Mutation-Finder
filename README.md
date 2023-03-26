# Kozak-nsp-Mutation-Finder
A mutation finder program adapted from Bloom Lab's SARS-CoV-2-RBD_MAP_HAARVI_sera, used to
collect data on single nucleotide polymorphisms (SNPs) in SARS-CoV-2 for a senior thesis.

How to Use:
- Install Anaconda or another python environment manager
- Run pip install -r requirements.txt
- Download mafft from https://mafft.cbrc.jp/alignment/software/
- In createmutfiles.py, fastasplitter.py, and masteragg.py, edit path names and parameters to fit project and gene of interest.
- Run fastasplitter.py if needed, to split files into more manageable portions and increase processing speed.
- Run createmutfiles.py to isolate mutations, then masteragg.py to isolate only mutations of interest (if needed).

Results are listed as a .csv file, which can be imported into spreadsheet programs like Microsoft Excel or Google Sheets, or used
in other projects.

Original project: https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP_HAARVI_sera