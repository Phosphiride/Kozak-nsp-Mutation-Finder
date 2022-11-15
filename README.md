# Kozak-nsp-Mutation-Finder
A mutation finder program adapted from Bloom Lab's SARS-CoV-2-RBD_MAP_HAARVI_sera

Original project: https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP_HAARVI_sera



alias py=python3
create virtual environment for this project
activate virtual enviornment
pip install -r requirement.txt


download mafft

upload datafile

scp -P [port numbers] test.txt login@ip:test.txt

ls to list files
mv to move a file (rename) mv sourcefile targetfile

cp (copy source/target file)

rm (delete)

sudo (admin perms, add to beginning of command)

tar xvf *filename (unzip tar file)

tar xvf ./allnuc0521.tar.xz

Transfer files:

1. sudo ls share (obtain list of files in share directory, use to find file you wish to transfer)
2. sudo scp -P [port numbers] share/[FILE NAME PLS REPLACE] login@ip:[FILE NAME PLS REPLACE] (setup transfer)
3. login on ssh server
4. ls (look at directory)
5. cd [destination] (change directory to file destination)
6. mv ../[FILE NAME PLS REPLACE] (move file under current directory, . is for same directory as above)
7. tar xvf [FILE NAME PLS REPLACE] (unzip tar file if needed)
8. rm -f -r allnuc0521 ( remove allnuc0521 directory to free space)

pyhthon3 -m venv .  ( to create virutal enviroment for this project , only run it once)

inputfile
