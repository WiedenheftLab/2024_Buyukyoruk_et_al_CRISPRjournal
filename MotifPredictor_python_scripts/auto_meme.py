#!/usr/bin/python

import os
import subprocess

workdir = os.path.abspath(os.getcwd())

filenames = []

for file in os.listdir(workdir):
    if file.endswith(".fasta"):
        filenames.append(file)

filenames.sort()

for i in range(len(filenames)):
    proc = subprocess.Popen("grep -c '>' " + filenames[i], shell=True, stdout=subprocess.PIPE, )
    seq_num = int(proc.communicate()[0].split('\n')[0])
    if seq_num >= 6:
        subtype = filenames[i].split("06102021_")[1].split("_leaders")[0]
        taxa = filenames[i].split("_")[9]
        out = subtype + "_" + taxa + "_motif.txt"
	if seq_num <= 6:
            os.system("python meme_cutoff.py -i " + filenames[i] + " -it 3 -o " + out)
            print("python meme_cutoff.py -i " + filenames[i] + " -it 3 -o " + out)
	else:
            os.system("python meme_cutoff.py -i " + filenames[i] + " -it 3 -o " + out)
            print("python meme_cutoff.py -i " + filenames[i] + " -it 3 -o " + out)

# os.system("fasta_splitter_by_taxa.py -i " + filenames[i] + " -l bac_arc_db_06102021_acc.txt")
