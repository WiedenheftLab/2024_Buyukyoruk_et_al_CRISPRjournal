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
    os.system("fasta_splitter_by_taxa.py -i " + filenames[i] + " -l bac_arc_db_06102021_acc.txt")
    os.system('mv ' + filenames[i] + ' originals')