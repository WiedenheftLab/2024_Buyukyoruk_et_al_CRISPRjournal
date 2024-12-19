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
    os.system("python CRISPRleader_remover.py -i " + filenames[i] + " -r bac_arc_db_06102021_CRISPRDetect_20210612_104054_summary_with_plasmid.txt -o " + filenames[i].split('.fasta')[0] + "_repeat_removed.fasta")
    os.system('mv ' + filenames[i] + ' originals')