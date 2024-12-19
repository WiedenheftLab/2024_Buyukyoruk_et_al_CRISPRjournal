#!/usr/bin/python

import os
import subprocess
from Bio import SeqIO
import re
import sys

workdir = os.path.abspath(os.getcwd())

filenames = []

for file in os.listdir(workdir):
    if file.endswith(".fasta"):
        filenames.append(file)

filenames.sort()

for i in range(len(filenames)):

    out = filenames[i].split('.fasta')[0] + "_currated.fasta"

    f = open(out, 'a')
    sys.stdout = f

    os.system("perl /home/muratb/Desktop/Local_Softwares/CRISPRDetect_3.0/CRISPRDetect3 -q -f " + filenames[i] + " -o CRISPRDetect3_search -array_quality_score_cutoff 3 -T 50 > log.txt")

    proc = subprocess.Popen("grep '>' CRISPRDetect3_search | cut -d'>' -f2 | cut -d' ' -f1", shell=True, stdout=subprocess.PIPE, )
    CRISPR = (proc.communicate()[0]).split('\n')[:-1]

    for record in SeqIO.parse(filenames[i], "fasta"):
        if record.id in CRISPR:
            continue
        else:
            print ">" + record.description
            print re.sub("(.{60})", "\\1\n", str(record.seq), 0, re.DOTALL)

    os.system('mv ' + filenames[i] + ' originals')

os.system('rm CRISPRDetect3_search*')
os.system('rm log.txt')