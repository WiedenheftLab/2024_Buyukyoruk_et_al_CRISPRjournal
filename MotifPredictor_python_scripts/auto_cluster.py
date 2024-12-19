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
    os.system("usearch -cluster_fast " + filenames[i] + " -id 0.6 -centroids " + filenames[i].split('.fasta')[0] + '_uclust_60.fasta -uc ' + filenames[i].split('.fasta')[0] + '_ucluster_60.uc')
    os.system('mv ' + filenames[i] + ' originals')
    os.system('mv *.uc originals')
