#!/usr/bin/python

import argparse
import re
import tqdm
import sys
import subprocess
from Bio import SeqIO
import os


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str, dest='filename', help='Specify a fastafile.\n')
parser.add_argument('-l', '--list', required=True, type=str, dest='taxa', help='Specify taxonomy file.\n')

results = parser.parse_args()
filename = results.filename
taxa = results.taxa

proc = subprocess.Popen("grep -c '>' " + filename, shell=True, stdout=subprocess.PIPE, )
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length)) as pbar:
    pbar.set_description('Splitting...')
    for record in SeqIO.parse(filename, "fasta"):
        pbar.update()
        acc = record.id.rsplit('_',2)[0]
        # print acc
        try:
            proc = subprocess.Popen("grep '" + acc + "' " + taxa, shell=True, stdout=subprocess.PIPE, )
            info = str(proc.communicate()[0].split('\n')[0])

            kingdom = info.split('\t')[2]

            if info.split('\t')[3] == "Proteobacteria":
                phylum = info.split('\t')[4].split('\n')[0].replace(" ","")
            else:
                phylum = info.split('\t')[3].replace(" ","")

            # phylum = info.split('\t')[3].replace(" ","")
            out = filename.split(".fasta")[0] + "_" + phylum + ".fasta"
            # print out
            # print phylum
            f = open(out, 'a')
            sys.stdout = f

            print ">" + record.description + " | " + kingdom + " | " + phylum
            print re.sub("(.{60})", "\\1\n", str(record.seq), 0, re.DOTALL)

        except:

            out = filename.split(".fasta")[0] + "_NA.fasta"
            # print out
            f = open(out, 'a')
            sys.stdout = f

            print ">" + record.description + " | NA | NA"
            print re.sub("(.{60})", "\\1\n", str(record.seq), 0, re.DOTALL)

os.system("mv " + filename + " originals")