#!/Users/muratbuyukyorukmsu/PycharmProjects/python/venv/bin/python

import argparse
import os
import tqdm
import sys
import time
import subprocess
from Bio import SeqIO
import re
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str, dest='filename', help='Specify a fastafile.\n')
parser.add_argument('-r', '--repeat', required=True, type=str, dest='repeat', help='Specify a fastafile with repeat only seq.\n')
parser.add_argument('-o', '--out', required=True, type=str, dest='out', help='Specify a output file.\n')

results = parser.parse_args()
filename = results.filename
repeat = results.repeat
out = results.out

repeat_id = []
repeat_seq = []

proc = subprocess.Popen("wc -l < " + repeat, shell=True, stdout=subprocess.PIPE, )
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length), desc = 'Reading...' ) as pbar:
    with open(repeat,'rU') as file:
        for line in file:
            pbar.update()
            if "Array_no" not in line:
                array = line.split('\t')[0]
                repeat_info = line.split('\t')[10]
                repeat_id.append(array)
                repeat_seq.append(repeat_info)

proc = subprocess.Popen("grep -c '>' " + filename, shell=True, stdout=subprocess.PIPE, )
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length), desc = 'Reading...' ) as pbar:
    f = open(out, 'a')
    sys.stdout = f
    for record in SeqIO.parse(filename, "fasta"):
        pbar.update()
        ind = repeat_id.index(record.id)
        repeat_found = repeat_seq[ind]
        print ">" + record.description
        print re.sub("(.{60})", "\\1\n", str(record.seq)[:-len(repeat_found)],0, re.DOTALL)


