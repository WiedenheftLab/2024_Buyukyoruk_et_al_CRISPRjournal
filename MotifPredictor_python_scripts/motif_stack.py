#!/usr/bin/python

import argparse
import sys
import os
import subprocess
import re
import textwrap
import time
import pdb
import pandas as pd
import math


try:
    import tqdm
except ImportError, e:
    print "tqdm module is not installed! Please install tqdm and try again."
    sys.exit()

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str, dest='filename',
                    help='Specify a original fasta file.\n')
parser.add_argument('-o', '--output', required=True, dest='out',
                    help='Specify a output fasta file name.\n')
parser.add_argument('-p', '--pval', required=True, dest='pval_expect',
                    help='Specify an expected pval.\n')

results = parser.parse_args()
filename = results.filename
out = results.out
pval_expect = float(results.pval_expect)

ori_out = sys.stdout

acc_list = []

os.system('rm tmp*')
#df_pre = pd.read_csv(filename)
#new = df_pre.columns.str.split('\t')
#df_pre = df_pre[df_pre.columns[0]].str.split('\t',expand=True)
#df_pre.columns = new.tolist()[0]
#min_score = df_pre['p-value'].astype(float).min()
min_score = pval_expect
min_score_conversion = -math.log10(min_score)

with open(filename,'rU') as file:
    for line in file:
        if "motif_id" not in line and line != "\n" and line[0] != "#":
            motif_id = line.split("\t")[0]
            acc = line.split("\t")[2]
            if "|" in acc:
                acc_fix = acc.split('|')[0]
                acc = acc_fix
            mid = 200 - ((int(line.split("\t")[3])+int(line.split("\t")[4]))/2)
            pval = float(line.split("\t")[7])
            score_conv = -math.log10(pval)*1/min_score_conversion
            
            f = open("tmp.txt", 'a')
            sys.stdout = f

            print acc + "\t" + str(mid) + '\t' + str(score_conv)
            # time.sleep(0.01)

sys.stdout = ori_out

print"Identifying unique accessions..."

acc_list = []
tmp = 'tmp.txt'
with open(tmp,'rU') as file:
    for line in file:
        acc = line.split('\t')[0]
        acc_list.append(acc)

res = [m for n, m in enumerate(acc_list) if m not in acc_list[:n]]

print"Merging same motif occurrances in a leader..."

out_file = out.split('.txt')[0] + '_' + str(motif_id) + '.txt'

os.system('> ' + out_file)

f = open(out_file, 'a')

with tqdm.tqdm(range(len(res)), desc = 'Processing:') as pbar:
    for k in range(len(res)):
        pbar.update()
        if "|" in res[k]:
            acc_search = res[k].split("|")[0]
            proc = subprocess.Popen("grep '" + acc_search + "|' " + tmp, shell=True, stdout=subprocess.PIPE, )
            info = (proc.communicate()[0]).replace(res[k],'').split("\n")[:-1]
        else:
            acc_search = res[k]
            proc = subprocess.Popen("grep '" + acc_search + "\t' " + tmp, shell=True, stdout=subprocess.PIPE, )
            info = (proc.communicate()[0]).replace(res[k],'').split("\n")[:-1]
        sys.stdout = f
        print acc_search + ''.join(info)

f.close()

sys.stdout = ori_out

print "Formatting..."

df = pd.read_csv(out_file,header=None)

length = len(df[0].str.split('\t', expand=True).columns)

add_line = []

for r in range((length+1)/2):
    if r != 0:
        add_line.append(out.split('.txt')[0] +"_" + str(motif_id) +'_' + str(r))
        add_line.append(out.split('.txt')[0] +"_" + str(motif_id) +'_' + str(r) + '_pvalue')
    if r == 0:
        add_line.append("Array")
merged_line = '\t'.join(add_line)
time.sleep(1)
os.system('echo "' + merged_line + '" | cat - ' + out_file + " > out" + str(motif_id) + ".txt")
os.system("mv out" + str(motif_id) + ".txt " + out_file)
time.sleep(1)
# os.system("rm tmp*")
