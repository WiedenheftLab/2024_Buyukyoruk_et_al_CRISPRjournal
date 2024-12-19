#!/usr/bin/python

import argparse
import os
import tqdm
import sys
import time
import subprocess
from Bio import SeqIO
import re
import pdb
import statistics
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str, dest='filename', help='Specify a fastafile.\n')
parser.add_argument('-it', '--iteration', required=True, type=int, dest='iter', help='Specify the number of iteration to perform to identify E-value cutoff.\n')
parser.add_argument('-o', '--out', required=True, type=str, dest='out', help='Specify a output file.\n')
parser.add_argument('-hs', '--hsfrac', required=False, default=0.5,type=float, dest='frac', help='Specify hsfact score.\n')


results = parser.parse_args()
filename = results.filename
iter = results.iter
out = results.out
frac = results.frac

orig_stdout = sys.stdout

os.mkdir('tmp')
os.chdir('tmp')

cutoff_list = []

stderr = None

meme_test = 'meme ../' + filename + ' -dna -oc test -mod anr -nmotifs 1 -minw 6 -maxw 40 -objfun se -allw -markov_order 1 -maxiter 1000 -p 120 -nostatus -shuf 2 -hsfrac ' + str(frac)
#print meme_test
if "_IE_" in filename:
    frac = 0.5
else:

    proc = subprocess.Popen(meme_test, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    if stderr != '':
        frac = float(str(stderr.split(" least ")[1].split()[0])[:-1])
    else:
        frac = 0.75


print "Using -hsfrac " + str(frac)

with tqdm.tqdm(range(iter), desc = 'Defining E-value cutoff for ' + filename ) as pbar:
    for i in range(iter):
        pbar.update()

        shuf_file ='shuffled_' + filename.split('.fasta')[0] + '_v' + str(i) + '.fasta'
        cmd = 'fasta-shuffle-letters -dna -copies 1 ../' + filename + ' ' + shuf_file
        os.system(cmd)

        bfile = shuf_file.split(".fasta")[0] + "_bg.txt"
        cmd = 'fasta-get-markov -nostatus -nosummary -dna -m 1 ' + shuf_file + ' ' + bfile
        os.system(cmd)

        # neg_file = shuf_file.split(".fasta")[0] + '_neg.fasta'
        # cmd = 'fasta-shuffle-letters -dna -copies 1000 ' + shuf_file + ' ' + neg_file
        # os.system(cmd)

        meme = 'meme ' + shuf_file + ' -dna -oc s' + str(i) + ' -mod anr -nmotifs 1 -minw 6 -maxw 40 -objfun se -allw -markov_order 1 -maxiter 1000 -p 120 -nostatus -shuf 2 -bfile ' + bfile + " -hsfrac " + str(frac) # + " -neg " + neg_file

        os.system(meme)

        os.chdir('s'+str(i))
        with open('meme.txt','rU') as file:
            for line in file:
                if 'E-value =' in line:
                    eval = line.split('E-value = ')[1].split('\n')[0]
                    if float(eval) < 0.1:
                        cutoff_list.append(float(eval))

        os.chdir('../')
    os.chdir('../')
time.sleep(1)
os.system('rm -r tmp/')

try:
    cutoff = min(cutoff_list) * 10
except:
    cutoff = 1

print "\nUsing " + str(cutoff) + " as E-value cutoff!!!\n"

bfile = filename.split(".fasta")[0] + "_bg.txt"
cmd = 'fasta-get-markov -dna -nostatus -nosummary -m 1 ' + filename + ' ' + bfile
os.system(cmd)

# neg_file = filename.split(".fasta")[0] + '_neg.fasta'
# cmd = 'fasta-shuffle-letters -dna -copies 1000 ' + filename + ' ' + neg_file
# os.system(cmd)

meme = 'meme ' + filename + ' -dna -oc ' + filename.split('.fasta')[0] + '_motif -mod anr -nmotifs 6 -minw 6 -maxw 40 -allw -objfun se -markov_order 1 -p 120 -maxiter 1000 -shuf 2 -nostatus -bfile ' + bfile + " -hsfrac " + str(frac) # + " -neg " + neg_file

print meme + '\n'

os.system(meme)

os.chdir(filename.split('.fasta')[0] + '_motif')

workdir = os.path.abspath(os.getcwd())

for files in os.listdir(workdir):
    if "meme.txt" in files:
        os.system("meme2meme meme.txt -bg ../" + bfile + ' > min_meme.txt')
    os.system('> ' + out)

    f = open(out, 'a')
    sys.stdout = f

current_array = []

stat = None

try:
    with open('min_meme.txt','rU') as file:
        for line in file:
            if line[0] != ' ' and "MOTIF" not in line and "MEME" not in line and "letter-probability matrix" not in line:
                print line.split('\n')[0]
            if "MOTIF" in line and "MEME" in line:
                current_array.append(line)
            if "letter-probability matrix" in line:
                eval_found = float(line.split("E= ")[1].split("\n")[0])
                if eval_found < cutoff:
                    stat = True
                    print current_array[0]
                    print line.split('\n')[0]
                else:
                    stat = False
                    current_array = []
            if line[0] == ' ' and "A" not in line and stat == True:
                print line.split('\n')[0]
            if line[0] == ' ' and "A" in line:
                print line.split('\n')[0]
            if "version" in line:
                print line.split('\n')[0]

    f.close()
    time.sleep(0.5)
    sys.stdout = orig_stdout

    os.system("meme2meme " + out + " -bg ../" + bfile + ' -numbers > ../' + out)
except:
    pass

os.chdir('../')

proc = subprocess.Popen("wc -l < " + out , shell=True, stdout=subprocess.PIPE, )
info = int(proc.communicate()[0])

print "\nParsing motifs found in " + out + "\n"

if info != 0:

    proc = subprocess.Popen("grep -c 'MOTIF' " + out, shell=True, stdout=subprocess.PIPE, )
    motif_count = int(proc.communicate()[0])

    for f in range(1,motif_count+1):

        fimo_test = "fimo --norc --motif " + str(f) + " --oc " + filename.split('.fasta')[0] + "_fimo_train " + out + " " + filename

        os.system(fimo_test)

        os.chdir(filename.split('.fasta')[0] + '_fimo_train')

        df = pd.read_csv("fimo.tsv",sep='\t')

        p_val = float(df["p-value"].max())

	p_val_expected = float(df["p-value"].min())

        os.chdir('../')

        fimo = "fimo --norc --thresh " + str(p_val) + " --motif " + str(f) + " --oc " + filename.split('.fasta')[0] + "_fimo_" + str(f) + ' ' + out + " database/Random_200bp_seq_20220708_142400.fasta"

        print '\n' + fimo + '\n'

        os.system(fimo)

        os.chdir(filename.split('.fasta')[0] + '_fimo_' + str(f))

        os.system("python ../motif_stack.py -i fimo.tsv -o " + out + " -p " + str(p_val_expected))

        os.system("mv " + out.split(".txt")[0] + "* ../")

        os.chdir('../')



