#!/software/mambaforge/envs/Murat_scripts/bin/python

import argparse
import sys
import os
import subprocess
import re
import textwrap
import distutils.spawn
import multiprocessing
from multiprocessing import Pool
import threading
import math
import time
import uuid
import multiprocessing.pool
from contextlib import closing
import pdb

orig_stdout = sys.stdout

timestr = time.strftime("%Y%m%d_%H%M%S")
#timestr="20230201_210508"
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str, dest='filename', help='Specify a original fasta file.\n')
parser.add_argument('-o', '--output', required=True, dest='out', help='Specify a output fasta file name.\n')

###### CRISPRDetect parameters ######

parser.add_argument('-f', '--file_type', required=False, type=str, dest='file_type', default='fasta',
                    help='Fasta "fasta" or Genbank "gbk" file.\n')
parser.add_argument('-t', '--thread', required=False, type=str, dest='thread', default='slow',
                    help='fast = 4/5 of all threads, slow = 1/5 of all threads, manual = used defined number.\n')
parser.add_argument('-wl', '--word_length', required=False, type=int, dest='word_length', default='11',
                    help='This is the default word length CRISPRDetect uses to find the putative CRISPRs. Any positive integer >=6 can be used.\n')
parser.add_argument('-wr', '--minimum_word_repeatation', required=False, type=int, dest='word_repeatation', default='2',
                    help='By default CRISPRDetect uses 2 repeating identical words to find putative CRISPRs.\n')
parser.add_argument('-g', '--max_gap_between_crisprs', required=False, type=int, dest='gap', default='125',
                    help='By default the maximum gap is set tp 125 nucleotides between the repeating identical seed words.\n')
parser.add_argument('-rlc', '--repeat_length_cutoff', required=False, type=int, dest='repeat_length_cutoff',
                    default='17',
                    help='After the intial processing, putative CRISPRs with repeat lengths less than this value will be rejected.\n')
parser.add_argument('-rl', '--minimum_repeat_length', required=False, type=int, dest='minimum_repeat_length',
                    default='23', help='Minimum length of repeats.\n')
parser.add_argument('-rn', '--minimum_no_of_repeats', required=False, type=int, dest='minimum_no_of_repeats',
                    default='2',
                    help='Predicted CRISPRs with number of repeats less than this value will be excluded.\n')
parser.add_argument('-s', '--array_quality_score_cutoff', required=False, type=float, dest='score', default='1.5',
                    help='Predicted CRISPRs with score less than this value will be excluded from the output file.(1.5 is for less stringent search)\n')
parser.add_argument('-lf', '--left_flank_length', required=False, type=int, dest='left_flank_length', default='500',
                    help="This is the default length of the 5' (upstream) region of the CRISPRs.\n")
parser.add_argument('-rf', '--right_flank_length', required=False, type=int, dest='right_flank_length', default='500',
                    help="This is the default length of the 3' (downstream) region of the CRISPRs.\n")
parser.add_argument('-c', '--annotate_cas_genes', required=False, type=int, dest='annotate_cas_genes', default='0',
                    help="By default this option is turned off. Use 1 to enable this option [Requires GeneMarkS package installed and gmhmmp.pl file in the PATH].\n")
parser.add_argument('-hse', '--hmm_seq_evalue_cutoff', required=False, type=float, dest='hmm_seq_evalue_cutoff',
                    default='0.00001', help="Specify sequence e-value cutoff for hmmsearch. Default value is 0.00001\n")
parser.add_argument('-hde', '--hmm_dom_evalue_cutoff', required=False, type=float, dest='hmm_dom_evalue_cutoff',
                    default='0.00001', help="Specify domain e-value cutoff for hmmsearch. Default value is 0.00001\n")
parser.add_argument('-w', '--wgs', required=False, type=int, dest='wgs', default='0',
                    help="Specify if all sequences in the input sequence file belongs to the same genome [useful with option '-annotate_cas_genes 1']. Default value is 0 [Zero].\n")
parser.add_argument('-p', '--check_plasmid', required=False, type=int, dest='check_plasmid', default='0',
                    help="By default this option is turned off. Use 1 to enable this option [Requires RefSeq Plasmid sequence DB].\n")

results = parser.parse_args()
filename = results.filename
out = results.out
file_type = results.file_type
thread = results.thread
word_length = results.word_length
word_repeatation = results.word_repeatation
gap = results.gap
repeat_length_cutoff = results.repeat_length_cutoff
minimum_repeat_length = results.minimum_repeat_length
minimum_no_of_repeats = results.minimum_no_of_repeats
score = results.score
left_flank_length = results.left_flank_length
right_flank_length = results.right_flank_length
annotate_cas_genes = results.annotate_cas_genes
hmm_seq_evalue_cutoff = results.hmm_seq_evalue_cutoff
hmm_dom_evalue_cutoff = results.hmm_dom_evalue_cutoff
wgs = results.wgs
check_plasmid = results.check_plasmid

working = "temp_" + timestr + "_" + str(uuid.uuid4().hex)
temp_dir = "tmp/tmp_" + str(uuid.uuid4().hex)
#temp_dir="tmp_ddfe72020d3c45c9be673970cc84e9c5"
print temp_dir

########################################################################################################################
########################################################################################################################

## Defining system ##
proc = subprocess.Popen("uname -s", shell=True, stdout=subprocess.PIPE, )
what_is_system = (proc.communicate()[0].split('\n')[0])
print what_is_system
## Define current user ##

proc = subprocess.Popen("whoami", shell=True, stdout=subprocess.PIPE, )
whoami = (proc.communicate()[0].split('\n')[0])
print whoami
## Defining number of threads to be used ##
no_of_threads_available = float(multiprocessing.cpu_count())

if thread == 'fast':
    thread_in_use = int(math.floor(no_of_threads_available * 0.8))

if thread == 'slow':
    thread_in_use = int(no_of_threads_available - math.floor(no_of_threads_available * 0.8))

if thread == 'manual':
    thread_in_use = int(raw_input('***Please specify number of threads: '))
    if thread_in_use > no_of_threads_available:
        print "Maximum number allowed is " + no_of_threads_available + 'threads.'

print thread_in_use

########### Checking dependencies #############

try:
    from Bio import SeqIO
except ImportError, e:
    print "SeqIO module is not installed! Please install SeqIO and try again."
    sys.exit()

try:
    from tqdm import tqdm
except ImportError, e:
    print "tqdm module is not installed! Please install tqdm and try again."
    sys.exit()

if (distutils.spawn.find_executable("macsyfinder")) == None:
    print("\nERROR: macsyfinder is not installed. Please make sure macsyfinder is defined in path.\n")
    sys.exit()

if (distutils.spawn.find_executable("hmmsearch")) == None:
    print("\nERROR: hmmsearch is not installed. Please make sure hmmsearch is defined in path.\n")
    sys.exit()

if (distutils.spawn.find_executable("prodigal")) == None:
    print("\nERROR: prodigal is not installed. Please make sure prodigal is defined in path.\n")
    sys.exit()

gmhmmp = None
if (distutils.spawn.find_executable("gmhmmp")) == None:
    print("WARNING! GeneMark is not installed or cannot find in the path. Cas gene prediction won't be available. Please install GeneMark software as described in CRISPRDetect3.")
    annotate_cas_genes = 0
elif (distutils.spawn.find_executable("gmhmmp")) != None:
    gmhmmp = True

if gmhmmp == True:
    proc = subprocess.Popen("ls /home/" + whoami + "/.gm_key 2>/dev/null", shell=True, stdout=subprocess.PIPE, )
    gm_key = (proc.communicate()[0].split('\n')[0])
    if gm_key == "":
        gm_key_found = None
        while gm_key_found ==None:
            gm_key_file = raw_input(
                "GeneMark Key is not found. Please specifiy the path to the key file or request one from 'http://exon.gatech.edu/Genemark/license_download.cgi'.")
            proc = subprocess.Popen("ls " + gm_key_file + ' 2>/dev/null', shell=True, stdout=subprocess.PIPE, )
            gm_key_found = (proc.communicate()[0].split('\n')[0])
            if gm_key_found != None or gm_key_found != "":
                break
        os.system('cp ' + gm_key_found + ' /home/' + whoami + '/.gm_key')
    else:
        pass

print "Found gmhmmp!"

# Search if profiles and definitions files are available #
path_2_CRISPR_CAS_search = ""

if what_is_system == 'Linux':
    proc = subprocess.Popen("locate -br '^CRISPR_CAS_search.py$'", shell=True, stdout=subprocess.PIPE, )
    path_2_CRISPR_CAS_search = (proc.communicate()[0].split('\n')[0])
    print path_2_CRISPR_CAS_search

if what_is_system == 'Darwin':
    proc = subprocess.Popen("mdfind -name 'CRISPR_CAS_search.py'", shell=True, stdout=subprocess.PIPE, )
    path_2_CRISPR_CAS_search = (proc.communicate()[0].split('\n')[0])
    print path_2_CRISPR_CAS_search

if what_is_system == 'Cygwin':
    print("\nERROR: Windows is not supported at the moment. Please use Mac or Linux based system.\n")
    sys.exit()

if path_2_CRISPR_CAS_search != "":
    pass
else:
    print(
        "\nERROR: CRISPR_CAS_search.py is not accessible. Please make sure CRISPR_CAS_search.pl is accessible by "
        "PATH.\n")
    sys.exit()

# Find Path to CRISPRDetectv3 #

path_2_CRISPRDetect = ""

if what_is_system == 'Linux':
    proc = subprocess.Popen("locate -br '^CRISPRDetect3$'", shell=True, stdout=subprocess.PIPE, )
    path_2_CRISPRDetect = (proc.communicate()[0].split('\n')[0])
    print path_2_CRISPRDetect
if what_is_system == 'Darwin':
    proc = subprocess.Popen("mdfind -name 'CRISPRDetect3'", shell=True, stdout=subprocess.PIPE, )
    path_2_CRISPRDetect = (proc.communicate()[0].split('\n')[0])
    print path_2_CRISPRDetect
if path_2_CRISPRDetect != "":
    pass

else:
    print("\nERROR: CRISPRDetect.pl is not installed. Please make sure CRISPRDetect.pl is installed and accessible.\n")
    sys.exit()


### Preparing for plasmid check database ###

def check_plasmid_db():
    path = '"' + path_2_CRISPRDetect + '"'
    check_plasmid_script = os.path.isfile('script_make_plasmid_db.pl')
    if check_plasmid_script == True:
        pass
    else:
        path_2_script_make_plasmid_db = ""
        if what_is_system == 'Linux':
            proc = subprocess.Popen("locate -br '^script_make_plasmid_db.pl$'", shell=True, stdout=subprocess.PIPE, )
            path_2_script_make_plasmid_db = (proc.communicate()[0].split('\n')[0])

        if what_is_system == 'Darwin':
            proc = subprocess.Popen("mdfind -name 'script_make_plasmid_db.pl'", shell=True, stdout=subprocess.PIPE, )
            path_2_script_make_plasmid_db = (proc.communicate()[0].split('\n')[0])

        if path_2_script_make_plasmid_db != "":
            os.system('perl ' + path_2_script_make_plasmid_db)

        else:
            print(
                "\nERROR: Couldn't locate script_make_plasmid_db.pl. Please make sure script_make_plasmid_db.pl is installed and accessible.\n")
            sys.exit()

# prep = []
#
# if check_plasmid == 1:
#     prep1 = multiprocessing.Process(target=check_plasmid_db)
#     prep.append(prep1)
#     prep1.start()
#
#     for prep1 in prep:
#         prep1.join(timeout=0)
#
# else:
#     pass


### Setting Multiprocessing options ###

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


class Pool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


########################################################################################################################
########################################################################################################################

manager = multiprocessing.Manager()
fetched_arrays_list = manager.list()
found_new_seq_list = manager.list()
found_old_seq_list = manager.list()
fetched_arrays_accurate_return_list = manager.list()
fetched_arrays_filtered_return_list = manager.list()
fetched_spacer_return_list = manager.list()
fetched_gff_return_list = manager.list()


###### Functions ######

### Look is there are duplicate accessions in file ###

def environment_setup():
    dir = path_2_CRISPR_CAS_search.split('/CRISPR_CAS_search.py')[0]
    print dir
    proc = subprocess.Popen("pwd", shell=True, stdout=subprocess.PIPE, )
    work_dir = (proc.communicate()[0].split('\n')[0])
    print work_dir
#    os.system('chmod -R 755 . && chmod 777 ' + dir)
#    os.system('chmod -R 755 . && chmod 777 ' + work_dir)
    os.chdir(dir)
    print "Creating temp_dir!"
    print temp_dir
    os.makedirs(temp_dir)
#    print "Adjusting permissions!"
#    os.system('chmod -R 755 . && chmod 777 ' + temp_dir)
    print "Created temp_dir!"
    check = os.path.isdir('profiles')
    if check == True:
        print "Found Profiles!"
    else:

        print(
            "\nERROR: 'profiles' directory is missing. Please make sure the profiles folder containing .HMM is available in the same folder with CRISPR_Cas_search.py.\n")
        sys.exit()

    check = os.path.isdir('definitions')
    if check == True:
        print "Found Definitions!"
    else:
        print(
            "\nERROR: 'definitions' directory is missing. Please make sure the definitions folder is available in the same folder with CRISPR_Cas_search.py.\n")
        sys.exit()

    check = os.path.isdir('storage')
    if check == True:
        print "Found Storage!"
    else:
        #os.makedirs('storage')
        os.system('chmod -R 755 . && chmod 777 storage/')

    os.chdir(work_dir)

    check_file = os.listdir(work_dir)
    if filename not in check_file:
        print "There is no file named " + filename + " in the current directory."
        sys.exit()
    else:
        print "Found file!"

    return dir, work_dir

### Check if file contains duplicate accessions ###

def remove_dub_acc(filename):
    print "Importing fasta..."

    temp = 'temp_' + str(uuid.uuid4().hex) + '.txt'

    os.system(
        """awk 'BEGIN{RS=">"}NR>1{sub("\\n","\\t"); gsub("\\n",""); print RS$0}' """ + filename + """| awk '!seen[$1]++' > """ + temp)

    proc = subprocess.Popen("grep -c '>' " + temp, shell=True, stdout=subprocess.PIPE, )
    length = int(proc.communicate()[0].split('\n')[0])

    proc = subprocess.Popen("grep -c '>' " + filename, shell=True,
                            stdout=subprocess.PIPE, )
    length_prev = int(proc.communicate()[0].split('\n')[0])

    if length_prev == length:
        print "No duplicates found..."

    if length_prev > length:
        print "Ignoring " + str(length_prev - length) + " duplicate Accessions..."
        os.system('mv ' + filename + ' ' + filename.split('.')[0] + '_original.' + filename.split('.')[1])
        os.system('>' + filename)
        with tqdm(range(length), desc="Removing duplicates...", leave=False) as pbar:
            with open(temp, 'rU') as file:
                for line in file:
                    if line[0] == '>':
                        pbar.update()
                        arr = line.split('\t')
                        f = open(filename, 'a')
                        sys.stdout = f
                        print arr[0]
                        print re.sub("(.{60})", "\\1\n", str(arr[1].split('\n')[0]), 0, re.DOTALL)

    os.system('rm ' + temp)
    sys.stdout = orig_stdout

### Check folder, Remove old files, Count number of sequences in file, split if >500 ###

def split_(filename):
#    working = "temp_" + timestr + "_" + str(uuid.uuid4().hex)
    os.mkdir(working)

    proc = subprocess.Popen("grep -c '>' " + filename, shell=True, stdout=subprocess.PIPE, )
    length = int(proc.communicate()[0].split('\n')[0])

    i = 1

    with tqdm(range(length)) as pbar:
        pbar.set_description('Reading...')
        for record in SeqIO.parse(filename, "fasta"):
            pbar.update()
            out_temp = 'temp_fasta_' + str(i) + '.fasta'
            f = open(working + "/" + out_temp, 'a')
            sys.stdout = f
            print ">" + record.description
            print re.sub("(.{60})", "\\1\n", str(record.seq), 0, re.DOTALL)
#	    time.sleep(5)
            f.close()
            i += 1

    os.chdir(working)

    sys.stdout = orig_stdout
#    time.sleep(5)
    return length

### Execute depending on number or files ###

def execute(filename):
    fasta_files = os.listdir(work_dir+"/" + working)
    try:
        fasta_files.remove('.*')

    except:
        pass

    fasta_files.sort()

#####sequential command#####
#    for file in fasta_files:
#         CRISPRDetect_run(file)

#####parallel command#####
    if no_of_seq <= 50:
        with closing(Pool(thread_in_use)) as p:
            p.map(CRISPRDetect_run, fasta_files)
            p.close()
            p.join()
    else:
        with tqdm(range(len(fasta_files))) as pbar:
            pbar.set_description('CRISPRDetect...')
            with closing(Pool(thread_in_use)) as p:
                for _ in p.imap_unordered(CRISPRDetect_run,fasta_files):
                    pbar.update()
                p.close()
                p.join()

### Prepare fasta file for processing ###

def fasta_formatter(filename):
    acc_list = []
    description_list = []
    seq_list = []

    proc = subprocess.Popen("grep -c '>' " + filename, shell=True, stdout=subprocess.PIPE, )
    length = int(proc.communicate()[0].split('\n')[0])

    if no_of_seq <= 50:
        pos_num = int(filename.split("temp_fasta_")[1].split('.')[0]) - 1

    else:
        pos_num = 1

    with tqdm(range(length), desc="Formatting fasta..." + filename,
              position=pos_num, leave=False) as pbar:
        fasta_formatted = filename.split('.')[0] + '_temp.' + filename.split('.')[1]
        os.system('> ' + fasta_formatted)
        for record in SeqIO.parse(filename, "fasta"):
            if no_of_seq <= 50:
                pbar.update()
            acc_list.append(record.id)
            description_list.append(record.description)
            seq_list.append(record.seq)
            f = open(fasta_formatted, 'a')
            sys.stdout = f
            print ">" + record.description.replace(";", "").replace(",", "").replace(":", "").replace("'", "").replace(
                "(", "").replace(")", "").replace("[", "").replace("]", "").replace("-", "").replace("~", "")
            print re.sub("(.{60})", "\\1\n", str(record.seq), 0, re.DOTALL)
        if no_of_seq > 50:
            pbar.update()
    os.system('mv ' + fasta_formatted + ' ' + filename)

    sys.stdout = orig_stdout
    return_list = []
    return_list.append(acc_list)
    return_list.append(description_list)
    return_list.append(seq_list)
    return return_list

### Background process on previous searches ###

def initiate_old_search():
    global fetched_arrays
    sys.stdout = orig_stdout
    filename_rem_list = []

    found_old_seq_dict = {}

    if file_type == 'fasta':
        proc = subprocess.Popen("grep '>' " + filename + " | cut -c 2- | cut -d' ' -f1", shell=True,
                                stdout=subprocess.PIPE, )
        acc_2_look = (proc.communicate()[0].split('\n'))

        acc_2_look.remove('')

        res = [k for n, k in enumerate(acc_2_look) if k not in acc_2_look[:n]]
        acc_2_look = res

    dir_list = os.listdir(dir + '/storage/')
    if any("_parameters.txt" in string for string in dir_list):

        filenames = os.listdir(dir + '/storage/')

        try:
            filenames.remove('.*')

        except:
            pass

        for i in range(len(filenames)):
            var = "parameters.txt" in str(filenames[i])
            if var == True:
                filename_remain = filenames[i]
                filename_rem_list.append(filename_remain)

        filename_rem_list.sort(reverse=True)

        if len(filename_rem_list) != 0:
            for i in range(len(filename_rem_list)):
                with open(dir + '/storage/' + filename_rem_list[i], 'rU') as file:
                    temp_arr = []
                    CRISPR_out_old = filename_rem_list[i].split("_parameters.txt")[0]
                    for line in file:
                        if "#" not in line and len(line.split()) != 0 and "," not in line:
                            prev_word_length = int(line.split('\t')[1])
                            prev_word_repeatation = int(line.split('\t')[2])
                            prev_gap = int(line.split('\t')[3])
                            prev_repeat_length_cutoff = int(line.split('\t')[4])
                            prev_minimum_repeat_length = int(line.split('\t')[5])
                            prev_minimum_no_of_repeats = int(line.split('\t')[6])
                            prev_score = float(line.split('\t')[7])
                            prev_left_flank_length = int(line.split('\t')[8])
                            prev_right_flank_length = int(line.split('\t')[9])
                            prev_annotate_cas_genes = int(line.split('\t')[10])
                            prev_hmm_seq_evalue_cutoff = float(line.split('\t')[11])
                            prev_hmm_dom_evalue_cutoff = float(line.split('\t')[12])
                            prev_wgs = int(line.split('\t')[13])
                            prev_check_plasmid = int(line.split('\t')[14].split('\n')[0])
                        if "," in line and len(line.split()) != 0:
                            prev_seq = line.split('\n')[0].split(',')
                            time.sleep(1)

                    ### Check if these settings appeared exactly in previous runs. ###

                    if prev_word_length == word_length and prev_word_repeatation == word_repeatation and prev_gap == gap and prev_repeat_length_cutoff == repeat_length_cutoff and prev_minimum_repeat_length == minimum_repeat_length and prev_minimum_no_of_repeats == minimum_no_of_repeats and prev_score == score and prev_left_flank_length == left_flank_length and prev_right_flank_length == right_flank_length and prev_annotate_cas_genes == annotate_cas_genes and prev_hmm_seq_evalue_cutoff == hmm_seq_evalue_cutoff and prev_hmm_dom_evalue_cutoff == hmm_dom_evalue_cutoff and prev_wgs == wgs and prev_check_plasmid == check_plasmid:

                        for l in range(len(acc_2_look)):
                            if acc_2_look[l] in prev_seq:
                                if acc_2_look[l] not in found_old_seq:
                                    found_old_seq.append(acc_2_look[l])
                                    temp_arr.append(acc_2_look[l])
                            else:
                                found_new_seq.append(acc_2_look[l])
                if len(found_old_seq) != 0:
                    found_old_seq_dict[CRISPR_out_old] = temp_arr

        for key in found_old_seq_dict:
            fetched_arrays = (get_old_results(key, found_old_seq_dict[key]))
            fetched_arrays_accurate_return.extend(fetched_arrays[0])
            fetched_arrays_filtered_return.extend(fetched_arrays[1])
            fetched_gff_return.extend(fetched_arrays[2])
            fetched_spacer_return.extend(fetched_arrays[3])
            fetched_arrays_list.append(fetched_arrays)

    found_new_seq_list.append(found_new_seq)
    found_old_seq_list.append(found_old_seq)
    fetched_arrays_accurate_return_list.append(fetched_arrays_accurate_return)
    fetched_arrays_filtered_return_list.append(fetched_arrays_filtered_return)
    fetched_spacer_return_list.append(fetched_spacer_return)
    fetched_gff_return_list.append(fetched_gff_return)
#    pdb.set_trace()
def get_old_results(CRISPR_out_old, found_old_seq):
    current_array = []
    arrays = []
    fetched_arrays = []
    fetched_arrays_accurate = []
    fetched_filtered_arrays = []
    fetched_gff = []
    fetched_spacer = []

    with open(dir + '/storage/' + CRISPR_out_old, 'rU') as data:
        for line_out in data:
            current_array.append(line_out.split('\n')[0])
            if '//' in line_out:
                array_info = '\n'.join(current_array)
                arrays.append(array_info)
                current_array = []
    time.sleep(0.1)
    for i in range(len(found_old_seq)):
        for l in range(len(arrays)):
            if found_old_seq[i] in arrays[l]:
                fetched_arrays_accurate.append(arrays[l])
    time.sleep(0.1)
    current_array = []
    arrays = []
    with open(dir + '/storage/' + CRISPR_out_old + '.fp', 'rU') as data:
        for line_out in data:
            current_array.append(line_out.split('\n')[0])
            if '//' in line_out:
                array_info = '\n'.join(current_array)
                arrays.append(array_info)
                current_array = []

    time.sleep(0.1)
    for i in range(len(found_old_seq)):
        for l in range(len(arrays)):
            if found_old_seq[i] in arrays[l]:
                fetched_filtered_arrays.append(arrays[l])
    time.sleep(0.1)
    for i in range(len(found_old_seq)):
        with open(dir + '/storage/' + CRISPR_out_old + ".gff", 'rU') as data:
            for line_out in data:
                if found_old_seq[i] in line_out:
                    fetched_gff.append(line_out.split('\n')[0])
    time.sleep(0.1)

    for i in range(len(found_old_seq)):
        for record in SeqIO.parse(dir + '/storage/' + CRISPR_out_old + '.spacers.fa', "fasta"):
            if found_old_seq[i] in record.id:
                fetched_spacer.append(str('>' + record.description + '~' + record.seq))

    fetched_arrays.append(fetched_arrays_accurate)
    fetched_arrays.append(fetched_filtered_arrays)
    fetched_arrays.append(fetched_gff)
    fetched_arrays.append(fetched_spacer)

    time.sleep(0.1)

    return fetched_arrays


### Run CRISPRDetect ###

def CRISPRDetect_run(filename):
    global fetched_arrays, found_new_seq, found_old_seq, fetched_arrays_accurate_return, fetched_arrays_filtered_return, fetched_gff_return, fetched_spacer_return

    CRISPR_out = filename.split('.fasta')[0] + '_CRISPRDetect_' + timestr

    ### Return Parameters & compare previous runs & decide if running CRISPRDetect ###

    sys.stdout = orig_stdout

    pos_num = int(filename.split("temp_fasta_")[1].split('.')[0]) - 1

    new_seq_file = filename.split('.fasta')[0] + '_new_seq.fasta'
    params = (filename.split('.fasta')[0] + '_CRISPRDetect_' + timestr + "_parameters.txt").split('/')[-1]

    procs = []

    list_2_parse = fasta_formatter(filename)
    acc_list = list_2_parse[0]
    description_list = list_2_parse[1]
    seq_list = list_2_parse[2]
#    pdb.set_trace()
    #### Remove duplicate entries coming from old searches ####

    if len(fetched_arrays_list) != 0:
        fetched_arrays = fetched_arrays_list[0]
        found_new_seq = found_new_seq_list[0]
        res = [k for n, k in enumerate(found_new_seq) if k not in found_new_seq[:n]]
        found_new_seq = res
        found_old_seq = found_old_seq_list[0]
        res = [k for n, k in enumerate(found_old_seq) if k not in found_old_seq[:n]]
        found_old_seq = res
        for i in range(len(found_old_seq)):
            if found_old_seq[i] in found_new_seq:
                ind = found_new_seq.index(found_old_seq[i])
                found_new_seq.remove(found_new_seq[ind])
        fetched_arrays_accurate_return = fetched_arrays_accurate_return_list[0]
        res = [k for n, k in enumerate(fetched_arrays_accurate_return) if k not in fetched_arrays_accurate_return[:n]]
        fetched_arrays_accurate_return = res
        fetched_arrays_filtered_return = fetched_arrays_filtered_return_list[0]
        res = [k for n, k in enumerate(fetched_arrays_filtered_return) if k not in fetched_arrays_filtered_return[:n]]
        fetched_arrays_filtered_return = res
        fetched_gff_return = fetched_gff_return_list[0]
        res = [k for n, k in enumerate(fetched_gff_return) if k not in fetched_gff_return[:n]]
        fetched_gff_return = res
        fetched_spacer_return = fetched_spacer_return_list[0]
        res = [k for n, k in enumerate(fetched_spacer_return) if k not in fetched_spacer_return[:n]]
        fetched_spacer_return = res

    ### Execute depending of found new and ols sequences ###

    if (len(found_new_seq) != 0 and len(found_old_seq) == 0) or (len(found_new_seq) == 0 and len(found_old_seq) == 0):
        # pdb.set_trace()
        p1 = multiprocessing.Process(target=CRISPRDetect_cmd, args=[filename, CRISPR_out])
        procs.append(p1)
        p1.start()

    if (len(found_old_seq) != 0 and len(found_new_seq) == 0):
        os.system('> ' + CRISPR_out)
        f = open(CRISPR_out, 'a')
        sys.stdout = f
        hits=[]
        hits=[s for s in fetched_arrays_accurate_return if acc_list[0] in s]
        print "\n".join(hits)
        sys.stdout = orig_stdout

        os.system('> ' + CRISPR_out + '.fp')
        f = open(CRISPR_out + '.fp', 'a')
        sys.stdout = f
        hits=[]
        hits=[s for s in fetched_arrays_filtered_return if acc_list[0] in s]
        print "\n".join(hits)
        sys.stdout = orig_stdout

        os.system('> ' + CRISPR_out + '.gff')
        f = open(CRISPR_out + '.gff', 'a')
        sys.stdout = f
        hits=[]
        hits=[s for s in fetched_gff_return if acc_list[0] in s]
        print "\n".join(hits)
        sys.stdout = orig_stdout

        os.system('> ' + CRISPR_out + '.spacers.fa')
        f = open(CRISPR_out + '.spacers.fa', 'a')
        sys.stdout = f
        hits=[]
        hits=[s for s in fetched_spacer_return if acc_list[0] in s]
        res = list(map(lambda st: str.replace(st, "~", "\n"), hits))
        hits= res
        print "\n".join(hits)
        sys.stdout = orig_stdout

    if len(found_new_seq) != 0 and len(found_old_seq) != 0:

        CRISPR_out_found = CRISPR_out + '_found'
        f = open(CRISPR_out_found, 'a')
        sys.stdout = f
        hits=[]
        hits=[s for s in fetched_arrays_accurate_return if acc_list[0] in s]
        print "\n".join(hits)
        sys.stdout = orig_stdout

        CRISPR_out_found = CRISPR_out + '_found.fp'
        f = open(CRISPR_out_found, 'a')
        sys.stdout = f
        hits=[]
        hits=[s for s in fetched_arrays_filtered_return if acc_list[0] in s]
        print "\n".join(hits)
        sys.stdout = orig_stdout

        CRISPR_out_found = CRISPR_out + '_found.gff'
        f = open(CRISPR_out_found, 'a')
        sys.stdout = f
        hits=[]
        hits=[s for s in fetched_gff_return if acc_list[0] in s]
        print "\n".join(hits)
        sys.stdout = orig_stdout

        CRISPR_out_found = CRISPR_out + '_found.spacers.fa'
        f = open(CRISPR_out_found, 'a')
        sys.stdout = f
        hits=[]
        hits=[s for s in fetched_spacer_return if acc_list[0] in s]
        res = list(map(lambda st: str.replace(st, "~", "\n"), hits))
        hits= res
        print "\n".join(hits)
        sys.stdout = orig_stdout

        os.system('> ' + new_seq_file)
        f = open(new_seq_file, 'a')
        sys.stdout = f
        for k in range(len(found_new_seq)):
            ind = acc_list.index(found_new_seq[k])
            print ">" + description_list[ind]
            print re.sub("(.{60})", "\\1\n", str(seq_list[ind]), 0, re.DOTALL)
        sys.stdout = orig_stdout
        input_filename = filename
        filename = new_seq_file

        p1 = multiprocessing.Process(target=CRISPRDetect_cmd, args=[filename, CRISPR_out])
        procs.append(p1)
        p1.start()

    #### Prgogress bar for CRISPRDetect ####

    hasOpened = False

    for p1 in procs:
        p1.join(timeout=0)

        while p1.is_alive():
            path = '"' + path_2_CRISPR_CAS_search + '"'
            check_log = os.path.isfile(CRISPR_out + '.log')
            if check_log == True:
                try:
                    proc = subprocess.Popen("grep 'Total putative CRISPRs to process:' " + CRISPR_out + '.log',
                                            shell=True, stdout=subprocess.PIPE, )
                    count = int((proc.communicate()[0].split('\n')[0]).split(': ')[1].split()[0])
                    if count == '' and not hasOpened:
                        print 'Searching for putative CRISPR arrays.'
                    elif count != '' and not hasOpened:
                        break
                except:
                    pass
            else:
                pass

        prog_old = 0

        if no_of_seq <= 50:
            with tqdm(range(count), desc='Processing Putative CRISPRs...' + filename,
                      position=pos_num, leave=False) as pbar:
                while p1.is_alive():
                    try:
                        proc = subprocess.Popen("grep -c 'Remaining' " + CRISPR_out + '.log', shell=True,
                                                stdout=subprocess.PIPE, )
                        prog = int(proc.communicate()[0].split('\n')[0])
                        if prog - prog_old == 0:
                            proc = subprocess.Popen("grep -c 'Total putative CRISPRs to process: 0 from 0 sequences' " + CRISPR_out + '.log', shell=True,
                                                    stdout=subprocess.PIPE, )
                            prog_fail = int(proc.communicate()[0].split('\n')[0])
                            if prog_fail == 1:
                                pbar.update()
                            else:
                                pass
                        else:
                            prog_old = prog
                            pbar.update()
                    except:
                        pass

    parameter_log(params, acc_list)
    CRISPRDetect_summary(CRISPR_out)

    ##### Tidy after CRISPRDetect ####

    check_matched = os.path.isfile(CRISPR_out + '_found')

    if check_matched == True:
        time.sleep(1)
        os.system('cat ' + CRISPR_out + '_found >>' + CRISPR_out)
        os.system('rm ' + CRISPR_out + "_found ")

        os.system('cat ' + CRISPR_out + '_found.fp >>' + CRISPR_out + ".fp")
        os.system('rm ' + CRISPR_out + "_found.fp ")

        os.system('cat ' + CRISPR_out + '_found.gff >>' + CRISPR_out + ".gff")
        os.system('rm ' + CRISPR_out + "_found.gff ")

        os.system('cat ' + CRISPR_out + '_found.spacers.fa >>' + CRISPR_out + ".spacers.fa")
        os.system('rm ' + CRISPR_out + "_found.spacers.fa ")

        check_newfile = os.path.isfile(new_seq_file)

        if check_newfile == True:
            os.system('rm ' + new_seq_file)



### CRISPRDetect Command ###

def CRISPRDetect_cmd(filename, CRISPR_out):
    time.sleep(1)
    os.system('perl ' + path_2_CRISPRDetect + ' -f ' + filename + ' -o ' + CRISPR_out + ' -word_length ' + str(word_length) + ' -minimum_word_repeatation ' + str(word_repeatation) + ' -max_gap_between_crisprs ' + str(gap) + ' -repeat_length_cutoff ' + str(repeat_length_cutoff) + ' -minimum_repeat_length ' + str(minimum_repeat_length) + ' -minimum_no_of_repeats ' + str(minimum_no_of_repeats) + ' -check_direction 1 -array_quality_score_cutoff ' + str(score) + ' -T ' + str(thread_in_use) + ' -annotate_cas_genes ' + str(annotate_cas_genes) + ' -hmm_seq_evalue_cutoff ' + str(hmm_seq_evalue_cutoff) + ' -hmm_dom_evalue_cutoff ' + str(hmm_dom_evalue_cutoff) + ' -wgs ' + str(wgs) + ' -check_plasmid ' + str(check_plasmid) + ' -tmp_dir ' +path_2_CRISPR_CAS_search.split('/CRISPR_CAS_search.py')[0] + '/' + temp_dir + '/ -ref_lib_file /home/muratb/Desktop/Local_Softwares/CRISPRDetect_3.0/Ref_lib_files/verified_repeats_with_family.txt > ' + CRISPR_out + '.log')

### Generate CRISPRDetect summary ###

def CRISPRDetect_summary(CRISPR_out):

    CRISPRDetect_summary = CRISPR_out + "_summary.txt"
#    print(CRISPR_out)
    os.system('> ' + CRISPRDetect_summary)
    f = open(CRISPRDetect_summary, 'a')
    sys.stdout = f

### need double check this part ###
    if "temp_fasta_" in CRISPR_out:
        if no_of_seq <= 50:
            proc = subprocess.Popen("grep -c '>' " + CRISPR_out, shell=True, stdout=subprocess.PIPE, )
            length = int(proc.communicate()[0].split('\n')[0])
            pos_num = int(CRISPR_out.split("temp_fasta_")[1].split('_CRISPRDetect')[0]) - 1
        else:
            length = no_of_seq
            pos_num = 2
    else:
        length = no_of_seq
        pos_num = 2
##############END#########

#    if no_of_seq <= 50:
#        proc = subprocess.Popen("grep -c '>' " + CRISPR_out, shell=True, stdout=subprocess.PIPE, )
#        length = int(proc.communicate()[0].split('\n')[0])
#        pos_num = int(CRISPR_out.split("temp_fasta_")[1].split('_CRISPRDetect')[0]) - 1
#    else:
#        length = no_of_seq
#        pos_num = 2

    print "Array_no\tArray_type\tAccession\tName\tStart\tStop\tStrand\tSubtpye\tRepeat_occurence\tRepeat_length\tRepeat_seq\tScore"

    with tqdm(range(length), desc='Creating summary file...' + CRISPR_out.split('_CRISPRDetect')[0],
              position=pos_num, leave=False) as pbar_sum:
        with open(CRISPR_out, 'rU') as file:
            for line in file:
                line_arr = line.split()
                if (len(line_arr) != 0):
                    if ">" in line:
                        if no_of_seq <= 50:
                            pbar_sum.update()
                    if line_arr[0] == "Array":
                        array_no = line_arr[1]
                    if "# Array family :" in line:
                        subtype = line.split("# Array family : ")[1].split('\n')[0]
                    if "# Summary:" in line:
                        arr = line.split("ID_START_STOP_DIR: ")[1].split(';')
                        name = arr[0]
                        fullname = name.split("-")[0]
                        if file_type == "fasta":
                            acc = name.split()[0]
                        low_bound = name.split("-")[-3]
                        high_bound = name.split("-")[-2]
                        ori = name.split("-")[-1]
                        dr_seq = arr[1].split(":")[1]
                        dr_len = arr[2].split(":")[1]
                        rep_occ = arr[3].split(":")[1]
                        score = arr[-2].split(":")[1]
                        print acc + '_Array_' + array_no + '\tCRISPR\t' + acc + '\t' + fullname + '\t' + low_bound + '\t' + high_bound + '\t' + ori + '\t' + subtype + '\t' + rep_occ + '\t' + dr_len + '\t' + dr_seq + '\t' + score
            if no_of_seq > 50:
                pbar_sum.update()

    sys.stdout = orig_stdout

### Generate parameter log ###

def parameter_log(params, acc_list):
    os.system('> ' + params)
    f = open(params, 'a')
    sys.stdout = f
    print '#filename\tword_length\tword_repeatation\tgap\trepeat_length_cutoff\tminimum_repeat_length\tminimum_no_of_repeats\tscore\tleft_flank_length\tright_flank_length\tannotate_cas_genes\thmm_seq_evalue_cutoff\thmm_dom_evalue_cutoff\twgs\tcheck_plasmid'
    print filename + '\t' + str(word_length) + '\t' + str(word_repeatation) + '\t' + str(
        gap) + '\t' + str(repeat_length_cutoff) + '\t' + str(
        minimum_repeat_length) + '\t' + str(minimum_no_of_repeats) + '\t' + str(
        score) + '\t' + str(left_flank_length) + '\t' + str(
        right_flank_length) + '\t' + str(annotate_cas_genes) + '\t' + str(
        hmm_seq_evalue_cutoff) + '\t' + str(hmm_dom_evalue_cutoff) + '\t' + str(
        wgs) + '\t' + str(check_plasmid)
    print (',' + ','.join(map(str, acc_list)))
    sys.stdout = orig_stdout
    # print ('\n'.join(map(str, acc_list)))

### Execute CRISPRDetect run ###

def post_CRISPRDetect():
    ### clean after CRISPRDetection ###

    os.system('find . -maxdepth 1 -type f -name "*' + timestr + '" -print0 | sort -zV | xargs -0 cat > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr)
    os.system('find . -maxdepth 1 -type f -name "*' + timestr + '.fp" -print0 | sort -zV | xargs -0 cat > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.fp')
    os.system('find . -maxdepth 1 -type f -name "*' + timestr + '.gff" -print0 | sort -zV | xargs -0 cat > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.gff')
    os.system('find . -maxdepth 1 -type f -name "*' + timestr + '.spacers.fa" -print0 | sort -zV | xargs -0 cat > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.spacers.fa')
    os.system('find . -maxdepth 1 -type f -name "*' + timestr + '.log" -print0 | sort -zV | xargs -0 cat > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.log')

    final_summary = filename.split('.')[0] + '_CRISPRDetect_' + timestr + '_parameters.txt'
    os.system('> ../' + final_summary)
    f = open("../" + final_summary, 'a')
    sys.stdout = f
    print '#filename\tword_length\tword_repeatation\tgap\trepeat_length_cutoff\tminimum_repeat_length\tminimum_no_of_repeats\tscore\tleft_flank_length\tright_flank_length\tannotate_cas_genes\thmm_seq_evalue_cutoff\thmm_dom_evalue_cutoff\twgs\tcheck_plasmid'
    print filename + '\t' + str(word_length) + '\t' + str(word_repeatation) + '\t' + str(
        gap) + '\t' + str(repeat_length_cutoff) + '\t' + str(
        minimum_repeat_length) + '\t' + str(minimum_no_of_repeats) + '\t' + str(
        score) + '\t' + str(left_flank_length) + '\t' + str(
        right_flank_length) + '\t' + str(annotate_cas_genes) + '\t' + str(
        hmm_seq_evalue_cutoff) + '\t' + str(hmm_dom_evalue_cutoff) + '\t' + str(
        wgs) + '\t' + str(check_plasmid)
    f.close()
    sys.stdout = orig_stdout

    os.system('find . -maxdepth 1 -type f -name "*' + timestr + '_parameters.txt" -print0 | sort -zV | xargs -0 tail -n 1 -q | awk "{print}" ORS="" >> ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '_parameters.txt')



    # os.system('cat *_CRISPRDetect_' + timestr + ' > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr)
    # os.system('cat *_CRISPRDetect_' + timestr + '.fp > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.fp')
    # os.system('cat *_CRISPRDetect_' + timestr + '.gff > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.gff')
    # os.system('cat *_CRISPRDetect_' + timestr + '.spacers.fa > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.spacers.fa')
    # os.system('cat *_CRISPRDetect_' + timestr + '_summary.txt > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '_summary.txt')
    # os.system('cat *_CRISPRDetect_' + timestr + '.log > ../' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.log')

    # filenames = os.listdir(work_dir+"/temp_"+timestr)
    # prev_acc_list = []
    # final_summary = filename.split('.')[0] + '_CRISPRDetect_' + timestr + '_parameters.txt'
    # os.system('> ../' + final_summary)
    # f = open("../"+final_summary, 'a')
    # sys.stdout = f
    # for i in range(len(filenames)):
    #     if "parameters.txt" in filenames[i] and timestr in filenames[i]:
    #         with open(filenames[i], 'rU') as file:
    #             for line in file:
    #                 if "#" not in line and len(line.split()) != 0 and "," not in line:
    #                     parameter = line.split('\n')[0]
    #                 if ',' in line and len(line.split()) != 0:
    #                     prev_acc = line.split('\n')[0].split(',')
    #                     prev_acc.remove('')
    #                     prev_acc_list.extend(prev_acc)
    # print '#filename\tword_length\tword_repeatation\tgap\trepeat_length_cutoff\tminimum_repeat_length\tminimum_no_of_repeats\tscore\tleft_flank_length\tright_flank_length\tannotate_cas_genes\thmm_seq_evalue_cutoff\thmm_dom_evalue_cutoff\twgs\tcheck_plasmid'
    # print parameter
    # print ',' + ','.join(prev_acc_list)
    # time.sleep(1)
    # f.close()
    # sys.stdout = orig_stdout
    os.chdir("../")

    CRISPRDetect_summary(filename.split('.')[0] + '_CRISPRDetect_' + timestr)

    time.sleep(0.5)
    os.system('cp ' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + ' ' + dir + '/storage/')
    os.system('cp ' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.fp ' + dir + '/storage/')
    os.system('cp ' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.gff ' + dir + '/storage/')
    os.system('cp ' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '.spacers.fa ' + dir + '/storage/')
    os.system('cp ' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '_summary.txt ' + dir + '/storage/')
    os.system('cp ' + filename.split('.')[0] + '_CRISPRDetect_' + timestr + '_parameters.txt ' + dir + '/storage/')

if __name__ == "__main__":

    found_new_seq = []
    found_old_seq = []
    fetched_arrays_accurate_return = []
    fetched_arrays_filtered_return = []
    fetched_gff_return = []
    fetched_spacer_return = []

    print "Starting environment setup!"

    env = environment_setup()
    dir = env[0]
    work_dir = env[1]

    print "Environment setup completed!"

    #initiate storage recall#

#    p_init = multiprocessing.Process(target=initiate_old_search)
#    p_init.start()
    remove_dub_acc(filename)
    no_of_seq = split_(filename)

    #wait until the storage is processed#

#    hasStarted = 0
#    while hasStarted == 0:
#        if p_init.is_alive():
#            pass
#        else:
#            hasStarted = 1
#            break

    #initiate CRISPRDetect#

    execute(filename)
    print "\nArranging..."
    post_CRISPRDetect()

    print "\nCleaning temporary files..."

    os.system("rm -r " + working)

    os.chdir(dir)
    os.system("rm -rf " + temp_dir)

    print "\nFinished!"

