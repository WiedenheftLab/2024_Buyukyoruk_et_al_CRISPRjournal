#!/usr/bin/python

import os

cwd = os.getcwd()

os.system('mkdir originals')

os.system('python ' + cwd + '/CRISPR_search_on_leader.py')
os.system('python ' + cwd + '/auto_divider.py')
os.system('python ' + cwd + '/auto_repeat_remover.py')
os.system('python ' + cwd + '/auto_cluster.py')
os.system('python ' + cwd + '/auto_meme.py')

