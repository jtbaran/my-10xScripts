# ########################################
# J. Baran-Gale
# University of Edinburgh
# V 1.0
# 15 Oct 2016
# Input: 10x Single cell barcode files (read-I1*.gz)
# Output: A density scatter plot detailing the overall quality of the reads
# Assumes all files are in a directory by sample name, the directory name is
# stripped from the input file path, and used as a prefix for the plot.
# ########################################

import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip, time, gc
from Bio import SeqIO
import os
import sys


print("BARCODE FILES:")
brc_dirs=[]
for i in range(1,len(sys.argv)):
    brc_dirs+=[sys.argv[i]]
    print(brc_dirs[i-1])
random.seed()
baseName=os.path.dirname(brc_dirs[0])

from collections import defaultdict

def encoding_map(ch):
    if ch=='A':return 0
    if ch=='G':return 1
    if ch=='C':return 2
    if ch=='T':return 3
    if ch=='N':return 0 # random.randint(0,3)

decoding_lst = ['A', 'G', 'C', 'T']

def encode(k):
    code = 0
    for ch in k:
        code *= 4
        code += encoding_map(ch)
    return code

######################################################


Qual=defaultdict(list)
for brc_dir in brc_dirs :
	with gzip.open(brc_dir) as f:
		for record in SeqIO.parse(f,"fastq"):
            		Qual[encode(record.seq)].append(np.mean(record.letter_annotations["phred_quality"]))

x=np.ones(len(Qual.keys()))
y=np.ones(len(Qual.keys()))
i=0
for k,v in Qual.items():
	x[i]=np.mean(v)
	y[i]=len(v)
	i+=1


np.random.seed(0)
xmin = 14
xmax = x.max()
ymin = y.min()
ymax = y.max()


plt.subplot()
plt.hexbin(x, y, yscale='log',bins='log', cmap=plt.get_cmap("plasma"))
plt.axis([xmin, xmax, ymin, ymax])
plt.ylabel('Cell read count')
plt.xlabel('Average read quality')
plt.title("Barcode quality")
cb = plt.colorbar()
cb.set_label('log10(N)')

plt.savefig(baseName +'qualDensity.pdf', bbox_inches='tight')
