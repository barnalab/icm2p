#!/usr/bin/python
# This script generates M2-seq correlated mutation counts matrix from ShapeMapper2 parsed mutation file

from Bio import SeqIO
from Bio.Seq import Seq
import os
import argparse
import numpy as np
from numpy import *
import string
from rdatkit.handler import RDATFile

parser = argparse.ArgumentParser( description='Generates M2-seq correlated mutation counts matrix')

parser.add_argument('--ref', type=str, help='FASTA reference')
parser.add_argument('--refname', type=str, help='FASTA record name')
parser.add_argument('--mutParse', type=argparse.FileType('r'), help='ShapeMapper2 mutation parser output file')
parser.add_argument('--outprefix',  type=str,  default='out', help='Output prefix')
parser.add_argument('--minlen', type=int, default=150, help='Minimum alignment length')
parser.add_argument('--minmut', type=int, default=4, help='Minimum number of mutations')
args = parser.parse_args()

refseq = str(SeqIO.to_dict(SeqIO.parse(args.ref,"fasta"))[args.refname].seq).upper()#.replace("U","T")

WTlen = len(refseq)
seqpos = np.arange( 1, WTlen+1 )

#Numpy arrays
data2d = np.zeros(([len(seqpos),len(seqpos)]),dtype=int32) #2D matrices; A/T/G/C separately
data2dA = np.zeros(([len(seqpos),len(seqpos)]),dtype=int32)
data2dT = np.zeros(([len(seqpos),len(seqpos)]),dtype=int32)
data2dG = np.zeros(([len(seqpos),len(seqpos)]),dtype=int32)
data2dC = np.zeros(([len(seqpos),len(seqpos)]),dtype=int32)

#Read ShapeMapper2 mutation parser output file and count correlated mutations
for line in args.mutParse:
    fields = line.strip().replace("\"","").split("\t")
    rname = fields[1] #Read name
    rstart = int(fields[2]) #Read start pos
    rend = int(fields[3]) #Read end pos
    if len(fields) < 10:
        continue
    if (rstart > 50):
        continue
    if (rend < (WTlen-50):
        continue
    if ( (rend-rstart+1) < args.minlen):
        continue
    muts = fields[9].split(' ')
    muts = [ muts[ (i*5):((i+1)*5) ] for i in range((len(muts)+5-1)//5 ) ]

    newseq = ['0']*len(refseq)
    subcount = 0 #Substitution count
    for y in muts:
        ref_start = int(y[0]) + 1
        ref_end = int(y[1])
        y_len = ref_end - ref_start
        y_seqlen = len(y[2])
        if ( y_len == y_seqlen ): #Substitutions only
            newseq[ref_start:ref_end] = y[4][1:(y_seqlen+1)]
            subcount += 1
    if(subcount>=args.minmut): #More than N substitutions
        start_pos = rstart + 1
        end_pos = rend + 1 

        useseq = newseq[rstart:(rend+1)]

        #Mutation indices
        mut_idx_init = [ (idx + start_pos - 1, val) for idx, val in enumerate(useseq) if val != "0"]
        mut_idx = [ pp for pp in mut_idx_init if int(pp[0]) < WTlen]
        mut_idx = transpose(array(mut_idx))
        mut_val = mut_idx[1,]
        mut_idx = array(mut_idx[0,], dtype=int32)

        #Count up matrix
        for ii in range(0,len(mut_idx)):
            data2d[mut_idx[ii], mut_idx] += 1
            if ( mut_val[ii] == "A"):
                data2dA[mut_idx[ii], mut_idx] += 1
            elif ( mut_val[ii] == "T"):
                data2dT[mut_idx[ii], mut_idx] += 1
            elif ( mut_val[ii] == "G"):
                data2dG[mut_idx[ii], mut_idx] += 1
            elif ( mut_val[ii] == "C"):
                data2dC[mut_idx[ii], mut_idx] += 1

np.savetxt(args.outprefix+'.raw',data2d,delimiter="\t",fmt='%i')
np.savetxt(args.outprefix+'.rawA',data2dA,delimiter="\t",fmt='%i')
np.savetxt(args.outprefix+'.rawT',data2dT,delimiter="\t",fmt='%i')
np.savetxt(args.outprefix+'.rawG',data2dG,delimiter="\t",fmt='%i')
np.savetxt(args.outprefix+'.rawC',data2dC,delimiter="\t",fmt='%i')
