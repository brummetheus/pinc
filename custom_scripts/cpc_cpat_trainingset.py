#! /usr/bin/env python3
# written by Matthias Schmal
# distributed under GNU GPLv3
# split data into training and test data in ratio of 4:1

# import required libraries and modules
import random
import argparse
import custom_functions as cf

# generate parser to read in different arguments
parser = argparse.ArgumentParser()
parser.add_argument('--coding', '-c', help='give sequence file in fasta format for coding sequences, \
                    wildcards possible', nargs='*')
parser.add_argument('--noncoding', '-n', help='give sequence file in fasta format for noncoding sequences, \
                    wildcards possible', nargs='*')
parser.add_argument('--split', '-s', help='give precentage to split data into training and test set', \
                    type=int)
parser.add_argument('--outname', '-o', help='give prefix for outputfiles')
args = parser.parse_args()

# intitialize different variables for later use
cod = {}
noncod = {}
train_cod = ()
train_noncod = ()
test_cod = ()
test_noncod = ()
not_allowed = ['N', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V']

# read in fasta files 
for filename in args.coding:
    cod.update(cf.fasta_reader(filename))
for filename in args.noncoding:
    noncod.update(cf.fasta_reader(filename))

# convert possible RNA sequences in DNA and remove every sequence with non-standard bases in it
cod = {k.split()[0]: v.replace('U', 'T') if 'U' in k else v for k,v in cod.items() \
       if all(n not in v for n in not_allowed)}
noncod = {k.split()[0]: v.replace('U', 'T') if 'U' in k else v for k,v in noncod.items() \
          if all(n not in v for n in not_allowed)}

# convert dictionaries in list for subsequent steps
cod = list(cod.items())
noncod = list(noncod.items())

# randomize sequences based on a given seed to ensure reproducability
random.seed(31415926)
random.shuffle(cod)
random.shuffle(noncod)


# generate training and test dataset
train_noncod = noncod[:int(args.split/100*len(noncod))]
train_cod = cod[:int(args.split/100*len(cod))]

test_cod = cod[len(train_cod):]
test_noncod = noncod[len(train_noncod):]

# split coding RNAs training set to have a coding RNA dataset with 25% the size of noncoding training set
train_cod_cpc = train_cod[:int(len(train_noncod)/4)]
train_cod = train_cod[len(train_cod_cpc):]


# train_noncod will be used for CPC2
cf.fasta_writer(dict(train_cod), args.outname+'_training_cod.fasta')
cf.fasta_writer(dict(train_cod_cpc), args.outname+'_training_cod_cpc.fasta')
cf.fasta_writer(dict(train_noncod), args.outname+'_training_noncod.fasta')

cf.fasta_writer(dict(test_cod), args.outname+'_test_cod.fasta')
cf.fasta_writer(dict(test_noncod), args.outname+'_test_noncod.fasta')
