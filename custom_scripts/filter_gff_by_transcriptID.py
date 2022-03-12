#! /usr/bin/env python3
# written by Matthias Schmal
# distributed under GNU GPLv3
# filter a gff3 file based on transcript IDs

# import required libraries and modules
import argparse

# generate parser to read in different arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='give gtf file for filtering')
parser.add_argument('--filter', '-f', help='give list of filter transcript IDs')
parser.add_argument('--output', '-o', help='give name of output file')
args = parser.parse_args()

# initialize arrays for later use
hits = []
input_header = []
output = []

# read in transcript IDs
with open(args.filter, 'r') as file:
    for line in file:
        hits.append(line.strip())

# filter gff3
with open(args.input, 'r') as file:
    tmp = hits.pop(0)
    last_hit = ""
    for line in file:
        if hits:
            if line[0] == '#':
                input_header.append(line.strip())
            else:
                if tmp in line.split('\t')[8].split('; ')[1]:
                    output.append(line)
                else:
                    if hits[0] in line.split('\t')[8].split('; ')[1]:
                        output.append(line)
                        tmp = hits.pop(0)
        else: 
            if tmp in line.split('\t')[8].split('; ')[1]:
                output.append(line)
            else:
                break

# write resulting gff3 to disk
with open(args.output, 'w+') as file:
    for line in input_header:
        file.write(line+'\n')
    for line in output:
        file.write(line)


