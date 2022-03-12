#! /usr/bin/env python3
# written by Matthias Schmal
# distributed under GNU GPLv3
# merge individual output files of htseq-count into one tsv, suitable to import into edgeR

# import required libraries and modules
import argparse
import csv

# generate parser to read in different arguments
parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='give htseq-count output files')
args = parser.parse_args()

combined = [["gene"], [args.files[0]]]
tmp = []
with open(args.files[0]) as file:
    for line in file:
        combined[0].append(line.split()[0])
        combined[1].append(line.split()[1])

    
for file_name in args.files[1:]:
    with open(file_name) as file:
        tmp.append(file_name)
        for line in file:
            tmp.append(line.split()[1])
    combined.append(tmp)
    tmp=[]

# transpose combined
combined_trans = [*zip(*combined[:])]


with open("count_combined.tsv", "w+") as file:
    csv_file = csv.writer(file, delimiter="\t")
    csv_file.writerows(combined_trans)


