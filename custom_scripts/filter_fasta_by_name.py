#! /usr/bin/env python3
# written by Matthias Schmal
# distributed under GNU GPLv3
# filter a fasta file based on fasta header given as filter

# import required libraries and modules
import argparse
import custom_functions as cf

# generate parser to read in different arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='give fasta file for to be filtered')
parser.add_argument('--filter', '-f', help='give list of sequence identifier; they need to be exact matches')
parser.add_argument('--output', '-o', help='give name of output file')
args = parser.parse_args()

fasta = cf.fasta_reader(args.input)
output = {}

# filter input fasta file based on fasta header in filter
with open(args.filter, 'r+') as file:
    for hit in file:
        output[">"+hit.strip()] = fasta[">"+hit.strip()]

cf.fasta_writer(output, args.output)
