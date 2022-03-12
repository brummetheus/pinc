#! /usr/bin/env python3
# written by Matthias Schmal
# distributed under GNU GPLv3
# trains and evaluates a CPAT model

# import required libraries and modules
import argparse
import subprocess
import custom_functions as cf

from sklearn.model_selection import KFold

# generate parser to read in different arguments
parser = argparse.ArgumentParser()
parser.add_argument('--coding', '-c', help='give sequence file in fasta format for coding sequences')
parser.add_argument('--noncoding', '-n', help='give sequence file in fasta format for noncoding sequences')
parser.add_argument('--cpc-coding', '-cpc', help='give sequence file in fasta format for cpc test set')
parser.add_argument('--test-coding', '-ct', help='give sequence file in fasta format for test sequences')
parser.add_argument('--test-noncoding', '-nt', help='give sequence file in fasta format for test sequences')
parser.add_argument('--cod-cod-ids', '-cc', help='give file with ids of predicted coding transcripts of coding RNAs')
parser.add_argument('--non-cod-ids', '-nc', help='give file with ids of predicted noncoding transcripts of coding RNAs')
parser.add_argument('--cod-non-ids', '-cn', help='give file with ids of predicted coding transcripts of noncoding RNAs')
parser.add_argument('--non-non-ids', '-nn', help='give file with ids of predicted noncoding transcripts of noncoding RNAs')
parser.add_argument('--outname', '-o', help='give prefix for filenames generated')
parser.add_argument('--weight', '-w', help='give weight for false-positive rate')
args = parser.parse_args()

# read in data from arguments
train_cod = cf.fasta_reader(args.coding)
train_cpc_cod = cf.fasta_reader(args.cpc_coding)
train_noncod = cf.fasta_reader(args.noncoding)

# initialize arrays for later use
cc_ids = []
nc_ids = []
cn_ids = []
nn_ids = []
cutoff = []

# read in trnascript IDs from CPC2 output
with open(args.cod_cod_ids) as file:
    for line in file:
        cc_ids.append(">"+line.strip())

with open(args.non_cod_ids) as file:
    for line in file:
        nc_ids.append(">"+line.strip())

with open(args.cod_non_ids) as file:
    for line in file:
        cn_ids.append(">"+line.strip())

with open(args.non_non_ids) as file:
    for line in file:
        nn_ids.append(">"+line.strip())

# generate tuples of fasta header and sequence
cc = [(x,train_cpc_cod[x]) for x in cc_ids]
nc = [(x,train_cpc_cod[x]) for x in nc_ids]
cn = [(x,train_noncod[x]) for x in cn_ids]
nn = [(x,train_noncod[x]) for x in nn_ids]
train_cod = list(train_cod.items())

# merge nn and nc to generate one dataset of predicted noncoding RNAs 
nn += nc

# write files for final training set
cf.list_fasta_writer(cc, args.outname + "_cc.fna")
cf.list_fasta_writer(nc, args.outname + "_nc.fna")
cf.list_fasta_writer(cn, args.outname + "_cn.fna")
cf.list_fasta_writer(nn, args.outname + "_nn.fna")


# create balanced training set
if (len(nn) >= len(train_cod)):
    nn = nn[:len(train_cod)]
else:
    train_cod = train_cod[:len(nn)]

# perform 10-fold cross-validation to determine optimal cutoff according to weighted Youden's Index
direc = args.outname+"_roc_curves"
subprocess.run(["mkdir", direc])
iter = 1
kf = KFold(n_splits=10, shuffle=True)
# train and evaluate for each split
for train, test in kf.split(nn):
    train_noncod_trans = [nn[i] for i in train]
    train_cod_trans = [train_cod[i] for i in train]
    test_noncod_trans = [nn[i] for i in test]
    test_cod_trans = [train_cod[i] for i in test]
    cf.list_fasta_writer(train_noncod_trans, "train_noncod_trans.fasta")
    cf.list_fasta_writer(test_noncod_trans, "test_noncod_trans.fasta")  
    cf.list_fasta_writer(train_cod_trans, "train_cod_trans.fasta")
    cf.list_fasta_writer(test_cod_trans, "test_cod_trans.fasta")

    with open("hexamer.tsv", "w+") as file:
        subprocess.run(["make_hexamer_tab.py", "-c", "train_cod_trans.fasta", "-n", "train_noncod_trans.fasta"], stdout=file)
        
    subprocess.run(["make_logitModel.py", "-x", "hexamer.tsv", "-c", "train_cod_trans.fasta", "-n", "train_noncod_trans.fasta", "-o", "train"])

    subprocess.run(["cpat.py", "-x", "hexamer.tsv", "-d", "train.logit.RData", "--top-orf=5", "--min-orf=24", "--antisense", "-g", "test_cod_trans.fasta", "-o", "cod_genes"])
    subprocess.run(["cpat.py", "-x", "hexamer.tsv", "-d", "train.logit.RData", "--top-orf=5", "--min-orf=24", "--antisense", "-g", "test_noncod_trans.fasta", "-o", "noncod_genes"])
    
    cutoff.append(cf.find_cutoff("cod_genes.ORF_prob.best.tsv", "noncod_genes.ORF_prob.best.tsv", float(args.weight), iter, direc))

    iter += 1

# write file with calculated cutoff as well as ROC auc-score for each iteration
with open(args.outname+"_cutoff_data.tsv", 'w+') as file:
    file.write("number\tcutoff\tauc_score\n")
    for i, split in enumerate(cutoff):
        file.write("Iteration"+str(i)+"\t"+str(split[0]) +"\t" + str(split[1]) + "\n")
    file.write("mean\t" + str(sum([x[0] for x in cutoff])/len(cutoff)) + "\t"+ str(sum([x[1] for x in cutoff])/len(cutoff)))

# do final training with complete dataset 
with open(args.outname + "_hexamer.tsv", "w+") as file:
    subprocess.run(["make_hexamer_tab.py", "-c", args.coding, "-n", args.noncoding], stdout=file)

subprocess.run(["make_logitModel.py", "-x", args.outname + "_hexamer.tsv", "-c", args.coding, "-n", args.noncoding, "-o", args.outname])

subprocess.run(["cpat.py", "-x", args.outname + "_hexamer.tsv", "-d", args.outname + ".logit.RData", "--top-orf=5", "--min-orf=0", "--antisense", "-g", args.test_coding, "-o", args.outname+"_coding"])
subprocess.run(["cpat.py", "-x", args.outname + "_hexamer.tsv", "-d", args.outname + ".logit.RData", "--top-orf=5", "--min-orf=0", "--antisense", "-g", args.test_noncoding, "-o", args.outname+"_noncoding"])
