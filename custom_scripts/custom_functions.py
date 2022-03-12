import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt
# from bisect import bisect_left, bisect_right
# from numpy import log


# reads a fasta file with multiple entries and returns a dictionary with the
# whole identifier as key and the sequence as value.
def fasta_reader(filename):
    with open(filename, 'r+') as file:
        fasta = {}
        key = None
        value = ""
        for line in file:
            if line[0] == ">":
                # just for first line
                if len(value) == 0:
                    key = line.strip()
                    continue
                else:
                    fasta[key] = value
                    key = line.strip()
                    value = ""

            elif line[0] == "\n":
                continue
            else:
                value = value+line.strip()
        fasta[key] = value
        return fasta


def fasta_writer(fasta_dict, filename):
    with open(filename, 'w+') as file:
        for k,v in fasta_dict.items():
            file.write(k+"\n")
            file.write(v+"\n")


# read "genes" from gff3 file derived from bp_genbank2gff3.pl script
# if strange_ID is TRUE: ID will be replaced with gene1, gene2, etc
def gff_gene_reader(filename, strange_ID):
    with open(filename, 'r+') as file:
        header = ""
        entries = []
        gene_counter = 1
        for line in file:
            if line[0] == '\n':
                continue
            if line[0] == '#':
                header += line
            else:
                tmp = line.split()
                if tmp[2] == 'gene':
                    if strange_ID:
                        id = tmp[8].split(';')
                        id[0]="ID=gene"+str(gene_counter)
                        tmp[8] = ';'.join(id)
                        gene_counter += 1
                    entries.append(tmp)
    return entries


# read in CPAT output and add "true label" column, true labels are given by 
# cpat_trainingset.py

def cpat_reader(filename, label):
    with open(filename, "r+") as file:
        # skip header
        file.readline()
        data = []
        for line in file:
            tmp = line.strip().split("\t")
            tmp[10] = float(tmp[10])
            tmp.append(int(label))
            data.append(tmp)
    return data


def list_fasta_writer(l, filename):
    with open(filename, "w+") as file:
        for entry in l:
            file.write(entry[0] + "\n" +
                       entry[1] + "\n")
    return


# find "best" cutoff for prediction of ncRNAs based on Youden's Index.
def find_cutoff(file_cod, file_noncod, w, n, d):

    cpat = cpat_reader(file_cod, 1) + cpat_reader(file_noncod, 0)
     
    # sort data according to coding probability
    cpat.sort(key=lambda entry: entry[10])

    # "transpose" data matrix
    cpat_trans = []
    for col in range(len(cpat[0])):
        cpat_trans.append([i[col] for i in cpat])
    
    fpr, tpr, thresh = roc_curve(cpat_trans[11], cpat_trans[10], drop_intermediate=False, pos_label=1)
    opt = np.argmax(2*(w*tpr-(1-w)*fpr))
    cutoff = (thresh[opt] + thresh[opt+1])/2
    auc = roc_auc_score(cpat_trans[11], cpat_trans[10])
    # plot roc_curve and cutoff vs performance 
    lw = 1.5
    fig, (ax1,ax2) = plt.subplots(1,2)
    fig.set_size_inches(14,7)
    fig.suptitle("Iteration " + str(n))
    ax1.plot(fpr, tpr, color="darkorange", lw=lw)
    ax1.plot([0,1], [0,1], color="navy", lw=lw, linestyle="--", label="ROC curve (AUC = %.03f" %auc)
    ax1.set_title("ROC curve")
    ax1.set_xlabel("False positive rate")
    ax1.set_ylabel("True positive rate")
    ax1.legend(loc="lower right")
    ax2.plot(thresh, tpr , color="darkorange", lw=lw, label="Sensitivity")
    ax2.plot(thresh, 1-fpr, color="navy", lw=lw, label="Specificity")
    ax2.set_title("Cutoff vs Performance")
    ax2.legend(loc="lower right")
    ax2.set_xlim(0,1)
    ax2.set_xlabel("Cutoff")
    ax2.set_ylabel("Performance")
    fig.savefig(d+"/Iteration_"+str(n)+".eps")
    return [cutoff, auc]
  
