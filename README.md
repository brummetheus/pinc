# Pinc
Prediction of lncRNA based on RNA-Seq data. 

*Pinc* is available in two modes: either as standalone nextflow pipeline, or as a docker image!

# Pinc as a Nextflow pipeline

Running *Pinc* as a Nextflow pipeline requires, that all used programs are installed and executable from $PATH. In addition all scripts in the folder *custom_scripts* also need to be in $PATH. 

## Prerequisite programs

Most tools can be installed by just following the basic installation instructions. Prerequisites for individual programs can be found on the respective websites. 

### Python3
available here: [Python3](https://www.python.org/downloads/)
### R
available here: [R](https://www.r-project.org/)

### Nextflow 
available here: [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

### fastp
available here: [fastp](https://github.com/OpenGene/fastp)

### HISAT2
available here: [HISAT2](https://github.com/DaehwanKimLab/hisat2)

### StringTie
available here: [StringTie](https://github.com/gpertea/stringtie)

### gffcompare
available here: [gffcompare](https://github.com/gpertea/gffcompare)

### gffread
available here: [gffread](https://github.com/gpertea/gffread)

### CPC2
available here: [CPC2](http://cpc2.gao-lab.org/download.php)

Recommended is the standalone version which supports **Python3** as this version was tested as well as the support for **Python2** ended. 

### CPAT
available here: [CPAT](https://cpat.readthedocs.io/en/latest/#installation)

**CPAT** can be installed via **pip** or **pip3** depending on the default **Python** version on the system: 

`pip install CPAT`

### HTSeq
available here: [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/install.html)

### edgeR
available here: [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

## Running Pinc as Nextflow script
In order to run *Pinc*, a [configuration file](https://github.com/brummetheus/pinc/edit/main/README.md#the-pinc-configuration-file) file needs to be given. This file contains all necessary information on where the data is located and configuration settings for the different programs used in the pipeline. A detailed description can be found later. 

After the setup of the configuration file is finished *Pinc* can be started with the following command:

`nextflow run pinc.nf -c pinc.config` 

# Pinc as a Docker container
Running *Pinc* as a Docker container requires only one program to be installed: [Docker](https://docs.docker.com/get-docker/). All used programs and tools come with the image and do not require any further time consuming installations. After setup for Docker is finished, *Pinc* can be downloaded with the following command:

`docker pull brummetheus/pinc`

After downloading the image successfully, *Pinc* can be run via this command: 

`docker run -v /path/to/project_Directory:/home/docker -w /home/docker pinc nextflow run pinc.nf -c pinc.config`

In order for the container to access files on the host filesystem, a directory needs to be mounted to the containers **/home/docker** directory. As a recommendation, the **project_Directory** should contain all necessary files for the analysis. An examplery directory might look like this: 
```
.
+-- pinc.nf
+-- pinc.config
+-- aux_files
|   +-- genome.fna
|   +-- genome_annotation.gff
|   +-- design.tsv
|   +-- mRNA.fna
+-- reads
|   +-- e_coli_a1.fastq.gz
|   +-- e_coli_a2.fastq.gz
|   +-- e_coli_t1.fastq.gz
|   +-- e_coli_t2.fastq.gz
+-- output_folder
```

## The Pinc configuration file
```
params {
    outdir = "$projectDir/path/to/output_folder"
    data_dir = "$projectDir/path/to/sequencing_reads.gz"
    genome = "$projectDir/path/to/genome.fna"
    annot = "$projectDir/path/to/genome_annotation.gff"
    design = "$projectDir/path/to/design.tsv"
    mRNA = "$projectDir/path/to/mRNA.fna"
    strand = "rf"
    singleEnds = false
    de_analysis = true
    weight = 0.5
    cpu_ind = 4
}

executor.$local.cpus = 12
```

`$projectDir` points to the directory, in which **Pinc.nf** is located. Therefore, relative paths can be used very efficiently. However, when using Pinc as a nextflow pipeline, absolute paths may also be given. 

For users of the containered version of *Pinc*: When mounting a directory to the container, all files required for *Pinc* need to be located either in the directory or in subdirectories. All paths in **pinc.config** need to be relative paths to `project_Dir`. 

#### outdir
Specifies the directory where all output files will be saved.
#### data_dir
Specifies the directory where the sequencing reads are stored. Either fastq files or gzipped files are recognized automatically. In case of paired-end reads *Pinc* will group files together, which have the same basename. E.g. **e_coli_a1_S1_L001_R1.fastq** and **e_coli_a1_S1_L001_R2.fastq** will be grouped together as the basename **e_coli_a1_S1_L001_R** is identical. 
#### genome
Specifies which genome file has to be used. 
#### annot
Specifies which genome annotation file has to be used. THe anntotation file has to be in gff3 format and chromosome names must correspond to the chromosome names in the genome file.
#### design.tsv
Specifies which design file has to be used. The design file is only needed, if a DE-analysis shall be performed. It is a plain csv-file which contains the sample name in the first column, and the group tag in the second column. Automatic file extensions, which are added to sample names, like **S1_L001** in the case of Illumina instruments may be ommitted. 
```
e_coli_a1   1
e_coli_a2   1
e_coli_t1   2
e_coli_t2   2
```
#### mRNA
Specifies the file which contains all mRNA sequences in fasta format. 
#### strand
Specifies the strandedness of the RNA-Seq library. For single-end read sequencing data options can be either **"f"**, **"r"** or **"u"** representing forward-strand, reverse-strand and unstranded libraries respectively. Similarly for paired-end sequencing options are **"fr"**, **"rf"** or **"un"**. If in doubt, this [summary](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/) will help to decide which option has to be used. 
#### singleEnds
Specifies if single-end reads are present or not. Set this to false, if paired-end sequencing was performed. E.g:
`singleEnd = false`
#### de_analysis
Specifies if a DE-analysis shall be performed. Reminder: For now DE-analysis can only be performed between two conditions. However, the number of replicates for each condition is unlimited. If one or both conditions is represented by only one sample, DE-analysis is restricted to calculations of foldchanges.
#### weight
Specifies the weight for the Weighted Youden's Index to optimize the cutoff used by CPAT. The weight can be between 0 and 1. A weight of 0.5 results in equal contributions of True-Positive-Rate (TPR) and False-Positive-Rate (FPR) to the Youden's Index. Increasing the weight will favor TPR over FPR resulting in more correct predicted ncRNAs at the price of more false-positives (higher contamination with coding RNAs). Decreasing the weight will increase the impact of FPR which leads to less contamination with coding RNAs at the price of potentially missing more ncRNAs. 
#### cpu_ind
Specifies how many cores are available for each individual process. 

#### executor.$local.cpus
Specifies how many core are available for *Pinc* in general. For example cpu_ind is set to 4 and executor.$local.cpus means that 3 processes can run in parallel. It is recommended to use multiple of cpu_ind to ensure smooth transitions between processes. 
