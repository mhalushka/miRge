miRge2.0
========

A Python program to perform comprehensive analysis of miRNA sequencing Data, including miRNA annotation, A-to-I analysis, novel miRNA detection, etc.


Documentation
-------------

* [Installation](#installation)
  * [Download libraries](#download-libraries)
  * [Install miRge2.0](#download-libraries)
* [How to use it](#how-to-use-it)
  * [Annoate mode](#annotate-mode)
  * [Predict mode](#predict-mode)
* [Citation](#citation)

Installation
------------

### Download libraries

miRge2.0 relies on a huge number of libraries like: <br />
a. Bowtie indexes of genome, hairping, mature miRNAs in miRBase, mature miRNAs in miRGeneDB, mRNA, rRNA, snoRNA, tRNA, other ncRNA and spike-in sequeces (optional) <br />
b. Sequeces of genome, mature miRNAs (including SNP information) in miRBase and miRGeneDB <br />
c. Corrdinates of repetitive elements and mature miRNAs in the genome and miRNA merge information in miRBase and miRGeneDB <br />

Libraries of six species including human, mouse, rat, zebrafish, nematode and fruitfly can be downloaded from: https://jh.box.com/s/wfe2wyhpmutxtu74np6pr1ulas0bgzbi.
In each species, there are three folders: index.Libs (libraries Part a), fasta.Libs (libraries Part b) and annotation.Libs (libraries Part c).

### Install miRge2.0

miRge2.0 is implemented as a Python program running on a Linux/Unix platform that requires pre-installation of Bowtie (http://bowtie-bio.sourceforge.net/index.shtml), SAMtools (http://samtools.sourceforge.net/) and RNAfold (http://www.tbi.univie.ac.at/RNA). <br />
It was built with Python (v2.7.*) programming language and Python-relate libraries, including cutadapt(v1.11), biopython(v1.68), numpy(v1.11.3), scipy(v0.17.0), matplotlib(v2.1.1), pandas(v0.21.0), sklearn(v0.18.1), reportlab(v3.3.0) and forgi(v0.20). <br />
The source code is hosted at: https://github.com/mhalushka/miRge-2.0. miRge2.0 is available in Bioconda, so the quick installation is to run:<br />
1. conda install mirge<br />
2. pip install forgi==0.20 (Since currently forgi is not available in conda)<br />

How to use it
-------------

miRge2.0 has two module: annotate and predict.

### Annotate mode

    Usage: miRge2.0 annotate [-h] [<args>]

    Example:
    miRge2.0 annotate -s seq_file.fastq -d miRBase -pb /home/yin/tools/bowtie-1.1.1/ -lib /home/yin/miRge.Libs/ -sp human -ad illumina -ai -gff -cpu 4

    optional arguments:
      -h, --help            show this help message and exit
      -s [sample <required> [sample <required> ...]]
                            Two options: 1. A file where each row represents one
                            sample name; 2. *.fastq, *fastq, ...
      -d <string required>  The miRNA database (default: miRBase. miRGeneDB is
                            optional)
      -pb <dir required>    The path to the system's bowtie binary
      -lib <dir required>   The path to the miRge libraries
      -sp <string required>
                            The species can be human, mouse, fruitfly, nematode,
                            rat and zebrafish (novle miRNA detection is confined
                            in human and mouse)
      -ad <string>          The adapter need to be removed which could be
                            illumina, ion or a defined sequence (default: none)
      -phred64              phred64 format(default: 64)
      -spikeIn              switch to annotate spike-ins if the bowtie index files
                            are loacted at the path of bowtie's index
                            files(default: off)
      -tcf                  switch to write trimmed and collapsed fasta file
                            (default: off)
      -di                   switch to calculate of isomirs entropy (default: off)
      -cpu <int>            the number of processors to use for trimming, qc, and
                            alignment (default: 1)
      -ai                   switch to calculate of A to I editing (default: off)
      -gff                  switch to output results in gff format (default: off)
      --version             show program's version number and exit


### Predict mode:

    Usage: miRge2.0.py predict [-h] [<args>]

    Example:
    miRge2.0.py predict -s seq_file.fastq -d miRBase -pb /home/yin/tools/bowtie-1.1.1/ -lib /home/yin/miRge.Libs/ -ps /usr/local/bin/ -pr /usr/local/bin/ -sp human -ad illumina -ai -gff -cpu 4

    optional arguments:
      -h, --help            show this help message and exit
      -s [sample <required> [sample <required> ...]]
                            Two options: 1. A file where each row represents one
                            sample name; 2. *.fastq, *fastq, ...
      -d <string required>  The miRNA database (default: miRBase. miRGeneDB is
                            optional)
      -pb <dir required>    The path to the system's bowtie binary
      -lib <dir required>   The path to the miRge libraries
      -sp <string required>
                            The species can be human, mouse, fruitfly, nematode,
                            rat and zebrafish (novle miRNA detection is confined
                            in human and mouse)
      -ps <dir required>    The path to the system's samtools binary
      -pr <dir required>    The path to the system's rnafold binary
      -ad <string>          The adapter need to be removed which could be
                            illumina, ion or a defined sequence (default: none)
      -phred64              phred64 format(default: 64)
      -spikeIn              switch to annotate spike-ins if the bowtie index files
                            are loacted at the path of bowtie's index
                            files(default: off)
      -tcf                  switch to write trimmed and collapsed fasta file
                            (default: off)
      -di                   switch to calculate of isomirs entropy (default: off)
      -cpu <int>            the number of processors to use for trimming, qc, and
                            alignment (default: 1)
      -ai                   switch to calculate of A to I editing (default: off)
      -gff                  switch to output results in gff format (default: off)
      -ws <file>            the file containing the overall samples to analysis
                            for novel miRNA prediction
                            (GodenStandardSampleList.txt)
      -minl <int>           the minimum length of the reatined reads (default: 16)
      -maxl <int>           the maximum length of the reatined reads (default: 25)
      -cc <int>             the maximum read count of the reatined reads (default:
                            2)
      -ml <int>             the maximum number of mapping loci for the retained
                            reads (default: 3)
      -sl <int>             the seed length when invoking Bowtie (default: 25)
      -olc <int>            the length of overlapped seqence when joining reads
                            into longer sequences based on the coordinate on the
                            genome (default: 14)
      -clc <int>            the maximum length of the clustered sequences
                            (default: 30)
      --version             show program's version number and exit

Citation
--------

Yin Lu, Alexander S. Baras, Marc K Halushka. miRge2.0: An updated tool to comprehensively analyze microRNA sequencing data. bioRxiv.org.