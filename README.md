miRge2.0
========

A Python program to perform comprehensive analysis of miRNA sequencing Data, including miRNA annotation, A-to-I analysis, novel miRNA detection, etc.


Documentation
-------------

* [Installation](#installation)
  * [Download libraries](#download-libraries)
  * [Install miRge2.0](#download-libraries)
  * [Troubleshooting installation](#troubleshooting-installation)
* [How to use it](#how-to-use-it)
  * [Annoate mode](#annotate-mode)
  * [Predict mode](#predict-mode)
* [Citation](#citation)

Installation
------------

### Download libraries

miRge2.0 relies on a huge number of libraries like: <br />
1) Bowtie indexes of genome, hairping, mature miRNAs in miRBase, mature miRNAs in miRGeneDB, mRNA, rRNA, snoRNA, tRNA, other ncRNA and spike-in sequeces (optional) <br />
2) Sequeces of genome, mature miRNAs (including SNP information) in miRBase and miRGeneDB <br />
3) Corrdinates of repetitive elements and mature miRNAs in the genome and miRNA merge information in miRBase and miRGeneDB <br />

Libraries of six species including ___human___, ___mouse___, ___rat___, ___zebrafish___, ___nematode___ and ___fruitfly___ can be downloaded separately by wget. `wget -O human.tar.gz https://jh.box.com/shared/static/rj7ufy5v15uw7ytsyyrsryw99u7ml82j.gz`; <br/>`wget -O mouse.tar.gz https://jh.box.com/shared/static/z2bcey8j9e9nxnvpmb4fm88zzq3da4m1.gz`; `wget -O rat.tar.gz https://jh.box.com/shared/static/mmztv42j8h7snk0eo80o7a7t30it4q9f.gz`; `wget -O zebrafish.tar.gz https://jh.box.com/shared/static/nwn7jzn5ekgm51k7jlk43a6h75aasgr1.gz`; `wget -O nematode.tar.gz https://jh.box.com/shared/static/boubl0g28b0xj952e8bhfmbds056gdvi.gz`; `wget -O fruitfly.tar.gz ***`.
After unpacking the downloaded ***.tar.gz files to the new directory of miRge.Libs, the folder named by species contains three folders: index.Libs (libraries Part 1), fasta.Libs (libraries Part 2) and annotation.Libs (libraries Part 3). The absolute directory of miRge.Libs is used as the value of parameter ‘-lib’ in miRge2.0. <br /><br />
If the users want to build the libraries for other species, they can use scripts miRge_bowtie_build.py and miRge_pckls_build.py which can be downloaded from https://github.com/mhalushka/miRge2.0_build, wherein miRge_bowtie_build.py is used to build bowtie index files and miRge_pckls_build.py is used to transform ***_genome_repeats.GTF and ***_genome.fa in oder to accelerate the speed of reading larg files into memory. 

### Install miRge2.0

miRge2.0 is implemented as a Python program running on a Linux/Unix platform that requires pre-installation of Bowtie (v1.1.1 or v1.1.2; http://bowtie-bio.sourceforge.net/index.shtml), SAMtools (v1.5; http://samtools.sourceforge.net/) and RNAfold (v2.3.5; http://www.tbi.univie.ac.at/RNA). <br />
It was built with Python (v2.7.*) programming language and Python-related libraries, including cutadapt(v1.11), biopython(v1.68), numpy(v1.11.3), scipy(v0.17.0), matplotlib(v2.1.1), pandas(v0.21.0), sklearn(v0.18.1), reportlab(v3.3.0) and forgi(v0.20). <br />
The source code is hosted at: https://github.com/mhalushka/miRge-2.0.<br />
miRge2.0 can be installed from the source code by pip:<br />
1) Download miRge2.0 source code from https://github.com/mhalushka/miRge-2.0 and unzip the zipped file folder.<br />
2) If the package of wheel is not installed, run `pip install wheel` to install it.<br />
3) Change the directory to miRge2.0's directory and run `python setup.py bdist_wheel` to build a wheel file for the subsequent installation via pip.<br />
4) Run `pip install ./dist/mirge-2.0-py2-none-any.whl` to install miRge2.0.<br />

miRge2.0 will be in Bioconda as well (in process), the installation is to run:<br />
1) `conda install mirge`<br />
2) `pip install forgi==0.20` (since currently forgi is not available in conda)<br />



### Troubleshooting installation<br />
miRge2.0 was tested on the specific version of required softwares and python packages. Please make sure the version is correct.<br />
1) If Bowtie, SAMtools or RNAfold have been already installed in the system, please run `which bowtie`, `which samtools` or `which RNAfold` to find their installation paths. If the versions are incorrect, please install them with right version. 
2) Running `pip freeze` to check th version of current python packages. If some python packages can't work, please mannually install them by running `pip install package==*.**`.<br />
3) If the required python pacakages of the specific version can't be installed by pip or imported by python, make sure the installed python is complied by 4-byte Unicode so that pip can install UCS4 wheels (supporting cp27mu not cp27m). Type python and enter following commands `import sys` `print  sys.maxunicode`. If output is 1114111 then it is UCS4 otherwise if output is 65535 then it is UCS2.<br />
   If it is UCS2, please re-compile already installed python with 4-bype Unicode from the source code by running: a) `./configure --enable-unicode=ucs4 --prefix=***` b) `make` c) `make install` <br />
   

How to use it
-------------

miRge2.0 has two module: annotate and predict.

### Annotate mode

    Usage: miRge2.0 annotate [-h] [<args>]

    Example:
    miRge2.0 annotate -s seq_file.fastq -d miRBase -pb /home/yin/tools/bowtie-1.1.1 -lib /home/yin/miRge.Libs -sp human -ad illumina -ai -gff -cpu 4

    optional arguments:
      -h, --help            show this help message and exit
      -s [sample <required> [sample <required> ...]]
                            two options: 1. A file where each row represents one
                            sample name; 2. *.fastq *.fastq ...
      -o <dir>              the directory of the outputs (default: current directory)
      -d <string required>  the miRNA database (default: miRBase. miRGeneDB is optional)
      -pb <dir required>    the path to the system's bowtie binary
      -lib <dir required>   the path to the miRge libraries
      -sp <string required>
                            the species can be human, mouse, fruitfly, nematode,
                            rat and zebrafish (novel miRNA detection is confined
                            in human and mouse)
      -ad <string>          the adapter need to be removed which could be
                            illumina, ion or a defined sequence (default: none)
      -phred64              phred64 format (default: 64)
      -spikeIn              switch to annotate spike-ins if the bowtie index files
                            are loacted at the path of bowtie's index
                            files (default: off)
      -tcf                  switch to write trimmed and collapsed fasta file
                            (default: off)
      -di                   switch to calculate of isomirs entropy (default: off)
      -cpu <int>            the number of processors to use for trimming, qc, and
                            alignment (default: 1)
      -ai                   switch to calculate of A to I editing (default: off)
      -gff                  switch to output results in gff format (default: off)
      --version             show program's version number and exit


### Predict mode:

    Usage: miRge2.0 predict [-h] [<args>]

    Example:
    miRge2.0.py predict -s seq_file.fastq -d miRBase -pb /home/yin/tools/bowtie-1.1.1 -lib /home/yin/miRge.Libs -ps /usr/local/bin -pr /usr/local/bin -sp human -ad illumina -ai -gff -cpu 4

    optional arguments:
      -h, --help            show this help message and exit
      -s [sample <required> [sample <required> ...]]
                            two options: 1. A file where each row represents one
                            sample name; 2. *.fastq *.fastq ...
      -o <dir>              the directory of the outputs (default: current directory)
      -d <string required>  the miRNA database (default: miRBase. miRGeneDB is optional)
      -pb <dir required>    the path to the system's bowtie binary
      -lib <dir required>   the path to the miRge libraries
      -sp <string required>
                            the species can be human, mouse, fruitfly, nematode,
                            rat and zebrafish (novel miRNA detection is confined
                            in human and mouse)
      -ps <dir required>    the path to the system's samtools binary
      -pr <dir required>    the path to the system's rnafold binary
      -ad <string>          the adapter need to be removed which could be
                            illumina, ion or a defined sequence (default: none)
      -phred64              phred64 format (default: 64)
      -spikeIn              switch to annotate spike-ins if the bowtie index files
                            are loacted at the path of bowtie's index
                            files (default: off)
      -tcf                  switch to write trimmed and collapsed fasta file
                            (default: off)
      -di                   switch to calculate of isomirs entropy (default: off)
      -cpu <int>            the number of processors to use for trimming, qc, and
                            alignment (default: 1)
      -ai                   switch to calculate of A to I editing (default: off)
      -gff                  switch to output results in gff format (default: off)
      -ws <file>            the file containing the overall samples to analysis
                            for novel miRNA prediction. No header, just a list of *.fastq
                            file names in a column. Names of files can be to your choosing
                            (e.g. filestochecknovel.txt)
      -minl <int>           the minimum length of the reatined reads for novel miRNA
                            detection (default: 16)
      -maxl <int>           the maximum length of the reatined reads for novel miRNA
                            detection (default: 25)
      -cc <int>             the maximum read count of the reatined reads for novel miRNA
                            detection (default: 2)
      -ml <int>             the maximum number of mapping loci for the retained
                            reads for novel miRNA detection (default: 3)
      -sl <int>             the seed length when invoking Bowtie for novel miRNA detection
                            (default: 25)
      -olc <int>            the length of overlapped seqence when joining reads
                            into longer sequences based on the coordinate on the
                            genome for novel miRNA detection (default: 14)
      -clc <int>            the maximum length of the clustered sequences for novel miRNA
                            detection (default: 30)
      --version             show program's version number and exit

Citation
--------

Yin Lu, Alexander S. Baras, Marc K Halushka. miRge2.0: An updated tool to comprehensively analyze microRNA sequencing data. bioRxiv.org.
