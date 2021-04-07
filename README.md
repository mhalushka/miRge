miRge2.0
========

A Python program to perform comprehensive analysis of miRNA sequencing Data, including miRNA annotation, A-to-I analysis, novel miRNA detection, isomiR analysis and tRF detection etc.

:warning:ATTENTION:
==========
**The latest software update of [miRge](https://github.com/mhalushka/miRge3.0) is now available with much more advanced features, including support to Python3.0 and improvised speed performance. Check the documentation of miRge3.0 [here](https://mirge3.readthedocs.io/).**


Documentation
-------------

* [Installation](#installation)
  * [Download libraries](#download-libraries)
  * [Install miRge2.0](#download-libraries)
  * [Troubleshooting installation](#troubleshooting-installation)
* [How to use it](#how-to-use-it)
  * [Annoate mode](#annotate-mode)
  * [Predict mode](#predict-mode)
* [Changelog](#changelog)
* [Citation](#citation)

Installation
------------

### Download libraries

miRge2.0 relies on a huge number of libraries like: <br />
1) Bowtie indexes of genome, hairping, mature miRNAs in miRBase, mature miRNAs in miRGeneDB, mRNA, rRNA, snoRNA, mature tRNA, primary tRNA, other ncRNA and spike-in sequences (optional) <br />
2) Sequences of genome, mature miRNAs (including SNP information) in miRBase and miRGeneDB <br />
3) Corrdinates of repetitive elements and mature miRNAs in the genome and miRNA merge information in miRBase and MirGeneDB <br />

Libraries of six species including ___human___, ___mouse___, ___rat___, ___zebrafish___, ___nematode___ and ___fruitfly___ can be downloaded separately. Our institution changed storage option and they can now be obtained using these links:

[Human](https://livejohnshopkins-my.sharepoint.com/:u:/g/personal/mhalush1_jh_edu/EUkL2ZbRpKRXJpemY32dObgBCyiiPbIiuzKCglrtBUIsFA?e=UIqTDC) <br />
[Mouse](https://livejohnshopkins-my.sharepoint.com/:u:/g/personal/mhalush1_jh_edu/EcoSDcx9DxNb295CfKEdEhYBiYNWnuZLAuYKBHYVjfjSfA?e=dXPGbS) <br />
[Rat](https://livejohnshopkins-my.sharepoint.com/:u:/g/personal/mhalush1_jh_edu/ET3zIO3qDY1bTwLqYZ50MnIBchdOkG8BKtF0Lf2epIeEHA?e=rBkXv3) <br />
[Zebrafish](https://livejohnshopkins-my.sharepoint.com/:u:/g/personal/mhalush1_jh_edu/EZ88aXFvLthRX6b-TaawAPoB1wxSV0cLL-EOpG9HsmQlSw?e=K2WmKe) <br />
[Nematode](https://livejohnshopkins-my.sharepoint.com/:u:/g/personal/mhalush1_jh_edu/Eczh_ADH8cRfMRBy9PIXS20Bv2HBBYlv3zSYcMoQihzvCQ?e=Rj1vhp) <br />
[Fruitfly](https://livejohnshopkins-my.sharepoint.com/:u:/g/personal/mhalush1_jh_edu/EUme-YXTWw5c7ywnOI5OLQ0B-_itElNcKQ3hFFFO3x8aBw?e=qwR9B0) <br />


After unpacking the downloaded ***.tar.gz files to the new directory of miRge.Libs, the folder named by species contains three folders: index.Libs (libraries Part 1), fasta.Libs (libraries Part 2) and annotation.Libs (libraries Part 3). The absolute directory of miRge.Libs is used as the value of parameter ‘-lib’ in miRge2.0. <br /><br />
If the users want to build the libraries for other species, they can use scripts miRge_bowtie_build.py and miRge_pckls_build.py which can be downloaded from https://github.com/mhalushka/miRge_build, wherein miRge_bowtie_build.py is used to build bowtie index files and miRge_pckls_build.py is used to transform ***_genome_repeats.GTF and ***_genome.fa in oder to accelerate the speed of reading larg files into memory. 

### Install miRge2.0

miRge2.0 is implemented as a Python program running on a Linux/Unix platform that requires pre-installation of Bowtie (v1.1.1 or v1.1.2; http://bowtie-bio.sourceforge.net/index.shtml), SAMtools (v1.5; http://samtools.sourceforge.net/) and RNAfold (v2.3.5; http://www.tbi.univie.ac.at/RNA). <br />
It was built with Python (v2.7.*) programming language and Python-related libraries, including cutadapt(v1.11 to v1.16), biopython(>= v1.68), numpy(>= v1.11.3), scipy(>= v0.17.0), matplotlib(>= v2.1.1), pandas(>= v0.21.0), sklearn(>= v0.18.1), reportlab(>= v3.3.0) and forgi(v0.20). <br />

PLEASE NOTE: miRge2.0 is currently incompatible with cutadapt v1.18.  Using v1.18 will give a "TypeError: __call__() takes exactly 3 arguments (2 given)" error.<br />

The source code is hosted at: https://github.com/mhalushka/miRge.<br />

miRge2.0 is available in bioconda now and can be installed via conda. You will need to add the bioconda channel as well as the other channels bioconda depends on. It is important to add them in this order so that the priority is set correctly (that is, bioconda is highest priority).<br />
    `conda config --add channels defaults`<br />
    `conda config --add channels conda-forge`<br />
    `conda config --add channels bioconda`<br />
Then install mirge by running:<br />
    `conda install mirge`<br />

miRge2.0 also can be installed from the source code by pip (THIS VERSION IS CURRENTLY UNDER REVISION. USE BIOCONDA FOR AN ACCURATE INSTALL UNTIL THIS MESSAGE DISAPPEARS):<br />
1) Download miRge2.0 source code from https://github.com/mhalushka/miRge and unzip the zipped file folder.<br />
2) If the package of wheel is not installed, run `pip install wheel` to install it.<br />
3) Change the directory to miRge2.0's directory and run `python setup.py bdist_wheel` to build a wheel file for the subsequent installation via pip.<br />
4) Run `pip install ./dist/mirge-2.0-py2-none-any.whl` to install miRge2.0.<br />


### Troubleshooting installation<br />
miRge2.0 was tested on the specific version of required softwares and python packages. Please make sure the version is correct.<br />
1) If Bowtie, SAMtools or RNAfold have been already installed in the system, please run `which bowtie`, `which samtools` or `which RNAfold` to find their installation paths. If the versions are incorrect, please install them with right version. 
2) Running `pip freeze` to check th version of current python packages. If some python packages can't work, please mannually install them by running `pip install package==*.**`.<br />
3) If the required python pacakages of the specific version can't be installed by pip or imported by python, make sure the installed python is complied by 4-byte Unicode so that pip can install UCS4 wheels (supporting cp27mu not cp27m). Type python and enter following commands `import sys` `print  sys.maxunicode`. If output is 1114111 then it is UCS4 otherwise if output is 65535 then it is UCS2.<br />
   If it is UCS2, please re-compile already installed python with 4-bype Unicode from the source code by running: a) `./configure --enable-unicode=ucs4 --prefix=***` b) `make` c) `make install` <br />
   

How to use it
-------------

miRge2.0 has two modules: annotate and predict.

### Annotate mode

    Usage: miRge2.0 annotate [-h] [<args>]

    Example:
    miRge2.0 annotate -s seq_file.fastq -d miRBase -pb /home/yin/tools/bowtie-1.1.1 -lib /home/yin/miRge.Libs -sp human -ad illumina -ai -gff -trf -cpu 4

    optional arguments:
      -h, --help            show this help message and exit
      -s [sample <required> [sample <required> ...]]
                            two options: 1. A file where each row represents one
                            sample name; 2. *.fastq *.fastq ...
      -o <dir>              the directory of the outputs (default: current directory)
      -d <string required>  the miRNA database (default: miRBase. MirGeneDB is optional)
      -pb <dir required>    the path to the system's bowtie binary
      -lib <dir required>   the path to the miRge libraries
      -sp <string required>
                            the species can be human, mouse, fruitfly, nematode,
                            rat and zebrafish (novel miRNA detection is confined
                            in human and mouse)
      -ex <float>           the threshold of the proportion of canonical reads for the
                            miRNAs to determine whether keeping them or not when counting.
                            Users can set it between 0 and 0.5 (default: 0.1)
      -ad <string>          the adapter need to be removed which could be
                            illumina, ion or a defined sequence (default: none)
      -phred64              phred64 format (default: 33)
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
      -trf                  switch to analyze tRNA fragment (default: off)
      --version             show program's version number and exit


### Predict mode:

    Usage: miRge2.0 predict [-h] [<args>]

    Example:
    miRge2.0 predict -s seq_file.fastq -d miRBase -pb /home/yin/tools/bowtie-1.1.1 -lib /home/yin/miRge.Libs -ps /usr/local/bin -pr /usr/local/bin -sp human -ad illumina -ai -gff -trf -cpu 4

    optional arguments:
      -h, --help            show this help message and exit
      -s [sample <required> [sample <required> ...]]
                            two options: 1. A file where each row represents one
                            sample name; 2. *.fastq *.fastq ...
      -o <dir>              the directory of the outputs (default: current directory)
      -d <string required>  the miRNA database (default: miRBase. MirGeneDB is optional)
      -pb <dir required>    the path to the system's bowtie binary
      -lib <dir required>   the path to the miRge libraries
      -sp <string required>
                            the species can be human, mouse, fruitfly, nematode,
                            rat and zebrafish (novel miRNA detection is confined
                            in human and mouse)
      -ps <dir required>    the path to the system's samtools binary
      -pr <dir required>    the path to the system's rnafold binary
      -ex <float>           the threshold of the proportion of canonical reads for the
                            miRNAs to determine whether keeping them or not when counting.
                            Users can set it between 0 and 0.5 (default: 0.1)
      -ad <string>          the adapter need to be removed which could be
                            illumina, ion or a defined sequence (default: none)
      -phred64              phred64 format (default: 33)
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
      -trf                  switch to analyze tRNA fragment (default: off)
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

Changelog
---------

**06-26-2018**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Add a `-trf` option, so that tRF detection can be performed.

**05-06-2018**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Rebuild the miRNA libraries from the newly released miRBase v22 and MirGeneDB v2.0.

**02-22-2018**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Add a `-ex` option. If the canonical reads count is less than 2 or the proportion of canonical reads is less that the set threshold, the miRNA will not be counted in the ouput files: miR.Counts.csv and miR.RPM.csv.

**07-15-2020**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Update the links to the libraries.

Citation
--------

Yin Lu, Alexander S. Baras, Marc K Halushka. miRge 2.0 for comprehensive analysis of microRNA sequencing data. <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2287-y" target="_blank">BMC Bioinformatics</a>. 2018.
