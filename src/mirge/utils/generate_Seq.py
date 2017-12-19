import os
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def generate_Seq(infTmp):
	# The arugment of this function is: unmapped_mirna_SRR944031_vs_genome_sorted_clusters_overlap14_trimmed.tsv
	outf1 = open(infTmp[:-4]+'.fa','w')
	outf2 = open(infTmp[:-4]+'_orig.fa','w')
	inf = open(infTmp,'r')
	line = inf.readline()
	line = inf.readline()
	while line != '':
		content = line.strip().split('\t')
		if content[2] == '1':
			outf1.write('>'+content[0]+":"+content[4]+":"+content[6]+"_"+content[7]+content[5]+'\n')
			outf2.write('>'+content[0]+":"+content[4]+":"+content[6]+"_"+content[7]+content[5]+'\n')
			outf1.write(content[8]+'\n')
			outf2.write(content[1]+'\n')
		else:
			pass
		line = inf.readline()
	inf.close()
	outf1.close()
	outf2.close()
