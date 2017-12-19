import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np

def preTrimClusteredSeq(CoordinateDic, file1, LenCutoff):
	inf1 = open(file1, 'r')
	outf = open(file1[:-4]+'_trimmed.tsv', 'w')
	line = inf1.readline()
	outf.write("\t".join(['miRClusterID\tOriginalSeq\tFlag\trepetitiveElementName']+line.split('\t')[1:]))
	line = inf1.readline()
	while line != "":
		# flag = "0" means the seq will be abandoned.
		# flag = "1" means the seq will be retained.
		flag = "0"
		repetitiveElementName = '*'
		content = line.strip().split('\t')
		# remove representative seqs with ength larger than 30
		if content[2] == '-':
			seq = Seq(content[5], generic_dna)
			orginalSeq = str(seq.reverse_complement())
		else:
			orginalSeq = content[5]
		if int(content[6]) <= int(LenCutoff):
			chr = content[1]
			start = int(content[3])
			end = int(content[4])
			# remove representative seqs with "AAAAAA" in the 3' and "TTTTTT" in the 5'
			if (re.search('A{6,}$', orginalSeq) is None) and (re.search('^T{6,}', orginalSeq) is None):
			#if (re.search('A{6,}$', orginalSeq) is None):
				if chr in CoordinateDic.keys():
					dist, idx =CoordinateDic[chr][0][0].query([(start, 0)], 2)
					if len(idx) == 1 and np.isfinite(dist[0][0]) and np.isfinite(dist[0][1]):
						candidate1 = CoordinateDic[chr][1][idx[0][0]]
						candidate2 = CoordinateDic[chr][1][idx[0][1]]
						list1 = range(candidate1[0],candidate1[1]+1)
						list2 = range(candidate2[0],candidate2[1]+1)
						listSeq = range(start, end+1)
						if len(set(list1) & set(listSeq)) == 0 and len(set(list2) & set(listSeq)) == 0:
							flag = "1"
						else:
							if len(set(list1) & set(listSeq)) > 0:
								repetitiveElementName = candidate1[2]
							else:
								repetitiveElementName = candidate2[2]
					else:
						candidate1 = CoordinateDic[chr][1][idx[0][0]]
						list1 = range(candidate1[0],candidate1[1]+1)
						listSeq = range(start, end+1)
						if len(set(list1) & set(listSeq)) == 0:
							flag = "1"
						else:
							repetitiveElementName = candidate1[2]
				else:
					flag = "1"
			else:
				pass
		else:
			pass
		outf.write('\t'.join([content[0], orginalSeq, flag, repetitiveElementName]+content[1:])+'\n')
		line = inf1.readline()
	inf1.close()
	outf.close()
