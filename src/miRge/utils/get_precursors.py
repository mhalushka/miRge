import os
import sys
import cPickle
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def headDashCount1(clusterSeq):
	# clusterSeq is: CTAATACTGCCTGGTAATGATGACGG-
	headDashCount = 0
	for symbol in clusterSeq:
		if symbol == '-':
			headDashCount = headDashCount + 1
		else:
			break
	return headDashCount

def tailDashCount1(clusterSeq):
	# clusterSeq is: CTAATACTGCCTGGTAATGATGACGG-
	tailDashCount = 0
	for symbol in clusterSeq[::-1]:
		if symbol == '-':
			tailDashCount = tailDashCount + 1
		else:
			break
	return tailDashCount

def get_precursors(infTmp, chrSeqDic):
	# The arguments of this function are as follows:
	# /home/yin/lib/human/hg38.pckl  unmapped_mirna_SRR944031_vs_representative_seq_bowtie2_modified_selected.tsv
	#with open(infTmp1, 'rb') as f:
	#	chrSeqDic = cPickle.load(f)
	clusterSeqType = 'stableClusterSeq'
	# To do, redefine the startPos and endPos....
	outf = open(infTmp[:-4]+"_precusor.fa","w")
	#selectRangeList = [(24, 74), (22, 72), (20, 70), (18, 68), (16, 66)]
	selectRangeList = [(20, 70)]
	precusorList = []
	with open(infTmp,'r') as inf:
		line = inf.readline()
		clusterSeqTypeLabel = line.strip().split('\t').index(clusterSeqType)
		headUnstableLengthLabel = line.strip().split('\t').index('headUnstableLength')
		tailUnstableLengthLabel = line.strip().split('\t').index('tailUnstableLength')
		alignedClusterSeqLabel = line.strip().split('\t').index('alignedClusterSeq')
		line = inf.readline()
		#i = 1
		while line != "":
			content = line.strip().split('\t')
			temp =content[5]
			neighborState = content[-1]
			alignedClusterSeq = content[alignedClusterSeqLabel]
			headDashCountTmp = headDashCount1(alignedClusterSeq)
			tailDashCountTmp = tailDashCount1(alignedClusterSeq)
			#if neighborState == 'Good':
			#	if i % 2 == 1:
			#		side = 'downstream'
			#	else:
			#		side = 'upstream'
			#	i = i + 1
			#else:
			#	side = 'both'
			side = 'both'
			if temp not in precusorList:
				chr = temp.split(':')[2]
				if chr in chrSeqDic.keys():
					startPos = int(temp.split(':')[3][:-1].split('_')[0].strip())
					endPos = int(temp.split(':')[3][:-1].split('_')[1].strip())
					strand = temp[-1]
					if clusterSeqType == 'stableClusterSeq':
						if strand == '+':
							startPos = startPos - headDashCountTmp + int(content[headUnstableLengthLabel])
							endPos = endPos + tailDashCountTmp - int(content[tailUnstableLengthLabel])
							#print '%s: %s'%(temp, str(Seq(chrSeqDic[chr][startPos-1:endPos], generic_dna).transcribe()))
						else:
							startPos = startPos - tailDashCountTmp + int(content[tailUnstableLengthLabel])
							endPos = endPos + headDashCountTmp - int(content[headUnstableLengthLabel])
							#startPos = startPos - headDashCountTmp + int(content[tailUnstableLengthLabel])
							#endPos = endPos + tailDashCountTmp - int(content[headUnstableLengthLabel])
							#print '%s: %s'%(temp, str(Seq(chrSeqDic[chr][startPos-1:endPos], generic_dna).reverse_complement().transcribe()))
					for index, item in enumerate(selectRangeList):
						unpstream = item[1] # upstream = 70
						downstream = item[0] # downstream = 20
						# :precusor1_1 means extend around 70 bp upstream and 20 bp downstream
						# :precusor1_2 means extend around 20 bp upstream and 70 bp downstream
						if side == 'both':
							if startPos-1-unpstream >= 0:
								#outf.write('>'+temp+':precusor'+str(index+1)+'_1\n')
								outf.write('>'+temp+':precusor_1\n')
								if strand == '+':
									rnaSeq = str(Seq(chrSeqDic[chr][startPos-1-unpstream:endPos+downstream], generic_dna).transcribe())
								else:
									rnaSeq = str(Seq(chrSeqDic[chr][startPos-1-unpstream:endPos+downstream], generic_dna).reverse_complement().transcribe())
								outf.write(rnaSeq+'\n')
							else:
								outf.write('>'+temp+':precusor_1\n')
								if strand == '+':
									rnaSeq = str(Seq(chrSeqDic[chr][:endPos+downstream], generic_dna).transcribe())
								else:
									rnaSeq = str(Seq(chrSeqDic[chr][:endPos+downstream], generic_dna).reverse_complement().transcribe())
								outf.write(rnaSeq+'\n')
								#pass
							if startPos-1-downstream >= 0:
								#outf.write('>'+temp+':precusor'+str(index+1)+'_2\n')
								outf.write('>'+temp+':precusor_2\n')
								if strand == '+':
									rnaSeq = str(Seq(chrSeqDic[chr][startPos-1-downstream:endPos+unpstream], generic_dna).transcribe())
								else:
									rnaSeq = str(Seq(chrSeqDic[chr][startPos-1-downstream:endPos+unpstream], generic_dna).reverse_complement().transcribe())
								outf.write(rnaSeq+'\n')
							else:
								outf.write('>'+temp+':precusor_2\n')
								if strand == '+':
									rnaSeq = str(Seq(chrSeqDic[chr][:endPos+unpstream], generic_dna).transcribe())
								else:
									rnaSeq = str(Seq(chrSeqDic[chr][:endPos+unpstream], generic_dna).reverse_complement().transcribe())
								outf.write(rnaSeq+'\n')
								#pass
						elif side == 'downstream':
							if startPos-1-downstream >= 0:
								#outf.write('>'+temp+':precusor'+str(index+1)+'_2\n')
								outf.write('>'+temp+':precusor_2\n')
								if strand == '+':
									rnaSeq = str(Seq(chrSeqDic[chr][startPos-1-downstream:endPos+unpstream], generic_dna).transcribe())
								else:
									rnaSeq = str(Seq(chrSeqDic[chr][startPos-1-downstream:endPos+unpstream], generic_dna).reverse_complement().transcribe())
								outf.write(rnaSeq+'\n')
							else:
								outf.write('>'+temp+':precusor_2\n')
								if strand == '+':
									rnaSeq = str(Seq(chrSeqDic[chr][:endPos+unpstream], generic_dna).transcribe())
								else:
									rnaSeq = str(Seq(chrSeqDic[chr][:endPos+unpstream], generic_dna).reverse_complement().transcribe())
								outf.write(rnaSeq+'\n')
								#pass
						else:
							if startPos-1-unpstream >= 0:
								#outf.write('>'+temp+':precusor'+str(index+1)+'_1\n')
								outf.write('>'+temp+':precusor_1\n')
								if strand == '+':
									rnaSeq = str(Seq(chrSeqDic[chr][startPos-1-unpstream:endPos+downstream], generic_dna).transcribe())
								else:
									rnaSeq = str(Seq(chrSeqDic[chr][startPos-1-unpstream:endPos+downstream], generic_dna).reverse_complement().transcribe())
								outf.write(rnaSeq+'\n')
							else:
								outf.write('>'+temp+':precusor_1\n')
								if strand == '+':
									rnaSeq = str(Seq(chrSeqDic[chr][:endPos+downstream], generic_dna).transcribe())
								else:
									rnaSeq = str(Seq(chrSeqDic[chr][:endPos+downstream], generic_dna).reverse_complement().transcribe())
								outf.write(rnaSeq+'\n')
								#pass
					precusorList.append(temp)
				else:
					pass
			else:
				pass
			line = inf.readline()
	outf.close()