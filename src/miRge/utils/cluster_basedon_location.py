#
#Here, overlapLenThreshold is the length of the overlapped sequence between the two adjoining sequnces according to their coordinates on the genome.
#
# 1)If overlapLenThreshold = 10, mir63645_2 and mir4190_2 are linked into one miRCluster_370_32 in mapped_mirna_SRR944031_vs_genome_sorted_clusters_10.tsv.
#However,
#                               CCACAACTAACCCATAATATGGACATGTTATA	miRCluster_370_32
#          TCCATATTATGGGTTAGTTGTGG	mir63645_2 (reverse); hsa-miR-3912-5p
# TATAACATGTCCATATTATGGG	        mir4190_2 (reverse); hsa-miR-3912-3p
#In this case, the overlapped length is 13.
#So the optimal overlapLenThreshold could be selected as 14.
#
# 2)If overlapLenThreshold = 14, mir1331_2,mir42299_3,mir27985_4,mir75999_2 are linked into one miRCluster_607_29 in mapped_mirna_SRR944031_vs_genome_sorted_clusters_14.tsv.
#   CTAGGTCGGTGCAAAAGTCATCACGGTTT	miRCluster_607_29
#   CTAGGTCGGTGCAAAAGTCATC	        mir1331_2 (forward); hsa-miR-548aa(isoMIR)
#        TCGGTGCAAAAGTCATCACGGT	    mir42299_3(forward); hsa-miR-548aw(isoMIR)
#         CGGTGCAAAAGTCATCACGGTT	mir27985_4 (forward); hsa-miR-548aw
#           GTGCAAAAGTCATCACGGTTT	mir75999_2 (forward); hsa-miR-548aw
# 3)mir1331_2 is wrongly clustered int to miRCluster_607_29. So another restrictive condition should be added, which is the lenthg of the unoverlapped sequence in the 5'.
#  The situation in 2) is unavoidable. Howvever, in order generate more continuous cluster sequence, unoverlappedLenThreshold shold not be set.
# 4)Direction information should be considered as well.
# In *.sam, if flag == "16", reads are mapped to the reverse strand.
# If flag == "0", reads are mapped to the forward strand.
#
import os

def cluster_basedon_location(infTmp, overlapLenThresholdTmp):
	# The arguments of this function are: unmapped_mirna_vs_genome_sorted.sam(Or mapped_mirna_vs_genome_sorted.sam) overlapLenThreshold(12, this value must be less than  16)
	overlapLenThreshold = int(overlapLenThresholdTmp)
	sampleName = '_'.join(os.path.basename(infTmp).split('_')[:-3])
	inf = open(infTmp,'r')
	chrList = []
	chrSeqDic = {}
	for line in inf:
		if line[0] != '@':
			contentTemp = line.strip().split('\t')
			chr = contentTemp[2]
			if 'chr' in chr:
				if chr not in chrList:
					chrList.append(chr)
					#in chrSeqDic[chr], the fist list store the information about "+" forward strand match.
					#and the second one store the information about "-" reverse strand match
					chrSeqDic.update({chr:[[],[]]})
				else:
					pass
				flag = contentTemp[1]
				startPos = int(contentTemp[3])
				seq = contentTemp[9]
				endPos = startPos+len(seq)-1
				mirID = contentTemp[0]
				
				if flag == "0" and len(chrSeqDic[chr][0]) == 0:
					chrSeqDic[chr][0].append([startPos, endPos, seq, [mirID]])
				elif flag == "16" and len(chrSeqDic[chr][1]) == 0:
					chrSeqDic[chr][1].append([startPos, endPos, seq, [mirID]])
				else:
					if flag == "0":
						temp = chrSeqDic[chr][0][-1]
						if startPos >= temp[0] and startPos <= temp[1] and temp[1]-startPos+1 >= overlapLenThreshold:
							mirIDList = temp[3]
							mirIDList.append(mirID)
							startPosNew = temp[0]
							if endPos <= temp[1]:
								endPosNew = temp[1]
								seqNew = temp[2]
							else:
								endPosNew = endPos
								seqNew = temp[2] + seq[temp[1]-startPos+1:]
							chrSeqDic[chr][0][-1] = [startPosNew, endPosNew, seqNew, mirIDList]
						else:
							chrSeqDic[chr][0].append([startPos, endPos, seq,[mirID]])
					elif flag == "16":
						temp = chrSeqDic[chr][1][-1]
						if startPos >= temp[0] and startPos <= temp[1]and temp[1]-startPos+1 >= overlapLenThreshold:
							mirIDList = temp[3]
							mirIDList.append(mirID)
							startPosNew = temp[0]
							if endPos <= temp[1]:
								endPosNew = temp[1]
								seqNew = temp[2]
							else:
								endPosNew = endPos
								seqNew = temp[2] + seq[temp[1]-startPos+1:]
							chrSeqDic[chr][1][-1] = [startPosNew, endPosNew, seqNew, mirIDList]
						else:
							chrSeqDic[chr][1].append([startPos, endPos, seq,[mirID]])
					else:
						pass
			else:
				pass
		else:
			pass
	inf.close()
	with open(infTmp[:-4]+'_clusters.tsv','w') as outf:
		i = 1
		outf.write('miRClusterID\tChr\tStrand\tStart\tEnd\tSequence\tSequenceLenght\tCoutOfReads\tCountOfMembers\tMembers\n')
		for chr in chrList:
			for content in chrSeqDic[chr][0]:
				readCount = 0
				for item in content[3]:
					readCount = readCount + int(item.split('_')[1])
				outf.write('\t'.join([sampleName+':'+'miRCluster_'+str(i)+'_'+str(len(content[2])),
							chr, "+", str(content[0]), str(content[1]), content[2],
							str(len(content[2])), str(readCount), str(len(content[3])),
							','.join(content[3])
							])+ '\n')
				i = i + 1
			for content in chrSeqDic[chr][1]:
				readCount = 0
				for item in content[3]:
					readCount = readCount + int(item.split('_')[1])
				outf.write('\t'.join([sampleName+':'+'miRCluster_'+str(i)+'_'+str(len(content[2])),
							chr, "-", str(content[0]), str(content[1]), content[2],
							str(len(content[2])), str(readCount), str(len(content[3])),
							','.join(content[3])
							])+ '\n')
				i = i + 1