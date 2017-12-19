import os

def quantReads(inputfile, seqDic, readLengthDic, sampleCount, sampleIndex, sampleList, trimmed_collapsed_fa, spikeIn):
	with open(inputfile, 'r') as inf:
		line = inf.readline()
		while line != '':
			line = inf.readline()
			seq = line.strip()
			try:
				seqDic[seq]['quant'][sampleIndex] = seqDic[seq]['quant'][sampleIndex] + 1
			except KeyError:
				if spikeIn:
					seqDic.update({seq:{'quant':[0]*sampleCount, 'annot':[0,'','','','','','','','',''], 'length':len(seq)}})
				else:
					seqDic.update({seq:{'quant':[0]*sampleCount, 'annot':[0,'','','','','','','',''], 'length':len(seq)}})
				seqDic[seq]['quant'][sampleIndex] = 1
			try:
				readLengthDic[len(seq)][sampleIndex] = readLengthDic[len(seq)][sampleIndex] + 1
			except KeyError:
				readLengthDic.update({len(seq):[0]*sampleCount})
				readLengthDic[len(seq)][sampleIndex] = 1
			line = inf.readline()
			line = inf.readline()
			line = inf.readline()
	# Output the collapsed sequences into a file if trimmed_collapsed_fa is True.
	if trimmed_collapsed_fa:
		dir_TMP = os.path.dirname(inputfile)
		with open(os.path.join(dir_TMP, os.path.splitext(sampleList[sampleIndex])[0]+'.trim.collapse.fa'), 'w') as outfTmp:
			t = 1
			#for seq in seqDic.keys():
			#	if seqDic[seq]['quant'][sampleIndex] > 0:
			#		seqName = 'seq'+str(t)+'_'+str(seqDic[seq]['quant'][sampleIndex])
			#		outfTmp.write('>'+seqName+'\n')
			#		outfTmp.write(seq+'\n')
			#		t = t + 1
			decorated = []
			for seq in seqDic.keys():
				if seqDic[seq]['quant'][sampleIndex] > 0:
					decorated.append((seqDic[seq]['quant'][sampleIndex], seq))
					t = t + 1
			decorated.sort(reverse=True)
			for i, tmp in enumerate(decorated):
				outfTmp.write('>seq'+str(i+1)+'_'+str(tmp[0])+'\n')
				outfTmp.write(tmp[1]+'\n')
	# Remove the '*.trim.fastq' to save space
	os.system('rm %s'%(inputfile))