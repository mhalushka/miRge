import os
import math
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment
from scipy import stats
import math
import sys
import subprocess

# add function to calculate A->I editing
def startAndEnd(seq1, seq2):
	# seq1 is the substring of seq2
	start = seq2.index(seq1)
	end = start+len(seq1)-1
	return (start, end)

def DashCount(seq):
	headDashCount = 0
	tailDashCount = 0
	for i in seq:
		if i == '-':
			headDashCount = headDashCount + 1
		else:
			break
	for i in range(len(seq)-1,-1,-1):
		if seq[i] == '-':
			tailDashCount = tailDashCount + 1
		else:
			break
	return (headDashCount, tailDashCount)

def judgeAllign(targetSeqTmpNew, seqTmp):
	mismatchCount = 0
	matchCount = 0
	mismatchLimit = 1
	headShift = 1
	tailShift = 3
	state = True
	headDash_target, tailDash_target = DashCount(targetSeqTmpNew)
	headDash_seqTmp, tailDash_seqTmp = DashCount(seqTmp)
	len1 = len(targetSeqTmpNew) - headDash_target - tailDash_target
	matchLimit = len1 - tailShift - mismatchLimit
	#len2 = len(seqTmp) - headDash_seqTmp - tailDash_seqTmp
	start_pos1 = headDash_target
	end_pos1 = len(targetSeqTmpNew)-headDash_target-1-tailShift
	start_pos2 = headDash_seqTmp
	end_pos2 = len(seqTmp)-tailDash_seqTmp-1
	if start_pos2-start_pos1 > headShift:
		state = False
	else:
		for pos in range(headDash_target, min([end_pos1, end_pos2])+1):
			if seqTmp[pos] == '-':
				pass
			elif targetSeqTmpNew[pos] != seqTmp[pos]:
				mismatchCount = mismatchCount + 1
			else:
				matchCount = matchCount + 1
		if start_pos2-start_pos1 == headShift:
			matchLimitNew = matchLimit - headShift
		else:
			matchLimitNew = matchLimit
		if mismatchCount > mismatchLimit:
			state = False
		if matchCount < matchLimitNew:
			state = False
	return state

def removeDash(seq):
	headDashCount, tailDashCount = DashCount(seq)
	return seq[headDashCount:len(seq)-tailDashCount]

def refineName(name):
	if '.fastq' in name:
		label = name.index('.fastq')
		newName = ''
		for i,item in enumerate(name):
			if i in range(label, label+len('.fastq')):
				pass
			else:
				newName = newName + item
	else:
		newName = name
	return newName

def checkSeqList(seqList, seqDic):
	# Check whether there are at least one sequences belonging to known miRNA. 
	state = False
	for seq in seqList:
		seqNew = removeDash(seq)
		try:
			if seqDic[seqNew]['annot'][1] != '':
				state = True
				break
		except KeyError:
			print 'seq %s and %s does not exit.'%(seqNew, seq)
	return state

def align2TargetSeq(targetSeq, seqList):
	align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
	alignSeqList = []
	stateList = []
	for i in range(len(seqList)):
		# Perform pairwise local alignment. And increase the penalty of gap opening and extending into -10 and -10,
		# so that there will be no gap in the alignned result.
		seq = seqList[i]
		alignTmp = pairwise2.align.localms(targetSeq, seq, 2, -1, -20, -20)
		targetSeqTmpNew = alignTmp[0][0]
		seqTmp = alignTmp[0][1]

		state = judgeAllign(targetSeqTmpNew, seqTmp)
		stateList.append(state)

		if len(alignSeqList) == 0:
			alignSeqList.append(targetSeqTmpNew)
			alignSeqList.append(seqTmp)
		else:
			if targetSeqTmpNew == alignSeqList[0]:
				alignSeqList.append(seqTmp)
			else:
				headAdd1 = targetSeqTmpNew.index(targetSeq)
				tailAdd1 = len(targetSeqTmpNew)-headAdd1-len(targetSeq)
				headAdd2 = alignSeqList[0].index(targetSeq)
				tailAdd2 = len(alignSeqList[0])-headAdd2-len(targetSeq)
				if headAdd1 >= headAdd2 and tailAdd1 >= tailAdd2:
					for j in range(len(alignSeqList)):
						alignSeqList[j] = (headAdd1 - headAdd2)*'-'+alignSeqList[j]+(tailAdd1 - tailAdd2)*'-'
					alignSeqList.append(seqTmp)
				elif headAdd1 >= headAdd2 and tailAdd1 < tailAdd2:
					for j in range(len(alignSeqList)):
						alignSeqList[j] = (headAdd1 - headAdd2)*'-'+alignSeqList[j]
					alignSeqList.append(seqTmp+(tailAdd2 - tailAdd1)*'-')
				elif headAdd1 < headAdd2 and tailAdd1 >= tailAdd2:
					for j in range(len(alignSeqList)):
						alignSeqList[j] = alignSeqList[j]+(tailAdd1 - tailAdd2)*'-'
					alignSeqList.append((headAdd2 - headAdd1)*'-'+seqTmp)
				elif headAdd1 < headAdd2 and tailAdd1 < tailAdd2:
					alignSeqList.append((headAdd2 - headAdd1)*'-'+seqTmp+(tailAdd2 - tailAdd1)*'-')
				else:
					pass
	return (alignSeqList, stateList)

def A2IEditing(targetSeq, seqList, countList, mirName, outfTmp, retainedSeqDic, startBase, endBase):
	p_mismatch = 0.001
	tailShift = 5
	alignSeqListTmp, stateListTmp = align2TargetSeq(targetSeq, seqList)
	startLabel, endLabel = startAndEnd(targetSeq, alignSeqListTmp[0])
	targetSeqAligned = alignSeqListTmp[0]
	targetSeqAlignedHeadCount = DashCount(targetSeqAligned)[0]
	#count = 0
	canonicalSeqCount = 0
	countSumTrue = 0
	seqCountTrue = 0
	a2IPositionCountDic = {}
	positionList = []
	alignSeqListKeptNew = []

	for j, seq in enumerate(alignSeqListTmp[1:]):
		if stateListTmp[j] and removeDash(seq) in retainedSeqDic:
		#if stateListTmp[j]:
			if seqList[j] in targetSeq:
				canonicalSeqCount = canonicalSeqCount + countList[j]
			alignSeqListKeptNew.append(seq)
			seqCountTrue = seqCountTrue + 1
			countSumTrue = countSumTrue + countList[j]
			for i in range(startLabel, endLabel+1-tailShift):
				try:
					if targetSeqAligned[i] == startBase and alignSeqListTmp[j+1][i] == endBase:
						positionNew = i+1-targetSeqAlignedHeadCount
						if positionNew not in positionList:
							positionList.append(positionNew)
							a2IPositionCountDic.update({positionNew:countList[j]})
						else:
							a2IPositionCountDic[positionNew] = a2IPositionCountDic[positionNew] + countList[j]
						#count = count + countList[j]
						#break
				except IndexError:
					pass
	#print 'Canonical_Seq of %s: %s'%(mirName, targetSeq)
	#print 'seqList size is: %d, %d'%(len(seqList), len(countList))
	outfTmp.write('Canonical_Seq of %s: %s\n'%(mirName, targetSeq))
	outfTmp.write('seqList size is: %d, %d\n'%(len(seqList), len(countList)))
	for index, itemTmp in enumerate(alignSeqListTmp):
		if index != 0:
			#print '\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])])
			outfTmp.write('\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])]))
			outfTmp.write('\n')
	#print '****************'
	#print 'retained seqList size is: %d'%(seqCountTrue)
	outfTmp.write('****************\n')
	outfTmp.write('retained seqList size is: %d\n'%(seqCountTrue))
	alignSeqListKept = []
	for index, itemTmp in enumerate(alignSeqListTmp):
		if index != 0:
			if stateListTmp[index-1]:
				alignSeqListKept.append(itemTmp)
				#print '\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])])
				outfTmp.write('\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])]))
				outfTmp.write('\n')
	outfTmp.write('****************\n')
	outfTmp.write('retained sequences after filering are:\n')
	alignSeqListKept2 = []
	for index, itemTmp in enumerate(alignSeqListTmp):
		if index != 0:
			if stateListTmp[index-1] and removeDash(itemTmp) in retainedSeqDic:
				alignSeqListKept2.append(itemTmp)
				#print '\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])])
				outfTmp.write('\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])]))
				outfTmp.write('\n')

	a2IPositionRatioDic = {}
	a2IPositionPvalueDic = {}
	for position in positionList:
		try:
			a2IPositionRatioDic.update({position:float(a2IPositionCountDic[position])/countSumTrue})
		except ZeroDivisionError:
			a2IPositionRatioDic.update({position:0})
		if countSumTrue - a2IPositionCountDic[position] >= 0:
			p_value = stats.binom.cdf(countSumTrue - a2IPositionCountDic[position], countSumTrue, 1-p_mismatch)
		else:
			p_value = 1.0
		a2IPositionPvalueDic.update({position:p_value})
	#try:
	#	aIValue = float(count)/countSumTrue
	#except ZeroDivisionError:
	#	aIValue = 0
	return (alignSeqListKeptNew, positionList, a2IPositionCountDic, a2IPositionRatioDic, a2IPositionPvalueDic, countSumTrue, seqCountTrue, canonicalSeqCount)
	#return (aIValue, countSumTrue, seqCountTrue)
	#return (float(count)/sum(countList), count, seqCountTrue)

def mismatchCountAnalysis(targetSeq, seqList, countList, retainedSeqDic):
	tailShift = 5
	alignSeqListTmp, stateListTmp = align2TargetSeq(targetSeq, seqList)
	#print 'Canonical_Seq of %s: %s'%(mirName, targetSeq)
	#print 'seqList size is: %d, %d'%(len(seqList), len(countList))
	#outfTmp.write('Canonical_Seq of %s: %s\n'%(mirName, targetSeq))
	#outfTmp.write('seqList size is: %d, %d\n'%(len(seqList), len(countList)))
	#for index, itemTmp in enumerate(alignSeqListTmp):
	#	if index != 0:
	#		#print '\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])])
	#		outfTmp.write('\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])]))
	#		outfTmp.write('\n')
	startLabel, endLabel = startAndEnd(targetSeq, alignSeqListTmp[0])
	targetSeqAligned = alignSeqListTmp[0]
	targetSeqAlignedHeadCount = DashCount(targetSeqAligned)[0]
	basePairMismachCountList = []
	for basePair in [('A', 'G'),('A', 'C'),('A', 'T'),('T', 'G'),('T', 'A'),('T', 'C'),('C', 'G'),('C', 'A'),('C', 'T'),('G', 'A'),('G', 'C'),('G', 'T')]:
		startBaseTmp = basePair[0]
		endBaseTmp = basePair[1]
		#startBaseTmp = 'A'
		#endBaseTmp = 'G'
		#count = 0
		countSum = 0
		seqCount = 0
		countSumTrue = 0
		seqCountTrue = 0
		countSumTrue2 = 0
		seqCountTrue2 = 0
		
		mismatchCountDic_raw = {}
		positionList_raw = [] 
		mismatchCountDic = {}
		positionList = []
		mismatchCountDic_filter = {}
		positionList_filter = []

		for j, seq in enumerate(alignSeqListTmp[1:]):
			seqCount = seqCount + 1
			countSum = countSum + countList[j]
			for i in range(startLabel, endLabel+1-tailShift):
				try:
					if targetSeqAligned[i] == startBaseTmp and alignSeqListTmp[j+1][i] == endBaseTmp:
						positionNew = i+1-targetSeqAlignedHeadCount
						if positionNew not in positionList_raw:
							positionList_raw.append(positionNew)
							mismatchCountDic_raw.update({positionNew:countList[j]})
						else:
							mismatchCountDic_raw[positionNew] = mismatchCountDic_raw[positionNew] + countList[j]
						#count = count + countList[j]
						#break
				except IndexError:
					pass
			if stateListTmp[j]:
				seqCountTrue = seqCountTrue + 1
				countSumTrue = countSumTrue + countList[j]
				for i in range(startLabel, endLabel+1-tailShift):
					try:
						if targetSeqAligned[i] == startBaseTmp and alignSeqListTmp[j+1][i] == endBaseTmp:
							positionNew = i+1-targetSeqAlignedHeadCount
							if positionNew not in positionList:
								positionList.append(positionNew)
								mismatchCountDic.update({positionNew:countList[j]})
							else:
								mismatchCountDic[positionNew] = mismatchCountDic[positionNew] + countList[j]
							#count = count + countList[j]
							#break
					except IndexError:
						pass
			if stateListTmp[j] and removeDash(seq) in retainedSeqDic:
				seqCountTrue2 = seqCountTrue2 + 1
				countSumTrue2 = countSumTrue2 + countList[j]
				for i in range(startLabel, endLabel+1-tailShift):
					try:
						if targetSeqAligned[i] == startBaseTmp and alignSeqListTmp[j+1][i] == endBaseTmp:
							positionNew = i+1-targetSeqAlignedHeadCount
							if positionNew not in positionList_filter:
								positionList_filter.append(positionNew)
								mismatchCountDic_filter.update({positionNew:countList[j]})
							else:
								mismatchCountDic_filter[positionNew] = mismatchCountDic_filter[positionNew] + countList[j]
							#count = count + countList[j]
							#break
					except IndexError:
						pass
		basePairMismachCount_raw = 0
		basePairMismachCount = 0
		basePairMismachCount_filter = 0
		for key in mismatchCountDic_raw.keys():
			basePairMismachCount_raw = basePairMismachCount_raw + mismatchCountDic_raw[key]
		for key in mismatchCountDic.keys():
			basePairMismachCount = basePairMismachCount + mismatchCountDic[key]
		for key in mismatchCountDic_filter.keys():
			basePairMismachCount_filter = basePairMismachCount_filter + mismatchCountDic_filter[key]
		basePairMismachCountList.append((basePairMismachCount_raw, basePairMismachCount, basePairMismachCount_filter))
	return basePairMismachCountList

def calcEntropy(inputList):
	sum1 = sum(inputList)
	entropy = 0
	for i in range(len(inputList)):
		if inputList[i] > 1:
			freq = float(inputList[i])/sum1
			entropy = entropy + -1*freq*math.log(freq, 2)
	return entropy

def writeDataToCSV(outputdir, annotNameList, sampleList, isomirDiff, a_to_i, logDic, seqDic, mirDic, mirNameSeqDic, mirMergedNameDic, bowtieBinary, genome_index, numCPU, phred64, removedMiRNAList, spikeIn, gff_output, isomiRContentDic, miRNA_database, paraent_dir):
	a2IPercentageCutoff = 2.0
	mappedFile = os.path.join(outputdir, 'mapped.csv')
	isomirFile = os.path.join(outputdir, 'isomirs.csv')
	isomirSampleFile = os.path.join(outputdir, 'isomirs.samples.csv')
	unmappedFile = os.path.join(outputdir, 'unmapped.csv')
	mirRPMFile = os.path.join(outputdir, 'miR.RPM.csv')
	mirCountsFile = os.path.join(outputdir, 'miR.Counts.csv')
	mirFile = os.path.join(outputdir, 'miR.Counts.csv')
	a2IEditingFileTmp1 = os.path.join(outputdir, 'a2IEditing.report.tmp1.csv')
	a2IEditingFileTmp2 = os.path.join(outputdir, 'a2IEditing.report.tmp2.csv')
	a2IEditingFileTmp3 = os.path.join(outputdir, 'a2IEditing.report.tmp3.csv')
	a2IEditingFile = os.path.join(outputdir, 'a2IEditing.report.csv')
	a2IEditingFileTrans = os.path.join(outputdir, 'a2IEditing.report.newform.csv')
	a2IEditingDetailFile = os.path.join(outputdir, 'a2IEditing.detail.txt')
	mismatchCountFile = os.path.join(outputdir, 'mismatchCount.csv')

	with open(mappedFile, 'w') as outf:
		outf.write('uniqueSequence,annotFlag,'+','.join(annotNameList)+','+','.join(sampleList)+'\n')
		isomirDic = {}
		for seqKey in seqDic.keys():
			entry = [seqKey]
			isomirVals = [seqKey]
			if seqDic[seqKey]['annot'][0] >0:
				isomir = seqDic[seqKey]['annot'][8]
				mirna = seqDic[seqKey]['annot'][1]
				if isomir != '' or mirna != '':
					if isomir != '':
						key = isomir
						key2 = 'isomirs'
					else:
						key = mirna
						key2 = 'mirnas'
					# get rid of the SNP annotation for isoMIRs
					if '.SNP' in key:
						key = key.split('.SNP')[0]
					if key not in isomirDic.keys():
						dicTmp = {'mirnas':{}}
						dicTmp.update({'isomirs':{}})
						isomirDic.update({key:dicTmp})
					isomirDic[key][key2][seqKey] = []
					for i in range(len(sampleList)):
						isomirDic[key][key2][seqKey].append(seqDic[seqKey]['quant'][i])
				if spikeIn:
					for i in range(10):
						entry.append(seqDic[seqKey]['annot'][i])
				else:
					for i in range(9):
						entry.append(seqDic[seqKey]['annot'][i])
				for i in range(len(sampleList)):
					readCount = seqDic[seqKey]['quant'][i]
					entry.append(readCount)
					#isomirDic[key][key2][seqKey].append(readCount)
				outf.write(','.join([str(item) for item in entry]))
				outf.write('\n')

	if gff_output:
		# For each sample, calculate the RPM value, 
		for i in range(len(sampleList)):
			#sampleName = '.'.join(sampleList[i].split('.')[:-1])
			sampleName = os.path.splitext(sampleList[i])[0]
			gffFileTmp = os.path.join(outputdir, sampleName+'_isomiRs.gff')
			with open(gffFileTmp, 'w') as outf:
				outf.write('# GFF3 adapted for miRNA sequencing data\n## VERSION 0.0.1\n## source-ontology: ')
				if miRNA_database == 'miRBase':
					source = 'miRBase21'
					outf.write('miRBase21\n')
				else:
					source = miRNA_database
					outf.write('%s\n'%(miRNA_database))
				outf.write('## COLDATA: %s\n'%(sampleName))
				for seqKey in isomiRContentDic.keys():
					#RPM = 1000000.0*seqDic[seqKey]['quant'][i]/logDic['quantStats'][i]['mirnaReadsFiltered']
					redCount = seqDic[seqKey]['quant'][i]
					#if RPM > 0:
					if redCount >= 1:
						#isomiRContentDic[seqKey]['expression'] = '%.2f'%(RPM)
						isomiRContentDic[seqKey]['expression'] = str(redCount)
						outf.write('\t'.join([isomiRContentDic[seqKey]['miRName'], source, isomiRContentDic[seqKey]['type'], isomiRContentDic[seqKey]['pre_start'], isomiRContentDic[seqKey]['pre_end'], '.', isomiRContentDic[seqKey]['strand'], '.']))
						outf.write('\t')
						outf.write(';'.join(['Read '+seqKey, ' UID '+isomiRContentDic[seqKey]['uid'], ' Name '+isomiRContentDic[seqKey]['miRName'], ' Parent '+isomiRContentDic[seqKey]['preMiRName'], ' Variant '+isomiRContentDic[seqKey]['variant'], ' Cigar '+isomiRContentDic[seqKey]['cigar'], ' Expression '+isomiRContentDic[seqKey]['expression'], ' Filter '+isomiRContentDic[seqKey]['filter']]))
						outf.write('\n')
	if isomirDiff:
		outf1 = open(isomirFile, 'w')
		outf2 = open(isomirSampleFile, 'w')
		outf1.write('miRNA,sequence')
		outf2.write('miRNA')
		for i in range(len(sampleList)):
			outf1.write(','+sampleList[i])
			outf2.write(','+sampleList[i]+' isomir+miRNA Entropy')
			outf2.write(','+sampleList[i]+' Canonical Sequence')
			outf2.write(','+sampleList[i]+' Canonical RPM')
			outf2.write(','+sampleList[i]+' Top Isomir RPM')
		outf1.write(',Entropy\n')
		outf2.write('\n')
		#print '&&&&&&&&&'
		#print len(isomirDic)
		#print isomirDic.keys()
		#print isomirDic.keys()[0]
		#print isomirDic[isomirDic.keys()[0]]
		for miRNA in isomirDic.keys():
			sampleIsomirs= {}
			samplemiRNAs = []
			for i in range(len(sampleList)):
				sampleIsomirs.update({i:[]})
				samplemiRNAs.append(0)
			for miRNASeq in isomirDic[miRNA]['mirnas'].keys():
				sampleArray = isomirDic[miRNA]['mirnas'][miRNASeq]
				for i in range(len(sampleArray)):
					samplemiRNAs[i] = samplemiRNAs[i] + sampleArray[i]
			for miRNASeq in isomirDic[miRNA]['isomirs'].keys():
				entry = [miRNA]
				entry.append(miRNASeq)
				sampleArray = isomirDic[miRNA]['isomirs'][miRNASeq]
				for i in range(len(sampleArray)):
					sampleIsomirs[i].append(sampleArray[i])
				entropy = calcEntropy(sampleArray)
				maxEntropy = math.log(len(sampleArray), 2)
				if maxEntropy == 0:
					entropy = "NA"
				else:
					entropy = str(entropy/maxEntropy)
				for i in range(len(sampleArray)):
					sampleArray[i] = str(sampleArray[i]*1000000.0/logDic['quantStats'][i]['mirnaReadsFiltered'])
				entry = entry + sampleArray
				entry.append(str(entropy))
				outf1.write(','.join(entry))
				outf1.write('\n')

			isomirOut = [miRNA]
			for sampleLane in range(len(sampleList)):
				rpmFactor = 1000000.0/logDic['quantStats'][sampleLane]['mirnaReadsFiltered']
				sampleEntropy = calcEntropy(sampleIsomirs[sampleLane])
				#print '!!!!!!!!!!'
				#print miRNA
				#print sampleIsomirs[sampleLane]
				#print sampleEntropy
				if len(sampleIsomirs[sampleLane]) > 0:
					topIsomir = max(sampleIsomirs[sampleLane])*rpmFactor
					isomirSum = sum(sampleIsomirs[sampleLane])*rpmFactor
					sampleIsomirs[sampleLane].append(samplemiRNAs[sampleLane])
					sampleEntropyWithmiRNA = calcEntropy(sampleIsomirs[sampleLane])
					miRNARPM = samplemiRNAs[sampleLane]*rpmFactor
					
					maxEntropy = len(sampleIsomirs[sampleLane])
					if maxEntropy > 1:
						sampleEntropyWithmiRNA = str(sampleEntropyWithmiRNA/(math.log(maxEntropy, 2)))
					else:
						sampleEntropyWithmiRNA = 'NA'
					isomirOut.append(sampleEntropyWithmiRNA)
					
					combined = miRNARPM + isomirSum
					if combined > 0:
						isomirOut.append(str(100.0*miRNARPM/combined))
					else:
						isomirOut.append('NA')
					
					isomirOut.append(str(miRNARPM))
					isomirOut.append(str(topIsomir))
					outf2.write(','.join(isomirOut))
					outf2.write('\n')
		outf1.close()
		outf2.close()

	with open(unmappedFile, 'w') as outf:
		outf.write('uniqueSequence,annotFlag,'+','.join(annotNameList))
		for i in range(len(sampleList)):
			outf.write(','+sampleList[i])
		outf.write('\n')
		for seqKey in seqDic.keys():
			if seqDic[seqKey]['annot'][0] == 0:
				outf.write(seqKey)
				if spikeIn:
					for i in range(10):
						outf.write(','+str(seqDic[seqKey]['annot'][i]))
				else:
					for i in range(9):
						outf.write(','+str(seqDic[seqKey]['annot'][i]))
				for i in range(len(sampleList)):
					outf.write(','+str(seqDic[seqKey]['quant'][i]))
				outf.write('\n')

	tmp = mirDic.keys()
	tmp.sort()
	with open(mirCountsFile, 'w') as outf:
		outf.write('miRNA')
		for i in range(len(sampleList)):
			outf.write(','+sampleList[i])
		outf.write('\n')
		outf.write('miRNAtotal')
		for i in range(len(sampleList)):
			outf.write(','+str(logDic['quantStats'][i]['mirnaReadsFiltered']))
		outf.write('\n')
		for mirKey in tmp:
			outf.write(mirKey)
			for i in range(len(sampleList)):
				outf.write(','+str(mirDic[mirKey]['quant'][i]))
			outf.write('\n')

	with open(mirRPMFile, 'w') as outf:
		outf.write('miRNA')
		for i in range(len(sampleList)):
			outf.write(','+sampleList[i])
		outf.write('\n')
		for mirKey in tmp:
			outf.write(mirKey)
			for i in range(len(sampleList)):
				if logDic['quantStats'][i]['mirnaReadsFiltered'] > 0:
					outf.write(','+str(1000000.0*mirDic[mirKey]['quant'][i]/logDic['quantStats'][i]['mirnaReadsFiltered']))
				else:
					outf.write(',0')
			outf.write('\n')

	if a_to_i:
		startBase = 'A'
		endBase = 'G'
		mirNameSeqDicTmp = {}
		for seqKey in seqDic.keys():
			if seqDic[seqKey]['annot'][1] != '' or seqDic[seqKey]['annot'][8] != '':
				for item in [seqDic[seqKey]['annot'][1], seqDic[seqKey]['annot'][8]]:
					if item != '':
						miRNameTmp = item
						break
				# Substitute the miRNA's name with the merged name, if it was merged.
				if miRNameTmp in mirMergedNameDic.keys():
					miRName = mirMergedNameDic[miRNameTmp]
				else:
					miRName = miRNameTmp
				# Keep reads if they are exact miRNA, even it's RPM < 1, or at least one of the samples whose RPM >= 1.
				if seqDic[seqKey]['annot'][1] != '':
					try:
						mirNameSeqDicTmp[miRName].append(seqKey)
					except KeyError:
						mirNameSeqDicTmp.update({miRName:[seqKey]})
				else:
					if any([1000000.0*seqDic[seqKey]['quant'][i]/logDic['quantStats'][i]['mirnaReadsFiltered'] >= 1 for i in range(len(sampleList))]):
						try:
							mirNameSeqDicTmp[miRName].append(seqKey)
						except KeyError:
							mirNameSeqDicTmp.update({miRName:[seqKey]})
					else:
						pass
		
		# In order to avoid cross-mapping, the reads whose last 2 nucluetides will be ignored when mapping against genome by Bowtie.
		# up to one mismatch and unique best hits.
		
		with open(os.path.join(outputdir, 'SeqToMap.fasta'), 'w') as outf:
			for mirName in mirNameSeqDicTmp.keys():
				for seqTmp in mirNameSeqDicTmp[mirName]:
					outf.write('>'+seqTmp+'\n'+seqTmp+'\n')
		
		#bwtCmd = bowtieBinary+"bowtie --threads "+numCPU+' '
		bwtCmd = os.path.join(bowtieBinary, 'bowtie')+' --threads '+numCPU+' '
		if(phred64):
			bwtCmd = bwtCmd + ' --phred64-quals '
		bwtCmdLines = bwtCmd+genome_index+' -n 1 -f -a -3 2 '+os.path.join(outputdir, 'SeqToMap.fasta')+' 1> '+os.path.join(outputdir, 'SeqToMap.sam')+' 2> '+os.path.join(outputdir, 'SeqToMap.log')
		os.system(bwtCmdLines)
		retainedSeqDic = {}
		retainedSeqContentDicTmp = {}
		with open(os.path.join(outputdir, 'SeqToMap.sam'), 'r') as inf:
			for line in inf:
				if '@' not in line:
					content = line[:-1].split('\t')
					try:
						retainedSeqContentDicTmp[content[0]].append(content[-1].count(':'))
					except KeyError:
						retainedSeqContentDicTmp.update({content[0]:[content[-1].count(':')]})
		
		os.system('rm %s %s %s'%(os.path.join(outputdir, 'SeqToMap.sam'), os.path.join(outputdir, 'SeqToMap.log'), os.path.join(outputdir, 'SeqToMap.fasta')))
		for seq_tmp in retainedSeqContentDicTmp.keys():
			kept = False
			content = retainedSeqContentDicTmp[seq_tmp]
			if len(content) == 1:
				kept = True
			else:
				minimum = min(content)
				if len([subitem for subitem in content if subitem == minimum]) == 1:
					kept = True
			if kept:
				retainedSeqDic.update({seq_tmp:True})
		
		#print 'mirNameSeqDicTmp length is: %d'%(len(mirNameSeqDicTmp))
		mirNamePositionValueDict = {}
		mirNamePositionList = []
		
		outfTmp = open(a2IEditingDetailFile, 'w')
		for mirName in mirNameSeqDicTmp.keys():
			targetSeq = mirNameSeqDic[mirName]
			seqList = mirNameSeqDicTmp[mirName]
			for i in range(len(sampleList)):
				#if mirDic[mirName]['quant'][i] > 0:
				if 1000000.0*mirDic[mirName]['quant'][i]/logDic['quantStats'][i]['mirnaReadsFiltered'] >= 1:
					countListTmp = []
					seqListTmp = []
					for seqTmp in seqList:
						if 1000000.0*seqDic[seqTmp]['quant'][i]/logDic['quantStats'][i]['mirnaReadsFiltered'] >= 1 or (seqDic[seqTmp]['annot'][1] != '' and seqDic[seqTmp]['quant'][i] > 0):
							seqListTmp.append(seqTmp)
							countListTmp.append(seqDic[seqTmp]['quant'][i])
					if len(countListTmp) >1:
						if min(countListTmp) > 0:
							alignSeqListKept, positionList, a2IPositionCountDic, a2IPositionRatioDic, a2IPositionPvalueDic, countTrue, seqCountTrue, canonicalSeqCount = A2IEditing(targetSeq, seqListTmp, countListTmp, mirName, outfTmp, retainedSeqDic, startBase, endBase)
							for position in positionList:
								mirNamePosition = ':'.join([mirName, str(position)])
								if mirNamePosition not in mirNamePositionValueDict.keys():
									mirNamePositionList.append(mirNamePosition)
									mirNamePositionValueDict.update({mirNamePosition:[[] for t in range(len(sampleList))]})
									mirNamePositionValueDict[mirNamePosition][i] = mirNamePositionValueDict[mirNamePosition][i] + [alignSeqListKept, str(sum(countListTmp)), str(len([item1 for item1 in countListTmp if item1 >0])), str(countTrue), str(seqCountTrue), str(canonicalSeqCount), str(a2IPositionCountDic[position]), a2IPositionRatioDic[position], a2IPositionPvalueDic[position]]
								else:
									mirNamePositionValueDict[mirNamePosition][i] = mirNamePositionValueDict[mirNamePosition][i] + [alignSeqListKept, str(sum(countListTmp)), str(len([item1 for item1 in countListTmp if item1 >0])), str(countTrue), str(seqCountTrue), str(canonicalSeqCount), str(a2IPositionCountDic[position]), a2IPositionRatioDic[position], a2IPositionPvalueDic[position]]
						else:
							pass
					else:
						pass
				else:
					pass
		outfTmp.close()
		
		basePairMismachCountList_total = []
		for i in range(len(sampleList)):
			basePairMismachCountListTmp = [[0, 0,0] for j in range(12)]
			for mirName in mirNameSeqDicTmp.keys():
				targetSeq = mirNameSeqDic[mirName]
				seqList = mirNameSeqDicTmp[mirName]
				if 1000000.0*mirDic[mirName]['quant'][i]/logDic['quantStats'][i]['mirnaReadsFiltered'] >= 1:
					countListTmp = []
					seqListTmp = []
					for seqTmp in seqList:
						if 1000000.0*seqDic[seqTmp]['quant'][i]/logDic['quantStats'][i]['mirnaReadsFiltered'] >= 1 or (seqDic[seqTmp]['annot'][1] != '' and seqDic[seqTmp]['quant'][i] > 0):
							seqListTmp.append(seqTmp)
							countListTmp.append(seqDic[seqTmp]['quant'][i])
					if len(countListTmp) >1:
						if min(countListTmp) > 0:
							basePairMismachCountList = mismatchCountAnalysis(targetSeq, seqListTmp, countListTmp, retainedSeqDic)
							for indexTmp, s in enumerate(basePairMismachCountList):
								basePairMismachCountListTmp[indexTmp][0] = basePairMismachCountListTmp[indexTmp][0] + s[0]
								basePairMismachCountListTmp[indexTmp][1] = basePairMismachCountListTmp[indexTmp][1] + s[1]
								basePairMismachCountListTmp[indexTmp][2] = basePairMismachCountListTmp[indexTmp][2] + s[2]
						else:
							pass
					else:
						pass
				else:
					pass
			basePairMismachCountList_total.append(basePairMismachCountListTmp)

		# P-value required after multiple testing correction. Here, Benjamini-Hochberg correction is performed.
		for i in range(len(sampleList)):
			#mirNamePositionPicked = []
			#mirNamePositionPickedDic = {}
			p_value_list = []
			for mirNamePosition in mirNamePositionList:
				if len(mirNamePositionValueDict[mirNamePosition][i]) != 0:
					p_value_list.append([mirNamePositionValueDict[mirNamePosition][i][8], mirNamePosition])
			p_value_list.sort()
			for index_tmp, p_value in enumerate(p_value_list):
				adjusted_p_value = p_value[0]*len(p_value_list)/(index_tmp+1)
				mirNamePositionValueDict[p_value[1]][i].append(adjusted_p_value)
		# Write the a2IEditing file.
		
		#print mirNameSeqDicTmp['hsa-miR-376a-3p']
		#if 'hsa-miR-376a-3p' in mirNamePositionValueDict.keys():
		#	print 'hsa-miR-376a-3p: %s'%(mirNamePositionValueDict['hsa-miR-376a-3p'])
		#else:
		#	print 'hsa-miR-376a-3p can not be found'
		
		with open(a2IEditingFileTmp1, 'w') as outf:
			#outf.write('miRNA:Position,miRNA sequence')
			outf.write('miRNA,A-to-I position in the miRNA,miRNA sequence')
			for i in range(len(sampleList)):
				name_new = refineName(sampleList[i])
				#outf.write(','+','.join([sampleList[i]+'.readCount', sampleList[i]+'.sequenceTypeCount', sampleList[i]+'.readCount.kept', sampleList[i]+'.sequenceTypeCount.kept', sampleList[i]+'.readCount.spot', sampleList[i]+'.AtoI.percentage', sampleList[i]+'.AtoI.raw.pValue', sampleList[i]+'.AtoI.adjusted.pValue']))
				outf.write(','+','.join([name_new+'.readCount', name_new+'.readCount.canonical', name_new+'.RPM.canonical', name_new+'.readCount.mismatch', name_new+'.RPM.mismatch', name_new+'.AtoI.percentage', name_new+'.AtoI.adjusted.pValue']))
			outf.write('\n')
			for mirNamePosition in mirNamePositionList:
				mirName = mirNamePosition.split(':')[0].strip()
				position = mirNamePosition.split(':')[1].strip()
				targetSeq = mirNameSeqDic[mirName]
				a2Icontent = mirNamePositionValueDict[mirNamePosition]
				#outf.write(','.join([mirNamePosition, targetSeq]))
				outf.write(','.join([mirName, position, targetSeq]))
				for i in range(len(sampleList)):
					outf.write(',')
					if len(a2Icontent[i]) != 0:
						# In the sequenceTypeCount.kept, at least one of the sequence is known miRNA.
						alignSeqKept = a2Icontent[i][0]
						try:
							if checkSeqList(alignSeqKept, seqDic):
								#outf.write(','.join(a2Icontent[i][1:6]))
								#outf.write(',%.2f%%'%(a2Icontent[i][6]*100))
								#outf.write(',%.2E'%(a2Icontent[i][7]))
								#outf.write(',%.2E'%(a2Icontent[i][8]))
								outf.write('%s'%(a2Icontent[i][3]))
								outf.write(',%s'%(a2Icontent[i][5]))
								outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][5])/logDic['quantStats'][i]['mirnaReadsFiltered']))
								outf.write(',%s'%(a2Icontent[i][6]))
								outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][6])/logDic['quantStats'][i]['mirnaReadsFiltered']))
								outf.write(',%.2f%%'%(a2Icontent[i][7]*100))
								if a2Icontent[i][9] <= 0.05:
									outf.write(',%.2E'%(a2Icontent[i][9]))
								else:
									outf.write(',NS')
							else:
								#outf.write(','.join(a2Icontent[i][1:6]))
								#outf.write(',NE,NE,NE')
								outf.write('%s'%(a2Icontent[i][3]))
								outf.write(',%s'%(a2Icontent[i][5]))
								outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][5])/logDic['quantStats'][i]['mirnaReadsFiltered']))
								outf.write(',%s'%(a2Icontent[i][6]))
								outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][6])/logDic['quantStats'][i]['mirnaReadsFiltered']))
								outf.write(',NE,NE')
						except KeyError:
							print 'Error happens at %s: %s'%(mirName, '\t'.join(alignSeqKept))
							outf.write('%s'%(a2Icontent[i][3]))
							outf.write(',%s'%(a2Icontent[i][5]))
							outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][5])/logDic['quantStats'][i]['mirnaReadsFiltered']))
							outf.write(',%s'%(a2Icontent[i][6]))
							outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][6])/logDic['quantStats'][i]['mirnaReadsFiltered']))
							outf.write(',NE,NE')
					else:
						outf.write(','.join(['NE','NE','NE','NE','NE','NE','NE']))
				outf.write('\n')
		# Remove some unimportant rows in *.a2IEditing.report.tmp1.csv
		outf = open(a2IEditingFileTmp2, 'w')
		with open(a2IEditingFileTmp1, 'r') as inf:
			line = inf.readline()
			outf.write(line)
			line = inf.readline()
			while line != '':
				content = line.strip().split(',')
				portion = ((len(content)-9)/7)+1
				keptState = False
				for i in range(portion):
					try:
						tmpVale = float(content[9+i*7])
						keptState = True
						break
					except ValueError:
						pass
				if keptState:
					outf.write(line)
				line = inf.readline()
		outf.close()
		# Sort the report file based on the first column.
		os.system("(head -n 1 %s && tail -n +2 %s | sort -t',' -k1,1 -k2,2n) > %s"%(a2IEditingFileTmp2, a2IEditingFileTmp2, a2IEditingFileTmp3))
		# Remove the miRNAs if:
		# 1) their canonical sequences' RPM < 1 across all of the samples
		# 2) miRNAs are located at repetitive element region.
		# 3) the correponding 1 nucluetide mismathed sequence can be alligned to no less than 1 locations on the genome with trimming the last two nucluetides at 3' in the setting of Bowtie.
		miRNAPositionDic = {}
		with open(os.path.join(outputdir, 'SeqToJudge.fasta'), 'w') as outf:
			with open(a2IEditingFileTmp3, 'r') as inf:
				line = inf.readline()
				line = inf.readline()
				while line != '':
					content = line.strip().split(',')
					miRNAPositionName = ':'.join([content[0], content[1]])
					mismathedSeq = ''
					for index, necleotide in enumerate(content[2]):
						if index != int(content[1])-1:
							mismathedSeq = mismathedSeq + necleotide
						else:
							mismathedSeq = mismathedSeq + endBase
					canonicalRPMList = []
					for i in range(((len(content)-9)/7)+1):
						try :
							canonicalRPMList.append(float(content[5+i*7]))
						except ValueError:
							canonicalRPMList.append(0)
					canonicalRPMState = not any([item >= 1 for item in canonicalRPMList])
					dicTmp = {}
					dicTmp.update({'mismathedSeq':mismathedSeq})
					dicTmp.update({'RPMRemoveState':canonicalRPMState})
					dicTmp.update({'ReElementRemoveState':content[0] in removedMiRNAList})
					dicTmp.update({'SeqRemoveState':False})
					miRNAPositionDic.update({miRNAPositionName:dicTmp})
					outf.write('>'+mismathedSeq+'\n'+mismathedSeq+'\n')
					line = inf.readline()
		#bwtCmd = bowtieBinary+"bowtie --threads "+numCPU+' '
		if(phred64):
			bwtCmd = bwtCmd + ' --phred64-quals '
		bwtCmdLines = bwtCmd+genome_index+' -n 0 -f -a -3 2 '+os.path.join(outputdir, 'SeqToJudge.fasta')+' 1> '+os.path.join(outputdir, 'SeqToJudge.sam')+' 2> '+os.path.join(outputdir, 'SeqToJudge.log')
		os.system(bwtCmdLines)
		
		removedSeqList = []
		with open(os.path.join(outputdir, 'SeqToJudge.sam'), 'r') as inf:
			for line in inf:
				if '@' not in line:
					content = line[:-1].split('\t')
					if content[0] not in removedSeqList:
						removedSeqList.append(content[0])
		
		os.system('rm %s %s %s'%(os.path.join(outputdir, 'SeqToJudge.fasta'), os.path.join(outputdir, 'SeqToJudge.sam'), os.path.join(outputdir, 'SeqToJudge.log')))
		for key in miRNAPositionDic.keys():
			if miRNAPositionDic[key]['mismathedSeq'] in removedSeqList:
				miRNAPositionDic[key]['SeqRemoveState'] = True
		
		with open(a2IEditingFile, 'w') as outf:
			with open(a2IEditingFileTmp3, 'r') as inf:
				line = inf.readline()
				outf.write(line)
				line = inf.readline()
				while line != '':
					miRNAPositionName = ':'.join([line.strip().split(',')[0], line.strip().split(',')[1]])
					if miRNAPositionDic[miRNAPositionName]['RPMRemoveState'] or miRNAPositionDic[miRNAPositionName]['ReElementRemoveState'] or miRNAPositionDic[miRNAPositionName]['SeqRemoveState']:
						pass
					else:
						outf.write(line)
					line = inf.readline()

		# Remove the *.a2IEditing.report.tmp1.csv *.a2IEditing.report.tmp2.csv *.a2IEditing.report.tmp3.csv
		os.system('rm %s %s %s'%(a2IEditingFileTmp1, a2IEditingFileTmp2, a2IEditingFileTmp3))
		with open(mismatchCountFile, 'w') as outf:
			f1 = ['>'.join(item)+'_raw' for item in [('A', 'G'),('A', 'C'),('A', 'T'),('T', 'G'),('T', 'A'),('T', 'C'),('C', 'G'),('C', 'A'),('C', 'T'),('G', 'A'),('G', 'C'),('G', 'T')]]
			f2 = ['>'.join(item)+'' for item in [('A', 'G'),('A', 'C'),('A', 'T'),('T', 'G'),('T', 'A'),('T', 'C'),('C', 'G'),('C', 'A'),('C', 'T'),('G', 'A'),('G', 'C'),('G', 'T')]]
			f3 = ['>'.join(item)+'_filtered' for item in [('A', 'G'),('A', 'C'),('A', 'T'),('T', 'G'),('T', 'A'),('T', 'C'),('C', 'G'),('C', 'A'),('C', 'T'),('G', 'A'),('G', 'C'),('G', 'T')]]
			outf.write('sample,'+','.join(f1)+','+','.join(f2)+','+','.join(f3)+'\n')
			for i in range(len(sampleList)):
				outf.write(sampleList[i]+',')
				s1 = [str(item[0]) for item in basePairMismachCountList_total[i]]
				s2 = [str(item[1]) for item in basePairMismachCountList_total[i]]
				s3 = [str(item[2]) for item in basePairMismachCountList_total[i]]
				outf.write(','.join(s1))
				outf.write(',')
				outf.write(','.join(s2))
				outf.write(',')
				outf.write(','.join(s3))
				outf.write('\n')
		
		# Remove the mismatchCountFile
		os.system('rm %s'%(mismatchCountFile))
		
		# Transform the format of a2IEditing.report.csv in order to plot heatmap.
		
		sampleDicTmp = {}
		sampleDicTmp.update({'SRR837842':'Colon 1'})
		sampleDicTmp.update({'SRR837839':'Colon 2'})
		sampleDicTmp.update({'SRR5127219':'Colon cell'})
		sampleDicTmp.update({'SRR1646473':'Colon cancer 1'})
		sampleDicTmp.update({'SRR1646493':'Colon cancer 2'})
		sampleDicTmp.update({'SRR1917324':'DKO1'})
		sampleDicTmp.update({'SRR1917336':'DLD1'})
		sampleDicTmp.update({'SRR1917329':'DKS8'})
		sampleDicTmp.update({'SRR567638':'Placenta 2'})
		
		with open(a2IEditingFileTrans, 'w') as outf:
			outf.write('miRNA:position,sample,A-to-I percentage,log2RPM\n')
			sampleListTmp = []
			miRNAListTmp = []
			miRNAContentListTmp = []
			miRNAContentDicTmp = {}
			with open(a2IEditingFile, "r") as inf:
				line = inf.readline()
				for item in line.strip().split(','):
					if '.AtoI.percentage' in item and item.split('.')[0] not in sampleListTmp:
						sampleListTmp.append(item.split('.')[0] )
				line = inf.readline()
				while line != '':
					content = line.strip().split(',')
					miRNANameNew = ":".join(content[:2])
					tmp = content[3:]
					contentTmpList = []
					for item in [tmp[i:i + 7] for i in xrange(0, len(tmp), 7)]:
						if item[6] == 'NE' or item[6] == 'NS':
							# (percentage, RPM)
							contentTmp = ('NA', 'NA')
						else:
							contentTmp = (item[5][:-1], str(math.log(float(item[4]),2)))
						contentTmpList.append(contentTmp)
					# The miRNA will be kept if the percentage of at least one sample is larger than 5.0%
					state = False
					for item in contentTmpList:
						try:
							if float(item[0]) >= a2IPercentageCutoff:
								state = True
								break
						except ValueError:
							pass
					if state:
						miRNAListTmp.append(miRNANameNew)
						miRNAContentListTmp.append(contentTmpList)
						miRNAContentDicTmp.update({miRNANameNew:contentTmpList})
					line = inf.readline()
			miRNAList2 = []
			for index, item in enumerate(miRNAListTmp):
				count = 0
				for item2 in miRNAContentListTmp[index]:
					if item2[0] != 'NA':
						count = count + 1
				miRNAList2.append((count, item))
			miRNAList2.sort(reverse=True)

			for index1, sample in enumerate(sampleListTmp):
				for tmp in miRNAList2:
					miRNA = tmp[1]
					if sample in sampleDicTmp.keys():
						outf.write(miRNA+','+sampleDicTmp[sample]+','+','.join(miRNAContentDicTmp[miRNA][index1])+'\n')
					else:
						outf.write(miRNA+','+sample+','+','.join(miRNAContentDicTmp[miRNA][index1])+'\n')

		if len(miRNAList2) <= 5:
			print 'The number of A-to-I editing sites for is less than 10 so that no heatmap is drawn.'
		else:
			#PopenTmp2 = subprocess.Popen(['which', 'miRge2.0.py'], stdout=subprocess.PIPE)
			#result2 = PopenTmp2.communicate()[0]
			RscriptDir = os.path.join(paraent_dir, 'rSciprts', 'A-to-I_plot.R')
			outA2Ipdf = os.path.join(outputdir, 'a-to-I.heatmap.pdf')
			os.system('Rscript %s %s %s'%(RscriptDir, a2IEditingFileTrans, outA2Ipdf))
