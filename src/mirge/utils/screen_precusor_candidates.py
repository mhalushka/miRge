import os
import sys
import cPickle
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mirge.classes.readPrecusor import ReadPrecusor

def chunks(arr, n):
	return [arr[i:i+n] for i in range(0, len(arr), n)]

def isset(v): 
	try: 
		type (eval(v)) 
	except: 
		return 0 
	else: 
		return 1

def selcteOptimalPrecusor_old(precusorFeatureList, precusorList, clusterSeqType):
	print precusorFeatureList
	print clusterSeqType
	print precusorList
	#To do, select the precusor based on mFE or precusorOrignialMFE. Right now, I choose precusorOrignialMFE
	if clusterSeqType == 'clusterSeq':
		bindingCountThreshold = 19
	elif clusterSeqType == 'stableClusterSeq':
		bindingCountThreshold = 19-1-3
	else:
		bindingCountThreshold = 19
	if len(precusorFeatureList) != len(precusorList):
		print len(precusorFeatureList)
		print len(precusorList)
		print precusorFeatureList
		print precusorList
	retainedIndex = []
	nullIndex = []
	for i in range(len(precusorFeatureList)):
		#print precusorFeatureList[i]
		if not (None in precusorFeatureList[i]):
			armTypeOverlenDic, hairpinCount, bindingCount, percent, mFE, precusorOrignialMFE = precusorFeatureList[i]
			#print hairpinCount, bindingCount, percent, armTypeOverlenDic
			if hairpinCount > 0 and bindingCount >= bindingCountThreshold and percent >= 0.6:
				if 'arm5' in armTypeOverlenDic.keys() and 'arm3' in armTypeOverlenDic.keys():
					if armTypeOverlenDic['arm5'][0] <= 5 and armTypeOverlenDic['arm3'][0] == 0:
						retainedIndex.append(i)
				elif 'arm5' in armTypeOverlenDic.keys():
					if armTypeOverlenDic['arm5'][0] <= 5:
						retainedIndex.append(i)
				elif 'arm3' in armTypeOverlenDic.keys():
					if armTypeOverlenDic['arm3'][0] == 0:
						retainedIndex.append(i)
				else:
					pass
		else:
			nullIndex.append(i)
		#print retainedIndex
	#print retainedIndex
	if len(retainedIndex) >= 1:
		#To do, select the precusor based on mFE or precusorOrignialMFE. Right now, I choose precusorOrignialMFE
		tmpList = [[precusorFeatureList[j][5], j] for j in retainedIndex]
		tmpList.sort()
		#print tmpList
		outputIndex = tmpList[0][1]
		selected_precusor = precusorList[outputIndex]
		#print outputIndex
	else:
		newrangeList = []
		for i in range(len(precusorFeatureList)):
			if i not in nullIndex:
				newrangeList.append(i)
		if len(newrangeList) >= 1:
			tmpList = [[precusorFeatureList[j][5], j] for j in newrangeList]
			tmpList.sort()
			outputIndex = tmpList[0][1]
			selected_precusor = precusorList[outputIndex]
		else:
			selected_precusor = None
	return selected_precusor


def selcteOptimalPrecusor(precusorFeatureList, precusorList, clusterSeqType):
	#print precusorFeatureList
	#print clusterSeqType
	#print precusorList
	#To do, select the precusor based on mFE or precusorOrignialMFE. Right now, I choose precusorOrignialMFE
	if clusterSeqType == 'clusterSeq':
		bindingCountThreshold = 19
	elif clusterSeqType == 'stableClusterSeq':
		bindingCountThreshold = 19-1-3
	else:
		bindingCountThreshold = 19
	if len(precusorFeatureList) != len(precusorList):
		print 'Errors happen at precusorFeatureList and precusorList because of different length'
		#print len(precusorFeatureList)
		#print len(precusorList)
		#print precusorFeatureList
		#print precusorList
		sys.exit(1)
	retainedIndex = []
	nullIndex = []
	for i in range(len(precusorFeatureList)):
		#print precusorFeatureList[i]
		if not (None in precusorFeatureList[i]):
			armTypeOverlenDic, hairpinCount, bindingCount, percent, mFE, precusorOrignialMFE = precusorFeatureList[i]
			#print hairpinCount, bindingCount, percent, armTypeOverlenDic
			if hairpinCount > 0 and bindingCount >= bindingCountThreshold and percent >= 0.6:
				if 'arm5' in armTypeOverlenDic.keys() and 'arm3' in armTypeOverlenDic.keys():
					if armTypeOverlenDic['arm5'][0] <= 5 and armTypeOverlenDic['arm3'][0] == 0:
						retainedIndex.append(i)
				elif 'arm5' in armTypeOverlenDic.keys():
					if armTypeOverlenDic['arm5'][0] <= 5:
						retainedIndex.append(i)
				elif 'arm3' in armTypeOverlenDic.keys():
					if armTypeOverlenDic['arm3'][0] == 0:
						retainedIndex.append(i)
				else:
					pass
		else:
			nullIndex.append(i)
		#print retainedIndex
	#print retainedIndex
	if len(retainedIndex) >= 1:
		#To do, select the precusor based on mFE or precusorOrignialMFE. Right now precusorOrignialMFE is chosen.
		tmpList = [[precusorFeatureList[j][5], j] for j in retainedIndex]
		tmpList.sort()
		#print tmpList
		outputIndex = tmpList[0][1]
		selected_precusor = precusorList[outputIndex]
		#print outputIndex
	else:
		newrangeList = []
		for i in range(len(precusorFeatureList)):
			if i not in nullIndex:
				newrangeList.append(i)
		if len(newrangeList) >= 1:
			tmpList = [[precusorFeatureList[j][5], j] for j in newrangeList]
			tmpList.sort()
			outputIndex = tmpList[0][1]
			selected_precusor = precusorList[outputIndex]
		else:
			outputIndex = None
			selected_precusor = None
	return (outputIndex, selected_precusor)

def screen_precusor_candidates(featureFile, precusorStrFile, prunedLength, rnafoldCmdTmp):
	closestDistance = 9
	fathestDistance = 44
	# The directory of featureFile, precusorFastaFile and precusorStrFile must be absolute directory.
	# miRNA will be tested by clusterSeq, majority seq and stable region of clusterSeq.
	dirTmp = os.path.dirname(featureFile)
	# Load the structure file into a dictionary
	wholePrecusorNameContentDic = {}
	with open(precusorStrFile, 'r') as inf2:
		totalContent = inf2.readlines()
	totalContentSplit = chunks(totalContent, 3)
	for item in totalContentSplit:
		precusorName = item[0].strip()[1:]
		precusorSeq = item[1].strip()
		precusorStr = item[2].strip().split()[0]
		try:
			if len(item[2].strip().split()[1]) >= 3:
				mfeValue = float(item[2].strip().split()[1][1:-1].strip())
			else:
				mfeValue = float(item[2].strip().split()[2][:-1].strip())
		except:
			mfeValue = 0.0
			print 'Errors of mfe values happens at: %s. Please check it in %s'%(clusterName, precusorStrFile)
		wholePrecusorNameContentDic.update({precusorName:[precusorName, precusorSeq, precusorStr, mfeValue]})
	#print len(wholePrecusorNameContentDic)
	#print wholePrecusorNameContentDic['mapped_nonMirna_ERR409900:miRCluster_30_28:chr1:630216_630243+:precusor_2']
	#print wholePrecusorNameContentDic['mapped_nonMirna_ERR409900:miRCluster_39_29:chr1:630492_630520+:precusor_2']
	#print wholePrecusorNameContentDic['mapped_nonMirna_ERR409900:miRCluster_28_25:chr1:630129_630153+:precusor_2']
	
	#for clusterSeqType in ['clusterSeq', 'majoritySeq', 'stableClusterSeq']:
	for clusterSeqType in ['stableClusterSeq']:
	#for clusterSeqType in ['clusterSeq', 'stableClusterSeq']:
		clusterNameList = []
		miRNAlist = []
		neighborStateList = []
		upstreamDistanceList = []
		downstreamDistanceList = []
		clusterNameCandidatePreDic = {}
		clusterNameMiRNASeqDicTmp = {}
		clusterNameMiRNASeqDic = {}
		with open(featureFile , 'r') as inf1:
			line1 = inf1.readline()
			content1 = line1.strip().split('\t')
			miRNALabel = content1.index(clusterSeqType)
			clusterNameLabel = content1.index('clusterName')
			neighborStateLabel = content1.index('neighborState')
			upstreamDistanceLabel = content1.index('upstreamDistance')
			downstreamDistanceLabel = content1.index('downstreamDistance')
			line1 = inf1.readline()
			while line1 != '':
				content1 = line1.strip().split('\t')
				clusterName = content1[clusterNameLabel]
				miRNASeq = str(Seq(content1[miRNALabel], generic_dna).transcribe())
				neighborState = content1[neighborStateLabel]
				if content1[upstreamDistanceLabel] == 'None':
					upstreamDistance = 10000
				else:
					upstreamDistance = int(content1[upstreamDistanceLabel])
				if content1[downstreamDistanceLabel] == 'None':
					downstreamDistance = 10000
				else:
					downstreamDistance = int(content1[downstreamDistanceLabel])
				clusterNameList.append(clusterName)
				miRNAlist.append(miRNASeq)
				neighborStateList.append(neighborState)
				upstreamDistanceList.append(upstreamDistance)
				downstreamDistanceList.append(downstreamDistance)
				clusterNameCandidatePreDic.update({clusterName:[]})
				clusterNameMiRNASeqDicTmp.update({clusterName:[]})
				line1 = inf1.readline()
		#print len(clusterNameList)
		#print len(miRNAlist)
		#print len(neighborStateList)
		#print len(upstreamDistanceList)
		#print len(downstreamDistanceList)
		#print len(clusterNameCandidatePreDic)
		#print len(clusterNameMiRNASeqDicTmp)
		#print '**************'
		
		for i in range(len(clusterNameList)):
			clusterName = clusterNameList[i]
			direction = clusterName[-1]
			upstreamDistance = upstreamDistanceList[i]
			downstreamDistance = downstreamDistanceList[i]
			neighborState = neighborStateList[i]
			# precursor (70 bp upstream and 20 bp downstream)
			upPrecusor = wholePrecusorNameContentDic[clusterName+':precusor_1']
			# precursor (20 bp upstream and 70 bp downstream)
			downPrecusor = wholePrecusorNameContentDic[clusterName+':precusor_2']
			if neighborState != 'Good':
				clusterNameCandidatePreDic[clusterName].append(upPrecusor)
				clusterNameCandidatePreDic[clusterName].append(downPrecusor)
				clusterNameMiRNASeqDicTmp[clusterName].append(miRNAlist[i])
			else:
				if (upstreamDistance > fathestDistance or upstreamDistance < closestDistance) and (downstreamDistance <= fathestDistance and downstreamDistance >= closestDistance):
					# The first one is the downPrecursor of this seq
					clusterNameCandidatePreDic[clusterName].append(downPrecusor)
					# The second one is the upPrecursor of the next seq 
					clusterNameCandidatePreDic[clusterName].append(wholePrecusorNameContentDic[clusterNameList[i+1]+':precusor_1'])
					clusterNameMiRNASeqDicTmp[clusterName].append(miRNAlist[i])
					clusterNameMiRNASeqDicTmp[clusterName].append(miRNAlist[i+1])
				elif (upstreamDistance <= fathestDistance and upstreamDistance >=  closestDistance) and (downstreamDistance <= fathestDistance and downstreamDistance >= closestDistance):
					# The first one is the downPrecursor of the previous seq
					clusterNameCandidatePreDic[clusterName].append(wholePrecusorNameContentDic[clusterNameList[i-1]+':precusor_2'])
					# The first one is the upPrecursor of this seq
					clusterNameCandidatePreDic[clusterName].append(upPrecusor)
					# The third one is the downPrecursor of this seq
					clusterNameCandidatePreDic[clusterName].append(downPrecusor)
					# The forth one is the upPrecursor of the next seq
					clusterNameCandidatePreDic[clusterName].append(wholePrecusorNameContentDic[clusterNameList[i+1]+':precusor_1'])
					clusterNameMiRNASeqDicTmp[clusterName].append(miRNAlist[i])
					clusterNameMiRNASeqDicTmp[clusterName].append(miRNAlist[i-1])
					clusterNameMiRNASeqDicTmp[clusterName].append(miRNAlist[i+1])
				elif (upstreamDistance <= fathestDistance and upstreamDistance >=  closestDistance) and (downstreamDistance > fathestDistance or downstreamDistance < closestDistance):
					# The first one is the downPrecursor of the previous seq
					clusterNameCandidatePreDic[clusterName].append(wholePrecusorNameContentDic[clusterNameList[i-1]+':precusor_2'])
					# The second one is the upPrecursor of this seq
					clusterNameCandidatePreDic[clusterName].append(upPrecusor)
					clusterNameMiRNASeqDicTmp[clusterName].append(miRNAlist[i])
					clusterNameMiRNASeqDicTmp[clusterName].append(miRNAlist[i-1])
				else:
					print "Exception happens at: %s when picking candidate precursors"%(clusterName)
		#print len(clusterNameCandidatePreDic)
		#print len(clusterNameMiRNASeqDicTmp)

		for i in range(len(clusterNameList)):
			clusterName = clusterNameList[i]
			direction = clusterName[-1]
			candidatePreList = clusterNameCandidatePreDic[clusterName]
			clusterSeqListTmp = clusterNameMiRNASeqDicTmp[clusterName]
			clusterSeq = clusterSeqListTmp[0]
			numberOfCandidate = len(candidatePreList)
			if numberOfCandidate == 2:
				clusterSeqList = [clusterSeqListTmp for i in range(numberOfCandidate)]
			else:
				clusterSeqList = []
				clusterSeqList.append([clusterSeqListTmp[0], clusterSeqListTmp[1]])
				clusterSeqList.append([clusterSeqListTmp[0], clusterSeqListTmp[1]])
				clusterSeqList.append([clusterSeqListTmp[0], clusterSeqListTmp[2]])
				clusterSeqList.append([clusterSeqListTmp[0], clusterSeqListTmp[2]])
			clusterNameMiRNASeqDic.update({clusterName:clusterSeqList})
		#print clusterNameCandidatePreDic['mapped_nonMirna_ERR409900:miRCluster_11_29:chr1:629373_629401+']
		#print clusterNameMiRNASeqDic['mapped_nonMirna_ERR409900:miRCluster_11_29:chr1:629373_629401+']
		#print clusterNameCandidatePreDic['mapped_nonMirna_ERR409900:miRCluster_2324_23:chr1:630375_630397-']
		#print clusterNameMiRNASeqDic['mapped_nonMirna_ERR409900:miRCluster_2324_23:chr1:630375_630397-']
		#print clusterNameCandidatePreDic['mapped_nonMirna_ERR409900:miRCluster_134_26:chr1:634803_634828+']
		#print clusterNameMiRNASeqDic['mapped_nonMirna_ERR409900:miRCluster_134_26:chr1:634803_634828+']
		#print clusterNameCandidatePreDic['mapped_nonMirna_ERR409900:miRCluster_138_24:chr1:634887_634910+']
		#print clusterNameMiRNASeqDicTmp['mapped_nonMirna_ERR409900:miRCluster_138_24:chr1:634887_634910+']
		#print clusterNameCandidatePreDic['mapped_nonMirna_ERR409900:miRCluster_136_25:chr1:634842_634866+']
		#print clusterNameMiRNASeqDic['mapped_nonMirna_ERR409900:miRCluster_136_25:chr1:634842_634866+']
		#print clusterNameCandidatePreDic['mapped_nonMirna_ERR409900:miRCluster_138_24:chr1:634887_634910+']
		#print clusterNameMiRNASeqDic['mapped_nonMirna_ERR409900:miRCluster_138_24:chr1:634887_634910+']
		#print len(clusterNameCandidatePreDic)
		#print len(clusterNameMiRNASeqDic)
		
		# selcteOptimalPrecusor
		clustermiRNADic = {}
		clusterPrecusorDic = {}
		for i in range(len(clusterNameList)):
			clusterName = clusterNameList[i]
			precusorList = clusterNameCandidatePreDic[clusterName]
			precusorFeatureList = []
			prunedPrecusorList = []
			for index, precusor in enumerate(precusorList):
				precusorName = precusor[0]
				precusorSeq = precusor[1]
				precusorStr = precusor[2]
				precusorOrignialMFE = precusor[3]
				clusterSeqListTmp2 = clusterNameMiRNASeqDic[clusterName][index]
				if len(clusterSeqListTmp2) == 1:
					precsorInstance = ReadPrecusor(rnafoldCmdTmp, precusorName, precusorSeq, precusorStr, precusorOrignialMFE, dirTmp, prunedLength, clusterSeqListTmp2[0])
				else:
					precsorInstance = ReadPrecusor(rnafoldCmdTmp, precusorName, precusorSeq, precusorStr, precusorOrignialMFE, dirTmp, prunedLength, clusterSeqListTmp2[0], clusterSeqListTmp2[1])
				prunedPrecusorSeq, prunedPrecusorStr, prunedPrecusorMFE = precsorInstance.getPrunedPrecusor()
				if len(clusterSeqListTmp2) == 1:
					prunedPrecursorInstance = ReadPrecusor(rnafoldCmdTmp, precusorName, prunedPrecusorSeq, prunedPrecusorStr, prunedPrecusorMFE, dirTmp, prunedLength, clusterSeqListTmp2[0])
				else:
					prunedPrecursorInstance = ReadPrecusor(rnafoldCmdTmp, precusorName, prunedPrecusorSeq, prunedPrecusorStr, prunedPrecusorMFE, dirTmp, prunedLength, clusterSeqListTmp2[0], clusterSeqListTmp2[1])
				
				#if clusterName == 'mapped_nonMirna_ERR409900:miRCluster_11_29:chr1:629373_629401+':
				#	print prunedPrecusorSeq
				#	print prunedPrecusorStr
				#	print prunedPrecusorMFE
				tmpList1 = []
				if not (prunedPrecusorStr is None)  and '(' in prunedPrecusorStr and ')' in prunedPrecusorStr:
					percent = prunedPrecursorInstance.percentageOfPairedInMiRNA()
					hairpinCount, bindingCount, interiorLoopCount, armType, apicalLoopSize, armTypeOverlenDic, stemLen, flag, mFE, pairState =  prunedPrecursorInstance.featuresInPrecusor()
					precusorFeatureList.append([armTypeOverlenDic, hairpinCount, bindingCount, percent, mFE, precusorOrignialMFE])
					tmpList1.append(precusorName)
					tmpList1.append(prunedPrecusorSeq)
					tmpList1.append(prunedPrecusorStr)
					tmpList1.append(mFE)
				else:
					precusorFeatureList.append([None, None, None, None, None, None])
					tmpList1.append(precusorName)
					tmpList1.append(prunedPrecusorSeq)
					tmpList1.append(prunedPrecusorStr)
					tmpList1.append(None)
				prunedPrecusorList.append(tmpList1)
			outputIndex, candiatePrecusor = selcteOptimalPrecusor(precusorFeatureList, prunedPrecusorList, clusterSeqType)
			if not (outputIndex is None):
				clustermiRNADic.update({clusterName:clusterNameMiRNASeqDic[clusterName][outputIndex]})
			if not (candiatePrecusor is None):
				clusterPrecusorDic.update({clusterName : candiatePrecusor})
		#print len(clustermiRNADic)
		#print len(clusterPrecusorDic)
		#print clusterPrecusorDic['mapped_nonMirna_ERR409900:miRCluster_11_29:chr1:629373_629401+']
		#print clustermiRNADic['mapped_nonMirna_ERR409900:miRCluster_11_29:chr1:629373_629401+']
		
		outf = open(featureFile[:-4]+'_updated_'+clusterSeqType+'_'+str(prunedLength)+'.tsv', 'w')
		with open(featureFile , 'r') as inf:
			line = inf.readline()
			content1 = line.strip().split('\t')
			clusterNameLabel = content1.index('clusterName')
			outf.write(line.strip()+'\t')
			outf.write('pruned_precusor_seq\tpruned_precusor_str\tarmType\tdistanceToloop\tpercentage_PairedInMiRNA\thairpin_count\tbinding_count\tinteriorLoopCount\tapicalLoop_size\tstem_length\tmFE\tcount_bindings_in_miRNA\tUGU_UGUG_motif\tpair_state\n')
			line = inf.readline()
			while line != '':
				clusterName = line.strip().split('\t')[clusterNameLabel]
				if clusterName in clusterPrecusorDic.keys():
					premiRNASeqContentList = clusterPrecusorDic[clusterName]
					premiRNASeqName = premiRNASeqContentList[0]
					premiRNASeq = premiRNASeqContentList[1]
					premiRNAStructure = premiRNASeqContentList[2]
					originalMFE = premiRNASeqContentList[3]
					miRlist = clustermiRNADic[clusterName]
					if len(miRlist) == 1:
						miRNASeq = miRlist[0]
						precsorInstance = ReadPrecusor(rnafoldCmdTmp, premiRNASeqName, premiRNASeq, premiRNAStructure, originalMFE, dirTmp, prunedLength, miRNASeq)
					else:
						miRNASeq1 = miRlist[0]
						miRNASeq2 = miRlist[1]
						precsorInstance = ReadPrecusor(rnafoldCmdTmp, premiRNASeqName, premiRNASeq, premiRNAStructure, originalMFE, dirTmp, prunedLength, miRNASeq1, miRNASeq2)
					#print clusterName
					#print miRlist
					#prunedPrecusorSeq, prunedPrecusorStr = precsorInstance.getPrunedPrecusor()
					#print corePrecusorSeq
					#print corePrecusorStr
					#print '**************'
					percent = precsorInstance.percentageOfPairedInMiRNA()
					hairpinCount, bindingCount, interiorLoopCount, armType, apicalLoopSize, armTypeOverlenDic, stemLen, flag, mFE, pairState =  precsorInstance.featuresInPrecusor()
					countBindingsInMiRNA = precsorInstance.bindingsInMiRNA()
					outf.write(line.strip()+'\t')
					outf.write('%s\t%s\t%s\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%.2f\t%d\t%s\t%s\n'%(premiRNASeq, premiRNAStructure, armType, armTypeOverlenDic[armType][1], percent, hairpinCount, bindingCount, interiorLoopCount, apicalLoopSize, stemLen, mFE, countBindingsInMiRNA, flag, pairState))
				else:
					outf.write(line.strip()+'\t')
					outf.write('\t'.join(['None']*14))
					outf.write('\n')
				line = inf.readline()
		outf.close()
