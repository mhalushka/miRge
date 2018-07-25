import os
import sys
import time
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
import re

def parseBowtieLog(bowtieLog):
	s1Tmp = 0
	s2Tmp = 0
	with open(bowtieLog, 'r') as inf:
		for line in inf:
			if '# reads processed:' in line:
				s1Tmp = int(line.strip().split(' ')[-1])
			if '# reads with at least one reported alignment:' in line:
				s2Tmp = int(line.strip().split(' ')[-2])
	return (s1Tmp, s2Tmp)

def parseAlignment(bowtieSam):
	alignmentResult = {}
	with open(bowtieSam, 'r') as inf:
		for line in inf:
			if '@' not in line:
				content = line.strip().split('\t')
				#if content[2] != '*':
				alignmentResult.update({content[0]:[content[2], content[1], content[3], content[5]]})
	return alignmentResult

def parseAlignment2(bowtieSam):
	# tRF-1 has 3' poly-Us.
	alignmentResult = {}
	with open(bowtieSam, 'r') as inf:
		for line in inf:
			if '@' not in line:
				content = line.strip().split('\t')
				if content[2] != '*' and content[0][-3:] == 'TTT':
					alignmentResult.update({content[0]:[content[2], content[1], content[3], content[5]]})
	return alignmentResult

def parseAlignment3(bowtieSam):
	alignmentResult = {}
	with open(bowtieSam, 'r') as inf:
		for line in inf:
			if '@' not in line:
				content = line.strip().split('\t')
				if content[2] != '*':
					try:
						alignmentResult[content[0]].append((content[2], content[1], content[3], content[5]))
					except KeyError:
						alignmentResult.update({content[0]:[(content[2], content[1], content[3], content[5])]})
	return alignmentResult

def parseAlignment4(bowtieSam):
	alignmentResult = {}
	with open(bowtieSam, 'r') as inf:
		for line in inf:
			if '@' not in line:
				content = line.strip().split('\t')
				if content[2] != '*':
					if content[0][-3:] == 'TTT':
						try:
							alignmentResult[content[0]].append((content[2], content[1], content[3], content[5]))
						except KeyError:
							alignmentResult.update({content[0]:[(content[2], content[1], content[3], content[5])]})
	return alignmentResult

def leftDashCount(seq):
	count = 0
	for item in seq:
		if item == '-':
			count = count + 1
		else:
			break
	return count

def rightDashCount(seq):
	count = 0
	for item in reversed(seq):
		if item == '-':
			count = count + 1
		else:
			break
	return count

def fillTerminal(preMiRSeq, mirnaLibSeq, mirnaSeq, readSeq, start, indexValue):
	#print 'Index:%d'%(indexValue)
	#print start
	#print preMiRSeq
	#print mirnaLibSeq
	#print mirnaSeq
	#print readSeq
	#print '$$$$$$'
	# fill terminals of preMiRSeq and mirnaLibSeq
	try:
		indexTmp = preMiRSeq.index(mirnaSeq)
		indexState = True
		if indexTmp+len(mirnaSeq)+6 > len(preMiRSeq):
			preMiRSeq = preMiRSeq+'-'*(indexTmp+len(mirnaSeq)+6-len(preMiRSeq))
		else:
			mirnaLibSeq = mirnaLibSeq+'-'*(len(preMiRSeq)-indexTmp-len(mirnaSeq)-6)
		if indexTmp < 2:
			preMiRSeq = '-'*(2-indexTmp)+preMiRSeq
		else:
			mirnaLibSeq = '-'*(indexTmp-2)+mirnaLibSeq
		preMiRSeqOut = preMiRSeq
		mirnaLibSeqOut = mirnaLibSeq
		# fill terminals of mirnaSeq
		indexTmp2 = preMiRSeqOut.index(mirnaSeq)
		mirnaSeq = mirnaSeq+'-'*(len(preMiRSeqOut)-indexTmp2-len(mirnaSeq))
		mirnaSeq = '-'*indexTmp2+mirnaSeq
		mirnaSeqOut = mirnaSeq
		# fill terminals of readSeq
		# this is the situation of miRNAs
		if indexValue == 0:
			i = 0
			headAdd = ''
			for item in mirnaLibSeqOut:
				if item == '-':
					headAdd = headAdd + '-'
				else:
					i = i + 1
					if i == start:
						break
					else:
						headAdd = headAdd + '-'
			readSeq = headAdd + readSeq
			readSeq = readSeq + '-'*(len(mirnaLibSeqOut)-len(readSeq))
		# this is the situration of isomiRs
		else:
			trimmedHeadSeq = readSeq[:1]
			trimmedTailSeq = readSeq[-2:]
			trimmedSeq = readSeq[1:-2]
			if len(readSeq)-1 > len(mirnaLibSeqOut)-leftDashCount(mirnaLibSeqOut)-start+1:
				countTmp = len(readSeq)-1-len(mirnaLibSeqOut)+leftDashCount(mirnaLibSeqOut)+start-1
				mirnaLibSeqOut = mirnaLibSeqOut + '-'*countTmp
				preMiRSeqOut = preMiRSeqOut + '-'*countTmp
				mirnaSeqOut = mirnaSeqOut + '-'*countTmp
			else:
				countTmp = len(mirnaLibSeqOut)-leftDashCount(mirnaLibSeqOut)-start+1-len(readSeq)+1
				readSeq = readSeq+'-'*countTmp
			if start-1+leftDashCount(mirnaLibSeqOut)<1:
				preMiRSeqOut = '-'*(1-start+1-leftDashCount(mirnaLibSeqOut))+preMiRSeqOut
				mirnaSeqOut = '-'*(1-start+1-leftDashCount(mirnaLibSeqOut))+mirnaSeqOut
				mirnaLibSeqOut = '-'*(1-start+1-leftDashCount(mirnaLibSeqOut))+mirnaLibSeqOut
			else:
				readSeq = '-'*(start-1+leftDashCount(mirnaLibSeqOut)-1)+readSeq
		readSeqOut = readSeq
		#print preMiRSeqOut
		#print mirnaLibSeqOut
		#print mirnaSeqOut
		#print readSeqOut
		#print '######'
		#if preMiRSeqOut[0] == '-':
		#	print preMiRSeqOut
		#	print mirnaLibSeqOut
		#	print mirnaSeqOut
		#	print readSeqOut
		#	print '@@@@@@'
		#if preMiRSeqOut[-1] == '-':
		#	print preMiRSeqOut
		#	print mirnaLibSeqOut
		#	print mirnaSeqOut
		#	print readSeqOut
		#	print '&&&&&&'
	except ValueError:
		indexState = False
		preMiRSeqOut = ''
		mirnaLibSeqOut = ''
		mirnaSeqOut = ''
		readSeqOut = ''
	return (preMiRSeqOut, mirnaLibSeqOut, mirnaSeqOut, readSeqOut, indexState)

def _add_cigar_char(counter, cache):
	if counter == 1:
		return cache
	else:
		return str(counter) + cache

def make_cigar(seq, mature):
	"""
	Function that will create CIGAR string from aligment
	between read and reference sequence.
	"""
	cigar = ""
	for pos in range(0,len(seq)):
		if seq[pos] == mature[pos]:
			cigar += "M"
		elif seq[pos] != mature[pos] and seq[pos] != "-" and mature[pos] != "-":
			cigar += seq[pos]
		elif seq[pos] == "-":
			cigar += "D"
		elif mature[pos] == "-":
			cigar += "I"

	cache = ""
	counter = 1
	short = ""
	for c in cigar:
		if c != cache and cache != "" and cache == "M":
			short += _add_cigar_char(counter, cache)
			counter = 1
			cache = c
		if c != "M":
			short += c
		if c == cache and c == "M":
			counter += 1
		cache = c

	if cache == "M":
		short += _add_cigar_char(counter, cache)
	return short

def make_id(seq, NT2CODE):
	start = 0
	idName = ""
	for i in range(0, len(seq) + 1, 3):
		if i == 0:
			continue
		trint = seq[start:i]
		try:
			idName += NT2CODE[trint]
		except KeyError:
			idName = '.'
			break
		start = i
	if len(seq) > i and idName != '.':
		dummy = "A" * (3 - (len(seq) - i))
		trint = seq[i:len(seq)]
		try:
			idName += NT2CODE["%s%s" % (trint, dummy)]
			idName += str(len(dummy))
		except KeyError:
			idName = '.'
	return idName

def analyzeAlignment(preMiRSeqOut, mirnaLibSeqOut, mirnaSeqOut, readSeqOut):
	if len(set([len(preMiRSeqOut), len(mirnaLibSeqOut), len(mirnaSeqOut), len(readSeqOut)])) != 1:
		print "Errors happen when aligning preMiRSeq, mirnaLibSeq, mirnaSeq, readSeq.\n"
		print preMiRSeqOut
		print mirnaLibSeqOut
		print mirnaSeqOut
		print readSeqOut
		sys.exit(1)
	l = len(preMiRSeqOut)
	variantList = []
	if mirnaSeqOut == readSeqOut:
		type = 'ref_miRNA'
		variantList.append('NA')
	else:
		type = 'isomiR'
		# 0-based location
		start1 = leftDashCount(mirnaSeqOut)
		end1 = len(mirnaSeqOut)-rightDashCount(mirnaSeqOut)
		# juge whether the read belongs to iso_snp
		iso_snp_status = False
		for i in range(start1, end1):
			if readSeqOut[i] != '-' and readSeqOut[i] != mirnaSeqOut[i]:
				iso_snp_status = True
				if i >=1 and i <= 6:
					iso_snp_pos = '_seed'
				elif i == 7:
					iso_snp_pos = '_central_offset'
				elif i >= 8 and i <= 11:
					iso_snp_pos = '_central'
				elif i >=12 and i <= 16:
					iso_snp_pos = 'central_supp'
				else:
					iso_snp_pos = ''
				break
		# juge whether the read belongs to iso_5p
		iso_5p_status = False
		iso_add_status = False
		if leftDashCount(readSeqOut) > leftDashCount(mirnaSeqOut):
			iso_5p_status = True
			iso_5p_pos = str(leftDashCount(mirnaSeqOut)-leftDashCount(readSeqOut))
		elif leftDashCount(readSeqOut) < leftDashCount(mirnaSeqOut):
			#iso_5p_status = True
			iso_5p_seq_start = leftDashCount(readSeqOut)
			iso_5p_seq = readSeqOut[iso_5p_seq_start:start1]
			refer_5p_seq = preMiRSeqOut[iso_5p_seq_start:start1]
			if iso_5p_seq != refer_5p_seq:
				#iso_add_status = True
				iso_5p_status = True
				iso_5p_pos = '+'+str(leftDashCount(mirnaSeqOut)-leftDashCount(readSeqOut))
				iso_snp_status = True
				iso_snp_pos = ''
			else:
				iso_5p_status = True
				iso_5p_pos = '+'+str(leftDashCount(mirnaSeqOut)-leftDashCount(readSeqOut))
		else:
			pass
		# judge whether the read belongs to iso_3p_status
		iso_3p_status = False
		if rightDashCount(readSeqOut) > rightDashCount(mirnaSeqOut):
			iso_3p_status = True
			iso_3p_pos = str(rightDashCount(mirnaSeqOut)-rightDashCount(readSeqOut))
		elif rightDashCount(readSeqOut) < rightDashCount(mirnaSeqOut):
			#iso_3p_status = True
			iso_3p_seq_end = len(readSeqOut)-rightDashCount(readSeqOut)-1
			iso_3p_seq = readSeqOut[end1:iso_3p_seq_end+1]
			refer_3p_seq = preMiRSeqOut[end1:iso_3p_seq_end+1]
			if iso_3p_seq != refer_3p_seq:
				iso_add_status = True
				iso_add_pos = '+'+str(rightDashCount(mirnaSeqOut)-rightDashCount(readSeqOut))
			else:
				iso_3p_status = True
				iso_3p_pos = '+'+str(rightDashCount(mirnaSeqOut)-rightDashCount(readSeqOut))
		else:
			pass
		if iso_snp_status:
			variantList.append('iso_snp'+iso_snp_pos)
		if iso_add_status:
			variantList.append('iso_add'+':'+iso_add_pos)
		if iso_5p_status:
			variantList.append('iso_5p'+':'+iso_5p_pos)
		if iso_3p_status:
			variantList.append('iso_3p'+':'+iso_3p_pos)
	variant = ','.join(variantList)
	if leftDashCount(preMiRSeqOut) > 0:
		if leftDashCount(readSeqOut) >= leftDashCount(preMiRSeqOut):
			pre_start = leftDashCount(readSeqOut)-leftDashCount(preMiRSeqOut)+1
		else:
			pre_start = 1-(leftDashCount(preMiRSeqOut)-leftDashCount(readSeqOut))
	else:
		pre_start = leftDashCount(readSeqOut)+1
	pre_end = pre_start + len(readSeqOut)-leftDashCount(readSeqOut)-rightDashCount(readSeqOut)-1
	
	pre_start = pre_start -1
	pre_end = pre_end -1
	start2 = leftDashCount(readSeqOut)
	end2 = len(readSeqOut)-rightDashCount(readSeqOut)-1
	readSeqOut_extract = readSeqOut[start2:end2+1]
	preMiRSeqOut_extract = preMiRSeqOut[start2:end2+1]
	cigarValue = make_cigar(readSeqOut_extract, preMiRSeqOut_extract)
	# 1-based location
	pre_start = pre_start + 1
	pre_end = pre_end + 1
	return (type, variant, pre_start, pre_end, cigarValue)

def updateAnnotDic(seqDic, alignmentResult, index):
	for seqKey in alignmentResult.keys():
		if alignmentResult[seqKey][0] != '*':
			seqDic[seqKey]['annot'][0] = 1
			seqDic[seqKey]['annot'][index] = alignmentResult[seqKey][0]

def updateAnnotDic2(seqDic, alignmentResult, index, subseqSeqDic):
	for seqKey in alignmentResult.keys():
		if alignmentResult[seqKey][0] != '*':
			for subseq in subseqSeqDic[seqKey]:
				seqDic[subseq]['annot'][0] = 1
				seqDic[subseq]['annot'][index] = alignmentResult[seqKey][0]

def inferPremiRName(miRNameCanonical, miRNamePreNameDic, miRNA_database):
	try:
		preName = miRNamePreNameDic[miRNameCanonical]
	except KeyError:
		if miRNA_database == 'miRBase':
			try:
				preName = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])]
			except KeyError:
				try:
					preName = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])+'-5p']
				except KeyError:
					try:
						preName = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])+'-3p']
					except KeyError:
						# In this case, the pre miRNA name can't be retrieved from the *.gff3 file meaning the miRNA's pre miRNA doesn't have coordinates.
						# This may happen in miRBase, so the pre miRNA name will be speculated.
						preName = miRNameCanonical.replace('-5p', '')
						preName = preName.replace('-3p', '')
						preName = preName.replace('miR', 'mir')
		else:
			# infer pre miRNA name.
			preName = miRNameCanonical.replace('_5p*', '')
			preName = preName.replace('_3p*', '')
			preName = preName.replace('_5p', '')
			preName = preName.replace('_3p', '')
			preName = preName + '_pre'
	return preName

def updateIsomiRDic(isomiRContentDic, alignmentResult, miRNamePreNameDic, indexValue, NT2CODE, miRNA_database):
	for seqKey in alignmentResult.keys():
		if alignmentResult[seqKey][0] != '*':
			if '.' in alignmentResult[seqKey][0]:
				miRNameCanonical = alignmentResult[seqKey][0].split('.')[0]
			else:
				miRNameCanonical = alignmentResult[seqKey][0]
			try:
				isomiRContentDic[seqKey]['miRName'] = alignmentResult[seqKey][0]
			except KeyError:
				isomiRContentDic.update({seqKey:{}})
				isomiRContentDic[seqKey]['miRName'] = alignmentResult[seqKey][0]

			preMiRName = inferPremiRName(miRNameCanonical, miRNamePreNameDic, miRNA_database)
			isomiRContentDic[seqKey]['preMiRName'] = preMiRName
			
			isomiRContentDic[seqKey]['start'] = alignmentResult[seqKey][2]
			if indexValue == 0:
				isomiRContentDic[seqKey]['cigar'] = alignmentResult[seqKey][3]
			else:
				isomiRContentDic[seqKey]['cigar'] = str(len(seqKey))+'M'

			isomiRContentDic[seqKey]['annot'] = 0
			isomiRContentDic[seqKey]['filter'] = 'Pass'
			isomiRContentDic[seqKey]['uid'] = make_id(seqKey, NT2CODE)

def updateIsomiRDic2(isomiRContentDic, hairpinNameSeqDic, mirnaLibNameSeqDic, indexValue, miRNamePreNameDic, miRNA_database):
	for seqKey in isomiRContentDic.keys():
		if isomiRContentDic[seqKey]['annot'] == 0:
			mirnaLibSeq = mirnaLibNameSeqDic[isomiRContentDic[seqKey]['miRName']]
			start = int(isomiRContentDic[seqKey]['start'])
			mirnaSeq = mirnaLibSeq[2:-6]
			readSeq = seqKey
			preMiRSeq = hairpinNameSeqDic[isomiRContentDic[seqKey]['preMiRName']]
			# if '.SNP*' exist in the name of mirnaSeq, the preMIRSeq should be modified accordingly.
			if '.SNP' in isomiRContentDic[seqKey]['miRName'] and '.SNPC' not in isomiRContentDic[seqKey]['miRName']:
				mirnaLibSeqCanonical = mirnaLibNameSeqDic[isomiRContentDic[seqKey]['miRName'].split('.')[0]+'.SNPC']
				mirnaSeqCanonical = mirnaLibSeqCanonical[2:-6]
				miRNameCanonical = isomiRContentDic[seqKey]['miRName'].split('.')[0]
				PreName = inferPremiRName(miRNameCanonical, miRNamePreNameDic, miRNA_database)
				preMiRSeqCanonical = hairpinNameSeqDic[PreName]
				try:
					#index_tmp = preMiRSeqCanonical.index(mirnaLibSeqCanonical)
					index_tmp = preMiRSeqCanonical.index(mirnaSeqCanonical)
				except ValueError:
					print 'Error happens.'
					print PreName
					print preMiRSeqCanonical
					print miRNameCanonical
					print mirnaLibSeqCanonical
					sys.exit(1)
				preMiRSeq = preMiRSeqCanonical[:index_tmp]+mirnaSeq+preMiRSeqCanonical[index_tmp+len(mirnaSeq):]
			# fill the head and tail part across preMiRSeq, mirnaLibSeq, mirnaSeq and readSeq via '-', so that the attribute of type and variant can be determined.
			preMiRSeqOut, mirnaLibSeqOut, mirnaSeqOut, readSeqOut, indexState = fillTerminal(preMiRSeq, mirnaLibSeq, mirnaSeq, readSeq, start, indexValue)
			if indexState:
				type, variant, pre_start, pre_end, cigarValue = analyzeAlignment(preMiRSeqOut, mirnaLibSeqOut, mirnaSeqOut, readSeqOut)
				isomiRContentDic[seqKey]['type'] = type
				isomiRContentDic[seqKey]['pre_start'] = str(pre_start)
				isomiRContentDic[seqKey]['pre_end'] = str(pre_end)
				isomiRContentDic[seqKey]['variant'] = str(variant)
				isomiRContentDic[seqKey]['strand'] = '+'
				isomiRContentDic[seqKey]['annot'] = 1
				isomiRContentDic[seqKey]['cigar'] = cigarValue
			else:
				del isomiRContentDic[seqKey]

def trfTypes(seq, tRNAName, start, trnaStruDic, pretrnaNameSeqDic):
	if 'pre_' not in tRNAName:
		tRNASeq = trnaStruDic[tRNAName]['seq']
		tRNAStru = trnaStruDic[tRNAName]['stru']
		tRNASeqLen = len(tRNASeq)
		# anticodonStart and anticodonEnd is 1-based position, so change it into 0-based
		anticodonStart = trnaStruDic[tRNAName]['anticodonStart']-1
		anticodonEnd = trnaStruDic[tRNAName]['anticodonEnd']-1
		if start == 0:
			if start+len(seq) == tRNASeqLen:
				trfType = 'tRF-whole'
			elif start+len(seq)-1 >= anticodonStart-2 and start+len(seq)-1 <= anticodonStart+1:
				trfType = "5'-half"
			else:
				trfType = "5'-tRF"
		else:
			if start+len(seq)-1 >= tRNASeqLen-1-2 and start+len(seq)-1 <= tRNASeqLen-1:
				if start >= anticodonStart-1 and start <= anticodonStart+2:
					trfType = "3'-half"
				else:
					trfType = "3'-tRF"
			else:
				trfType = 'i-tRF'
	else:
		tRNASeq = pretrnaNameSeqDic[tRNAName]
		trfType = 'tRF-1'
	return trfType

def updatetrfContentDic(trfContentDic, alignmentResult, trnaStruDic, pretrnaNameSeqDic, NT2CODE):
	# This function is invalid now.
	for seqKey in alignmentResult.keys():
		#if alignmentResult[seqKey][0] != '*':
		trfContentDic[seqKey]['uid'] = make_id(seqKey, NT2CODE)
		for item in alignmentResult[seqKey]:
			tRNAName = item[0]
			if tRNAName != '*':
				# item[2] is 1-based, so change it into 0-based
				startTmp = int(item[2])-1
				# add tRNAName as a key.
				trfContentDic[seqKey][tRNAName] = {}
				# determine the tRF type based on the postion of the tRF for specific tRNA or pre-tRNA-trailer
				trfContentDic[seqKey][tRNAName]['tRFType'] = trfTypes(seqKey, tRNAName, startTmp, trnaStruDic, pretrnaNameSeqDic)
				trfContentDic[seqKey][tRNAName]['start'] = startTmp
				if 'pre_' not in tRNAName:
					trfContentDic[seqKey][tRNAName]['end'] = startTmp+len(seqKey)-1
				else:
					trfContentDic[seqKey][tRNAName]['end'] = startTmp+len(seqKey)-4
				trfContentDic[seqKey][tRNAName]['cigar'] = 'undifined'

def updatetrfContentDic2(trfContentDic, alignmentResultNew, trnaStruDic, pretrnaNameSeqDic, NT2CODE, seqDic, sampleList):
	# Deal with the alignment of mature tRNA
	for seqKey in alignmentResultNew.keys():
		readCountList = []
		for i in range(len(sampleList)):
			readCountList.append(int(seqDic[seqKey]['quant'][i]))
		trfContentDic.update({seqKey:{'count':readCountList}})
		trfContentDic[seqKey]['uid'] = make_id(seqKey, NT2CODE)
		for item in alignmentResultNew[seqKey]:
			tRNAName = item[0]
			if tRNAName != '*':
				# item[2] is 1-based, so change it into 0-based
				startTmp = int(item[2])-1
				# add tRNAName as a key.
				trfContentDic[seqKey][tRNAName] = {}
				# determine the tRF type based on the postion of the tRF for specific tRNA or pre-tRNA-trailer
				trfContentDic[seqKey][tRNAName]['tRFType'] = trfTypes(seqKey, tRNAName, startTmp, trnaStruDic, pretrnaNameSeqDic)
				trfContentDic[seqKey][tRNAName]['start'] = startTmp
				trfContentDic[seqKey][tRNAName]['end'] = startTmp+len(seqKey)-1
				trfContentDic[seqKey][tRNAName]['cigar'] = 'undifined'

def updatetrfContentDic3(trfContentDic, alignmentResultNew, trnaStruDic, pretrnaNameSeqDic, NT2CODE, seqDic, sampleList, subseqSeqDic):
	# Deal with the alignment of precusor tRNA
	for seqKey in alignmentResultNew.keys():
		for readSeq in subseqSeqDic[seqKey]:
			readCountList = []
			for i in range(len(sampleList)):
				readCountList.append(int(seqDic[readSeq]['quant'][i]))
			trfContentDic.update({readSeq:{'count':readCountList}})
			trfContentDic[readSeq]['uid'] = make_id(readSeq, NT2CODE)
			for item in alignmentResultNew[seqKey]:
				tRNAName = item[0]
				if tRNAName != '*':
					startTmp = int(item[2])-1
					trfContentDic[readSeq][tRNAName] = {}
					trfContentDic[readSeq][tRNAName]['tRFType'] = trfTypes(readSeq, tRNAName, startTmp, trnaStruDic, pretrnaNameSeqDic)
					trfContentDic[readSeq][tRNAName]['start'] = startTmp
					numberOfTailT2 = 0
					for t in readSeq[::-1]:
						if t == 'T':
							numberOfTailT2 += 1
						else:
							break
					trfContentDic[readSeq][tRNAName]['end'] = startTmp+len(readSeq)-1-numberOfTailT2
					trfContentDic[readSeq][tRNAName]['cigar'] = 'undifined'

def writeSeqToAnnot(lengthFilter, seqDic, outputdir):
	with open(os.path.join(outputdir, 'SeqToAnnot.fasta'), 'w') as outf:
		for seq in seqDic.keys():
			if seqDic[seq]['annot'][0] == 0:
				if lengthFilter < 0:
					if seqDic[seq]['length'] < 0-lengthFilter:
						outf.write('>'+seq+'\n'+seq+'\n')
				elif lengthFilter > 0:
					if  seqDic[seq]['length'] > lengthFilter:
						outf.write('>'+seq+'\n'+seq+'\n')
				else:
					outf.write('>'+seq+'\n'+seq+'\n')

def writetRFToAnnot(tRFsLabelDic, outputdir):
	with open(os.path.join(outputdir, 'tRFToAnnot.fasta'), 'w') as outf:
		for seqKey in tRFsLabelDic.keys():
			if tRFsLabelDic[seqKey] == 0:
				outf.write('>'+seqKey+'\n'+seqKey+'\n')

def updatetRFsLabelDic(tRFsLabelDic, alignmentResult):
	for seq in alignmentResult.keys():
		tRFsLabelDic[seq] = 1

def runAnnotationPipeline(bowtieBinary, seqDic, numCPU, phred64, annotNameList, outputdir, logDic, file_mirna, file_hairpin, file_mature_tRNA, file_pre_tRNA, file_snoRNA, file_rRNA, file_ncrna_others, file_mrna, spikeIn, file_spikeIn, gff_output, miRNamePreNameDic, isomiRContentDic, miRNA_database, trf_output, trnaStruDic, trfContentDic, sampleList):
	bwtCmd = os.path.join(bowtieBinary, 'bowtie')+' --threads '+numCPU+' '
	bwtCmd2 = os.path.join(bowtieBinary, 'bowtie-inspect ')
	#bwtCmd = bowtieBinary+"bowtie --threads "+numCPU+' '
	#bwtCmd = "bowtie --threads "+numCPU+' '
	if(phred64):
		bwtCmd = bwtCmd + ' --phred64-quals '
	if spikeIn:
		lengthFilters = [-26, 25, 0, 0, 0, 0, 0, 0, 0, 0]
		rnaLibrary = ['miRNA', 'hairpin', 'mature tRNA', 'precusor tRNA', 'snoRNA', 'rRNA', 'ncrna others', 'mRNA', 'isomiR', 'spikeIn']
		bwtCmdLines = [
		bwtCmd+file_mirna+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_hairpin+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mature_tRNA+' -v 1 -f -a --best --strata --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		'Multiple processings will be perfromed',
		bwtCmd+file_snoRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_rRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_ncrna_others+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mrna+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mirna+' -5 1 -3 2 -v 2 -f --norc --best -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_spikeIn+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log')]
	else:
		lengthFilters = [-26, 25, 0, 0, 0, 0, 0, 0, 0]
		rnaLibrary = ['miRNA', 'hairpin', 'mature tRNA', 'precusor tRNA', 'snoRNA', 'rRNA', 'ncrna others', 'mRNA', 'isomiR']
		bwtCmdLines = [
		bwtCmd+file_mirna+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_hairpin+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mature_tRNA+' -v 1 -f -a --best --strata --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		'Multiple processings will be perfromed',
		bwtCmd+file_snoRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_rRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_ncrna_others+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mrna+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mirna+' -5 1 -3 2 -v 2 -f --norc --best -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log')]

	# -- ALIGNMENT 1 -- length < 26, up to 0 mismatch to miRNA
	# -- ALIGNMENT 2 -- length > 25, up to 1 mismatch to hairpin
	# -- ALIGNMENT 3 -- any length, up to 1 mismatch to mature tRNA
	# -- ALIGNMENT 4 -- any length, up to 0 mismatch with special 3 prime considerations to primary tRNA 3' trailer 
	# -- ALIGNMENT 5, 6, 7 -- any length, up to 1 mismatch to snoRNA, rRNA and ncrna others
	# -- ALIGNMENT 8 -- any length, up to 0 mismatch to mRNA
	# -- ALIGNMENT 9 -- any length, up to 2 mismatches with special 5 and 3 prime considerations to miRNA (isoMiR)
	# -- ALIGNMENT (10) -- any length, up to 0 mismatches to spike-in 
	if gff_output:
		os.system(bwtCmd2+file_hairpin+' > '+os.path.join(outputdir, 'hairpinTmp.fa'))
		os.system(bwtCmd2+file_mirna+' > '+os.path.join(outputdir, 'mirnaTmp.fa'))
		hairpinNameSeqDic = {}
		mirnaLibNameSeqDic = {}
		for record in SeqIO.parse(os.path.join(outputdir, 'hairpinTmp.fa'), 'fasta'):
			hairpinNameSeqDic.update({record.id:str(record.seq)})
		for record in SeqIO.parse(os.path.join(outputdir, 'mirnaTmp.fa'), 'fasta'):
			mirnaLibNameSeqDic.update({record.id:str(record.seq)})
		os.system('rm %s'%(os.path.join(outputdir, 'hairpinTmp.fa')))
		os.system('rm %s'%(os.path.join(outputdir, 'mirnaTmp.fa')))
	NT2CODE = {'AAA': '@', 'AAC': 'f', 'AAG': 'c', 'AAT': 'o', 'ACA': 'l', 'ACC': 'a', 'ACG': 'd', 'ACT': 's',
                   'AGA': 'm', 'AGC': 'k', 'AGG': 'h', 'AGT': 'w', 'ATA': 'g', 'ATC': 'e', 'ATG': 'b', 'ATT': 'p',
                   'CAA': 'v', 'CAC': 't', 'CAG': 'D', 'CAT': 'n', 'CCA': 'x', 'CCC': '#', 'CCG': 'y', 'CCT': 'i',
                   'CGA': 'C', 'CGC': 'E', 'CGG': 'G', 'CGT': 'S', 'CTA': 'r', 'CTC': 'j', 'CTG': 'q', 'CTT': 'H',
                   'GAA': 'T', 'GAC': '8', 'GAG': '4', 'GAT': 'F', 'GCA': 'V', 'GCC': 'X', 'GCG': 'Z', 'GCT': '6',
                   'GGA': 'K', 'GGC': 'M', 'GGG': '$', 'GGT': 'A', 'GTA': 'W', 'GTC': 'Y', 'GTG': '3', 'GTT': '5',
                   'TAA': 'L', 'TAC': 'N', 'TAG': 'J', 'TAT': 'z', 'TCA': 'U', 'TCC': '9', 'TCG': 'P', 'TCT': '0',
                   'TGA': '7', 'TGC': 'I', 'TGG': 'u', 'TGT': 'B', 'TTA': 'Q', 'TTC': 'O', 'TTG': 'R', 'TTT': '%'}

	if trf_output:
		os.system(bwtCmd2+file_pre_tRNA+' > '+os.path.join(outputdir, 'pretrnaTmp.fa'))
		pretrnaNameSeqDic = {}
		for record in SeqIO.parse(os.path.join(outputdir, 'pretrnaTmp.fa'), 'fasta'):
			pretrnaNameSeqDic.update({record.id:str(record.seq)})
		os.system('rm %s'%(os.path.join(outputdir, 'pretrnaTmp.fa')))

	for i in range(len(lengthFilters)):
		writeSeqToAnnot(lengthFilters[i], seqDic, outputdir)
		# run aligment
		#print "Starting Annotation-%s"%(annotNameList[i])
		dicRmp = {}
		time_tmp1 = time.time()
		if i != 3:
			alignmentStatus = os.system(bwtCmdLines[i])
			time_tmp2 = time.time()
			dicRmp.update({'cpuTime': time_tmp2-time_tmp1})
			# parse alignment completed without error
			if alignmentStatus == 0:
				s1, s2 = parseBowtieLog(os.path.join(outputdir, 'SeqToAnnot.log'))
				dicRmp.update({'readsProcessed': s1})
				dicRmp.update({'readsAligned':s2})
				alignmentResult = parseAlignment(os.path.join(outputdir, 'SeqToAnnot.sam'))
				updateAnnotDic(seqDic, alignmentResult, i+1)
				if gff_output:
					if i == 0 or i == 8:
						updateIsomiRDic(isomiRContentDic, alignmentResult, miRNamePreNameDic, i, NT2CODE, miRNA_database)
						updateIsomiRDic2(isomiRContentDic, hairpinNameSeqDic, mirnaLibNameSeqDic, i, miRNamePreNameDic, miRNA_database)
				if i == 2 and trf_output:
					# update mature tRNA reads into trfContentDic
					alignmentResultNew = parseAlignment3(os.path.join(outputdir, 'SeqToAnnot.sam'))
					updatetrfContentDic2(trfContentDic, alignmentResultNew, trnaStruDic, pretrnaNameSeqDic, NT2CODE, seqDic, sampleList)
			else:
				print "Alignment to library %s exited with none-zero status.\n"%(rnaLibrary[i])
				sys.exit(1)
		else:
			# Step1: Keep the reads that ends with ploy-Ts (>= 3 Ts at 5'-end).
			outfTmp = open(os.path.join(outputdir, 'SeqToAnnot2.fasta'), 'w')
			subseqSeqDic = {}
			with open(os.path.join(outputdir, 'SeqToAnnot.fasta'), 'r') as infTmp:
				for line in infTmp:
					if line[0] != '>':
						readseq = line.strip()
						if not (re.search('T{3,}$', readseq) is None):
							numberOfTailT = 0
							for t in readseq[::-1]:
								if t == 'T':
									numberOfTailT += 1
								else:
									break
							readSubseq = readseq[:-numberOfTailT]
							if len(readSubseq) >= 11:
								outfTmp.write('>'+readSubseq+'\n'+readSubseq+'\n')
								try:
									subseqSeqDic[readSubseq].append(readseq)
								except KeyError:
									subseqSeqDic.update({readSubseq:[readseq]})
			outfTmp.close()
			# Step2: Align the reads that ends with ploy-Ts (>= 3 Ts at 5'-end) to the precusor tRNA libary.
			bwtCmdTmp = os.path.join(bowtieBinary, 'bowtie')+' --threads '+numCPU+' '+file_pre_tRNA+' -v 0 -f  -a --best --strata --norc -S '+os.path.join(outputdir, 'SeqToAnnot2.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log')
			alignmentStatus = os.system(bwtCmdTmp)
			time_tmp2 = time.time()
			dicRmp.update({'cpuTime': time_tmp2-time_tmp1})
			if alignmentStatus == 0:
				s1, s2 = parseBowtieLog(os.path.join(outputdir, 'SeqToAnnot.log'))
				dicRmp.update({'readsProcessed': s1})
				dicRmp.update({'readsAligned':s2})
				alignmentResult = parseAlignment(os.path.join(outputdir, 'SeqToAnnot.sam'))
				updateAnnotDic2(seqDic, alignmentResult, i+1, subseqSeqDic)
				if trf_output:
					# update precusor tRNA reads into trfContentDic
					alignmentResultNew = parseAlignment3(os.path.join(outputdir, 'SeqToAnnot.sam'))
					updatetrfContentDic3(trfContentDic, alignmentResultNew, trnaStruDic, pretrnaNameSeqDic, NT2CODE, seqDic, sampleList, subseqSeqDic)
			else:
				print "Alignment to library %s exited with none-zero status.\n"%(rnaLibrary[i])
				sys.exit(1)
		logDic['annotStats'].append(dicRmp)
	# Remove the intermediate files
	os.system('rm %s %s %s %s'%(os.path.join(outputdir, 'SeqToAnnot.fasta'), os.path.join(outputdir, 'SeqToAnnot.sam'), os.path.join(outputdir, 'SeqToAnnot.log'), os.path.join(outputdir, 'SeqToAnnot2.fasta')))