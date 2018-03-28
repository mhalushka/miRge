import os
import sys
import time
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2

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
	indexTmp = preMiRSeq.index(mirnaSeq)
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
	return (preMiRSeqOut, mirnaLibSeqOut, mirnaSeqOut, readSeqOut)

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
					iso_snp_pos = '_central_supp'
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
	
	pre_start = pre_start -3
	pre_end = pre_end -3
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
	#for seqKey in seqDic.keys():
		#if seqKey in alignmentResult.keys():
		if alignmentResult[seqKey][0] != '*':
			seqDic[seqKey]['annot'][0] = 1
			seqDic[seqKey]['annot'][index] = alignmentResult[seqKey][0]
		#else:
		#	seqDic[seqKey]['annot'][index] = alignmentResult[seqKey][0]

def updateIsomiRDic(isomiRContentDic, alignmentResult, miRNamePreNameDic, indexValue, NT2CODE):
	for seqKey in alignmentResult.keys():
		if alignmentResult[seqKey][0] != '*':
			if '.' in alignmentResult[seqKey][0]:
				miRNameCanonical = alignmentResult[seqKey][0].split('.')[0]
			else:
				miRNameCanonical = alignmentResult[seqKey][0]
			try:
				isomiRContentDic[seqKey]['miRName'] = alignmentResult[seqKey][0]
				#isomiRContentDic[seqKey]['preMiRName'] = miRNamePreNameDic[alignmentResult[seqKey][0]]
				try:
					isomiRContentDic[seqKey]['preMiRName'] = miRNamePreNameDic[miRNameCanonical]
				except KeyError:
					isomiRContentDic[seqKey]['preMiRName'] = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])]
				isomiRContentDic[seqKey]['start'] = alignmentResult[seqKey][2]
				if indexValue == 0:
					isomiRContentDic[seqKey]['cigar'] = alignmentResult[seqKey][3]
				else:
					isomiRContentDic[seqKey]['cigar'] = str(len(seqKey))+'M'
				isomiRContentDic[seqKey]['annot'] = 0
				isomiRContentDic[seqKey]['filter'] = 'Pass'
				isomiRContentDic[seqKey]['uid'] = make_id(seqKey, NT2CODE)
			except KeyError:
				isomiRContentDic.update({seqKey:{}})
				isomiRContentDic[seqKey]['miRName'] = alignmentResult[seqKey][0]
				#isomiRContentDic[seqKey]['preMiRName'] = miRNamePreNameDic[alignmentResult[seqKey][0]]
				try:
					isomiRContentDic[seqKey]['preMiRName'] = miRNamePreNameDic[miRNameCanonical]
				except KeyError:
					try:
						isomiRContentDic[seqKey]['preMiRName'] = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])]
					except KeyError:
						try:
							isomiRContentDic[seqKey]['preMiRName'] = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])+'-5p']
						except KeyError:
							isomiRContentDic[seqKey]['preMiRName'] = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])+'-3p']
				isomiRContentDic[seqKey]['start'] = alignmentResult[seqKey][2]
				if indexValue == 0:
					isomiRContentDic[seqKey]['cigar'] = alignmentResult[seqKey][3]
				else:
					isomiRContentDic[seqKey]['cigar'] = str(len(seqKey))+'M'
				isomiRContentDic[seqKey]['annot'] = 0
				isomiRContentDic[seqKey]['filter'] = 'Pass'
				isomiRContentDic[seqKey]['uid'] = make_id(seqKey, NT2CODE)

def updateIsomiRDic2(isomiRContentDic, hairpinNameSeqDic, mirnaLibNameSeqDic, indexValue, miRNamePreNameDic):
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
				try:
					PreName = miRNamePreNameDic[miRNameCanonical]
				except KeyError:
					try:
						PreName = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])]
					except KeyError:
						try:
							PreName = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])+'-5p']
						except KeyError:
							PreName = miRNamePreNameDic['-'.join(miRNameCanonical.split('-')[:-1])+'-3p']
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
			preMiRSeqOut, mirnaLibSeqOut, mirnaSeqOut, readSeqOut = fillTerminal(preMiRSeq, mirnaLibSeq, mirnaSeq, readSeq, start, indexValue)
			type, variant, pre_start, pre_end, cigarValue = analyzeAlignment(preMiRSeqOut, mirnaLibSeqOut, mirnaSeqOut, readSeqOut)
			isomiRContentDic[seqKey]['type'] = type
			isomiRContentDic[seqKey]['pre_start'] = str(pre_start)
			isomiRContentDic[seqKey]['pre_end'] = str(pre_end)
			isomiRContentDic[seqKey]['variant'] = str(variant)
			isomiRContentDic[seqKey]['strand'] = '+'
			isomiRContentDic[seqKey]['annot'] = 1
			isomiRContentDic[seqKey]['cigar'] = cigarValue
			#if seqKey == 'TTCTCCCACTACCAGGCTCCCA':
			#	print 'Found:'
			#	print preMiRSeqOut
			#	print mirnaLibSeqOut
			#	print mirnaSeqOut
			#	print readSeqOut
			#	print type
			#	print indexValue
			#	print 'Found done'
		

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

def runAnnotationPipeline(bowtieBinary, seqDic, numCPU, phred64, annotNameList, outputdir, logDic, file_mirna, file_hairpin, file_tRNA, file_snoRNA, file_rRNA, file_ncrna_others, file_mrna, spikeIn, file_spikeIn, gff_output, miRNamePreNameDic, isomiRContentDic):
	bwtCmd = os.path.join(bowtieBinary, 'bowtie')+' --threads '+numCPU+' '
	bwtCmd2 = os.path.join(bowtieBinary, 'bowtie-inspect ')
	#bwtCmd = bowtieBinary+"bowtie --threads "+numCPU+' '
	#bwtCmd = "bowtie --threads "+numCPU+' '
	if(phred64):
		bwtCmd = bwtCmd + ' --phred64-quals '
	if spikeIn:
		lengthFilters = [-26, 25, 0, 0, 0, 0, 0, 0, 0]
		bwtCmdLines = [
		bwtCmd+file_mirna+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_hairpin+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_tRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_snoRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_rRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_ncrna_others+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mrna+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mirna+' -5 1 -3 2 -v 2 -f --norc --best -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_spikeIn+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log')]
	else:
		lengthFilters = [-26, 25, 0, 0, 0, 0, 0, 0]
		bwtCmdLines = [
		bwtCmd+file_mirna+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_hairpin+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_tRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_snoRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_rRNA+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_ncrna_others+' -n 1 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mrna+' -n 0 -f --norc -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log'),
		bwtCmd+file_mirna+' -5 1 -3 2 -v 2 -f --norc --best -S '+os.path.join(outputdir, 'SeqToAnnot.fasta')+' 1> '+os.path.join(outputdir, 'SeqToAnnot.sam')+' 2> '+os.path.join(outputdir, 'SeqToAnnot.log')]

	# -- ALIGNMENT 1 -- length < 26, up to 0 mismatch to miRNA
	# -- ALIGNMENT 2 -- length > 25, up to 1 mismatch to hairpin
	# -- ALIGNMENT 3, 4, 5, 6 -- any length, up to 1 mismatch tRNA, snoRNA, rRNA and ncrna_others
	# -- ALIGNMENT 7 -- any length, up to 0 mismatch to mRNA
	# -- ALIGNMENT 8 -- any length, up to 2 mismatches with special 5 vs 3 prime considerations to miRNA (isoMiR)
	# -- ALIGNMENT (9) -- any length, up to 0 mismatches to spike-in 
	if gff_output:
		#if i == 0 or i == 7:
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

	for i in range(len(lengthFilters)):
		writeSeqToAnnot(lengthFilters[i], seqDic, outputdir)
		# run aligment
		#print "Starting Annotation-%s"%(annotNameList[i])
		dicRmp = {}
		time_tmp1 = time.time()
		#print bwtCmdLines[i]
		alignmentStatus = os.system(bwtCmdLines[i])
		#print alignmentStatus
		#alignmentStatus = alignmentStatus >> 8
		time_tmp2 = time.time()
		dicRmp.update({'cpuTime': time_tmp2-time_tmp1})
		
		# parse alignment completed without error
		if alignmentStatus == 0:
			s1, s2 = parseBowtieLog(os.path.join(outputdir, 'SeqToAnnot.log'))
			#print s1
			#print s2
			dicRmp.update({'readsProcessed': s1})
			dicRmp.update({'readsAligned':s2})
			alignmentResult = parseAlignment(os.path.join(outputdir, 'SeqToAnnot.sam'))
			updateAnnotDic(seqDic, alignmentResult, i+1)
			if gff_output:
				if i == 0 or i == 7:
					updateIsomiRDic(isomiRContentDic, alignmentResult, miRNamePreNameDic, i, NT2CODE)
					updateIsomiRDic2(isomiRContentDic, hairpinNameSeqDic, mirnaLibNameSeqDic, i, miRNamePreNameDic)
		else:
			print "Alignment exited with none-zero status.\n"
			sys.exit(1)
		logDic['annotStats'].append(dicRmp)
		#os.system('mv %s %s'%(os.path.join(outputdir, 'SeqToAnnot.fasta'), os.path.join(outputdir, 'SeqToAnnot_'+str(i)+'.fasta')))
		#os.system('mv %s %s'%(os.path.join(outputdir, 'SeqToAnnot.sam'), os.path.join(outputdir, 'SeqToAnnot_'+str(i)+'.sam')))
		#os.system('mv %s %s'%(os.path.join(outputdir, 'SeqToAnnot.log'), os.path.join(outputdir, 'SeqToAnnot_'+str(i)+'.log')))
	# Remove the intermediate files
	os.system('rm %s %s %s'%(os.path.join(outputdir, 'SeqToAnnot.fasta'), os.path.join(outputdir, 'SeqToAnnot.sam'), os.path.join(outputdir, 'SeqToAnnot.log')))
