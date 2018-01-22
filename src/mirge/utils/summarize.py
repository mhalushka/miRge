import os

def summarize(seqDic, sampleList, logDic, mirDic, file_mirna, outputdir, spikeIn, bowtieBinary):
	mirList = []
	os.system(os.path.join(bowtieBinary, 'bowtie-inspect')+' -n '+file_mirna+' > '+os.path.join(outputdir, 'mirnaList.tmp'))
	with open(os.path.join(outputdir, 'mirnaList.tmp'), 'r') as inf:
		for line in inf:
			mirList.append(line.strip())
	# mirList stores the name of miRNAs
	os.system('rm %s'%(os.path.join(outputdir, 'mirnaList.tmp')))
	for i in range(len(mirList)):
		mirDic.update({mirList[i]:{}})
		mirDic[mirList[i]].update({'quant': [0]*len(sampleList)})
		mirDic[mirList[i]].update({'iscan': [0]*len(sampleList)})
	for i in range(len(sampleList)):
		logDic['quantStats'][i].update({'trimmedUniq': 0})
		logDic['quantStats'][i].update({'mirnaReads': 0})
		logDic['quantStats'][i].update({'hairpinReads': 0})
		logDic['quantStats'][i].update({'trnaReads': 0})
		logDic['quantStats'][i].update({'snornaReads': 0})
		logDic['quantStats'][i].update({'rrnaReads': 0})
		logDic['quantStats'][i].update({'ncrnaOthersReads': 0})
		logDic['quantStats'][i].update({'mrnaReads': 0})
		logDic['quantStats'][i].update({'remReads': 0})
		if spikeIn:
			logDic['quantStats'][i].update({'spikeInReads': 0})

	for seqKey in seqDic.keys():
		for i in range(len(sampleList)):
			if seqDic[seqKey]['quant'][i] != 0:
				logDic['quantStats'][i]['trimmedUniq'] = logDic['quantStats'][i]['trimmedUniq']+1
				if seqDic[seqKey]['annot'][1] != '' or seqDic[seqKey]['annot'][8] != '':
					logDic['quantStats'][i]['mirnaReads'] = logDic['quantStats'][i]['mirnaReads'] + seqDic[seqKey]['quant'][i]
					if seqDic[seqKey]['annot'][1] != '':
						mirDic[seqDic[seqKey]['annot'][1]]['quant'][i] = mirDic[seqDic[seqKey]['annot'][1]]['quant'][i] + seqDic[seqKey]['quant'][i]
						mirDic[seqDic[seqKey]['annot'][1]]['iscan'][i] = mirDic[seqDic[seqKey]['annot'][1]]['iscan'][i] + seqDic[seqKey]['quant'][i]
					else:
						mirDic[seqDic[seqKey]['annot'][8]]['quant'][i] = mirDic[seqDic[seqKey]['annot'][8]]['quant'][i] + seqDic[seqKey]['quant'][i]
				elif seqDic[seqKey]['annot'][2] != '':
					logDic['quantStats'][i]['hairpinReads'] = logDic['quantStats'][i]['hairpinReads'] + seqDic[seqKey]['quant'][i]
				elif seqDic[seqKey]['annot'][3] != '':
					logDic['quantStats'][i]['trnaReads'] = logDic['quantStats'][i]['trnaReads'] + seqDic[seqKey]['quant'][i]
				elif seqDic[seqKey]['annot'][4] != '':
					logDic['quantStats'][i]['snornaReads'] = logDic['quantStats'][i]['snornaReads'] + seqDic[seqKey]['quant'][i]
				elif seqDic[seqKey]['annot'][5] != '':
					logDic['quantStats'][i]['rrnaReads'] = logDic['quantStats'][i]['rrnaReads'] + seqDic[seqKey]['quant'][i]
				elif seqDic[seqKey]['annot'][6] != '':
					logDic['quantStats'][i]['ncrnaOthersReads'] = logDic['quantStats'][i]['ncrnaOthersReads'] + seqDic[seqKey]['quant'][i]
				elif seqDic[seqKey]['annot'][7] != '':
					logDic['quantStats'][i]['mrnaReads'] = logDic['quantStats'][i]['mrnaReads'] + seqDic[seqKey]['quant'][i]
				else:
					if spikeIn:
						if seqDic[seqKey]['annot'][9] != '':
							logDic['quantStats'][i]['spikeInReads'] = logDic['quantStats'][i]['spikeInReads'] + seqDic[seqKey]['quant'][i]
						else:
							logDic['quantStats'][i]['remReads'] = logDic['quantStats'][i]['remReads'] + seqDic[seqKey]['quant'][i]
					else:
						logDic['quantStats'][i]['remReads'] = logDic['quantStats'][i]['remReads'] + seqDic[seqKey]['quant'][i]
