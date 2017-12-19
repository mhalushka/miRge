import os

def miRNAmerge(mergeLibFile, sampleList, mirDic, mirna_fa_tmp, mirNameSeqDic):
	with open(mirna_fa_tmp, 'r') as inf:
		line = inf.readline()
		while line != '':
			mirName = line.strip()[1:]
			line = inf.readline()
			sequence = line.strip()
			mirNameSeqDic.update({mirName:sequence})
			line = inf.readline()

	if os.path.isfile(mergeLibFile):
		delMirNameList = []
		with open(mergeLibFile, 'r') as inf:
			line = inf.readline()
			while line != '':
				contentTmp = line.strip().split(',')
				for i in range(1, len(contentTmp)):
					for j in range(len(sampleList)):
						try:
							if mirDic[contentTmp[i]]['quant'][j] > 0:
								try:
									mirDic[contentTmp[0]]['quant'][j] = mirDic[contentTmp[0]]['quant'][j] + mirDic[contentTmp[i]]['quant'][j]
									mirDic[contentTmp[0]]['iscan'][j] = mirDic[contentTmp[0]]['iscan'][j] + mirDic[contentTmp[i]]['iscan'][j]
								except KeyError:
									mirDic.update({contentTmp[0]:{}})
									mirDic[contentTmp[0]].update({'quant': [0]*len(sampleList)})
									mirDic[contentTmp[0]].update({'iscan': [0]*len(sampleList)})
									mirDic[contentTmp[0]]['quant'][j] = mirDic[contentTmp[i]]['quant'][j]
									mirDic[contentTmp[0]]['iscan'][j] = mirDic[contentTmp[i]]['iscan'][j]
							#del mirDic[contentTmp[i]]
							delMirNameList.append(contentTmp[i])
						except KeyError:
							pass
				line = inf.readline()
		for item in set(delMirNameList):
			try:
				del mirDic[item]
			except KeyError:
				pass 
	else:
		print "Cannot find merges file, skipping merge step.\n"