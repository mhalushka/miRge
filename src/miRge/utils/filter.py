import sys

def filter(mirDic, sampleList, logDic):
	iscanFilter = 2
	for mirKey in mirDic.keys():
		for i in range(len(sampleList)):
			if mirDic[mirKey]['iscan'][i] < iscanFilter:
				mirDic[mirKey]['quant'][i] = 0
	for mirKey in mirDic.keys():
		for i in range(len(sampleList)):
			if mirDic[mirKey]['quant'][i] > 0:
				try:
					logDic['quantStats'][i]['mirnaReadsFiltered'] = logDic['quantStats'][i]['mirnaReadsFiltered'] + mirDic[mirKey]['quant'][i]
					logDic['quantStats'][i]['mirnaUniqFiltered'] = logDic['quantStats'][i]['mirnaUniqFiltered'] + 1
				except KeyError:
					logDic['quantStats'][i].update({'mirnaReadsFiltered': 0})
					logDic['quantStats'][i].update({'mirnaUniqFiltered': 0})
					logDic['quantStats'][i]['mirnaReadsFiltered'] = logDic['quantStats'][i]['mirnaReadsFiltered'] + mirDic[mirKey]['quant'][i]
					logDic['quantStats'][i]['mirnaUniqFiltered'] = logDic['quantStats'][i]['mirnaUniqFiltered'] + 1
	for i in range(len(sampleList)):
		if logDic['quantStats'][i]['mirnaReadsFiltered'] == 0:
			print 'No miRNA reads were found in sample %s. Please check your files and provided arguments.\n'%(sampleFiles[i])
			sys.exit(1)