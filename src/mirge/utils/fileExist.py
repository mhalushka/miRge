def fileExist(fileName, sampleList):
	state = False
	sampleListNew = [item.split('.')[0] for item in sampleList]
	for item in sampleListNew:
		if item+'.fa' in fileName or item+'_raw.fa' in fileName:
			state = True
			break
	return state