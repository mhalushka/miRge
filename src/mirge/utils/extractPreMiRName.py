def allDigit(list):
	out = True
	for item in list:
		if not item.isdigit():
			out = False
			break
	return out

def pickOptimal(list):
	listTmp = [item.split('-')[-1] for item in list]
	if allDigit(listTmp):
		newList = [[int(listTmp[i]), list[i]] for i in range(len(list))]
	else:
		newList = [[list[i]] for i in range(len(list))]
	newList.sort()
	optimal = newList[0][-1]
	return optimal

def extractPreMiRName(miRNA_coordinate):
	miRNAMirCorDic = {}
	miRNAMirList = []
	stemLoopNameDic = {}
	stemLoopAliasNameDic = {}
	with open(miRNA_coordinate, 'r') as inf:
		for line in inf:
			if line[0] != '#':
				content = line.strip().split('\t')
				if content[2] == 'miRNA':
					name = content[-1].split('Name=')[1].split(';')[0].strip()
					chr = content[0]
					strand = content[6]
					start = int(content[3])-1
					end = int(content[4])-1
					stemLoopName = content[-1].split('Derives_from=')[1].split(';')[0].strip()
					if name not in miRNAMirList:
						miRNAMirList.append(name)
						miRNAMirCorDic.update({name:{stemLoopName:(chr,strand,start,end)}})
					else:
						miRNAMirCorDic[name].update({stemLoopName:(chr,strand,start,end)})
				if content[2] == 'miRNA_primary_transcript':
					stemLoopName = content[-1].split('Alias=')[1].split(';')[0].strip()
					stemLoopName2 = content[-1].split('Name=')[1].split(';')[0].strip()
					#if stemLoopName not in stemLoopNameDic.keys():
					stemLoopNameDic.update({stemLoopName:stemLoopName2})
					stemLoopAliasNameDic.update({stemLoopName2:stemLoopName})
	miRNamePreNameDic = {}
	for item in miRNAMirList:
		if len(miRNAMirCorDic[item]) > 1:
			stemLoopName = pickOptimal([stemLoopNameDic[subitem] for subitem in miRNAMirCorDic[item].keys()])
		else:
			stemLoopName = stemLoopNameDic[miRNAMirCorDic[item].keys()[0]]
		miRNamePreNameDic.update({item:stemLoopName})
	return miRNamePreNameDic