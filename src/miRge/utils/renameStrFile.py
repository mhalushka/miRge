import os


def renameStrFile(fastaFile, strFile, strFileNew):
	fastaNameList = []
	with open(fastaFile, 'r') as inf:
		line = inf.readline()
		while line != '':
			fastaNameList.append(line)
			line = inf.readline()
			line = inf.readline()
	outf = open(strFileNew, 'w')
	i = 0
	with open(strFile, 'r') as inf:
		line = inf.readline()
		while line != '':
			outf.write(fastaNameList[i])
			line = inf.readline()
			outf.write(line)
			line = inf.readline()
			outf.write(line)
			line = inf.readline()
			i = i + 1
	outf.close()