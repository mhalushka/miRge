# This version of convert2Fasta.py only extract the seq with the count larg[0]er than 1 and the length between 16 and 25.
# The raw fasta files will be also generated as well.
# Meanwhile, all the sequences are collapsed.
# Here, it's better to use the absolute path of the file "mapped.csv" and the outputDir

import os
import sys

def convert2Fasta(infTmp, minLengthTmp, maxLengthTmp, countCutoffTmp, outputDirTmp, species, speciesNameDic, spikeIn):
	# convert2Fasta should have 5 fixed parameters as follows:
	# mapped.csv(or unmapped.csv) minLength(16) maxLength(25) countCutoff(2) outputDir
	if spikeIn:
		startP = 11
	else:
		startP = 10
	minLength = int(minLengthTmp)
	maxLength = int(maxLengthTmp)
	countCutoff = int(countCutoffTmp)
	abbrName = speciesNameDic[species].lower()
	inf = open(infTmp,"r")
	outf = open(os.path.join(outputDirTmp, os.path.splitext(infTmp)[0].split('/')[-1]+"_mirna.fa"), "w")
	outf1 = open(os.path.join(outputDirTmp, os.path.splitext(infTmp)[0].split('/')[-1]+"_mirna_raw.fa"), "w")
	if os.path.splitext(infTmp)[0].split('/')[-1] == 'mapped':
		outf2 = open(os.path.join(outputDirTmp, os.path.splitext(infTmp)[0].split('/')[-1]+"_nonMirna.fa"), "w")
		outf3 = open(os.path.join(outputDirTmp, os.path.splitext(infTmp)[0].split('/')[-1]+"_nonMirna_raw.fa"), "w")
		contentList2 = []
		contentList3 = []
		s = 1
		t = 1
	line = inf.readline()
	line1 = line.strip()
	if line1[-1] == ",":
		line1 = line1[:-1]
	else:
		pass
	fileList = []
	for item1 in line1.split(",")[startP:]:
		fileList.append(item1.strip().split(".")[0])
	line = inf.readline()
	i = 1
	j = 1
	contentList = []
	contentList1 = []
	while line != "":
		newline = line.strip()
		if newline[-1] == ",":
			newline = newline[:-1]
		else:
			pass
		temp = newline.split(",")
		count1 = 0
		for item in temp[startP:]:
			count1 = count1 + int(item.strip())
		if temp[1].strip() == "1":
			if count1 >= 1:
				if abbrName in temp[2] or abbrName in temp[9]:
					name = "mir"+str(j)+"_"+str(count1)
					outf1.write(">")
					outf1.write(name)
					outf1.write("\n")
					outf1.write(temp[0].strip())
					outf1.write("\n")
					tempList = ["mir"+str(j)]
					tempList.append(temp[0].strip())
					for subitem in temp[startP:]:
						tempList.append(subitem.strip())
					contentList1.append(tempList)
					j = j + 1
				else:
					if abbrName not in temp[3]:
						if spikeIn :
							if temp[10] == '':
								name = "mir"+str(t)+"_"+str(count1)
								outf3.write(">")
								outf3.write(name)
								outf3.write("\n")
								outf3.write(temp[0].strip())
								outf3.write("\n")
								tempList = ["mir"+str(t)]
								tempList.append(temp[0].strip())
								for subitem in temp[startP:]:
									tempList.append(subitem.strip())
								contentList3.append(tempList)
								t = t + 1
							else:
								pass
						else:
							name = "mir"+str(t)+"_"+str(count1)
							outf3.write(">")
							outf3.write(name)
							outf3.write("\n")
							outf3.write(temp[0].strip())
							outf3.write("\n")
							tempList = ["mir"+str(t)]
							tempList.append(temp[0].strip())
							for subitem in temp[startP:]:
								tempList.append(subitem.strip())
							contentList3.append(tempList)
							t = t + 1
					else:
						pass
			else:
				pass
		elif temp[1].strip() == "0":
			if count1 >= 1:
				name = "mir"+str(j)+"_"+str(count1)
				outf1.write(">")
				outf1.write(name)
				outf1.write("\n")
				outf1.write(temp[0].strip())
				outf1.write("\n")
				tempList = ["mir"+str(j)]
				tempList.append(temp[0].strip())
				for subitem in temp[startP:]:
					tempList.append(subitem.strip())
				contentList1.append(tempList)
				j = j + 1
			else:
				pass
		else:
			pass
		if len(temp[0].strip()) >= minLength and len(temp[0].strip()) <= maxLength:
			count=0
			for item in temp[startP:]:
				count = count + int(item.strip())
			if temp[1].strip() == "1":
				if count >= countCutoff:
					if abbrName in temp[2] or abbrName in temp[9]:
						name = "mir"+str(i)+"_"+str(count)
						outf.write(">")
						outf.write(name)
						outf.write("\n")
						outf.write(temp[0].strip())
						outf.write("\n")
						tempList = ["mir"+str(i)]
						tempList.append(temp[0].strip())
						for subitem in temp[startP:]:
							tempList.append(subitem.strip())
						contentList.append(tempList)
						i = i + 1
					else:
						if abbrName not in temp[3]:
							if spikeIn :
								if temp[10] == '':
									name = "mir"+str(s)+"_"+str(count)
									outf2.write(">")
									outf2.write(name)
									outf2.write("\n")
									outf2.write(temp[0].strip())
									outf2.write("\n")
									tempList = ["mir"+str(s)]
									tempList.append(temp[0].strip())
									for subitem in temp[startP:]:
										tempList.append(subitem.strip())
									contentList2.append(tempList)
									s = s + 1
								else:
									pass
							else:
								name = "mir"+str(s)+"_"+str(count)
								outf2.write(">")
								outf2.write(name)
								outf2.write("\n")
								outf2.write(temp[0].strip())
								outf2.write("\n")
								tempList = ["mir"+str(s)]
								tempList.append(temp[0].strip())
								for subitem in temp[startP:]:
									tempList.append(subitem.strip())
								contentList2.append(tempList)
								s = s + 1
						else:
							pass
				else:
					pass
			elif temp[1].strip() == "0":
				if count >= countCutoff:
					name = "mir"+str(i)+"_"+str(count)
					outf.write(">")
					outf.write(name)
					outf.write("\n")
					outf.write(temp[0].strip())
					outf.write("\n")
					tempList = ["mir"+str(i)]
					tempList.append(temp[0].strip())
					for subitem in temp[startP:]:
						tempList.append(subitem.strip())
					contentList.append(tempList)
					i = i + 1
				else:
					pass
			else:
				pass
		else:
			pass
		line = inf.readline()
	inf.close()
	outf.close()
	outf1.close()
	for k in range(len(fileList)):
		outf = open(os.path.join(outputDirTmp, os.path.splitext(infTmp)[0].split('/')[-1]+"_mirna_"+fileList[k]+".fa"), "w")
		for item in contentList:
			count = int(item[2+k])
			if count >= countCutoff:
				outf.write(">"+item[0]+"_"+item[2+k]+"\n")
				outf.write(item[1]+"\n")
		outf.close()
	for k in range(len(fileList)):
		outf = open(os.path.join(outputDirTmp, os.path.splitext(infTmp)[0].split('/')[-1]+"_mirna_"+fileList[k]+"_raw.fa"), "w")
		for item in contentList1:
			count = int(item[2+k])
			if count >= 1:
				outf.write(">"+item[0]+"_"+item[2+k]+"\n")
				outf.write(item[1]+"\n")
		outf.close()
		
	if os.path.splitext(infTmp)[0].split('/')[-1] == 'mapped':
		for k in range(len(fileList)):
			outf = open(os.path.join(outputDirTmp, os.path.splitext(infTmp)[0].split('/')[-1]+"_nonMirna_"+fileList[k]+".fa"), "w")
			for item in contentList2:
				count = int(item[2+k])
				if count >= countCutoff:
					outf.write(">"+item[0]+"_"+item[2+k]+"\n")
					outf.write(item[1]+"\n")
			outf.close()
		for k in range(len(fileList)):
			outf = open(os.path.join(outputDirTmp, os.path.splitext(infTmp)[0].split('/')[-1]+"_nonMirna_"+fileList[k]+"_raw.fa"), "w")
			for item in contentList3:
				count = int(item[2+k])
				if count >= 1:
					outf.write(">"+item[0]+"_"+item[2+k]+"\n")
					outf.write(item[1]+"\n")
			outf.close()
