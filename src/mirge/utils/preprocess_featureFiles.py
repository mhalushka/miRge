import os

def preprocess_featureFiles(infTmp, feature_namelist_file):
	# infTmp is unmapped_mirna_JH-29_Pros_SMC_vs_representative_seq_modified_selected_sorted_features_updated_stableClusterSeq_15.tsv
	readCountLimit = 10
	seqCountLimit = 3
	tmpName = '_'.join(os.path.basename(infTmp).split('_')[:-9])+'_dataset_'+os.path.basename(infTmp).split('_')[-1].split('.')[0]
	with open(os.path.join(os.path.dirname(infTmp), tmpName+'.csv'), 'w') as outf:
		with open(infTmp, 'r') as inf:
			line = inf.readline()
		content = line.strip().split('\t')
		retainedIndexList = [0, 1, 5] + range(10, len(content))
		retainedColNameList = [content[i] for i in range(len(content)) if i in retainedIndexList]
		readCountLimitLabel = content.index('readCountSum')
		seqCountLimitLabel = content.index('seqCount')

		with open(infTmp, 'r') as inf:
			line = inf.readline()
			content = line.strip().split('\t')
			outf.write(','.join([content[i] for i in range(len(content)) if i in retainedIndexList]))
			outf.write('\n')
			line = inf.readline()
			while line != '':
				content = line.strip().split('\t')
				if int(content[readCountLimitLabel]) >= readCountLimit and int(content[seqCountLimitLabel]) >= seqCountLimit and ('N' not in content[seqCountLimitLabel:]):
					outf.write(','.join([content[i] for i in range(len(content)) if i in retainedIndexList]))
					outf.write('\n')
				line = inf.readline()
			line = inf.readline()

	outf = open(os.path.join(os.path.dirname(infTmp), tmpName+'_refined_tmp.csv'), 'w')
	#outf = open(os.path.splitext(infTmp)[0]+'_refined_tmp.csv', 'w')
	with open(os.path.join(os.path.dirname(infTmp), tmpName+'.csv'), 'r') as inf:
		line = inf.readline()
		outf.write(line)
		content = line.strip().split(',')
		pruned_precusor_strLabel = content.index('pruned_precusor_str')
		line = inf.readline()
		while line != '':
			content = line.strip().split(',')
			if content[pruned_precusor_strLabel] != 'None':
				outf.write(line)
			else:
				pass
			line = inf.readline()
	outf.close()

	nameList = []
	with open(feature_namelist_file, 'r') as inf:
		for line in inf:
			if line.strip() not in nameList:
				nameList.append(line.strip())

	outf = open(os.path.join(os.path.dirname(infTmp), tmpName+'_refined_features.csv'), 'w')
	with open(os.path.join(os.path.dirname(infTmp), tmpName+'_refined_tmp.csv'), 'r') as inf:
		retainedIndexList = []
		line = inf.readline()
		content = line.strip().split(',')
		reatinedContent = []
		for i, name in enumerate(content):
			if name in nameList:
				retainedIndexList.append(i)
				reatinedContent.append(name)
		outf.write(','.join(reatinedContent))
		outf.write('\n')
		line = inf.readline()
		while line != '':
			content = line.strip().split(',')
			reatinedContent = []
			for j in range(len(content)):
				if j in retainedIndexList:
					reatinedContent.append(content[j])
			outf.write(','.join(reatinedContent))
			outf.write('\n')
			line = inf.readline()
	outf.close()