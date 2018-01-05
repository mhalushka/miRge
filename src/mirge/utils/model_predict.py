import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('agg')
from sklearn.preprocessing import OneHotEncoder
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, chi2, f_classif
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.metrics import make_scorer, confusion_matrix, matthews_corrcoef, roc_auc_score, roc_curve, auc
#from sklearn.model_selection import cross_validation
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import time
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.externals import joblib

def model_predict(infTmp, modelFile):
	tmpName = '_'.join(os.path.basename(infTmp).split('_')[:-9])+'_dataset_'+os.path.basename(infTmp).split('_')[-1].split('.')[0]
	inputFile = os.path.join(os.path.dirname(infTmp), tmpName+'_refined_features.csv')
	cutoff_picked = 0.8
	sc, clf, selectFeatureNameList = joblib.load(modelFile)
	data = pd.read_csv(inputFile)
	dirTmp = os.path.dirname(inputFile)
	sampleName = '_'.join(os.path.basename(inputFile).split('_')[2:-5])
	data_raw = data.iloc[:,:].values
	featureList_raw = data.columns.values.tolist()[:3]
	featureList_raw.insert(0, 'probability')
	featureList_raw.insert(0, 'predictedValue')
	# Encode categorical features using one-hot-coding after deleting the 2nd and 3rd column in data
	data.drop(data.columns[[0,1,2]], axis=1, inplace=True)
	data = pd.get_dummies(data)
	# Check the dummy varibles, such as 'pair_state_Yes', 'pair_state_No', etc.
	# If some dummy varibles are missing, they will be added with value of 0 across all of the samples.
	totalfeatureListTmp = data.columns.values.tolist()
	missedFeatureList = []
	for item in selectFeatureNameList:
		if item[2] not in totalfeatureListTmp:
			missedFeatureList.append(item[2])
	for item in missedFeatureList:
		data[item] = np.array([0]*(data.shape[0]))
	# Subselect data with the selected features and reorder the feature sequence according to the feature sequence in selectFeatureNameList
	totalfeatureList = data.columns.values.tolist()
	subIndexList = [totalfeatureList.index(item[2]) for item in selectFeatureNameList]
	
	data_x = data.iloc[:,subIndexList].values
	x_test_std_new = sc.transform(data_x)
	#featureLabel_new = subIndexList
	x_test_std_new2 = x_test_std_new
	
	#featureList = data.columns.values.tolist()
	#data_x = data.values
	#x_test_std_new = sc.transform(data_x)
	# Reorder the feature sequence according to the feature sequence in x_train_std_filter_select and x_test_std_filter_select
	# which is also selectFeatureNameList
	#featureLabel_new = [featureList.index(item) for item in selectFeatureNameList]
	#x_test_std_new2 = x_test_std_new[:, featureLabel_new]
	# Predict the x_test_std_new2 using the predictive model
	predictedValue = clf.predict(x_test_std_new2)
	probability = [max(value) for value in clf.predict_proba(x_test_std_new2)]
	# Output the predicted result of the test data set into a file.
	array_predicted = np.column_stack((predictedValue, probability, data_raw[:,:3]))
	data_predicted_test = pd.DataFrame(array_predicted, columns=featureList_raw)
	data_predicted_test.to_csv(os.path.join(dirTmp, 'predicted_'+sampleName+'_detailed.csv'), index=False)
	outf = open(os.path.join(dirTmp, sampleName+'_novel_miRNAs_miRge2.0.csv'), 'w')
	contentList_tmp = []
	with open(os.path.join(dirTmp, 'predicted_'+sampleName+'_detailed.csv'), 'r') as inf:
		line = inf.readline()
		outf.write(line)
		line = inf.readline()
		while line != '':
			content = line.strip().split(',')
			if content[0] == '1' and float(content[1]) >= cutoff_picked:
				contentList_tmp.append([float(content[1]), line])
			line = inf.readline()
	contentList_tmp.sort(reverse=True)
	for item in contentList_tmp:
		outf.write(item[1])
	outf.close()
