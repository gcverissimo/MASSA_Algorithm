import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt
from matplotlib import colors as colrs

def split_train_test_sets(file, training_percent, test_percent):
	training, test = train_test_split(file, test_size=test_percent, train_size=training_percent, stratify=file['Cluster_General'])
	index_of_test_molecules = [i for i in test.index]
	dict_training = {'set':'training'}
	dict_test = {'set': 'test'}
	training = training.assign(**dict_training)
	test = test.assign(**dict_test)
	dataset = pd.concat([training, test])
	return dataset, index_of_test_molecules

def freq_clusters(file, directoryFileOutput, extension_type, Cluster_Domain, barplot_xfont_size):
	Total_count = []
	Training_count = []
	Test_count = []
	for index, row in file.iterrows():
		if row['set'] == 'training':
			Training_count.append(row[Cluster_Domain])
			Total_count.append(row[Cluster_Domain])
		if row['set'] == 'test':
			Test_count.append(row[Cluster_Domain])
			Total_count.append(row[Cluster_Domain])
	min_value = min(Total_count)
	max_value = max(Total_count)
	
	clusters_total = {}
	clusters_training = {}
	clusters_test = {}

	for i in range(min_value, (max_value+1)):
		if i <= 9:
			clusters_total[('cluster 0'+str(i))] = ((Total_count.count(i))*100)/len(Total_count)
			clusters_training[('cluster 0'+str(i))] = ((Training_count.count(i))*100)/len(Training_count)
			clusters_test[('cluster 0'+str(i))] = ((Test_count.count(i))*100)/len(Test_count)
		else:
			clusters_total[('cluster '+str(i))] = ((Total_count.count(i))*100)/len(Total_count)
			clusters_training[('cluster '+str(i))] = ((Training_count.count(i))*100)/len(Training_count)
			clusters_test[('cluster '+str(i))] = ((Test_count.count(i))*100)/len(Test_count)

	list_total = sorted(list(clusters_total.items()), key=lambda x:x[0])
	list_training = sorted(list(clusters_training.items()), key=lambda x:x[0])
	list_test = sorted(list(clusters_test.items()), key=lambda x:x[0])

	label_total = [list_total[i][0] for i in range(0,len(list_total))]
	percentages_total = [list_total[i][1] for i in range(0,len(list_total))]
	percentages_training = [list_total[i][1] for i in range(0,len(list_training))]
	percentages_test = [list_test[i][1] for i in range(0,len(list_test))]

	Total_ar = np.arange(len(list_total))
	Training_ar = np.arange(len(list_training))
	Test_ar = np.arange(len(list_test))
	
	plt.rcParams.update({'font.size':12})
	largura = 0.25
	fig,ax = plt.subplots(figsize = (10,6))
	group1 = ax.bar(Total_ar-largura, percentages_total, largura, label = 'Total', color='black', edgecolor = 'black')
	group2 = ax.bar(Training_ar, percentages_training, largura, label = 'Training set', color='dimgrey', edgecolor = 'black')
	group3 = ax.bar(Test_ar+largura, percentages_test, largura, label = 'Test set', color= 'whitesmoke', edgecolor = 'black')
	ax.legend()
	ax.set_ylabel('Frequency (%)')
	ax.set_ylim([0,100])
	ax.set_xticks(Total_ar)
	ax.set_xticklabels(label_total, fontsize=barplot_xfont_size)

	if Cluster_Domain == 'Cluster_Biological':
		ax.set_title('Distribution of dataset: Biological Clustering')
		plt.savefig(directoryFileOutput+'\\Images\\Distribution\\Biological.'+extension_type)
	elif Cluster_Domain == 'Cluster_Structural':
		ax.set_title('Distribution of dataset: Structural Clustering')
		plt.savefig(directoryFileOutput+'\\Images\\Distribution\\Structural.'+extension_type)
	elif Cluster_Domain == 'Cluster_Physicochemical':
		ax.set_title('Distribution of dataset: Physicochemical Clustering')
		plt.savefig(directoryFileOutput+'\\Images\\Distribution\\Physicochemical.'+extension_type)
	elif Cluster_Domain == 'Cluster_General':
		ax.set_title('Distribution of dataset: General Clustering')
		plt.savefig(directoryFileOutput+'\\Images\\Distribution\\General.'+extension_type)
	plt.close()
	return clusters_total, clusters_training, clusters_test