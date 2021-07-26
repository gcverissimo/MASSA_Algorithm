import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler

def Normalizer_or_matrix(file, biological_activity):
	# Biological activity:
	if len(biological_activity) == 1:
		bio_array = (file[biological_activity].to_numpy(copy=True)).reshape(-1,1)
	else:
		bio_array = file[biological_activity].to_numpy(copy=True)
	scaler = MinMaxScaler()
	bio_array_normalized = scaler.fit_transform(bio_array)
	
	# Physicochemical descriptors:
	descriptors_physicochemical = ['NumHAcceptors', 'NumHDonors', 'ExactMolWt', 'NumRotatableBonds', 'FractionCSP3', 'TPSA', 'LogP_WildmanCrippen']
	PhCh_array = file[descriptors_physicochemical].to_numpy(copy=True)
	scaler = MinMaxScaler()
	PhCh_array_normalized = scaler.fit_transform(PhCh_array)

	# Fingerprint:
	FP_array = np.array(list(file['AtomPairs']))
	return bio_array_normalized, PhCh_array_normalized, FP_array

def PCA_maker(file, nPCS, svd_parameter):
	if (nPCS >= file.shape[1]):
		pca_props = file
	elif (file.shape[1] <= 3):
		pca_props = file
	else:
		pca = PCA(n_components=nPCS, svd_solver=svd_parameter)
		pca_props = pca.fit_transform(file)
	return pca_props

def organize_HCA(file, HCAdict, ident):
	if ident == 'all':
		file['Cluster_General'] = pd.Series(HCAdict)
		return file
	elif ident == 'bio':
		file['Cluster_Biological'] = pd.Series(HCAdict)
		return file
	elif ident == 'PhCh':
		file['Cluster_Physicochemical'] = pd.Series(HCAdict)
		return file
	else:
		file['Cluster_Structural'] = pd.Series(HCAdict)
		return file

def organize_for_kmodes(file):
	kmodes_columns = file[['Cluster_Biological', 'Cluster_Physicochemical', 'Cluster_Structural']]
	dict_cluster_bio = {}
	dict_cluster_PhCh = {}
	dict_cluster_FP = {}
	for i, a, b, c in zip(kmodes_columns.index, kmodes_columns.Cluster_Biological, kmodes_columns.Cluster_Physicochemical, kmodes_columns.Cluster_Structural):
		dict_cluster_bio[i] = 'cluster '+str(a)
		dict_cluster_PhCh[i] = 'cluster '+str(b)
		dict_cluster_FP[i] = 'cluster '+str(c)

	kmodes_new_columns = pd.DataFrame()
	kmodes_new_columns['Cluster_Biological'] = pd.Series(dict_cluster_bio)
	kmodes_new_columns['Cluster_Physicochemical'] = pd.Series(dict_cluster_PhCh)
	kmodes_new_columns['Cluster_Structural'] = pd.Series(dict_cluster_FP)
	kmodes_matrix = np.array(kmodes_new_columns)
	return kmodes_matrix