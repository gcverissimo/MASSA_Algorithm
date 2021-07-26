import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colrs
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering
from kmodes import kmodes

def numerical_number_of_clusters(distances_list):
	shorted_list = distances_list[0:21]
	x1, y1 = 1, shorted_list[0]
	x2, y2 = len(shorted_list)-1, shorted_list[len(shorted_list)-1]
	distances = []
	for i in range(len(shorted_list)):
		x0 = i+1
		y0 = shorted_list[i]
		numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
		denominator = ((y2 - y1)**2 + (x2 - x1)**2)**(1/2)
		distances.append(numerator/denominator)
	return (distances.index(max(distances)) + 1)

def categorical_number_of_clusters(file):
	# Costs = Clustering cost, defined as the sum distance of all points to their respective cluster centroids.
	listClusCosts = []
	for n in range(2,21):
		km = kmodes.KModes(n_clusters=n, verbose=0, n_init=20, init='Cao', max_iter=300)
		km.fit(X=file)
		listClusCosts.append(km.cost_)
	x1, y1 = 2, listClusCosts[0]
	x2, y2 = 20, listClusCosts[len(listClusCosts)-1]

	distances = []
	for i in range(len(listClusCosts)):
		x0 = i+2
		y0 = listClusCosts[i]
		numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
		denominator = ((y2 - y1)**2 + (x2 - x1)**2)**(1/2)
		distances.append(numerator/denominator)
	return (distances.index(max(distances)) + 2)

def optimal_threshold(file, ident, directoryFileOutput, extension_type):
	clustering = AgglomerativeClustering(compute_distances=True, linkage='complete',affinity='euclidean').fit(file)
	distances_list = (sorted(list(clustering.distances_)))
	distances_list.reverse()

	number_of_clusters = numerical_number_of_clusters(distances_list)
	CutOff = distances_list[number_of_clusters - 2]

	if ident == 'bio':
		name_euc_dist_file = directoryFileOutput+'\\Images\\Distance_HCA_Images\\Dist_HCA_biological.'+extension_type
		title_distances = 'Plot of Euclidean Distances from Biological HCA'
		print('Biological:', number_of_clusters)
		plot_euclidean_distance(distances_list, name_euc_dist_file, title_distances)
	elif ident == 'PhCh':
		name_euc_dist_file = directoryFileOutput+'\\Images\\Distance_HCA_Images\\Dist_HCA_physicochemical.'+extension_type
		title_distances = 'Plot of Euclidean Distances from Physicochemical HCA'
		print('Physicochemical:', number_of_clusters)
		plot_euclidean_distance(distances_list, name_euc_dist_file, title_distances)
	else:
		name_euc_dist_file = directoryFileOutput+'\\Images\\Distance_HCA_Images\\Dist_HCA_structural.'+extension_type
		title_distances = 'Plot of Euclidean Distances from Structural HCA'
		print('Structural:', number_of_clusters)
		plot_euclidean_distance(distances_list, name_euc_dist_file, title_distances)

	return CutOff, number_of_clusters

def plot_euclidean_distance(distances_list,name_euc_dist_file,title_distances):
	plt.plot(list(range(1, len(distances_list)+1)), distances_list, '-o', color = 'black', markersize=4, linewidth=1, markeredgewidth=1)
	plt.ylabel('Euclidean distance')
	original_ticks = list((plt.xticks())[0])
	increase_rate = original_ticks[2]
	lista = list(range(1, len(distances_list)+1,int(increase_rate/2)))
	plt.xticks(lista, lista, fontsize=8)
	plt.title(title_distances)
	plt.savefig(name_euc_dist_file)
	plt.close()

def HCA_clusters(file, labels, ident, directoryFileOutput, extension_type):
	# Get the threshold value (CutOff) and the number of clusters:
	CutOff, number_of_clusters = optimal_threshold(file, ident, directoryFileOutput, extension_type)
	linkage = sch.linkage(file, method = 'complete', optimal_ordering=False)

	# Get the cluster labels:
	leaves_cluster = sch.fcluster(linkage, t=number_of_clusters, criterion='maxclust')
	cluster_dict = dict(zip(labels, list(leaves_cluster)))

	return leaves_cluster, cluster_dict, linkage, CutOff

def HCA_plot(linkage, labels, leaves_cluster, CutOff, ident, directoryFileOutput, extension_type, dendrogram_xfont_size, index_of_test_molecules):
	# Working with colors:
	hexRGB = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#691605', '#b3f084', '#ffee00', '#00ff6e', '#0400ff', '#c4c4c0', '#ff00fb', '#290e0e', '#12290e', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#691605', '#b3f084', '#ffee00', '#00ff6e', '#0400ff', '#c4c4c0', '#ff00fb', '#290e0e', '#12290e']
	sch.set_link_color_palette(hexRGB)

	# Create an arrow for molecules in test set:
	labels_arrow = labels.copy()
	for i in range(0,len(labels_arrow)):
		if labels_arrow[i] in index_of_test_molecules:
			labels_arrow[i] = ('→ '+labels_arrow[i])

	# Initial setting of HCA dendrograms:
	plt.figure(figsize=(20, 10), dpi= 300)
	dendrogram = sch.dendrogram(linkage, color_threshold=CutOff, above_threshold_color='k', labels=labels_arrow, get_leaves=True)

	# Get the color for each cluster:
	color_list = list(dendrogram['leaves_color_list'])

	# Set graphs font size for y axis and adjust the plot in image space:
	plt.yticks(fontsize=12)
	plt.xticks(fontsize=dendrogram_xfont_size)
	plt.ylabel('Euclidean distance',fontsize=12)
	plt.subplots_adjust(left=0.06 ,right=0.96)

	# Set new scale for y axis (which is 2x the previous scale):
	b = [float(i) for i in plt.yticks()[0]]
	scale = b[1]-b[0]
	new_scale = scale/2
	new_arange = np.arange(min(b)-scale/20,max(b)+new_scale,new_scale)
	plt.yticks(new_arange)

	# Assign the legend to each cluster:
	set__leaves_cluster = sorted(list(set(leaves_cluster)))
	set_leaves_cluster = []
	for i in set__leaves_cluster:
		if i < 10:
			set_leaves_cluster.append('cluster 0'+str(i))
		else:
			set_leaves_cluster.append('cluster '+str(i))
	plt.legend(set_leaves_cluster)

	# Add ticks to x-axis:
	plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)

	# Apply cluster color to x-labels:
	zip_label_color = list(zip(plt.gca().get_xticklabels(), color_list))
	for i,a in zip_label_color:
		i.set_color(a)

	# Defines the title of the graph and saves the figure:
	if ident == 'bio':
		plt.title('HCA of biological activity')
		plt.savefig(directoryFileOutput+'\\Images\\Dendrogram_HCA_Images\\HCA_biological.'+extension_type)
		plt.close()
	elif ident == 'PhCh':
		plt.title('HCA of physicochemical properties')
		plt.savefig(directoryFileOutput+'\\Images\\Dendrogram_HCA_Images\\HCA_physicochemical.'+extension_type)
		plt.close()
	else:
		plt.title('HCA of AtomPairs fingerprint')
		plt.savefig(directoryFileOutput+'\\Images\\Dendrogram_HCA_Images\\HCA_structural.'+extension_type)
		plt.close()

def KModes(file, names):
	n_clusters = categorical_number_of_clusters(file)
	print('General:', n_clusters)
	km_test = kmodes.KModes(n_clusters=n_clusters, n_init=20, init='Cao', max_iter=300)
	clusters = km_test.fit_predict(file)
	clusters_list = list(clusters)
	HCA_General_dict = {}
	for i, a in zip(names, clusters_list):
		HCA_General_dict[i] = int(a+1)
	return HCA_General_dict