import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colrs
import scipy.cluster.hierarchy as sch
from kmodes import kmodes
from sklearn.cluster import MiniBatchKMeans as KMeans
from multiprocessing import cpu_count as mpi_count
from scipy.spatial.distance import cdist


def numerical_number_of_clusters(distances_list):
    # Elbow method using the Euclidean distance:
    # Create a list of distances from the first 21 clusters.
    shorted_list = distances_list[0:21]
    # Define the x-axis coordinates of the first and last points on the line.
    x1, y1 = 1, shorted_list[0]
    # Define the y-axis coordinates of the first and last points on the line.
    x2, y2 = len(shorted_list)-1, shorted_list[len(shorted_list)-1]
    distances = []
    for i, label in enumerate(shorted_list):
        # Coordinate on the x-axis of the point to calculate the distance to the line.
        x0 = i+1
        # Coordinate on the y-axis of the point to calculate the distance to the line.
        y0 = label
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = ((y2 - y1)**2 + (x2 - x1)**2)**(1/2)
        # Distance between a point and a line.
        distances.append(numerator/denominator)
    return distances.index(max(distances)) + 1


def categorical_number_of_clusters(file, directory_fileoutput, extension_type):
    # Elbow method using the clustering cost:
    list_clus_costs = []
    maxcluster = 31
    ran_curve = [x for x in range(2, maxcluster, 2)]
    for n in range(2, maxcluster):  # Calculate the cost for each n cluster.
        km = kmodes.KModes(n_clusters=n, verbose=0,
                           n_init=20, init='Cao',
                           n_jobs=-1, max_iter=300)
        km.fit(X=file)
        # Costs = Clustering cost, defined as the sum distance of all points
        # to their respective cluster centroids.
        list_clus_costs.append(km.cost_)
    # Define the x-axis coordinates of the first and last points on the line.
    x1, y1 = 2, list_clus_costs[0]
    # Define the y-axis coordinates of the first and last points on the line.
    x2, y2 = 20, list_clus_costs[18]

    distances = []
    for i in range(19):
        # Coordinate on the x-axis of the point to calculate the distance to the line.
        x0 = i+2
        # Coordinate on the y-axis of the point to calculate the distance to the line.
        y0 = list_clus_costs[i]
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = ((y2 - y1)**2 + (x2 - x1)**2)**(1/2)
        # Distance between a point and a line.
        distances.append(numerator/denominator)
    number_of_clusters = (distances.index(max(distances)) + 2)
    ident = "all"
    plot_distortion(
        list_clus_costs, directory_fileoutput,
        ident, extension_type, number_of_clusters,
        ran_curve)
    return number_of_clusters


def kmeans_number_of_clusters(file, ident, directory_fileoutput, extension_type):
    # Elbow method using the clustering cost:
    list_clus_costs = []
    maxcluster = 31
    ran_curve = [x for x in range(2, maxcluster, 2)]
    for n in range(2, maxcluster):  # Calculate the cost for each n cluster.
        km = KMeans(n_clusters=n, verbose=0,
                    init='k-means++', max_iter=300,
                    batch_size=256*int(mpi_count()),
                    compute_labels=True, random_state=2025)
        km.fit(X=file)
        # Costs = Clustering cost, defined as the sum distance of all points
        # to their respective cluster centroids.
        list_clus_costs.append(
            np.sum(np.min(cdist(file, km.cluster_centers_, 'euclidean'), axis=1)**2) / file.shape[0])
    # Define the x-axis coordinates of the first and last points on the line.
    x1, y1 = 2, list_clus_costs[0]
    # Define the y-axis coordinates of the first and last points on the line.
    x2, y2 = 20, list_clus_costs[18]

    distances = []
    for i in range(19):
        # Coordinate on the x-axis of the point to calculate the distance to the line.
        x0 = i+2
        # Coordinate on the y-axis of the point to calculate the distance to the line.
        y0 = list_clus_costs[i]
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = ((y2 - y1)**2 + (x2 - x1)**2)**(1/2)
        # Distance between a point and a line.
        distances.append(numerator/denominator)

    number_of_clusters = (distances.index(max(distances)) + 2)
    plot_distortion(
        list_clus_costs, directory_fileoutput,
        ident, extension_type, number_of_clusters,
        ran_curve)
    return number_of_clusters


def optimal_threshold(file, ident, directory_fileoutput, extension_type, linkage):
    # Calculate the maximum distance between clusters.
    distances_list = sch.maxdists(linkage)
    # Extract the calculated Euclidean distance between clusters.
    distances_list = (sorted(list(distances_list)))
    distances_list.reverse()  # Sort the list in descending order.

    # Calculate the number of clusters for the HCA.
    number_of_clusters = numerical_number_of_clusters(distances_list)
    # Get the cutoff threshold for the dendrogram.
    cutoff = distances_list[number_of_clusters - 2]

    if ident == 'bio':
        # Define the name, directory and extension of the Euclidean distance graph.
        name_euc_dist_file = directory_fileoutput + \
            'Images/Distance_Images/Dist_HCA_biological.'+extension_type
        title_distances = 'Plot of Euclidean Distances from Biological HCA'
        print('Number of clusters:')
        # Print the number of clusters for the biological domain.
        print('Biological:', number_of_clusters)
        # Plot the graph of Euclidean distances.
        plot_euclidean_distance(
            distances_list, name_euc_dist_file, title_distances)
    elif ident == 'PhCh':
        # Define the name, directory and extension of the Euclidean distance graph.
        name_euc_dist_file = directory_fileoutput + \
            'Images/Distance_Images/Dist_HCA_physicochemical.'+extension_type
        title_distances = 'Plot of Euclidean Distances from Physicochemical HCA'
        # Print the number of clusters for the physicochemical domain.
        print('Physicochemical:', number_of_clusters)
        # Plot the graph of Euclidean distances.
        plot_euclidean_distance(
            distances_list, name_euc_dist_file, title_distances)
    else:
        # Define the name, directory and extension of the Euclidean distance graph.
        name_euc_dist_file = directory_fileoutput + \
            'Images/Distance_Images/Dist_HCA_structural.'+extension_type
        title_distances = 'Plot of Euclidean Distances from Structural HCA'
        # Print the number of clusters for the structural domain.
        print('Structural:', number_of_clusters)
        # Plot the graph of Euclidean distances.
        plot_euclidean_distance(
            distances_list, name_euc_dist_file, title_distances)

    return cutoff, number_of_clusters


# Plot the graph of Euclidean distances.
def plot_euclidean_distance(distances_list, name_euc_dist_file, title_distances):
    plt.plot(
        list(range(1, len(distances_list)+1)),
        distances_list, '-o',
        color='black',
        markersize=4, linewidth=1,
        markeredgewidth=1
    )  # Initial setting of graph values and parameters.
    plt.ylabel('Euclidean distance')  # Set the y-axis label.
    original_ticks = list((plt.xticks())[0])  # Get the x-axis ticks.
    increase_rate = original_ticks[2]  # Get the x-axis increment.
    # Create a list to define the new xticks that correspond
    # to double the points compared to the previous one.
    list_x = list(range(1, len(distances_list)+1, int(increase_rate/2)))
    # Define the new x-ticks and set the font size.
    plt.xticks(list_x, list_x, fontsize=8)
    plt.title(title_distances)  # Define the graph title.
    plt.savefig(name_euc_dist_file)  # Save figure.
    plt.close()  # Close plot.


def hca_clusters(file, labels, ident, directory_fileoutput, extension_type, linkage_method):
    # Get the threshold value (cutoff) and the number of clusters:
    linkage = sch.linkage(file, method=linkage_method, optimal_ordering=False)
    cutoff, number_of_clusters = optimal_threshold(
        file, ident, directory_fileoutput, extension_type, linkage)

    # Get the cluster labels:
    leaves_cluster = sch.fcluster(
        linkage, t=number_of_clusters, criterion='maxclust')
    cluster_dict = dict(zip(labels, list(leaves_cluster)))

    return leaves_cluster, cluster_dict, linkage, cutoff


def hca_plot(
        linkage,
        labels,
        leaves_cluster,
        cutoff,
        ident,
        directory_fileoutput,
        extension_type,
        dendrogram_xfont_size,
        test_or_val_molecules):

    # Adding arrows for molecules in test and validation sets:
    if len(test_or_val_molecules) == 1:
        index_of_test_molecules = test_or_val_molecules[0]
        labels_arrow = labels.copy()
        for i, label in enumerate(labels_arrow):
            if label in index_of_test_molecules:
                labels_arrow[i] = '→ ' + label
    else:
        index_of_test_molecules = test_or_val_molecules[0]
        index_of_val_molecules = test_or_val_molecules[1]
        labels_arrow = labels.copy()
        for i, label in enumerate(labels_arrow):
            if label in index_of_test_molecules:
                labels_arrow[i] = '→ ' + label
            elif label in index_of_val_molecules:
                labels_arrow[i] = '→→ ' + label

    # Working with colors:
    # Definition of which colors (colorblind friendly) will be used for each cluster (1-19).
    hex_rgb = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
               '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
               '#bcbd22', '#17becf', '#691605', '#b3f084',
               '#ffee00', '#00ff6e', '#0400ff', '#c4c4c0',
               '#ff00fb', '#290e0e', '#12290e', '#1f77b4',
               '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22',
               '#17becf', '#691605', '#b3f084', '#ffee00',
               '#00ff6e', '#0400ff', '#c4c4c0', '#ff00fb',
               '#290e0e', '#12290e']
    # Set hex_rgb as the dendrogram colors.
    sch.set_link_color_palette(hex_rgb)

    # Initial setting of HCA dendrograms:
    plt.figure(figsize=(20, 10), dpi=300)  # Set the image space.
    dendrogram = sch.dendrogram(linkage, color_threshold=cutoff,
                                above_threshold_color='k',
                                labels=labels_arrow, get_leaves=True)

    # Set graphs font size for both axis and adjust the plot in image space:
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=dendrogram_xfont_size)
    plt.ylabel('Euclidean distance', fontsize=12)
    plt.subplots_adjust(left=0.06, right=0.96)

    # Set new scale for y axis (which is 2x the number of points in the previous scale):
    b = [float(i) for i in plt.yticks()[0]]
    scale = b[1]-b[0]
    new_scale = scale/2
    new_arange = np.arange(min(b)-scale/20, max(b)+new_scale, new_scale)
    plt.yticks(new_arange)

    # Add ticks to x-axis:
    plt.tick_params(axis='x', which='both', bottom=True,
                    top=False, labelbottom=True)

    # Assign the wrong legend to each cluster:
    set__leaves_cluster = sorted(list(set(leaves_cluster)))
    legend_value = plt.legend(set__leaves_cluster)

    # Get the color for each molecule:
    color_list = list(dendrogram['leaves_color_list'])
    ivl_list = list(dendrogram['ivl'])
    dictionary_of_clusters = dict(zip(labels, leaves_cluster))
    for i in range(0, len(ivl_list)):
        if '→→ ' in ivl_list[i]:
            ivl_list[i] = ivl_list[i].replace('→→ ', '')
        elif '→ ' in ivl_list[i]:
            ivl_list[i] = ivl_list[i].replace('→ ', '')
    dictionary_of_colors = dict(zip(ivl_list, color_list))

    # Get the right cluster colors:
    list_cluster_color = []
    for i in dictionary_of_clusters.keys():
        list_cluster_color.append(
            [i, dictionary_of_clusters[i], dictionary_of_colors[i]])

    list_cluster_color = sorted(list_cluster_color, key=lambda x: x[1])
    lcc1 = {}

    if ('k' in color_list) or ('#000000' in color_list):
        lcc1['#000000'] = []

    for i in range(0, len(list_cluster_color)):
        if list_cluster_color[i][2] == 'k':
            lcc1['#000000'].append(list_cluster_color[i][1])
        else:
            lcc1[list_cluster_color[i][2]] = []
            lcc1[list_cluster_color[i][2]].append(list_cluster_color[i][1])

    right_legend = []
    for i in range(0, len(legend_value.get_lines())):
        right_legend.append(colrs.to_hex(
            legend_value.get_lines()[i].get_color()))

    # Assign the right legend to each cluster:
    set_leaves_cluster = []
    for i in right_legend:
        if (len(lcc1[i])) > 1:
            c_name = ''
            for z in range(0, len(lcc1[i])):
                a = lcc1[i][z]
                if z == len(lcc1[i])-1:
                    if a < 10:
                        c_name += 'cluster 0' + str(a) + ' (outliers)'
                    else:
                        c_name += 'cluster ' + str(a) + ' (outliers)'
                else:
                    if a < 10:
                        c_name += 'cluster 0' + str(a) + ', '
                    else:
                        c_name += 'cluster ' + str(a) + ', '
            set_leaves_cluster.append(c_name)
        else:
            a = int(list(set(lcc1[i]))[0])
            if a < 10:
                c_name = 'cluster 0' + str(a)
                set_leaves_cluster.append(c_name)
            else:
                c_name = 'cluster ' + str(a)
                set_leaves_cluster.append(c_name)
    legend_value = plt.legend(set_leaves_cluster)

    # Apply cluster color to x-labels:
    zip_label_color = list(zip(plt.gca().get_xticklabels(), color_list))
    for i, a in zip_label_color:
        i.set_color(a)

    # Define the title of the graph and saves the figure:
    if ident == 'bio':
        plt.title('HCA of biological activity')
        plt.savefig(directory_fileoutput +
                    'Images/Dendrogram_HCA_Images/HCA_biological.'+extension_type)
        plt.close()
    elif ident == 'PhCh':
        plt.title('HCA of physicochemical properties')
        plt.savefig(directory_fileoutput +
                    'Images/Dendrogram_HCA_Images/HCA_physicochemical.'+extension_type)
        plt.close()
    else:
        plt.title('HCA of AtomPairs fingerprint')
        plt.savefig(directory_fileoutput +
                    'Images/Dendrogram_HCA_Images/HCA_structural.'+extension_type)
        plt.close()


def kmodes_clusters(file, names, directory_fileoutput, extension_type):
    # Calculate the number of clusters for the Kmodes:
    n_clusters = categorical_number_of_clusters(
        file, directory_fileoutput, extension_type)

    # Run Kmodes with the calculated number of clusters:
    # Parameter setting.
    km_test = kmodes.KModes(n_clusters=n_clusters,
                            n_init=20, init='Cao', max_iter=300)
    # Run Kmodes for the dataset and return an array with cluster labels for each molecule.
    clusters = km_test.fit_predict(file)
    clusters_list = list(clusters)  # Transform np.array into a list.

    # Create a dictionary with the pair "molecule name": "General cluster ID":
    HCA_General_dict = {}
    for i, a in zip(names, clusters_list):
        HCA_General_dict[i] = int(a+1)
    return HCA_General_dict


def plot_distortion(
        list_clus_costs, directory_fileoutput,
        ident, extension_type, number_of_clusters,
        ran_curve):
    plt.plot(
        list(range(2, len(list_clus_costs)+2)),
        list_clus_costs, '-o',
        color='black',
        markersize=4, linewidth=1,
        markeredgewidth=1
    )  # Initial setting of graph values and parameters.
    if ident == 'bio':
        # Define the name, directory and extension of the Distortion graph.
        name_euc_dist_file = directory_fileoutput + \
            'Images/Distance_Images/Distortion_biological.'+extension_type
        title_distances = 'Plot of Distortion from Biological Space'
        print('Number of clusters:')
        # Print the number of clusters for the biological domain.
        print('Biological:', number_of_clusters)
        plt.ylabel('KMeans Distortion')  # Set the y-axis label.
    elif ident == 'PhCh':
        # Define the name, directory and extension of the Distortion graph.
        name_euc_dist_file = directory_fileoutput + \
            'Images/Distance_Images/Distortion_physicochemical.'+extension_type
        title_distances = 'Plot of Distortion from Physicochemical Space'
        # Print the number of clusters for the physicochemical domain.
        print('Physicochemical:', number_of_clusters)
        plt.ylabel('KMeans Distortion')  # Set the y-axis label.
    elif ident == 'FP':
        # Define the name, directory and extension of the Distortion graph.
        name_euc_dist_file = directory_fileoutput + \
            'Images/Distance_Images/Distortion_structural.'+extension_type
        title_distances = 'Plot of Distortion from Structural Space'
        # Print the number of clusters for the structural domain.
        print('Structural:', number_of_clusters)
        plt.ylabel('KMeans Distortion')  # Set the y-axis label.
    else:
        # Define the name, directory and extension of the Distortion graph.
        name_euc_dist_file = directory_fileoutput + \
            'Images/Distance_Images/Cost_general.'+extension_type
        title_distances = 'Cost Plot from General (KModes)'
        # Print the number of clusters for the structural domain.
        print('General:', number_of_clusters)
        plt.ylabel('KModes Cost')  # Set the y-axis label.
    plt.title(title_distances)  # Define the graph title.
    plt.xticks(ran_curve)  # Define the x-axis ticks.
    plt.savefig(name_euc_dist_file)  # Save figure.
    plt.close()  # Close plot.


def kmeans_number_of_clusters(file, ident, directory_fileoutput, extension_type):
    # Elbow method using the clustering cost:
    list_clus_costs = []
    maxcluster = 31
    ran_curve = [x for x in range(2, maxcluster, 2)]
    for n in range(2, maxcluster):  # Calculate the cost for each n cluster.
        km = KMeans(n_clusters=n, verbose=0,
                    init='k-means++', max_iter=300,
                    batch_size=256*int(mpi_count()),
                    compute_labels=True, random_state=2025)
        km.fit(X=file)
        # Costs = Clustering cost, defined as the sum distance of all points
        # to their respective cluster centroids.
        list_clus_costs.append(
            np.sum(np.min(cdist(file, km.cluster_centers_, 'euclidean'), axis=1)**2) / file.shape[0])
    # Define the x-axis coordinates of the first and last points on the line.
    x1, y1 = 2, list_clus_costs[0]
    # Define the y-axis coordinates of the first and last points on the line.
    x2, y2 = 20, list_clus_costs[18]

    distances = []
    for i in range(19):
        # Coordinate on the x-axis of the point to calculate the distance to the line.
        x0 = i+2
        # Coordinate on the y-axis of the point to calculate the distance to the line.
        y0 = list_clus_costs[i]
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = ((y2 - y1)**2 + (x2 - x1)**2)**(1/2)
        # Distance between a point and a line.
        distances.append(numerator/denominator)

    number_of_clusters = (distances.index(max(distances)) + 2)
    plot_distortion(
        list_clus_costs, directory_fileoutput,
        ident, extension_type, number_of_clusters,
        ran_curve)
    return number_of_clusters


def kmeans_clusters(file, names, ident, directory_fileoutput, extension_type):
    # Calculate the number of clusters for the Kmodes:
    n_clusters = kmeans_number_of_clusters(
        file, ident, directory_fileoutput, extension_type)

    # Run Kmodes with the calculated number of clusters:
    # Parameter setting.
    km = KMeans(n_clusters=n_clusters, verbose=0,
                init='k-means++', max_iter=300,
                batch_size=256*int(mpi_count()),
                compute_labels=True, random_state=2025)
    # Run Kmodes for the dataset and return an array with cluster labels for each molecule.
    clusters = km.fit_predict(file)
    clusters_list = list(clusters)  # Transform np.array into a list.

    # Create a dictionary with the pair "molecule name": "General cluster ID":
    kmeans_cluster_dict = {}
    for i, a in zip(names, clusters_list):
        kmeans_cluster_dict[i] = int(a+1)
    return kmeans_cluster_dict
