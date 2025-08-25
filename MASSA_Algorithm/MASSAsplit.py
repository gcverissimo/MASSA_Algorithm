import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt


def split_train_test_sets(file, training_percent, splitting_strategy):
    """This function splits the molecules into subsets (training/test/validation).

    Args:
        file (pd.DataFrame): The dataset of molecules
        training_percent (float): The percentage of molecules in the training set.
                                  It will also be used to calculate percentages of molecules
                                  in the validation and test sets.
        splitting_strategy (str): Training-test split or Training-test-validation split.

    Returns:
        dataset: The dataset with the labels for training-test-validation.
        test_or_val_molecules: An identification of which molecules are on the validation set.

    """
    if splitting_strategy == "tt":
        # Split the set of molecules into training set and test set,
        # preserving the cluster diversity of the total set according
        # to the stratification of the general clustering.

        test_percent = round(float(1 - training_percent), 3)
        training, test = train_test_split(
            file,
            test_size=test_percent,
            stratify=file['Cluster_General'])
        # Extraction of the name of the molecules included in the test set.
        test_molecules = [i for i in test.index]
        dict_training = {'set': 'training'}
        dict_test = {'set': 'test'}
        # Assign the training label to the molecules in the training set.
        training = training.assign(**dict_training)
        # Assign the test label to the molecules in the tests set.
        test = test.assign(**dict_test)
        # Group all molecules into a single dataframe.
        dataset = pd.concat([training, test])
        # Define the output sets for the algorithm.
        test_or_val_molecules = [test_molecules]
    else:
        # Split the set of molecules into training, test and validation sets,
        # preserving the cluster diversity of the total set according
        # to the stratification of the general clustering.
        test_percent = round(float(1 - training_percent), 3)/2
        val_percent = round(
            float(test_percent/(training_percent+test_percent)), 3)
        temp_training, test = train_test_split(
            file, test_size=test_percent, stratify=file['Cluster_General'])
        training, validation = train_test_split(
            temp_training, test_size=val_percent, stratify=temp_training['Cluster_General'])
        # Extraction of the name of the molecules included in the validation and test sets.
        validation_molecules = [i for i in validation.index]
        test_molecules = [i for i in test.index]
        dict_training = {'set': 'training'}
        dict_test = {'set': 'test'}
        dict_val = {'set': 'validation'}
        # Assign the training label to the molecules in the training set.
        training = training.assign(**dict_training)
        # Assign the test label to the molecules in the test set.
        test = test.assign(**dict_test)
        # Assign the validation label to the molecules in the validation set.
        validation = validation.assign(**dict_val)
        # Group all molecules into a single dataframe.
        dataset = pd.concat([training, test, validation])
        # Define the output sets for the algorithm.
        test_or_val_molecules = [test_molecules, validation_molecules]
    return dataset, test_or_val_molecules


def freq_clusters(
        file, dir_fileoutput,
        extension_type, cluster_domain,
        barplot_xfont_size, splitting_strategy):
    """_summary_

    Args:
        file (pd.DataFrame): Dataframe with the molecules.
        dir_fileoutput (str): Path to the output directory.
        extension_type (str): Extension type of output images.
        cluster_domain (str): Type of domain used (General, Bio, PhCh, or FP)
        barplot_xfont_size (int): Font size on X-axis.
        splitting_strategy (str): Algorithm of splitting.

    Returns:
        clusters_total (dict): Percentage of molecules in each cluster for the whole dataset.
        clusters_training (dict): Percentage of molecules in each cluster for the training subset.
        clusters_test (dict): Percentage of molecules in each cluster for the test subset.
        clusters_validation (dict): Percentage of molecules in each cluster for the validation subset.
    """
    if splitting_strategy == "tt":
        # Creating lists with the cluster IDs of each molecule:
        total_count = []
        training_count = []
        test_count = []

        for index, row in file.iterrows():
            if row['set'] == 'training':
                training_count.append(row[cluster_domain])
                total_count.append(row[cluster_domain])
            if row['set'] == 'test':
                test_count.append(row[cluster_domain])
                total_count.append(row[cluster_domain])
        min_value = min(total_count)
        max_value = max(total_count)

        # Frequency calculations for each cluster in each set saved in dictionaries:
        clusters_total = {}
        clusters_training = {}
        clusters_test = {}
        for i in range(min_value, (max_value+1)):
            if i <= 9:
                clusters_total[('cluster 0'+str(i))
                               ] = ((total_count.count(i))*100)/len(total_count)
                clusters_training[(
                    'cluster 0'+str(i))] = ((training_count.count(i))*100)/len(training_count)
                clusters_test[('cluster 0'+str(i))
                              ] = ((test_count.count(i))*100)/len(test_count)
            else:
                clusters_total[('cluster '+str(i))
                               ] = ((total_count.count(i))*100)/len(total_count)
                clusters_training[(
                    'cluster '+str(i))] = ((training_count.count(i))*100)/len(training_count)
                clusters_test[('cluster '+str(i))
                              ] = ((test_count.count(i))*100)/len(test_count)

        # Transform the dictionaries into a list of tuples and sort the list based
        # on the first element of the tuple (the cluster ID):
        list_total = sorted(list(clusters_total.items()), key=lambda x: x[0])
        list_training = sorted(
            list(clusters_training.items()), key=lambda x: x[0])
        list_test = sorted(list(clusters_test.items()), key=lambda x: x[0])

        # Separate the cluster IDs from their frequencies:
        label_total = [list_total[i][0] for i in range(0, len(list_total))]
        percentages_total = [list_total[i][1]
                             for i in range(0, len(list_total))]
        percentages_training = [list_training[i][1]
                                for i in range(0, len(list_training))]
        percentages_test = [list_test[i][1] for i in range(0, len(list_test))]

        # Create an "np.arange" for each set:
        total_ar = np.arange(len(list_total))
        training_ar = np.arange(len(list_training))
        test_ar = np.arange(len(list_test))

        # Bar Plot (matplotlib)
        # Set the font size in all fields of the chart.
        plt.rcParams.update({'font.size': 12})
        bar_width = 0.2
        fig, ax = plt.subplots(figsize=(10, 6))  # Set the figure size.
        group1 = ax.bar(total_ar-bar_width, percentages_total, bar_width, label='Total',
                        color='black', edgecolor='black')  # Bar plot for the total set.
        group2 = ax.bar(training_ar, percentages_training, bar_width, label='Training set',
                        color='dimgrey', edgecolor='black')  # Bar plot for the training set.
        group3 = ax.bar(test_ar+bar_width, percentages_test, bar_width, label='Test set',
                        color='lightgrey', edgecolor='black')  # Bar plot for the test set.
        ax.legend()  # Set the plot legend.
        ax.set_ylabel('Frequency (%)')  # Set the y-axis label.
        ax.set_ylim([0, 100])  # Set the y-axis range.
        ax.set_xticks(total_ar)  # Set the x-axis increment.
        # Set the font size and label of x-axis.
        ax.set_xticklabels(label_total, fontsize=barplot_xfont_size)

        if cluster_domain == 'Cluster_Biological':
            ax.set_title('Distribution of dataset: Biological Clustering')
            plt.savefig(dir_fileoutput +
                        '/Images/Distribution/Biological.'+extension_type)
        elif cluster_domain == 'Cluster_Structural':
            ax.set_title('Distribution of dataset: Structural Clustering')
            plt.savefig(dir_fileoutput +
                        '/Images/Distribution/Structural.'+extension_type)
        elif cluster_domain == 'Cluster_Physicochemical':
            ax.set_title('Distribution of dataset: Physicochemical Clustering')
            plt.savefig(dir_fileoutput +
                        '/Images/Distribution/Physicochemical.'+extension_type)
        elif cluster_domain == 'Cluster_General':
            ax.set_title('Distribution of dataset: General Clustering')
            plt.savefig(dir_fileoutput +
                        '/Images/Distribution/General.'+extension_type)
        plt.close()
        clusters_validation = None
    else:
        # Creating lists with the cluster IDs of each molecule:
        total_count = []
        training_count = []
        test_count = []
        validation_count = []

        for index, row in file.iterrows():
            if row['set'] == 'training':
                training_count.append(row[cluster_domain])
                total_count.append(row[cluster_domain])
            if row['set'] == 'test':
                test_count.append(row[cluster_domain])
                total_count.append(row[cluster_domain])
            if row['set'] == 'validation':
                validation_count.append(row[cluster_domain])
                total_count.append(row[cluster_domain])
        min_value = min(total_count)
        max_value = max(total_count)

        # Frequency calculations for each cluster in each set saved in dictionaries:
        clusters_total = {}
        clusters_training = {}
        clusters_test = {}
        clusters_validation = {}

        for i in range(min_value, (max_value+1)):
            if i <= 9:
                clusters_total[('cluster 0'+str(i))
                               ] = ((total_count.count(i))*100)/len(total_count)
                clusters_training[(
                    'cluster 0'+str(i))] = ((training_count.count(i))*100)/len(training_count)
                clusters_test[('cluster 0'+str(i))
                              ] = ((test_count.count(i))*100)/len(test_count)
                clusters_validation[('cluster 0'+str(i))
                                    ] = ((validation_count.count(i))*100)/len(validation_count)
            else:
                clusters_total[('cluster '+str(i))
                               ] = ((total_count.count(i))*100)/len(total_count)
                clusters_training[(
                    'cluster '+str(i))] = ((training_count.count(i))*100)/len(training_count)
                clusters_test[('cluster '+str(i))
                              ] = ((test_count.count(i))*100)/len(test_count)
                clusters_validation[('cluster '+str(i))
                                    ] = ((validation_count.count(i))*100)/len(validation_count)

        # Transform the dictionaries into a list of tuples and
        # sort the list based on the first element of the tuple (the cluster ID):
        list_total = sorted(list(clusters_total.items()), key=lambda x: x[0])
        list_training = sorted(
            list(clusters_training.items()), key=lambda x: x[0])
        list_test = sorted(list(clusters_test.items()), key=lambda x: x[0])
        list_validation = sorted(
            list(clusters_validation.items()), key=lambda x: x[0])

        # Separate the cluster IDs from their frequencies:
        label_total = [list_total[i][0] for i in range(0, len(list_total))]
        percentages_total = [list_total[i][1]
                             for i in range(0, len(list_total))]
        percentages_training = [list_training[i][1]
                                for i in range(0, len(list_training))]
        percentages_test = [list_test[i][1] for i in range(0, len(list_test))]
        percentages_validation = [list_validation[i][1]
                                  for i in range(0, len(list_validation))]

        # Create an "np.arange" for each set:
        total_ar = np.arange(len(list_total))
        training_ar = np.arange(len(list_training))
        test_ar = np.arange(len(list_test))
        validation_ar = np.arange(len(list_validation))

        # Bar Plot (matplotlib)
        # Set the font size in all fields of the chart.
        plt.rcParams.update({'font.size': 12})
        bar_width = 0.2
        fig, ax = plt.subplots(figsize=(10, 6))  # Set the figure size.
        group1 = ax.bar(total_ar-(bar_width*1.5), percentages_total, bar_width,
                        label='Total',
                        color='black', edgecolor='black')  # Bar plot for the total set.
        group2 = ax.bar(total_ar-(bar_width*1.5)+bar_width, percentages_training, bar_width,
                        label='Training set',
                        color='dimgrey', edgecolor='black')  # Bar plot for the training set.
        group3 = ax.bar(total_ar-(bar_width*1.5)+bar_width*2, percentages_test, bar_width,
                        label='Test set',
                        color='lightgrey', edgecolor='black')  # Bar plot for the test set.
        group4 = ax.bar(total_ar-(bar_width*1.5)+bar_width*3, percentages_validation, bar_width,
                        label='Validation set',
                        color='whitesmoke', edgecolor='black')  # Bar plot for the validation set.
        ax.legend()  # Set the plot legend.
        ax.set_ylabel('Frequency (%)')  # Set the y-axis label.
        ax.set_ylim([0, 100])  # Set the y-axis range.
        ax.set_xticks(total_ar)  # Set the x-axis increment.
        # Set the font size and label of x-axis.
        ax.set_xticklabels(label_total, fontsize=barplot_xfont_size)

        if cluster_domain == 'Cluster_Biological':
            ax.set_title('Distribution of dataset: Biological Clustering')
            plt.savefig(dir_fileoutput +
                        '/Images/Distribution/Biological.'+extension_type)
        elif cluster_domain == 'Cluster_Structural':
            ax.set_title('Distribution of dataset: Structural Clustering')
            plt.savefig(dir_fileoutput +
                        '/Images/Distribution/Structural.'+extension_type)
        elif cluster_domain == 'Cluster_Physicochemical':
            ax.set_title('Distribution of dataset: Physicochemical Clustering')
            plt.savefig(dir_fileoutput +
                        '/Images/Distribution/Physicochemical.'+extension_type)
        elif cluster_domain == 'Cluster_General':
            ax.set_title('Distribution of dataset: General Clustering')
            plt.savefig(dir_fileoutput +
                        '/Images/Distribution/General.'+extension_type)
        plt.close()
    return clusters_total, clusters_training, clusters_test, clusters_validation


def log_of_distributions(
    distribution_type,
    clusters_total,
    clusters_training,
    clusters_test,
    clusters_validation,
    splitting_strategy,
    writelog
):
    if splitting_strategy == "tt":
        writelog.write('['+distribution_type+']\nTotal: '+str(clusters_total)+'\n' +
                       'Training: '+str(clusters_training)+'\n'+'Test: '+str(clusters_test)+'\n\n')
    else:
        writelog.write('['+distribution_type+']\nTotal: '+str(clusters_total)+'\n' +
                       'Training: '+str(clusters_training)+'\n'+'Test: '+str(clusters_test)+'\n' +
                       'Validation: '+str(clusters_validation)+'\n''\n\n')
