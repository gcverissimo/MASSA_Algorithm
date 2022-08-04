from MASSA_Algorithm import MASSAlogos
from MASSA_Algorithm import MASSAargs
from MASSA_Algorithm import MASSAod
from MASSA_Algorithm import MASSAopen_files
from MASSA_Algorithm import MASSAextraction
from MASSA_Algorithm import MASSAdescriptors
from MASSA_Algorithm import MASSApreparation
from MASSA_Algorithm import MASSAcluster
from MASSA_Algorithm import MASSAsplit
from MASSA_Algorithm import MASSAmoloutput

def returns_zero(total, test):
	# It evaluates if the distribution is not adequate = the iterated cluster has a percentage greater than 0.5% in the complete data set, but a percentage lower than 0.5% in the test set.
	definer = False
	for i in total.keys():
		if (total[i] > 0.5) and (test[i] <= 0.5):
			definer = True # Definer = True (Distribution was not done properly).
	return definer

def main(): # Main subroutine, allows the program to run directly from the command line after installed via pip.
	## Initializing from the command line:
	MASSAlogos.initial_print() # Print the program logo.
	FileInput, FileOutput, directoryFileOutput, extension_type, dendrogram_Xfont_size, barplot_Xfont_size, training_percent, test_percent, numberBioAct, BioActAsArgs, nPCS, svd_parameter, linkage_method, flag_dendrogram = MASSAargs.capture_args() # It captures command line arguments.
	print('Initializing, wait...\n')


	## Create log.txt file in directory:
	ArqLog = directoryFileOutput+'/log.txt'
	WriteLog = open(ArqLog, 'w')
	
	
	## Initial file management:
	MASSAod.output_directory(directoryFileOutput) # It creates the output directories.
	mols = MASSAopen_files.read_molecules(FileInput, WriteLog) # Read molecules.
	sdf_property_names = MASSAopen_files.get_sdf_property_names(mols) # Extracting the property names from the ".sdf" input file.
	molsH = MASSAopen_files.hydrogen_add(mols) # Structure 3D management - It adds hydrogens keeping 3D coordenates.


	## Extraction properties from ".sdf":
	names, dataframe = MASSAextraction.name_extraction(molsH) # It extracts the names of the molecules and creates a name:molecule dictionary and a dataframe.
	biological_activity = MASSAextraction.the_biological_handler(sdf_property_names, numberBioAct, BioActAsArgs) # It defines a list of what biological activities are being extracted.
	dataframe = MASSAextraction.list_activities(dataframe, biological_activity) # It adds the biological activities to the dataframe.


	## Get fingeprint and other descriptors:
	dataframe = MASSAdescriptors.physicochemical_descriptors(dataframe) # Get physicochemical descriptors.
	dataframe = MASSAdescriptors.atompairs_fingerprint(dataframe) # Get AtomPairs fingerprint.


	## Normalizes physicochemical and biological properties and creates matrices for the three domains:
	bio_matrix, PhCh_matrix, FP_matrix = MASSApreparation.normalizer_or_matrix(dataframe, biological_activity)


	## PCA:
	bio_PCA = MASSApreparation.pca_maker(bio_matrix, nPCS, svd_parameter) # PCA for the biological domain.
	PhCh_PCA = MASSApreparation.pca_maker(PhCh_matrix, nPCS, svd_parameter) # PCA for the physicochemical domain.
	FP_PCA = MASSApreparation.pca_maker(FP_matrix, nPCS, svd_parameter) # PCA for the structural domain.


	## First clustering (HCA):
	leaves_cluster_bio, bioHCA, linkage_bio, CutOff_bio = MASSAcluster.hca_clusters(bio_PCA, names, 'bio', directoryFileOutput, extension_type, linkage_method) # It performs HCA clustering without generating the dendrogram for the biological domain.
	leaves_cluster_phch, phchHCA, linkage_phch, CutOff_phch = MASSAcluster.hca_clusters(PhCh_PCA, names, 'PhCh', directoryFileOutput, extension_type, linkage_method) # It performs HCA clustering without generating the dendrogram for the physicochemical domain.
	leaves_cluster_fp, fpHCA, linkage_fp, CutOff_fp = MASSAcluster.hca_clusters(FP_PCA, names, 'FP', directoryFileOutput, extension_type, linkage_method) # It performs HCA clustering without generating the dendrogram for the structural domain.

	dataframe = MASSApreparation.organize_df_clusterization(dataframe, bioHCA, 'bio') # It adds the biological cluster identification to the spreadsheet.
	dataframe = MASSApreparation.organize_df_clusterization(dataframe, phchHCA, 'PhCh') # It adds the physicochemical cluster identification to the spreadsheet.
	dataframe = MASSApreparation.organize_df_clusterization(dataframe, fpHCA, 'FP') # It adds the structural cluster identification to the spreadsheet.


	## Second clustering (Kmodes):
	matrix_for_kmodes = MASSApreparation.organize_for_kmodes(dataframe) # It creates a matrix with cluster identifications for each of the three domains, in order to prepare for Kmodes.
	allHCA = MASSAcluster.kmodes_clusters(matrix_for_kmodes, names) # It performs Kmodes clustering for the general domain.
	dataframe = MASSApreparation.organize_df_clusterization(dataframe, allHCA, 'all') # It adds the general cluster identification to the spreadsheet.


	## Split into training, test:
	dataframe, test_molecules = MASSAsplit.split_train_test_sets(dataframe, training_percent, test_percent)


	## Bar plot of frequencies (Calculates the percentages of molecules in each cluster for each dataset and generates a bar graph for each domain):
	bio_total, bio_training, bio_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Biological', barplot_Xfont_size) # Biological Bar Plot
	PhCh_total, PhCh_training, PhCh_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Physicochemical', barplot_Xfont_size) # Physicochemical Bar Plot
	FP_total, FP_training, FP_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Structural', barplot_Xfont_size) # Structural Bar Plot
	all_total, all_training, all_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_General', barplot_Xfont_size) # General Bar Plot


	## Verifying percentages:
	bio_ok = returns_zero(bio_total, bio_test) # For biological domain.
	PhCh_ok = returns_zero(PhCh_total, PhCh_test) # For physicochemical domain.
	FP_ok = returns_zero(FP_total, FP_test) # For structural domain.
	ok = [bio_ok, PhCh_ok, FP_ok]
	max_iters = 0
		
	# Redo the distribution in case of errors (up to 10 times):
	while (True in ok) and (max_iters < 10):
		## Split into training, test:
		dataframe, test_molecules = MASSAsplit.split_train_test_sets(dataframe, training_percent, test_percent)

		## Bar plot of frequencies (Calculates the percentages of molecules in each cluster for each dataset and generates a bar graph for each domain):
		bio_total, bio_training, bio_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Biological', barplot_Xfont_size) # Biological Bar Plot
		PhCh_total, PhCh_training, PhCh_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Physicochemical', barplot_Xfont_size) # Physicochemical Bar Plot
		FP_total, FP_training, FP_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Structural', barplot_Xfont_size) # Structural Bar Plot
		all_total, all_training, all_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_General', barplot_Xfont_size) # General Bar Plot

		## Verifying percentages:
		bio_ok = returns_zero(bio_total, bio_test) # For biological domain.
		PhCh_ok = returns_zero(PhCh_total, PhCh_test) # For physicochemical domain.
		FP_ok = returns_zero(FP_total, FP_test) # For structural domain.
		ok = [bio_ok, PhCh_ok, FP_ok]
		max_iters += 1


	## Write distribution information to log file:
	bio_distr = 'Biological Distribution'
	MASSAsplit.log_of_distributions(bio_distr, bio_total, bio_training, bio_test, WriteLog)
	PhCh_distr = 'Physicochemical Distribution'
	MASSAsplit.log_of_distributions(PhCh_distr, PhCh_total, PhCh_training, PhCh_test, WriteLog)
	FP_distr = 'Structural (FP) Distribution'
	MASSAsplit.log_of_distributions(FP_distr, FP_total, FP_training, FP_test, WriteLog)
	all_distr = 'General Distribution'
	MASSAsplit.log_of_distributions(all_distr, all_total, all_training, all_test, WriteLog)
	WriteLog.close()


	## Plot HCAs:
	if flag_dendrogram == True:
		print('\nGenerating dendrogram images. Please wait...')
		MASSAcluster.hca_plot(linkage_bio, names, leaves_cluster_bio, CutOff_bio, 'bio', directoryFileOutput, extension_type, dendrogram_Xfont_size, test_molecules) #Bio_Plot: Plot the HCA dendrogram
		MASSAcluster.hca_plot(linkage_phch, names, leaves_cluster_phch, CutOff_phch, 'PhCh', directoryFileOutput, extension_type, dendrogram_Xfont_size, test_molecules) #PhCh_Plot: Plot the HCA dendrogram
		MASSAcluster.hca_plot(linkage_fp, names, leaves_cluster_fp, CutOff_fp, 'FP', directoryFileOutput, extension_type, dendrogram_Xfont_size, test_molecules) #FP_Plot: Plot the HCA dendrogram


	## Output management:
	MASSAmoloutput.output_mols(dataframe, FileOutput) # It adds, for each molecule, the values of the calculated properties, the identifications of each cluster and which set the molecule belongs to.
	print('Completed')