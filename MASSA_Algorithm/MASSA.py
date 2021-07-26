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
	total_values = list(total.values())
	test_values = list(test.values())
	definer = False
	for i in range(0, len(total_values)):
		if (total_values[i] > 0.5) and (test_values[i] <= 0.5):
			definer = True
	return definer

def main():
	## Initializing from the command line:
	MASSAlogos.initial_print()
	FileInput, FileOutput, directoryFileOutput, extension_type, dendrogram_Xfont_size, barplot_Xfont_size, training_percent, test_percent, numberBioAct, BioActAsArgs, nPCS, svd_parameter = MASSAargs.capture_args()
	MASSAlogos.print_initializing()


	## Initial file management:
	MASSAod.output_directory(directoryFileOutput) # Create the output directories.
	mols = MASSAopen_files.read_molecules(FileInput) # Read molecules.
	MASSAopen_files.error_check(mols) # Alert of RDKit fail to read.
	sdf_property_names = MASSAopen_files.get_sdf_property_names(mols) # Capture the properties of .sdf.
	molsH = MASSAopen_files.hydrogen_add(mols) # Structure 3D management - Adds hydrogens keeping 3D coordenates.


	## Extraction properties from ".sdf":
	names, dataframe = MASSAextraction.name_extraction(molsH)
	biological_activity = MASSAextraction.the_biological_handler(sdf_property_names, numberBioAct, BioActAsArgs)
	dataframe = MASSAextraction.list_activities(dataframe, biological_activity)


	## Get fingeprint and other descriptors:
	dataframe = MASSAdescriptors.physicochemical_descriptors(dataframe) # Get physicochemical descriptors.
	dataframe = MASSAdescriptors.AtomPairs_fingerprint(dataframe) # Get AtomPairs fingerprint.


	## Normalizes physicochemical and biological properties and creates matrices of the three domains:
	bio_matrix, PhCh_matrix, FP_matrix = MASSApreparation.Normalizer_or_matrix(dataframe, biological_activity)


	## PCA:
	bio_PCA = MASSApreparation.PCA_maker(bio_matrix, nPCS, svd_parameter)
	PhCh_PCA = MASSApreparation.PCA_maker(PhCh_matrix, nPCS, svd_parameter)
	FP_PCA = MASSApreparation.PCA_maker(FP_matrix, nPCS, svd_parameter)


	## First clustering (HCA):
	leaves_cluster_bio, bioHCA, linkage_bio, CutOff_bio = MASSAcluster.HCA_clusters(bio_PCA, names, 'bio', directoryFileOutput, extension_type)
	leaves_cluster_phch, phchHCA, linkage_phch, CutOff_phch = MASSAcluster.HCA_clusters(PhCh_PCA, names, 'PhCh', directoryFileOutput, extension_type)
	leaves_cluster_fp, fpHCA, linkage_fp, CutOff_fp = MASSAcluster.HCA_clusters(FP_PCA, names, 'FP', directoryFileOutput, extension_type)

	dataframe = MASSApreparation.organize_HCA(dataframe, bioHCA, 'bio')
	dataframe = MASSApreparation.organize_HCA(dataframe, phchHCA, 'PhCh')
	dataframe = MASSApreparation.organize_HCA(dataframe, fpHCA, 'FP')


	## Second clustering (Kmodes):
	matrix_for_kmodes = MASSApreparation.organize_for_kmodes(dataframe)
	allHCA = MASSAcluster.KModes(matrix_for_kmodes, names)
	dataframe = MASSApreparation.organize_HCA(dataframe, allHCA, 'all')


	## Split into training, test:
	dataframe, index_of_test_molecules = MASSAsplit.split_train_test_sets(dataframe, training_percent, test_percent)


	## Bar plot of frequencies:
	all_total, all_training, all_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_General', barplot_Xfont_size)
	bio_total, bio_training, bio_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Biological', barplot_Xfont_size)
	PhCh_total, PhCh_training, PhCh_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Physicochemical', barplot_Xfont_size)
	FP_total, FP_training, FP_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Structural', barplot_Xfont_size)


	## Verifying percentages:
	bio_ok = returns_zero(bio_total, bio_test)
	PhCh_ok = returns_zero(PhCh_total, PhCh_test)
	FP_ok = returns_zero(FP_total, FP_test)
	ok = [bio_ok, PhCh_ok, FP_ok]
	max_iters = 0
	
	while (True in ok) and (max_iters < 20):
		## Split into training, test:
		dataframe = MASSAsplit.split_train_test_sets(dataframe, training_percent, test_percent)

		## Bar plot of frequencies:
		all_total, all_training, all_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_General', barplot_Xfont_size)
		bio_total, bio_training, bio_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Biological', barplot_Xfont_size)
		PhCh_total, PhCh_training, PhCh_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Physicochemical', barplot_Xfont_size)
		FP_total, FP_training, FP_test = MASSAsplit.freq_clusters(dataframe, directoryFileOutput, extension_type, 'Cluster_Structural', barplot_Xfont_size)

		bio_ok = returns_zero(bio_total, bio_test)
		PhCh_ok = returns_zero(PhCh_total, PhCh_test)
		FP_ok = returns_zero(FP_total, FP_test)
		ok = [bio_ok, PhCh_ok, FP_ok]
		max_iters += 1


	## Plot HCAs:
	MASSAcluster.HCA_plot(linkage_bio, names, leaves_cluster_bio, CutOff_bio, 'bio', directoryFileOutput, extension_type, dendrogram_Xfont_size, index_of_test_molecules) #BioPlot
	MASSAcluster.HCA_plot(linkage_phch, names, leaves_cluster_phch, CutOff_phch, 'PhCh', directoryFileOutput, extension_type, dendrogram_Xfont_size, index_of_test_molecules) #PhchPlot
	MASSAcluster.HCA_plot(linkage_fp, names, leaves_cluster_fp, CutOff_fp, 'FP', directoryFileOutput, extension_type, dendrogram_Xfont_size, index_of_test_molecules) #FpPlot


	## Output management:
	MASSAmoloutput.outputMOLS(dataframe, FileOutput)