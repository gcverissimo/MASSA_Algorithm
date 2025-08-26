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


def returns_zero(total, test, validation):
    """
    Description:
        It evaluates whether the distribution is inadequate. That is,
        when the iterated cluster has a percentage greater than 0.5%
        in the complete dataset but less than 0.5% in the test set.

    Args:
        total (dict):
        Percentage of molecules in each cluster for the whole dataset.

        test (dict):
        Percentage of molecules in each cluster for the test subset.

        validation (dict):
        Percentage of molecules in each cluster for the validation subset.

    Returns:
        definer (boolean): Returns True if the distribution is
                           inadequate, and False if it is adequate.
    """
    definer = False
    for i in total.keys():
        if (total[i] > 0.5) and (test[i] <= 0.5):
            # Definer = True (Distribution was not done properly).
            definer = True
        elif validation != None:
            if (total[i] > 0.5) and (validation[i] <= 0.5):
                definer = True
    return definer


def main():
    """
    Description:
        Main subroutine, allows the program to run directly from
        the command line when installed via pip.
    """
    # Initializing from the command line.
    # Print the program logo:
    MASSAlogos.initial_print()

    # It captures command line arguments.
    (file_input,
     file_output,
     directoryfile_output,
     extension_type,
     dendrogram_xfont_size,
     barplot_xfont_size,
     training_percent,
     number_bioact,
     bioactivities_as_args,
     npcs,
     svd_parameter,
     linkage_method,
     flag_dendrogram,
     splitting_strategy,
     drop_errors) = MASSAargs.capture_args()
    print('Initializing, wait...\n')

    # Create log.txt file in directory:
    arqlog = directoryfile_output+'/log.txt'
    writelog = open(arqlog, 'w')

    # Initial file management:
    # It creates the output directories.
    MASSAod.output_directory(directoryfile_output)
    # Read molecules.
    mols = MASSAopen_files.read_molecules(file_input, writelog, drop_errors)
    # Extracting the property names from the ".sdf" input file.
    sdf_property_names = MASSAopen_files.get_sdf_property_names(mols)
    # Structure 3D management - It adds hydrogens keeping 3D coordenates.
    mols_h = MASSAopen_files.hydrogen_add(mols)

    # Extraction properties from ".sdf":
    # It extracts the names of the molecules and creates a name:molecule dictionary and a dataframe.
    names, dataframe = MASSAextraction.name_extraction(mols_h)
    # It defines a list of what biological activities are being extracted.
    biological_activity = MASSAextraction.the_biological_handler(
        sdf_property_names, number_bioact, bioactivities_as_args)
    # It adds the biological activities to the dataframe.
    dataframe = MASSAextraction.list_activities(dataframe, biological_activity)

    # Get fingeprint and other descriptors:
    # Get physicochemical descriptors.
    dataframe = MASSAdescriptors.physicochemical_descriptors(dataframe)
    # Get AtomPairs fingerprint.
    dataframe = MASSAdescriptors.atompairs_fingerprint(dataframe)

    # Normalizes physicochemical and biological properties.
    # Creates matrices for the three domains (struc, phys-chem, bio):
    bio_matrix, phch_matrix, fp_matrix = MASSApreparation.normalizer_or_matrix(
        dataframe, biological_activity)

    # PCA:
    # PCA for the biological domain.
    bio_pca = MASSApreparation.pca_maker(bio_matrix, npcs, svd_parameter)
    # PCA for the physicochemical domain.
    phch_pca = MASSApreparation.pca_maker(phch_matrix, npcs, svd_parameter)
    # PCA for the structural domain.
    fp_pca = MASSApreparation.pca_maker(fp_matrix, npcs, svd_parameter)

    # First clustering (HCA):
    # It performs HCA clustering without generating the dendrogram for the biological domain.
    leaves_cluster_bio, bio_hca, linkage_bio, cutoff_bio = MASSAcluster.hca_clusters(
        bio_pca, names, 'bio', directoryfile_output, extension_type, linkage_method)
    # It performs HCA clustering without generating the dendrogram for the physicochemical domain.
    leaves_cluster_phch, phch_hca, linkage_phch, cutoff_phch = MASSAcluster.hca_clusters(
        phch_pca, names, 'PhCh', directoryfile_output, extension_type, linkage_method)
    # It performs HCA clustering without generating the dendrogram for the structural domain.
    leaves_cluster_fp, fp_hca, linkage_fp, cutoff_fp = MASSAcluster.hca_clusters(
        fp_pca, names, 'FP', directoryfile_output, extension_type, linkage_method)

    # It adds the biological cluster identification to the spreadsheet.
    dataframe = MASSApreparation.organize_df_clusterization(
        dataframe, bio_hca, 'bio')
    # It adds the physicochemical cluster identification to the spreadsheet.
    dataframe = MASSApreparation.organize_df_clusterization(
        dataframe, phch_hca, 'PhCh')
    # It adds the structural cluster identification to the spreadsheet.
    dataframe = MASSApreparation.organize_df_clusterization(
        dataframe, fp_hca, 'FP')

    # Second clustering (Kmodes):
    # It creates a matrix with cluster identifications for each of the three domains for Kmodes.
    matrix_for_kmodes = MASSApreparation.organize_for_kmodes(dataframe)
    # It performs Kmodes clustering for the general domain.
    all_hca = MASSAcluster.kmodes_clusters(matrix_for_kmodes, names)
    # It adds the general cluster identification to the spreadsheet.
    dataframe = MASSApreparation.organize_df_clusterization(
        dataframe, all_hca, 'all')

    # Split into training, test:
    dataframe, test_or_val_molecules = MASSAsplit.split_train_test_sets(
        dataframe, training_percent, splitting_strategy)

    # Bar plot of frequencies
    # (Calculates the percentages of molecules in each cluster for each
    #  dataset and generates a bar graph for each domain):

    # Biological Bar Plot
    bio_total, bio_training, bio_test, bio_val = MASSAsplit.freq_clusters(
        dataframe,
        directoryfile_output,
        extension_type,
        'Cluster_Biological',
        barplot_xfont_size,
        splitting_strategy)
    # Physicochemical Bar Plot
    phch_total, phch_training, phch_test, phch_val = MASSAsplit.freq_clusters(
        dataframe,
        directoryfile_output,
        extension_type,
        'Cluster_Physicochemical',
        barplot_xfont_size,
        splitting_strategy)
    # Structural Bar Plot
    fp_total, fp_training, fp_test, fp_val = MASSAsplit.freq_clusters(
        dataframe,
        directoryfile_output,
        extension_type,
        'Cluster_Structural',
        barplot_xfont_size,
        splitting_strategy)
    # General Bar Plot
    all_total, all_training, all_test, all_val = MASSAsplit.freq_clusters(
        dataframe,
        directoryfile_output,
        extension_type,
        'Cluster_General',
        barplot_xfont_size,
        splitting_strategy)

    # Verifying percentages:
    # For biological domain:
    bio_ok = returns_zero(bio_total, bio_test, bio_val)
    # For physicochemical domain:
    phch_ok = returns_zero(phch_total, phch_test, phch_val)
    # For structural domain:
    fp_ok = returns_zero(fp_total, fp_test, fp_val)
    ok = [bio_ok, phch_ok, fp_ok]
    max_iters = 0

    # Redo the distribution in case of errors (up to 10 times):
    while (True in ok) and (max_iters < 10):
        # Split into training, test:
        dataframe, test_or_val_molecules = MASSAsplit.split_train_test_sets(
            dataframe, training_percent, splitting_strategy)

        # Bar plot of frequencies (Calculates the percentages of molecules in each
        # cluster for each dataset and generates a bar graph for each domain):

        # Biological Bar Plot:
        bio_total, bio_training, bio_test, bio_val = MASSAsplit.freq_clusters(
            dataframe, directoryfile_output, extension_type,
            'Cluster_Biological', barplot_xfont_size, splitting_strategy
        )
        # Physicochemical Bar Plot
        phch_total, phch_training, phch_test, phch_val = MASSAsplit.freq_clusters(
            dataframe, directoryfile_output, extension_type,
            'Cluster_Physicochemical', barplot_xfont_size, splitting_strategy
        )
        # Structural Bar Plot
        fp_total, fp_training, fp_test, fp_val = MASSAsplit.freq_clusters(
            dataframe, directoryfile_output, extension_type,
            'Cluster_Structural', barplot_xfont_size, splitting_strategy
        )
        # General Bar Plot
        all_total, all_training, all_test, all_val = MASSAsplit.freq_clusters(
            dataframe, directoryfile_output, extension_type,
            'Cluster_General', barplot_xfont_size, splitting_strategy
        )

        # Verifying percentages:

        # For biological domain:
        bio_ok = returns_zero(bio_total, bio_test, bio_val)
        # For physicochemical domain:
        phch_ok = returns_zero(phch_total, phch_test, phch_val)
        # For structural domain:
        fp_ok = returns_zero(fp_total, fp_test, fp_val)
        ok = [bio_ok, phch_ok, fp_ok]
        max_iters += 1

    # Write distribution information to log file:
    bio_distr = 'Biological Distribution'
    MASSAsplit.log_of_distributions(
        bio_distr, bio_total, bio_training,
        bio_test, bio_val, splitting_strategy,
        writelog
    )
    phch_distr = 'Physicochemical Distribution'
    MASSAsplit.log_of_distributions(
        phch_distr, phch_total, phch_training,
        phch_test, phch_val, splitting_strategy,
        writelog
    )
    fp_distr = 'Structural (FP) Distribution'
    MASSAsplit.log_of_distributions(
        fp_distr, fp_total, fp_training,
        fp_test, fp_val, splitting_strategy,
        writelog
    )
    all_distr = 'General Distribution'
    MASSAsplit.log_of_distributions(
        all_distr, all_total, all_training,
        all_test, all_val, splitting_strategy,
        writelog
    )
    writelog.close()

    # Plot HCAs:
    if flag_dendrogram == True:
        print('\nGenerating dendrogram images. Please wait...')
        MASSAcluster.hca_plot(
            linkage_bio, names, leaves_cluster_bio,
            cutoff_bio, 'bio', directoryfile_output,
            extension_type, dendrogram_xfont_size, test_or_val_molecules
        )  # Bio_Plot: Plot the HCA dendrogram
        MASSAcluster.hca_plot(
            linkage_phch, names, leaves_cluster_phch,
            cutoff_phch, 'PhCh', directoryfile_output,
            extension_type, dendrogram_xfont_size, test_or_val_molecules
        )  # phch_Plot: Plot the HCA dendrogram
        MASSAcluster.hca_plot(
            linkage_fp, names, leaves_cluster_fp,
            cutoff_fp, 'FP', directoryfile_output,
            extension_type, dendrogram_xfont_size, test_or_val_molecules
        )  # fp_Plot: Plot the HCA dendrogram

    # Output management:
    # It adds, for each molecule, the values of the calculated properties, the identifications of each cluster and which set the molecule belongs to.
    MASSAmoloutput.output_mols(dataframe, file_output)
    print('Completed')
