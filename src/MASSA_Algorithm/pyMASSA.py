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
import time


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


def py_massa(
    file_input=str(),
    working_directory=str(),
    job_name="MASSAoutput.sdf",
    extension_type="svg",
    dendrogram_xfont_size=5,
    barplot_xfont_size=9,
    training_percent=round(float(0.8), 3),
    bioactivities_as_args=None,
    npcs=0.85,
    svd_parameter="full",
    linkage_method="complete",
    flag_dendrogram=True,
    splitting_strategy="tt",
    drop_errors=True,
    large_library="auto"
):
    # Initializing from Python.
    # Print the program logo:
    MASSAlogos.python_logo()
    start_time = time.time()

    # Defining the working_directory alias:
    working_directory = f"{working_directory}/"
    directoryfile_output = working_directory

    # It creates the output directories:
    MASSAod.output_directory(directoryfile_output, flag_dendrogram)

    # Test nPCs:
    if npcs > 1:
        npcs = int(npcs)
    elif (npcs > 0) and (npcs <= 1):
        npcs = float(npcs)
    else:
        raise ValueError(
            f"ERROR: The defined number of main components is invalid. npcs = {str(npcs)}")

    # Autodetect if the biological activities are set.
    if isinstance(bioactivities_as_args, str):
        number_bioact = 1
        bioactivities_as_args = [bioactivities_as_args]
    elif isinstance(bioactivities_as_args, list):
        number_bioact = len(bioactivities_as_args)
    elif bioactivities_as_args == None:
        raise ValueError(
            "The biological activities are not specified. bioactivities_as_args == None")

    # Create log.txt file in directory:
    log_path = directoryfile_output+'log.txt'
    write_log = open(log_path, 'w')

    if isinstance(file_input, str):
        # Initial file management:
        # Read molecules.
        mols = MASSAopen_files.read_molecules(
            file_input, write_log, drop_errors)
    else:
        mols = file_input

    # Test if large_library is set:
    if isinstance(large_library, str):
        if large_library.lower() == 'true':
            large_library = True
        elif large_library.lower() == 'false':
            large_library = False
        elif large_library.lower() == 'auto':
            if int(len(mols)) > 9999:
                large_library = True
            else:
                large_library = False

    # If auto detection is enable to select the initial cluster algorithm.
    if large_library == 'auto':
        if int(len(mols)) > 9999:
            large_library = True
        else:
            large_library = False

    if large_library == True:
        # Deactivates the HCA dendrogram.
        flag_dendrogram = False

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

    if large_library == True:
        # First clustering (MiniBatch-KMeans):
        # It performs KMeans clustering for the biological domain.
        bio_hca = MASSAcluster.kmeans_clusters(
            bio_pca, names, 'bio', directoryfile_output, extension_type)
        # It performs KMeans clustering for the physicochemical domain.
        phch_hca = MASSAcluster.kmeans_clusters(
            phch_pca, names, 'PhCh', directoryfile_output, extension_type)
        # It performs KMeans clustering for the structural domain.
        fp_hca = MASSAcluster.kmeans_clusters(
            fp_pca, names, 'FP', directoryfile_output, extension_type)
    else:
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
    all_hca = MASSAcluster.kmodes_clusters(matrix_for_kmodes, names,
                                           directoryfile_output, extension_type)
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
        write_log
    )
    phch_distr = 'Physicochemical Distribution'
    MASSAsplit.log_of_distributions(
        phch_distr, phch_total, phch_training,
        phch_test, phch_val, splitting_strategy,
        write_log
    )
    fp_distr = 'Structural (FP) Distribution'
    MASSAsplit.log_of_distributions(
        fp_distr, fp_total, fp_training,
        fp_test, fp_val, splitting_strategy,
        write_log
    )
    all_distr = 'General Distribution'
    MASSAsplit.log_of_distributions(
        all_distr, all_total, all_training,
        all_test, all_val, splitting_strategy,
        write_log
    )

    # Output management:
    # It adds, for each molecule, the values of the calculated properties, the identifications of each cluster and which set the molecule belongs to.
    out_name = working_directory + job_name
    out_mols = MASSAmoloutput.output_file_mols(dataframe, out_name)

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

    end_time = time.time()
    elapsed_time = end_time - start_time
    write_log.write(f'Completed! Elapsed time = {elapsed_time} s.')
    write_log.close()
    print(f'Completed! Elapsed time = {elapsed_time} s.')
    return out_mols
