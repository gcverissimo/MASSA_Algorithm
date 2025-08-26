import re
import os
import argparse
import numpy as np


def capture_args():
    """Capture arguments from command-line

    Returns:
        args: Return the arguments to run the algorithm.
        For more information read the --help|-h.
    """
    # Capturing arguments from command line:
    parser = argparse.ArgumentParser(
        prog='MASSA_Algorithm',
        epilog='For more information visit: https://github.com/gcverissimo/MASSA_Algorithm',
        description=(
            'Help for MASSA Algorithm: Molecular data set sampling  '
            '—  Training-Test Separation'
        ),
        usage=(
            '%(prog)s [-h] [-i <input_file_name>] [-o <output_file_name] '
            '[-y <splitting_strategy>] [-p <percentage_of_molecules_in_training_set>] '
            '[-b <number_of_biological_activities_in_input_file>] [-s <activity1>,<activity2>,<activityN>]'
            '[-n <number_of_principal_components>] [-v <svd_parameter>] [-l <linkage_method] '
            '[-t <extension_for_images>] [-d <dendrogram_font_size_on_Xaxis>] '
            '[-x <frequency_bar_plot_font_size_on_Xaxis>] [-f <toggle_dendrogram>] '
            '[-e <drop_molecules_with_errors]'
        )
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', type=str, help=(
                          'Input the file containing the dataset or its file path. '
                          'Recommended to use an .sdf file to avoid errors.'
                          ), required=True)
    required.add_argument('-o', '--output', type=str, help=(
                          'Enter the output file name or file path. '
                          'Recommended to use an .sdf file to avoid errors. '
                          'Image files will be saved to a '
                          'folder within the same directory as the output file.'
                          ), required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-y', '--splitting_strategy', type=str, choices=['tt', 'tt-val'],
                          metavar='SPLITTING_STRATEGY', default='tt', help=(
                          'Defines the splitting strategy, either into training and test '
                          'sets or into training, test, and validation sets. Options = '
                          '\'tt\', \'tt-val\'. tt, means splitting into training and test sets. '
                          'tt-val, means splitting into training, test, and validation sets. '
                          'Default = \'tt\'.'
                          ))
    optional.add_argument('-p', '--percentage_of_training', type=float,
                          choices=[round(i, 2) for i in list(
                              np.arange(0.00, 1.00, 0.01))],
                          metavar='PERCENTAGE_OF_TRAINING',
                          default=0.8, help=(
                              'Percentage of molecules in training set. '
                              'Must be a number from 0 to 1. Default = 0.8.'
                          ))
    optional.add_argument('-b', '--number_of_biological', type=int, default=1,
                          help=(
                              'Number of biological activities that will '
                              'be used to separate the set into training and test. Default = 1'
                          ))
    optional.add_argument('-s', '--the_biological_activities', type=str, default=None,
                          help=(
                              'Enter a list with the names of biological activities '
                              'separated by commas and no spaces. Example: pIC50,pMIC.')
                          )
    optional.add_argument('-n', '--number_of_PCs', type=float,
                          metavar='NUMBER_OF_PRINCIPAL_COMPONENTS', default=0.85, help=(
                              'Defines the number of principal components to reduce the '
                              'dimensionality of variables related to biological, physicochemical '
                              'and structural domains. If the value is a decimal between 0 and 1, '
                              'the number of principal components is what justifies for '
                              '(<input number> * 100) percent of the variance. If the value '
                              'is greater than 1, the number of PCs will be exactly the input '
                              'integer, but PAY ATTENTION: 1st) If the number of PCs is an integer '
                              'and equal to or greater than the number of physicochemical '
                              'properties (7), the PCA step will be bypassed for this domain. '
                              '2nd) The same for the biological domain. 3rd) If the number of '
                              'biological activities is less than 3, the PCA step will be '
                              'bypassed for this domain. Default = 0.85.'
                          ))
    optional.add_argument('-v', '--svd_solver_for_PCA', type=str, default='full',
                          metavar='SVD_SOLVER_FOR_PCA',
                          help=('See the sklearn.decomposition.PCA topic on '
                                'https://scikit-learn.org/ '
                                'for more info. Default = full.'))
    optional.add_argument('-l', '--linkage_method', type=str,
                          choices=['complete', 'single', 'ward', 'average',
                                   'weighted', 'centroid', 'median'], metavar='LINKAGE_METHOD',
                          default='complete', help=(
                              'The linkage criterion to use. The algorithm '
                              'will merge the pairs of cluster that minimize '
                              'this criterion. Options = complete, single, ward, '
                              'average, weighted, centroid, median. Default = complete.'
                          ))
    optional.add_argument('-t', '--image_type', type=str, default='png',
                          help=(
                              'Extension of the image files that will be generated. '
                              'Suggested = png or svg. Default = png.'
                          ))
    optional.add_argument('-d', '--dendrogram_xfont_size', type=int, default=5,
                          help=(
                              'Sets the font size on the x-axis '
                              'of the dendrogram (molecule labels). Default = 5.'
                          ))
    optional.add_argument('-x', '--barplot_xfont_size', type=int, default=9,
                          help=(
                              'Sets the font size on the x-axis '
                              'of the bar plot (cluster labels). Default = 9.'
                          ))

    optional.add_argument('-f', '--dendrogram_plot', type=str, choices=['true', 'false'],
                          metavar='ENABLE_DENDROGRAM_PLOT', default='true',
                          help=('Defines whether or not dendrogram images will '
                                'be generated. Options = true (dendrogram will be generated), '
                                'false (dendrogram will not be generated). Default = true.'
                                )
                          )
    optional.add_argument('-e', '--drop_errors', type=str, choices=['true', 'false'],
                          metavar='DROP_MOLECULES_WITH_ERRORS', default='false', help=(
                          'Ignore chemistry errors, saving only molecules without any errors. '
                          'Default = \'false\'.'
                          ))
    args = parser.parse_args()

    # Definition of arguments:
    # It captures the input name and directory+name:
    file_input = args.input
    # It captures the output name or directory+name:
    undefined_file_output = args.output
    # It captures the file extension for images:
    extension_type = args.image_type
    # It captures the font size for the x-axis of the dendrogram:
    dendrogram_xfont_size = int(args.dendrogram_xfont_size)
    # It captures the font size for the x-axis of the bar plot:
    barplot_xfont_size = int(args.barplot_xfont_size)
    # It captures the percentage for training set:
    training_percent = round(float(args.percentage_of_training), 3)
    # It captures the number of biological activities:
    nbio_activities = int(args.number_of_biological)
    # It captures the input number of PCs:
    ninput_pcs = args.number_of_PCs
    # It captures the svd_solver parameter:
    svd_parameter = args.svd_solver_for_PCA
    # It captures the linkage method:
    linkage_method = args.linkage_method
    # It captures the flag for the dendrogram plot.
    flag_dendrogram = str(args.dendrogram_plot).lower()
    # It captures the selected splitting strategy.
    flag_strategy = str(args.splitting_strategy).lower()
    # It captures whether or not to ignore chemistry errors.
    drop = str(args.drop_errors).lower()

    if args.the_biological_activities != None:
        # It creates a list with the names of the biological activities to be searched for.
        bioactivities_as_args = args.the_biological_activities.split(',')
    else:
        bioactivities_as_args = None

    # If a directory path is passed along with the filename it copies only the directory path.
    if ('\\' in undefined_file_output) or ('/' in undefined_file_output):
        outputsplitted = re.split('\\\\|/', undefined_file_output)
        outputsplitted[-1] = ''
        dir_fileoutput = '/'.join(outputsplitted)
        del outputsplitted
    else:  # If only the filename is passed, the directory path is taken by os.getcwd().
        dir_fileoutput = str(os.getcwd()) + '/'

    if ninput_pcs > 1:
        npcs = int(ninput_pcs)
    elif (ninput_pcs > 0) and (ninput_pcs <= 1):
        npcs = float(ninput_pcs)
    else:
        print(' ERROR: The defined number of main components is invalid.')
        exit()

    if flag_dendrogram == 'true':
        flag_dendrogram = True
    elif flag_dendrogram == 'false':
        flag_dendrogram = False
    else:
        print(
            ' ERROR: ENABLE_DENDROGRAM_PLOT parameter only accepts \'true\' or \'false\'.')
        exit()

    if drop == 'true':
        drop = True
    elif drop == 'false':
        drop = False
    else:
        print(
            ' ERROR: DROP_MOLECULES_WITH_ERRORS parameter only accepts \'true\' or \'false\'.')
        exit()

    return (
        file_input,
        undefined_file_output,
        dir_fileoutput,
        extension_type,
        dendrogram_xfont_size,
        barplot_xfont_size,
        training_percent,
        nbio_activities,
        bioactivities_as_args,
        npcs,
        svd_parameter,
        linkage_method,
        flag_dendrogram,
        flag_strategy,
        drop
    )
