import os


def output_directory(directory_file_output, flag_dendrogram):
    """
    Creates the output directories.
    """
    try:
        os.makedirs(directory_file_output+'Images')
    except FileExistsError:
        pass
    try:
        os.makedirs(directory_file_output+'Images/Distance_Images')
    except FileExistsError:
        pass
    if flag_dendrogram == True:
        try:
            os.makedirs(directory_file_output+'Images/Dendrogram_HCA_Images')
        except FileExistsError:
            pass
    try:
        os.makedirs(directory_file_output+'Images/Distribution')
    except FileExistsError:
        pass
