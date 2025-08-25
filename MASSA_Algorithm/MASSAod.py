import os


def output_directory(directory_file_output):
    """
    Creates the output directories.
    """
    try:
        os.makedirs(directory_file_output+'Images')
    except FileExistsError:
        pass
    try:
        os.makedirs(directory_file_output+'Images/Distance_HCA_Images')
    except FileExistsError:
        pass
    try:
        os.makedirs(directory_file_output+'Images/Dendrogram_HCA_Images')
    except FileExistsError:
        pass
    try:
        os.makedirs(directory_file_output+'Images/Distribution')
    except FileExistsError:
        pass
