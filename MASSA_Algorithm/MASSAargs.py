import re
import os
import argparse
import numpy as np

def capture_args():
	## Capturing arguments from command line:
	parser = argparse.ArgumentParser(prog='MASSA_Algorithm', epilog='For more information visit: https://github.com/gcverissimo/MASSA_Algorithm', description='Help for MASSA Algorithm: Molecular data set sampling  —  Training-Test Separation', usage='%(prog)s [-h] [-i <input_file_name>] [-o <output_file_name] [-t <extension_for_images>] [-p <percentage_of_molecules_in_training_set>] [-b <number_of_biological_activities_in_input_file>], [-s <activity1>,<activity2>,<activityN>], [-n <number_of_principal_components>], [-v <svd_parameter>], [-d <dendrogram_font_size_on_Xaxis>], [-x <frequency_bar_plot_font_size_on_Xaxis>]')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i','--input', type=str, help='Input the file containing the dataset or its file path. Recommended to use an .sdf file to avoid errors.', required=True)
	required.add_argument('-o','--output', type=str, help='Enter the output file name or file path. Recommended to use an .sdf file to avoid errors. Image files will be saved to a folder within the same directory as the output file.', required=True)
	optional = parser.add_argument_group('optional arguments')
	optional.add_argument('-t','--image_type', type=str, default='png', help='Extension of the image files that will be generated. Suggested = png or svg. Default = png.')
	optional.add_argument('-d','--dendrogram_Xfont_size', type=int, default=5, help='Sets the font size on the x-axis of the dendrogram (molecule labels). Default = 5.')
	optional.add_argument('-x','--barplot_Xfont_size', type=int, default=11, help='Sets the font size on the x-axis of the bar plot (cluster labels). Default = 12.')
	optional.add_argument('-p','--percentage_of_training', type=float, choices=[round(i,2) for i in list(np.arange(0.00, 1.00, 0.01))], metavar='PERCENTAGE_OF_TRAINING', default=0.8, help='Percentage of molecules in training set. Must be a number from 0 to 1. Default = 0.8.')
	optional.add_argument('-b','--number_of_biological', type=int, default=1, help='Number of biological activities that will be used to separate the set into training and test. Default = 1')
	optional.add_argument('-s','--the_biological_activities', type=str, default=None, help='Enter a list with the names of biological activities separated by commas and no spaces. Example: pIC50,pMIC.')
	optional.add_argument('-n','--number_of_PCs', type=float, metavar='NUMBER_OF_PRINCIPAL_COMPONENTS', default=0.85, help='Defines the number of principal components to reduce the dimensionality of variables related to biological, physicochemical and structural domains. If the value is a decimal between 0 and 1, the number of principal components is what justifies for (<input number> * 100) percent of the variance. If the value is greater than 1, the number of PCs will be exactly the input integer, but PAY ATTENTION: 1st) If the number of PCs is an integer and equal to or greater than the number of physicochemical properties (7), the PCA step will be bypassed for this domain. 2nd) The same for the biological domain. 3rd) If the number of biological activities is less than 3, the PCA step will be bypassed for this domain. Default = 0.85.')
	optional.add_argument('-v','--svd_solver_for_PCA', type=str, default='auto', help='See the sklearn.decomposition.PCA topic on https://scikit-learn.org/ for more info. Default = auto.')
	args = parser.parse_args()

	## Definition of arguments:
	FileInput = args.input # Captures the input name and directory+name.
	undefinedFileOutput = args.output # Captures the output name or directory+name.
	extension_type = args.image_type # Captures the file extension for images.
	dendrogram_Xfont_size = args.dendrogram_Xfont_size # Captures the font size for the x-axis of the dendrogram.
	barplot_Xfont_size = args.barplot_Xfont_size # Captures the font size for the x-axis of the bar plot.
	training_percent = args.percentage_of_training # Captures the percentage for training set.
	test_percent = 1 - training_percent # Calculates the percentage of molecules in the test set.
	numberBioAct = args.number_of_biological # Captures the number of biological activities.
	inputnumberPCs = args.number_of_PCs # Captures the input number of PCs.
	svd_parameter = args.svd_solver_for_PCA # Captures the svd_solver parameter.

	if args.the_biological_activities != None:
		BioActAsArgs = args.the_biological_activities.split(',')
	else:
		BioActAsArgs = None

	if ('\\' in undefinedFileOutput) or ('/' in undefinedFileOutput):
		outputsplitted = re.split('\\\\|/', undefinedFileOutput)
		outputsplitted[-1] = ''
		directoryFileOutput = '\\'.join(outputsplitted)
		del outputsplitted
	else:
		directoryFileOutput = str(os.getcwd()) + '\\'
	
	if (inputnumberPCs > 1):
		nPCS = int(inputnumberPCs)
	elif (inputnumberPCs > 0) and (inputnumberPCs <= 1):
		nPCS = float(inputnumberPCs)
	else:
		print(' ERROR: The defined number of main components is invalid.')
		exit()

	return FileInput, undefinedFileOutput, directoryFileOutput, extension_type, dendrogram_Xfont_size, barplot_Xfont_size, training_percent, test_percent, numberBioAct, BioActAsArgs, nPCS, svd_parameter