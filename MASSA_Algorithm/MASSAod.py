import os

def output_directory(directoryFileOutput):
	try:
		os.makedirs(directoryFileOutput+'Images')
	except FileExistsError:
		pass
	try:
		os.makedirs(directoryFileOutput+'Images/Distance_HCA_Images')
	except FileExistsError:
		pass
	try:
		os.makedirs(directoryFileOutput+'Images/Dendrogram_HCA_Images')
	except FileExistsError:
		pass
	try:
		os.makedirs(directoryFileOutput+'Images/Distribution')
	except FileExistsError:
		pass