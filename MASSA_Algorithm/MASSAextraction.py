import pandas as pd
from rdkit.Chem import AllChem

def name_extraction(file):
	names = []
	dict_names = {}
	for i in file:
		try:
			names.append(i.GetProp('_Name')) # Get name of molecules to a list
			dict_names[i.GetProp('_Name')] = i # Get name of molecules to a dict
		except:
			names.append(None)
	if None in names: # If any molecule is unnamed it returns an error.
		print('ERROR: Couldn\'t get any name from a molecule. Adding names before running the script is essential. Check which molecule is in the ordered list below: \n', names)
		print('Closing...')
		exit()
	else:
		dataframe = pd.DataFrame({'molecules': pd.Series(dict_names)})
	return names, dataframe

def the_biological_handler(sdf_property_names, numberBioAct, BioActAsArgs):
	if BioActAsArgs == None: # If no biological activity name is passed as a command line argument.
		if (numberBioAct == len(sdf_property_names)): # If the number of biological activities equals the number of properties in the file (not counting the name), the only available properties will automatically be identified as the biological activities.
			name_biological_list = sdf_property_names.copy()
			return name_biological_list
		elif (numberBioAct > len(sdf_property_names)): # If the number of biological activities is greater than the number of properties in the file (not counting the name), returns an error.
			print('ERROR: The number of biological activities defined is greater than the number of properties in the .sdf file.')
			exit()
		else: # If no biological activity name was passed as a command line argument and the number of biological activities is less than the number of properties: Prompts user iteration to define the names of biological activities.
			name_biological_list = []
			print('Found the following properties in the input file: ' + str(sdf_property_names) + '.')
			for i in range(0,numberBioAct):
				a = (input('Enter the name of the ' + str(i+1) + 'st/nd/rd/th biological activity:  '))
				max_iterations = 0
				while (a not in sdf_property_names) and (max_iterations < 3): # If any property is typed incorrectly, the user will have 3 attempts to type correctly.
					print('ERROR: The biological property name entered does not match any of the properties found. Try again or exit the program from the command line using Ctrl+C.')
					a = (input('Enter the name of the ' + str(i+1) + 'st/nd/rd/th biological activity:  '))
					max_iterations += 1
				if a not in sdf_property_names: # If even after 3 attempts the name of the activity is still incorrect, the program will close.
					print('ERROR: The name of the biological activity was not typed correctly.')
					exit()
				name_biological_list.append(a)
			return name_biological_list

	else: # If biological activity names are passed as a command-line argument.
		if (numberBioAct == len(BioActAsArgs)): # The number of biological activities defined by the command line must be equal to the number of biological activities passed by the name.
			name_biological_list = BioActAsArgs.copy()
			for i in name_biological_list:
				if i not in sdf_property_names: # If the properties entered by the command line are typed incorrectly the program will close.
					print('ERROR: The name of the biological activity was not typed correctly.')
					exit()
			return name_biological_list
		else:
			print('ERROR: The number of biological activities passed as list \'-s\' is different from the number of biological activities passed as int \'b\'.')
			exit()

def list_activities(file, name_biological_activity):
	for bioACname in name_biological_activity: # For each biological activity:
		bio_dict = {}
		for i in file['molecules']: # For each molecule in the dataframe:
			try:
				bio_name = i.GetProp('_Name') # Capture the name of the molecule.
				bio_value = (float(i.GetProp(bioACname))) # Capture the value of the biological activity.
				bio_dict[bio_name] = bio_value # Add the "molecule name":"value" pair to the dictionary.
			except: # If any molecule does not have the biological activity in question it will add NoneType in place of the value.
				bio_name = i.GetProp('_Name')
				bio_value = None
				bio_dict[bio_name] = bio_value
		if(set(list(bio_dict.values()))) == {None}: # If no molecule has the biological property in question, it closes the program reporting the problem found.
			print('ERROR: The chosen property does not contain data or the name of the biological activity was not typed correctly.')
			exit()
		elif None in bio_dict.values(): # If any molecule has no biological activity it will return an error informing this and will print the biological activity in question and the dictionary "molecule name":"value".
			print('ERROR: Some molecule has no biological activity value.')
			print('Biological activity:', bioACname)
			print('Molecules:', bio_dict)
			exit()
		file[bioACname] = pd.Series(bio_dict)
	return file