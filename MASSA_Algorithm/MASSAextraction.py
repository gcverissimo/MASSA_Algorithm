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
	if None in names:
		print(' ERROR: Couldn\'t get any name from a molecule. Adding names before running the script is essential. Check which molecule is in the ordered list below: \n', names)
		print('Closing...')
		exit()
	else:
		dataframe = pd.DataFrame({'molecules': pd.Series(dict_names)})
	return names, dataframe

def the_biological_handler(sdf_property_names, numberBioAct, BioActAsArgs):
	if BioActAsArgs == None:
		if (numberBioAct == len(sdf_property_names)):
			name_biological_list = sdf_property_names.copy()
			return name_biological_list
		elif (numberBioAct > len(sdf_property_names)):
			print(' ERROR: The number of biological activities defined is greater than the number of properties in the .sdf file.')
			exit()
		else:
			name_biological_list = []
			print(' Found the following properties in the \".sdf\" file: ' + str(sdf_property_names) + '.')
			for i in range(0,numberBioAct):
				a = (input(' Enter the name of the ' + str(i+1) + 'st/nd/rd/th biological activity:  '))
				while a not in sdf_property_names:
					print(' ERROR: The biological property name entered does not match any of the properties found. Try again or exit the program from the command line using Ctrl+C.')
					a = (input(' Enter the name of the ' + str(i+1) + 'st/nd/rd/th biological activity:  '))
				name_biological_list.append(a)
			return name_biological_list
	else:
		if (numberBioAct == len(BioActAsArgs)):
			name_biological_list = BioActAsArgs.copy()
			return name_biological_list
		else:
			print(' ERROR: The number of biological activities passed as list \'-s\' is different from the number of biological activities passed as int \'b\'.')

def list_activities(file, name_biological_activity):
	for bioACname in name_biological_activity:
		bio_dict = {}
		for i in file['molecules']:
			try:
				bio_name = i.GetProp('_Name')
				bio_value = (float(i.GetProp(bioACname)))
				bio_dict[bio_name] = bio_value
			except:
				bio_name = i.GetProp('_Name')
				bio_value = None
				bio_dict[bio_name] = bio_value
		if (set(list(bio_dict.values()))) == {None}:
			print(' ERROR: The name of the biological activity was not typed correctly or the chosen property does not contain data.')
			exit()
		file[bioACname] = pd.Series(bio_dict)
	return file