from rdkit.Chem import AllChem

def read_molecules(file):
	mols = AllChem.SDMolSupplier(file) # Open the input file.
	return mols

def error_check(file):
	for i in file:
		if i == None: # If any molecule is not able to be read: 
			print('ERROR: Unable to read a molecule') # Print an error.
			exit()
		else:
			pass

def hydrogen_add(file):
	molsH = []
	for i in file:
		molsH.append(AllChem.AddHs(i, addCoords = True)) # Add hydrogens keeping 3D coordenates.
	return molsH

def get_sdf_property_names(file): # Extracting the property names from the ".sdf" input file.
	sdf_property_names = []
	for i in file:
		for a in list(i.GetPropNames()):
			sdf_property_names.append(a)
	sdf_property_names = list(set(sdf_property_names))
	if 'Name' in sdf_property_names: # Drop the molecule name from the property list.
		sdf_property_names.remove('Name')
	return sdf_property_names