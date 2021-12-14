from rdkit.Chem import AllChem

def read_molecules(file, WriteLog):
	mols_supplier = AllChem.SDMolSupplier(file, sanitize=False) # Open the input file.
	mols = []
	for mol in mols_supplier:
		try:
			mols.append(AllChem.rdmolops.SanitizeMol(mol))
		except Exception as e:
			mols.append(None)
			mol_of_error = str(mol.GetProp('_Name'))
			error = str(e)
			precomposed = '[Error while reading molecule: \"' + error +  '\". Molecule name: \"' +  mol_of_error + '\".]'
			composed = precomposed+'\n'
			print(precomposed)
			WriteLog.write(composed)
	if None in mols:
		print('ERROR: Error while reading molecules. Check the log.txt file for more information. \nERROR: Closing...') # Print an error.
		exit()
	mols_supplier_sanitized = AllChem.SDMolSupplier(file, sanitize=True)
	return mols_supplier_sanitized

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