from rdkit.Chem import AllChem

def read_molecules(file):
	mols = AllChem.SDMolSupplier(file)
	return mols

def error_check(file):
	for i in file:
		if i == None:
			print(' ERROR: Unable to read a molecule')
			exit()
		else:
			pass

def hydrogen_add(file):
	molsH = []
	for i in file:
		molsH.append(AllChem.AddHs(i, addCoords = True))
	return molsH

def get_sdf_property_names(file):
	sdf_property_names = []
	for i in file:
		for a in list(i.GetPropNames()):
			sdf_property_names.append(a)
	sdf_property_names = list(set(sdf_property_names))
	if 'Name' in sdf_property_names:
		sdf_property_names.remove('Name')
	return sdf_property_names