from rdkit.Chem import AllChem

def output_mols(dataframe, FileOutput): # Add, for each molecule, the values of the calculated properties, the identifications of each cluster and which set the molecule belongs to.
	# Creation of dictionaries with the name of the molecule:
	dict_molecules = dict(dataframe['molecules']) # name:molecule
	dict_set = dict(dataframe['set']) # name:set
	dict_NumHAcceptors = dict(dataframe['NumHAcceptors']) # name:HBA
	dict_NumHDonors = dict(dataframe['NumHDonors']) # name:HBD
	dict_ExactMolWt = dict(dataframe['ExactMolWt']) # name:Molecular Weight
	dict_NumRotatableBonds = dict(dataframe['NumRotatableBonds']) # name:nRotBounds
	dict_FractionCSP3 = dict(dataframe['FractionCSP3']) # name:Fsp3
	dict_TPSA = dict(dataframe['TPSA']) # name:TPSA
	dict_LogP_WildmanCrippen = dict(dataframe['LogP_WildmanCrippen']) # name:logP
	dict_Cluster_General = dict(dataframe['Cluster_General']) # name:general_cluster
	dict_Cluster_Biological = dict(dataframe['Cluster_Biological']) # name:biological_cluster
	dict_Cluster_Physicochemical = dict(dataframe['Cluster_Physicochemical']) # name:physicochemical_cluster
	dict_Cluster_Structural = dict(dataframe['Cluster_Structural']) # name:strucutral_cluster

	# Add the properties to molecules from the name:molecule dictionary:
	for i in dict_molecules.keys():
		dict_molecules[i].SetProp('Set', str(dict_set[i]))
		dict_molecules[i].SetProp('HBA count', str(dict_NumHAcceptors[i]))
		dict_molecules[i].SetProp('HBD count', str(dict_NumHDonors[i]))
		dict_molecules[i].SetProp('ExactMolWt', str(dict_ExactMolWt[i]))
		dict_molecules[i].SetProp('NumRotatableBounds', str(dict_NumRotatableBonds[i]))
		dict_molecules[i].SetProp('FractionCsp3', str(dict_FractionCSP3[i]))
		dict_molecules[i].SetProp('TPSA', str(dict_TPSA[i]))
		dict_molecules[i].SetProp('LogP_WildmanCrippen', str(dict_LogP_WildmanCrippen[i]))
		dict_molecules[i].SetProp('General_Cluster', str(dict_Cluster_General[i]))
		dict_molecules[i].SetProp('Biological_Cluster', str(dict_Cluster_Biological[i]))
		dict_molecules[i].SetProp('Physicochemical_Cluster', str(dict_Cluster_Physicochemical[i]))
		dict_molecules[i].SetProp('Structural_Cluster', str(dict_Cluster_Structural[i]))

	# Create a dict_values with just the molecules (the name of the molecule is already part of molecule object):
	output_mols = dict_molecules.values()

	# Create output ".sdf" file:
	with AllChem.SDWriter(FileOutput) as w:
		for m in output_mols:
			w.write(m)