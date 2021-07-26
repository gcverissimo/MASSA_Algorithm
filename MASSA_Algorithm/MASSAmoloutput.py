from rdkit.Chem import AllChem

def outputMOLS(dataframe, FileOutput):
	dict_molecules = dict(dataframe['molecules'])
	dict_set = dict(dataframe['set'])
	dict_NumHAcceptors = dict(dataframe['NumHAcceptors'])
	dict_NumHDonors = dict(dataframe['NumHDonors'])
	dict_ExactMolWt = dict(dataframe['ExactMolWt'])
	dict_NumRotatableBonds = dict(dataframe['NumRotatableBonds'])
	dict_FractionCSP3 = dict(dataframe['FractionCSP3'])
	dict_TPSA = dict(dataframe['TPSA'])
	dict_LogP_WildmanCrippen = dict(dataframe['LogP_WildmanCrippen'])
	dict_Cluster_General = dict(dataframe['Cluster_General'])
	dict_Cluster_Biological = dict(dataframe['Cluster_Biological'])
	dict_Cluster_Physicochemical = dict(dataframe['Cluster_Physicochemical'])
	dict_Cluster_Structural = dict(dataframe['Cluster_Structural'])

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

	output_mols = dict_molecules.values()

	with AllChem.SDWriter(FileOutput) as w:
		for m in output_mols:
			w.write(m)