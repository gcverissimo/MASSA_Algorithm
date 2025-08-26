from rdkit.Chem import AllChem


def output_mols(dataframe, file_output):
    """
    For each molecule, adds the values of the calculated properties,
    the identifications of each cluster and which set the molecule belongs to.
    Also, it creates the final output .sdf file.

    Args:
        dataframe (pd.DataFrame): The dataframe of molecules.
        file_output (str): The path to the output file.
    """
    # Creation of dictionaries with the name of the molecule:
    dict_molecules = dict(
        dataframe['molecules'])  # name:molecule
    dict_set = dict(
        dataframe['set'])  # name:set
    dic_hba = dict(
        dataframe['NumHAcceptors'])  # name:HBA
    dic_hbd = dict(
        dataframe['NumHDonors'])  # name:HBD
    dict_exactmolwt = dict(
        dataframe['ExactMolWt'])  # name:Molecular Weight
    dict_rotbonds = dict(
        dataframe['NumRotatableBonds'])  # name:nRotBonds
    dict_fsp3 = dict(
        dataframe['FractionCSP3'])  # name:Fsp3
    dict_tpsa = dict(
        dataframe['TPSA'])  # name:TPSA
    dict_logp_wc = dict(
        dataframe['LogP_WildmanCrippen'])  # name:logP
    dict_cluster_general = dict(
        dataframe['Cluster_General'])  # name:general_cluster
    dict_cluster_biological = dict(
        dataframe['Cluster_Biological'])  # name:biological_cluster
    dict_cluster_physicochemical = dict(
        dataframe['Cluster_Physicochemical'])  # name:physicochemical_cluster
    dict_cluster_structural = dict(
        dataframe['Cluster_Structural'])  # name:strucutral_cluster

    # Add the properties to molecules from the name:molecule dictionary:
    for i in dict_molecules.keys():
        dict_molecules[i].SetProp('Set', str(dict_set[i]))
        dict_molecules[i].SetProp('HBA count', str(dic_hba[i]))
        dict_molecules[i].SetProp('HBD count', str(dic_hbd[i]))
        dict_molecules[i].SetProp('ExactMolWt', str(dict_exactmolwt[i]))
        dict_molecules[i].SetProp(
            'NumRotatableBonds', str(dict_rotbonds[i]))
        dict_molecules[i].SetProp('FractionCsp3', str(dict_fsp3[i]))
        dict_molecules[i].SetProp('TPSA', str(dict_tpsa[i]))
        dict_molecules[i].SetProp(
            'LogP_WildmanCrippen', str(dict_logp_wc[i]))
        dict_molecules[i].SetProp(
            'General_Cluster', str(dict_cluster_general[i]))
        dict_molecules[i].SetProp(
            'Biological_Cluster', str(dict_cluster_biological[i]))
        dict_molecules[i].SetProp('Physicochemical_Cluster', str(
            dict_cluster_physicochemical[i]))
        dict_molecules[i].SetProp(
            'Structural_Cluster', str(dict_cluster_structural[i]))

    # Create a dict_values with just the molecules
    # (the name of the molecule is already part of molecule object):
    output_mols = dict_molecules.values()

    # Create output ".sdf" file:
    with AllChem.SDWriter(file_output) as w:
        for m in output_mols:
            w.write(m)
