import pandas as pd
from rdkit.Chem import AllChem


def find_smiles_column(columns):
    """
    Find the name of the "SMILES" column.

    Args:
        columns (list): receives a list with the name of pandas columns.

    Raises:
        ValueError: if there is no column named "SMILES" raises an error.

    Returns:
        column_name (str): returns the name of the "SMILES" column.
    """
    smiles_columns = []
    smile_column = ""

    for column in columns:
        if "smiles" == column.lower():
            smiles_columns.append(column)
    if len(smiles_columns) > 1:
        """
        Warning: If there are more than one "SMILES" (case insensitive) column
        the first one will be chosen.
        """

        print("Warning: There is more than one"
              " column named 'smiles' (case insensitive)."
              " The first one was chosen.")
        smile_column = smiles_columns[0]
        return smile_column
    elif len(smiles_columns) == 1:
        smile_column = smiles_columns[0]
        return smile_column
    else:
        error = "Error: There is no column named 'SMILES'"
        raise ValueError(error)


def find_name_column(columns):
    """
    Find the name of the molecule name column.

    Args:
        columns (list): receives a list with the name of pandas columns.

    Returns:
        column_name (str): returns the name of the molecule name column.
    """

    for column in columns:
        if "molecule name" == column.lower():
            print((f"Found '{column}' column."
                   " Using it to name the molecules."))
            name_column = column
            return name_column
        elif "name" == column.lower():
            print((f"Found '{column}' column."
                   " Using it to name the molecules."))
            name_column = column
            return name_column
        elif "molecule chembl id" == column.lower():
            print((f"Found '{column}' column."
                   " Using it to name the molecules."))
            name_column = column
            return name_column
        elif "id" == column.lower():
            print((f"Found '{column}' column."
                   " Using it to name the molecules."))
            name_column = column
            return name_column
        else:
            print(("Unable to find a column with molecule names."
                   " Using the index to name the molecules."))
            name_column = "index"
            return name_column


def read_df_smiles(df, WriteLog):
    """Reads an dataframe and returns a sanitized molecule supplier.

    Args:
        df_supplier (pd.DataFrame): pd.DataFrame object with SMILES.
        WriteLog (file object): Log file object to write error messages.

    Returns:
        AllChem.SDMolSupplier like: Sanitized molecule supplier.
    """
    columns = list(df.columns)

    smi_name_column = find_smiles_column(columns)
    name_name_column = find_name_column(columns)

    # name_name_column = "index"
    columns.remove(smi_name_column)
    columns.remove(name_name_column)

    mols = []
    for index, row in df.iterrows():
        try:
            mol = AllChem.MolFromSmiles(row[smi_name_column])
            AllChem.rdmolops.SanitizeMol(mol)
            if name_name_column == "index":
                mol.SetProp("_Name", str(index))
            else:
                mol.SetProp("_Name", str(row[name_name_column]))
            for c in columns:
                mol.SetProp(c, str(row[c]))
            mols.append(mol)
        except Exception as e:
            mols.append(None)
            mol_of_error = str(mol.GetProp('_Name'))
            error = str(e)
            precomposed = (f"[Error while reading molecule: {error}. "
                           f"Molecule name: {mol_of_error}].")
            composed = precomposed+'\n'
            print(precomposed)
            WriteLog.write(composed)
    if None in mols:
        print("ERROR: Error while reading molecules. "
              "Check the log.txt file for more information."
              "\nERROR: Closing...")
        exit()
    return mols


def read_EXCEL_smiles(file, WriteLog):
    df_supplier = pd.read_excel(file)
    mols_sanitized = read_df_smiles(df_supplier, WriteLog)
    return mols_sanitized


def read_CSV_smiles(file, WriteLog):
    df_supplier = pd.read_csv(file, sep=",")
    mols_sanitized = read_df_smiles(df_supplier)
    return mols_sanitized


def read_SDF(file, WriteLog):
    """Reads an SDF file and returns a sanitized molecule supplier.

    Args:
        file (str): Path to the SDF file.
        WriteLog (file object): Log file object to write error messages.

    Returns:
        AllChem.SDMolSupplier: Sanitized molecule supplier.
    """
    # Open the input file.
    mols_supplier = AllChem.SDMolSupplier(file, sanitize=False)
    mols = []
    for mol in mols_supplier:
        try:
            mols.append(AllChem.rdmolops.SanitizeMol(mol))
        except Exception as e:
            mols.append(None)
            mol_of_error = str(mol.GetProp('_Name'))
            error = str(e)
            precomposed = (f"[Error while reading molecule: {error}. "
                           f"Molecule name: {mol_of_error}].")
            composed = precomposed+'\n'
            print(precomposed)
            WriteLog.write(composed)
    if None in mols:
        print("ERROR: Error while reading molecules. "
              "Check the log.txt file for more information."
              "\nERROR: Closing...")
        exit()
    mols_sanitized = AllChem.SDMolSupplier(file, sanitize=True)
    return mols_sanitized


def remove_comment_line(data):
    # Remove lines starting with '#'
    return '\n'.join(
        line for line in data.splitlines() if not line.strip().startswith('#')
        )


def mol2_properties(data):
    data_split = data.split("@")
    properties = {}
    for d in data_split:
        if "MOL_PROPERTY" in d:
            list_values = d.split("\n")
            p_name = list_values[1]
            p_value = list_values[3]
            properties[p_name] = p_value
    return properties


def read_MOL2(file, WriteLog):
    """Reads an MOL2 file from Discovery Studio
       and returns a sanitized molecule supplier.

    Args:
        file (str): Path to the MOL2 file.
        WriteLog (file object): Log file object to write error messages.

    Returns:
        AllChem.SDMolSupplier like: Sanitized molecule supplier.
    """
    with open(file, 'r') as f:
        mol_data = f.read()
    # Split the file contents on the delimiter
    # Remove comment lines
    mol_data = remove_comment_line(mol_data)

    # Alert of MOL2 limitations.
    print("""
ATTENTION: .mol2 file detected.
MOL2 files have limited support for storing molecular properties.
Therefore, the biological properties in this file may not be \
detected or even stored.
We do not recommend using a .mol2 file.
However, we provide support for .mol2 files generated by \
Discovery Studio Visualizer.
The biological properties in these files follow the format below:

@<SCITEGIC>MOL_PROPERTY
Name_of_property (e.g., pIC50 value)
Type_of_property (e.g., SciTegic.value.DoubleValue)
Value_of_the_property (e.g., 8.100000)
          """)

    # Split the file contents based on a delimiter for each molecule.
    mol_blocks = mol_data.split("@<TRIPOS>MOLECULE")

    # List to store molecules.
    mols = []

    # Iterate over each block and reassemble the molecule data
    for mol_block in mol_blocks:
        if mol_block.strip():  # To avoid empty blocks
            mol_block = "@<TRIPOS>MOLECULE" + mol_block
            properties = mol2_properties(mol_block)
            mol = AllChem.MolFromMol2Block(mol_block)
            try:
                AllChem.rdmolops.SanitizeMol(mol)
                for p, p_value in properties.items():
                    mol.SetProp(p, str(p_value))
                mols.append(mol)
            except Exception as e:
                mols.append(None)
                mol_of_error = str(mol.GetProp('_Name'))
                error = str(e)
                precomposed = (f"[Error while reading molecule: {error}. "
                               f"Molecule name: {mol_of_error}].")
                composed = precomposed+'\n'
                print(precomposed)
                WriteLog.write(composed)
    if None in mols:
        print("ERROR: Error while reading molecules. "
              "Check the log.txt file for more information."
              "\nERROR: Closing...")
        exit()
    return mols


def read_molecules(file, WriteLog):
    if file.endswith(".mol2"):
        molecules = read_MOL2(file, WriteLog)
    elif file.endswith(".sdf") or file.endswith(".mol"):
        molecules = read_SDF(file, WriteLog)
    elif file.endswith(".xlsx") or file.endswith(".xls"):
        molecules = read_EXCEL_smiles(file, WriteLog)
    elif file.endswith(".csv"):
        molecules = read_CSV_smiles(file, WriteLog)
    return molecules


def hydrogen_add(file):
    molsH = []
    for i in file:
        i = AllChem.AddHs(i, addCoords=True)
        molsH.append(i)
        # Add hydrogens keeping 3D coordenates.
    return molsH


def get_sdf_property_names(file):
    # Extracting the property names from the ".sdf" input file.
    sdf_property_names = []
    for i in file:
        for a in list(i.GetPropNames()):
            sdf_property_names.append(a)
    sdf_property_names = list(set(sdf_property_names))
    if 'Name' in sdf_property_names:
        # Drop the molecule name from the property list.
        sdf_property_names.remove('Name')
    return sdf_property_names
