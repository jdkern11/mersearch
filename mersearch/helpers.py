from rdkit import Chem

def smiles_to_mol_file_format(smiles: str) -> str:
    """Converts smiles string to mol file format

    Args:
        smiles (str): smiles to convert.
    Returns:
        String in mol file format.
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToMolBlock(mol)
