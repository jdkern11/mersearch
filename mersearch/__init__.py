from typing import List

import pandas as pd

from mersearch.helpers import smiles_to_mol_file_format
from mersearch.vendors import search_vwr

def search(smiles_list: List[str]) -> pd.DataFrame:
    """Search for molecules from a variety of vendors and return vendor information

    Args:
        smiles_list: List of smiles strings to collect information on
    Returns:
        pandas dataframe with smiles column and vendors column
    """
    data = [{'smiles': smiles} for smiles in smiles_list]
    mol_files = {smiles: smiles_to_mol_file_format(smiles) for smiles in smiles_list}
