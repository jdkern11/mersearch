from enum import Enum
from typing import List, Set, Dict

import pandas as pd
from selenium.webdriver.remote.webdriver import WebDriver

from mersearch.helpers import smiles_to_mol_file_format
from mersearch.vendors import vwr


class VendorsEnum(Enum):
    """Up to date enumeration of vendors the search is implemented for"""

    vwr = "vwr"


def exact_search(
    driver: WebDriver,
    smiles_list: List[str],
    vendors: Set[VendorsEnum] = {VendorsEnum.vwr},
):  # TODO add type hint
    """Search for molecules from a variety of vendors and return vendor information

    Args:
        driver: web driver to search for smiles using
        smiles_list: List of smiles strings to collect information on
        vendors: Set of vendors to run search using
    Returns:
        pandas dataframe with smiles column and vendors column
    """
    data = {}
    mol_files = {smiles: smiles_to_mol_file_format(smiles) for smiles in smiles_list}
    if VendorsEnum.vwr in vendors:
        data["vwr_data"] = vwr.exact_search(driver, mol_files)
    return data
