from pathlib import Path

import pytest
from rdkit import Chem

from mersearch import helpers


@pytest.fixture
def mol_file_loc():
    return Path(__file__).parent.absolute() / "test_data" / "mol_files"


def test_smiles_to_mol_file_format(mol_file_loc):
    """Testing our mol file format output against those from https://datascience.unm.edu/tomcat/biocomp/convert"""
    files = [f for f in mol_file_loc.iterdir() if f.suffix == ".mol"]
    header_len = 2
    for f in files:
        smiles = f.name.split(".")[0]
        mol_from_smiles = helpers.smiles_to_mol_file_format(smiles)
        m1 = Chem.rdmolfiles.MolFromMolBlock(mol_from_smiles)
        with open(f, "r") as fin:
            mol_from_file = "".join(fin.readlines())
        m2 = Chem.rdmolfiles.MolFromMolBlock(mol_from_file)
        m1 = Chem.MolToSmiles(m1)
        m2 = Chem.MolToSmiles(m2)
        assert Chem.CanonSmiles(m1) == Chem.CanonSmiles(m2)
