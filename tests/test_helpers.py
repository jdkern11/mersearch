from pathlib import Path

import pytest

from mersearch import helpers

@pytest.fixture
def mol_file_loc():
    return Path(__file__).parent.absolute() / 'test_data' / 'mol_files'

def test_smiles_to_mol_file_format(mol_file_loc):
    """Testing our mol file format output against those from https://datascience.unm.edu/tomcat/biocomp/convert"""
    files = [f for f in mol_file_loc.iterdir() if f.suffix == '.mol']
    header_len = 2
    for f in files:
        smiles = f.name.split('.')[0]
        mol_from_smiles = helpers.smiles_to_mol_file_format(smiles)
        with open(f, 'r') as fin:
            mol_from_file = ''.join(fin.readlines())
        mol_from_smiles = mol_from_smiles.split('\n')
        mol_from_file = mol_from_file.split('\n')
        num_atoms = 0
        num_bonds = 0
        for i in range(header_len + 1, len(mol_from_smiles) - 2):
            s_line = ' '.join(mol_from_smiles[i].strip().split()).split()
            f_line = ' '.join(mol_from_file[i].strip().split()).split()
            if i == header_len + 1:
                num_atoms = int(s_line[0])
                num_bonds = int(s_line[1])
                assert num_atoms == int(f_line[0])
                assert num_bonds == int(f_line[1])
            elif i < header_len + 1 + num_atoms + 1:
                element = s_line[3]
                assert element == f_line[3]
            else:
                print(s_line)
                print(f_line)
                s_bond = s_line[0:4]
                f_bond = f_line[0:4]
                s_bond.sort()
                f_bond.sort()
                assert s_bond == f_bond
