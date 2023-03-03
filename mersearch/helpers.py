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

def scroll_shim(passed_in_driver, object):
    """Allows firefox driver to scroll to elements without throwing an error"""
    x = object.location['x']
    y = object.location['y']
    scroll_by_coord = 'window.scrollTo(%s,%s);' % (
        x,
        y
    )
    scroll_nav_out_of_way = 'window.scrollBy(0, -120);'
    passed_in_driver.execute_script(scroll_by_coord)
    passed_in_driver.execute_script(scroll_nav_out_of_way)
