from rdkit import Chem
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import (
    StaleElementReferenceException,
    ElementNotInteractableException,
)


def smiles_to_mol_file_format(smiles: str) -> str:
    """Converts smiles string to mol file format

    Args:
        smiles (str): smiles to convert.
    Returns:
        String in mol file format.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol,explicitOnly=True)
    return Chem.MolToMolBlock(mol)


def scroll(driver, object):
    scroll_shim(driver, object)
    actions = ActionChains(driver)
    # scroll to the object
    try:
        actions.move_to_element(object).perform()
    except ElementNotInteractableException:
        pass


def scroll_shim(driver, object):
    """Allows firefox driver to scroll to elements without throwing an error"""
    x = object.location["x"]
    y = object.location["y"]
    scroll_by_coord = "window.scrollTo(%s,%s);" % (x, y)
    scroll_nav_out_of_way = "window.scrollBy(0, -120);"
    driver.execute_script(scroll_by_coord)
    driver.execute_script(scroll_nav_out_of_way)


def click_until_stale(object):
    try:
        while True:
            object.click()
    except StaleElementReferenceException:
        pass
