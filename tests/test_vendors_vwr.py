import time

import pytest
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.support import expected_conditions as EC

from mersearch.vendors import vwr
from mersearch.helpers import smiles_to_mol_file_format
from .conf import configuration, webdriver


@pytest.fixture
def driver():
    driver = webdriver(**configuration)
    yield driver
    driver.close()


def test_login(driver):
    vwr.login(driver)
    try:
        WebDriverWait(driver=driver, timeout=30).until(
            EC.presence_of_element_located(
                (By.XPATH, '//span[@title="Import Molfile"]')
            )
        )
    except NoSuchElementException:
        pytest.fail("Login failed")


def test_search_for_mol(driver):
    vwr.login(driver)
    # this is a valid molecule to search for
    with pytest.raises(ValueError):
        vwr.search_for_mol(
            driver,
            smiles_to_mol_file_format("O=C1CCCCCCCCCCCCCCO1"),
            search_type="Bad search",
        )
    vwr.search_for_mol(driver, smiles_to_mol_file_format("O=C1CCCCCCCCCCCCCCO1"))
    try:
        WebDriverWait(driver=driver, timeout=30).until(
            EC.presence_of_element_located((By.CLASS_NAME, "hitlist-table"))
        )
    except NoSuchElementException:
        pytest.fail("Search failed")


def test_search_for_substructure(driver):
    vwr.login(driver)
    # this is a valid molecule to search for
    vwr.search_for_mol(
        driver, smiles_to_mol_file_format("O=C1CCCO1"), search_type="substructure"
    )
    try:
        WebDriverWait(driver=driver, timeout=30).until(
            EC.presence_of_element_located((By.CLASS_NAME, "hitlist-table"))
        )
    except NoSuchElementException:
        pytest.fail("Search failed")


def test_extract_mol_data(driver):
    vwr.login(driver)
    # this is a valid molecule to search for
    vwr.search_for_mol(driver, smiles_to_mol_file_format("O=C1CCCCCCCCCCCCCCO1"))
    assert vwr.click_mol(driver)
    data = vwr.extract_mol_data(driver)
    assert "property_data" in data
    assert "name_data" in data
    assert "supplier_data" in data
    assert "price_data" in data


def test_exact_search(driver):
    smiles_list = ["O=C1CCCCCCCCCCCCCCO1", "COC(=O)C1(C)CCCC2(C)C3CC(=O)OCC3CCC12"]
    mol_files = {smiles: smiles_to_mol_file_format(smiles) for smiles in smiles_list}
    data = vwr.exact_search(driver, mol_files)
    for row in data:
        if row["smiles"] == "O=C1CCCCCCCCCCCCCCO1":
            assert len(row["data"]) == 4
        else:
            assert len(row["data"]) == 0


def test_substructure_search(driver):
    data = vwr.substructure_search(
        driver, smiles_to_mol_file_format("O=C1CCCO1"), pages=2, add_sleep=False
    )
    assert len(data) == 2
