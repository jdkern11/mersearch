from typing import Dict, List

from bs4 import BeautifulSoup
from selenium.webdriver.remote.webdriver import WebDriver
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import (
    NoSuchElementException, 
    StaleElementReferenceException
)

from mersearch.helpers import scroll_shim

def search(driver: WebDriver, mol_files: Dict[str, str]) -> Dict[str, str]:
    """Search VWR International website for specific molecule

    Args:
        driver: web driver to search for molecule using.
        mol_files: dictionary with key of smiles string and values of mol files to 
            search for.
    Returns:
        Dictionary with data parsed from the VWR search. If no data found, dictionary
        will be empty.
    """
    login(driver)
    for smiles, mol_file in mol_files.items():
        search_for_mol(driver, mol_file)
        extract_mol_data(driver)

def login(driver: WebDriver):
    """Logs in to vwr by navigating to https://us.vwr.com/store/search/searchMol.jsp,
    then continuing to https://vwr.emolecules.com/index.php from a button on the first
    page

    Args:
        driver: WebDriver used for search
    """
    # navigate to this site to get auto logged in so search can be performed
    url = "https://us.vwr.com/store/search/searchMol.jsp"
    driver.get(url)
    # navigate to actual search site
    button_link = WebDriverWait(
        driver=driver, timeout=30
    ).until(EC.presence_of_element_located((By.ID, "emolPunchout")))
    button = button_link.find_element(By.TAG_NAME, 'input')
    scroll_shim(driver, button)
    actions = ActionChains(driver)
    # scroll to the button
    actions.move_to_element(button_link).perform()
    # navigate to https://vwr.emolecules.com/index.php
    # click button until button is actually registered
    try:
        while True:
            button.click()
    except StaleElementReferenceException:
        pass

def search_for_mol(driver: WebDriver, mol_file: str):
    """Search vwr for the molecule represented by the mol file.

    Navigates to results page for the passed mol_file

    Args:
        driver: WebDriver used for search. It should already be at 
            https://vwr.emolecules.com/index.php
        mol_file: molecule to search for.
    """
    # find button to load mol_file and click it
    WebDriverWait(driver=driver, timeout=30).until(EC.presence_of_element_located(
        (By.XPATH, '//span[@title="Import Molfile"]')
    )).click()
    # find textarea to add mol_file data to
    chemwriter = WebDriverWait(
        driver=driver, timeout=30
    ).until(EC.presence_of_element_located((By.CLASS_NAME, "chemwriter")))
    mol_file_pane = chemwriter.find_element(By.CLASS_NAME, 'content')
    mol_file_textarea = mol_file_pane.find_element(By.TAG_NAME, 'textarea')
    mol_file_textarea.clear()
    mol_file_textarea.send_keys(smiles_to_mol_file_format(mol_file))
    # find button to write molfile data to chemwriter
    chemwriter.find_element(By.TAG_NAME, 'button').click()
    # run exact structure search
    WebDriverWait(
        driver=driver, timeout=30
    ).until(EC.presence_of_element_located((By.NAME, "ex_button"))).click()
