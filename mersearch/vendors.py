from typing import Dict, List

from bs4 import BeautifulSoup
from selenium.webdriver.remote.webdriver import WebDriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.support.wait import WebDriverWait

def search_vwr(driver: WebDriver, smiles_list: List[str]) -> Dict[str, str]:
    """Search VWR International website for specific molecule

    Args:
        driver: web driver to search for molecule using.
        smiles_list: list of smiles string to search for
    Returns:
        Dictionary with data parsed from the VWR search. If no data found, dictionary
        will be empty.
    """
    login_to_vwr(driver)

def login_to_vwr(driver: WebDriver):
    """Logs in to vwr by navigating to https://us.vwr.com/store/search/searchMol.jsp,
    then continuing to https://vwr.emolecules.com/index.php from a button on the first
    page

    Args:
        driver: WebDriver used for search
    """
