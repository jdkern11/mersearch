import json
from datetime import datetime
from time import sleep
from typing import Dict, List, Union, Optional
from copy import deepcopy
import logging

from numpy.random import rand
from bs4 import BeautifulSoup
from selenium.webdriver.remote.webdriver import WebDriver
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException, TimeoutException

from mersearch.helpers import scroll, click_until_stale

logger = logging.getLogger(__name__)


def exact_search(
    driver: WebDriver, mol_files: Dict[str, str], add_sleep: bool = True
) -> List[Dict]:
    """Search VWR International website for specific molecule

    Args:
        driver: web driver to search for molecule using.
        mol_files: dictionary with key of smiles string and values of mol files to
            search for.
    Returns:
        List of dictionaries with keys "smiles" and "data", with smiles pointing to
        the smiles string and data pointing to the data parse from the vwr search, or
        an empty dictionary if data missing.
    """
    data = []
    login(driver)
    for smiles, mol_file in mol_files.items():
        logger.info(f"Searching for {smiles} in vwr")
        search_for_mol(driver, mol_file)
        if click_mol(driver):
            data.append({"smiles": smiles, "data": extract_mol_data(driver)})
        else:
            logger.info("Smiles failed")
            data.append({"smiles": smiles, "data": {}})
        if add_sleep:
            sleep(rand())
        driver.get("https://vwr.emolecules.com/index.php")
        if add_sleep:
            sleep(rand())
    return data


def substructure_search(
    driver: WebDriver,
    mol_file: str,
    pages: int = -1,
    add_sleep: bool = True,
    filename: Optional[str] = None,
    substructure: Optional[str] = None,
) -> List[Dict]:
    """Search VWR International website for substructures

    Args:
        driver: web driver to search for molecule using.
        mol_file: mol file to search for.
        pages: how many pages to search. If < 1, all pages will be searched.
        add_sleep: whether or not to add fake sleeping.
        filename: Optional name of file to save data to in json format.
            If not provided, data is not saved.
        substructure: Optional smiles string representing the substructure searched.
    Returns:
        List of dictionaries with keys "smiles" and "data", with smiles pointing to
        the smiles string and data pointing to the data parse from the vwr search, or
        an empty dictionary if data missing.

    """
    if substructure is not None:
        logger.info(f"Searching for substructure {substructure}")
    if filename is not None:
        with open(filename, "a+") as f:
            f.write("[\n")
    data = []
    page = 1
    login(driver)
    search_for_mol(driver, mol_file, "substructure")
    row = 0
    while pages < 1 or page <= pages:
        page_data = []
        logger.info(f"Searching page {page} in vwr")
        while click_mol(driver, row=f"row_{row}.0"):
            try:
                page_data.append(extract_mol_data(driver))
            except Exception as e:
                logger.warning(e)
            if add_sleep:
                sleep(rand())
            driver.back()
            if add_sleep:
                sleep(rand())
            row += 1
        now = str(datetime.now())
        logger.info(f"Page finished reading at time {now}")
        data.append(
            {
                "data": deepcopy(page_data),
                "time": now,
                "substructure": substructure,
                "page": page,
            }
        )
        if filename is not None:
            with open(filename, "a+") as f:
                if page > 1:
                    f.write(",\n")
                json.dump(data[-1], f, indent=4)

        if add_sleep:
            sleep(rand())
        try:
            next_button = WebDriverWait(driver=driver, timeout=30).until(
                EC.presence_of_element_located(
                    (By.XPATH, '//a[contains(text(),"Next")]')
                )
            )
        except (NoSuchElementException, TimeoutException) as e:
            if filename is not None:
                with open(filename, "a+") as f:
                    f.write("\n]")
            return data
        scroll(driver, next_button)
        click_until_stale(next_button)
        page += 1
        if add_sleep:
            sleep(rand())
    if filename is not None:
        with open(filename, "a+") as f:
            f.write("\n]")
    return data


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
    button_link = WebDriverWait(driver=driver, timeout=30).until(
        EC.presence_of_element_located((By.ID, "emolPunchout"))
    )
    # navigate to https://vwr.emolecules.com/index.php
    button = button_link.find_element(By.TAG_NAME, "input")
    scroll(driver, button)
    click_until_stale(button)


def search_for_mol(
    driver: WebDriver, mol_file: str, search_type: str = "exact", add_sleep: bool = True
):
    """Search vwr for the molecule represented by the mol file.

    Navigates to results page for the passed mol_file

    Args:
        driver: WebDriver used for search. It should already be at
            https://vwr.emolecules.com/index.php
        mol_file: molecule to search for.
        search_type: what type of search to run. Can be 'exact' or 'substructure'.
    """
    if search_type not in {"exact", "substructure"}:
        raise ValueError("Only 'exact' or 'substructure' search allowed")
    # find button to load mol_file and click it
    WebDriverWait(driver=driver, timeout=30).until(
        EC.presence_of_element_located((By.XPATH, '//span[@title="Import Molfile"]'))
    ).click()
    # find textarea to add mol_file data to
    chemwriter = WebDriverWait(driver=driver, timeout=30).until(
        EC.presence_of_element_located((By.CLASS_NAME, "chemwriter"))
    )
    mol_file_pane = chemwriter.find_element(By.CLASS_NAME, "content")
    mol_file_textarea = mol_file_pane.find_element(By.TAG_NAME, "textarea")
    if add_sleep:
        sleep(rand())
    mol_file_textarea.clear()
    mol_file_textarea.send_keys(mol_file)
    if add_sleep:
        sleep(rand())
    # find button to write molfile data to chemwriter
    chemwriter.find_element(By.TAG_NAME, "button").click()
    # run exact structure search
    if search_type == "exact":
        WebDriverWait(driver=driver, timeout=30).until(
            EC.presence_of_element_located((By.NAME, "ex_button"))
        ).click()
    elif search_type == "substructure":
        WebDriverWait(driver=driver, timeout=30).until(
            EC.presence_of_element_located((By.NAME, "ss_button"))
        ).click()


def click_mol(driver: WebDriver, row: str = "row_0.0", add_sleep: bool = True) -> bool:
    """Clicks on the mol returned by vwr

    Args:
        driver: WebDriver used for mol search. It should already be searching for the
        mol
        row: the row to click on in the vwr table. Default is row_0.0
    """
    try:
        button = (
            WebDriverWait(driver=driver, timeout=30)
            .until(EC.presence_of_element_located((By.CLASS_NAME, "hitlist-table")))
            .find_element(By.ID, row)
            .find_element(By.CLASS_NAME, "compound_actions_left")
            .find_element(By.TAG_NAME, "img")
        )
        scroll(driver, button)
        click_until_stale(button)
    except (NoSuchElementException, TimeoutException) as e:
        return False
    return True


def extract_mol_data(
    driver: WebDriver, add_sleep: bool = True
) -> Dict[str, List[Dict[str, str]]]:
    """Extract properties, name, supplier and price data for the molecule searched for
    via search_for_mol

    Args:
        driver: WebDriver used for mol search. It should already have run click_mol
    Returns:
        A dictionary with the keys "property_data", "name_data", "supplier_data",
        and "price_data". Each key points to a list of of dictionaries that contain
        the specified data. If search was unsuccessful, an empty dictionary is returned.

        The property_data list contains dictionaries with keys "property" and "value".
        The name_data list contains dictionaries with keys "name_type" and "name".
        The supplier_data list contains dictionaries with keys "supplier" and "supplier_id".
        The price_data list contains dictionaries with keys "supplier", "supplier_id",
        "amount", "unit", and "price".
    """
    soup = BeautifulSoup(driver.page_source, "lxml")

    # find property data
    table = soup.find("table", attrs={"id": "properties_table"})
    table_data = table.tbody.find_all("tr")
    properties_data = []
    # ignore properties row by starting at 1
    for tr in table_data[1:]:
        td = tr.find_all("td")
        try:
            properties_data.append(
                {"property": td[0].string.strip(), "value": td[1].string.strip()}
            )
        except Exception as e:
            logger.warning(e)

    # find name data
    name_data = []
    table = soup.find("div", attrs={"id": "name_table"})
    table = table.find("table", attrs={"class": "data_table"})
    table_data = table.tbody.find_all("tr")
    # ignore known names row
    for tr in table_data[1:]:
        td = tr.find_all("td")
        try:
            name_data.append(
                {
                    "name_type": td[0].string.strip().strip(":"),
                    "name": td[1].string.strip(),
                }
            )
        except Exception as e:
            logger.warning(e)

    # find supplier data
    supplier_data = []
    table = soup.find("div", attrs={"id": "supplier_table"})
    table = table.find("table", attrs={"class": "data_table"})
    table_data = table.tbody.find_all("tr")
    # ignore source, compound id row
    for tr in table_data[1:]:
        td = tr.find_all("td")
        try:
            supplier_data.append(
                {
                    "supplier": td[0].string.strip().strip(":"),
                    "supplier_id": td[1].string.strip(),
                }
            )
        except Exception as e:
            logger.warning(e)

    # find and click more info button to see prices
    if add_sleep:
        sleep(rand())
    WebDriverWait(driver=driver, timeout=30).until(
        EC.presence_of_element_located((By.ID, "add_item_0"))
    ).click()

    soup = BeautifulSoup(driver.page_source, "lxml")
    while soup.find("table", attrs={"class": "bbpricetable"}) is None:
        soup = BeautifulSoup(driver.page_source, "lxml")

    # get costs
    table = soup.find("table", attrs={"class": "bbpricetable"})
    table_data = table.tbody.find_all("tr")
    price_data = []
    for tr in table_data:
        td = tr.find_all("td")
        # Some compounds are prohibited by law or company policy
        if td[0].string == "Prohibited Compound":
            price_data.append("Prohibited Compound")
            break
        # large headers so skip
        if len(td) == 1:
            continue
        # column headers so skip
        elif len(td) == 7:
            continue
        elif len(td) == 8:
            try:
                row = {
                    "supplier": td[0].text.strip(),
                    "supplier_id": td[1].text.split("Name")[0].strip(),
                }
                row["amount"] = td[5].text.strip()
                row["units"] = td[6].text.strip()
                row["price"] = td[7].text.strip()
                price_data.append(row.copy())
            except Exception as e:
                logger.warning(e)
        elif len(td) == 5:
            try:
                row["amount"] = td[2].text.strip()
                row["unit"] = td[3].text.strip()
                row["price"] = td[4].text.strip()
                price_data.append(row.copy())
            except Exception as e:
                logger.warning(e)

    return {
        "property_data": deepcopy(properties_data),
        "name_data": deepcopy(name_data),
        "supplier_data": deepcopy(supplier_data),
        "price_data": deepcopy(price_data),
    }
