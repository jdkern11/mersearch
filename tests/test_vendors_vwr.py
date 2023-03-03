import pytest
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.support import expected_conditions as EC

from mersearch.vendors import vwr
from .conf import configuration, webdriver

@pytest.fixture
def driver():
    driver = webdriver(**configuration)
    yield driver
    driver.quit()

def test_login(driver):
    vwr.login(driver)
    try:
        WebDriverWait(driver=driver, timeout=30).until(EC.presence_of_element_located(
            (By.XPATH, '//span[@title="Import Molfile"]')
        ))
    except NoSuchElementException:
        pytest.fail("Login failed")

