from pathlib import Path

from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options

profile_path = str(
    Path("~")
    / "Library"
    / "Application Support"
    / "Firefox"
    / "Profiles"
    / "pwofczxw.default"
)
PROXY_HOST = "12.12.12.123"
PROXY_PORT = "1234"
options = Options()
options.set_preference("profile", profile_path)
options.set_preference("network.proxy.type", 1)
options.set_preference("network.proxy.http", PROXY_HOST)
options.set_preference("network.proxy.http_port", int(PROXY_PORT))
options.set_preference("dom.webdriver.enabled", False)
options.set_preference("useAutomationExtension", False)
options.set_preference("excludeSwitches", "enable-automation")
options.add_argument("--headless")
service = Service(r"/usr/local/bin/geckodriver")

configuration = {"service": service, "options": options}
webdriver = webdriver.Firefox
