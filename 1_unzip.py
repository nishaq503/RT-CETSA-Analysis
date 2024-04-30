import zipfile
from pathlib import Path

path = Path(".data")
path.mkdir(exist_ok=True)

with zipfile.ZipFile("Data for Nick.zip", "r") as zip_ref:
    zip_ref.extractall(".data/")
