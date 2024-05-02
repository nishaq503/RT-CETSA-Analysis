# %%
import string
from itertools import product
from pathlib import Path

from bfio import BioReader
from filepattern import FilePattern

from lib import get_plate_params
from lib import extract_intensity
from lib import index_to_battleship


IMAGE_PATH = Path("/home/jovyan/work/RT-CETSA-Analysis/.data/output/1.ome.tif")

with BioReader(IMAGE_PATH) as br:
    image = br[:]

params = get_plate_params(image)

# Some wells seem overexposed, so we try to grab all light in each grid to compensate
R = min(params.X[1] - params.X[0], params.Y[1] - params.Y[0]) // 2

fp = FilePattern(str(IMAGE_PATH.parent), "{index:d+}.ome.tif")

temperature_range = [37, 90]

# Write a csv file
with open(".data/output/plate.csv", "w") as fw:

    # Write the headers
    headers = ["Temperature"]
    for x, y in product(range(len(params.X)), range(len(params.Y))):
        headers.append(index_to_battleship(x, y, params.size))
    fw.write(",".join(headers) + "\n")

    for index, files in fp(pydantic_output=False):
        for file in files:

            # Get the temperature
            temp = temperature_range[0] + (index["index"] - 1) / (len(fp) - 1) * (
                temperature_range[1] - temperature_range[0]
            )
            plate = [f"{temp:.1f}"]
            with BioReader(file) as br:
                image = br.read()

                for x, y in product(range(len(params.X)), range(len(params.Y))):
                    plate.append(extract_intensity(image, params.X[x], params.Y[y], R))

            fw.write(",".join([str(p) for p in plate]) + "\n")
