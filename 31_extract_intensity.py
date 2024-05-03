# %%
import string
from itertools import product
from pathlib import Path

import os

from bfio import BioReader
from filepattern import FilePattern

from lib import get_plate_params
from lib import extract_intensity
from lib import index_to_battleship


IMAGE_PATH = Path(".data/output/1.ome.tif")

with BioReader(IMAGE_PATH) as br:
    image = br[:]

params = get_plate_params(image)

print(params.X)
print(params.Y)

# Some wells seem overexposed, so we try to grab all light in each grid to compensate
# NOTE take half distance of the smallest distance between wells
# and take that as center of a square around the well.
# this is close to what the matlab code is doing
R = min(params.X[1] - params.X[0], params.Y[1] - params.Y[0]) // 2
R2 = ((params.X[1] - params.X[0]) // 2, (params.Y[1] - params.Y[0]) // 2) 
print(R)
print(R2)

fp = FilePattern(str(IMAGE_PATH.parent), "{index:d+}.ome.tif")

temperature_range = [37, 90]

# Write a csv file
with open(".data/output/plate.csv", "w") as fw:

    headers = []
    for y, x in product(range(len(params.Y)),range(len(params.X))):
        headers.append(index_to_battleship(x, y, params.size))
    fw.write(",".join(headers) + "\n")

    for index, files in sorted(fp(pydantic_output=False), key=lambda x: x[0]['index']):
        with BioReader(files[0]) as br:
            image = br.read()
            time_step = [
                extract_intensity(image, params.X[x], params.Y[y], R)
                for y, x in product(range(len(params.Y)), range(len(params.X)))
            ]
            fw.write(",".join([str(p) for p in time_step]) + "\n")                        

