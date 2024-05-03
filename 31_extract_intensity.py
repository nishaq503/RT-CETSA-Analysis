# %%
import string
from itertools import product
from pathlib import Path
import numpy as np

import os

from bfio import BioReader
from filepattern import FilePattern

from lib import get_plate_params
from lib import extract_intensity
from lib import index_to_battleship

import csv


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

# This get us very close from the original values.
R += 1
fp = FilePattern(str(IMAGE_PATH.parent), "{index:d+}.ome.tif")

temperature_range = [37, 90]
temperatures = []

with open("data/cleaned_expt1.csv", newline='\n') as csvfile:
    results = csv.reader(csvfile, delimiter='\n')
    for row in list(results)[1:]:
        val = row[0].split(",")[0]
        temperatures.append(val)

print(temperatures)
# Write a csv file
with open(".data/output/plate.csv", "w") as fw:

    headers = ["Temperature"]
    for y, x in product(range(len(params.Y)),range(len(params.X))):
        headers.append(index_to_battleship(x, y, params.size))
    fw.write(",".join(headers) + "\n")

    for index, files in sorted(fp(pydantic_output=False), key=lambda x: x[0]['index']):
        with BioReader(files[0]) as br:
            image = br.read()
            time_step = [temperatures[index['index']-1]] + [
                extract_intensity(image, params.X[x], params.Y[y], R)
                for y, x in product(range(len(params.Y)), range(len(params.X)))
            ]
            fw.write(",".join([str(p) for p in time_step]) + "\n")                        

