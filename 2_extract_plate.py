# %%

from pathlib import Path

from bfio import BioWriter
from filepattern import FilePattern
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from skimage.filters import threshold_otsu
from skimage.transform import rotate
from tifffile import imread
import skimage
import itertools

from lib import get_plate_params, get_wells
from dotenv import load_dotenv, find_dotenv
import os
env_file_path = find_dotenv()
load_dotenv(env_file_path, override=True)

def draw_region(img, cx,cy,radius):
    max = img.max()
    img[cy-2:cy+2, cx-2:cx+2] = max
    circle = skimage.draw.circle_perimeter(cy,cx, radius, shape=img.shape)
    img[circle] = max

DATA_DIR = os.getenv("DATA_DIR")
if not DATA_DIR:
    raise Exception("create an .env file with DATA_DIR=/path/to/data")

IMAGE_PATH = Path(
    f"{DATA_DIR}/20210318 LDHA compound plates/20210318 LDHA compound plate 1 6K cells/1.tif"
)

image = imread(IMAGE_PATH)

# NOTE Compute for first plate only
params = get_plate_params(image)

fp = FilePattern(str(IMAGE_PATH.parent), "{index:d+}.tif")

out_dir= Path(".data/output")
out_dir_check = Path(".data/output_check")
out_dir_check.mkdir(exist_ok=True)

for index, files in fp(pydantic_output=False):
    for file in files:
        image = imread(file)

        with BioWriter(out_dir.joinpath(file.name.replace(".tif", ".ome.tif"))) as bw:
            image_rotated = rotate(image, params.rotate, preserve_range=True)[
                params.bbox[0] : params.bbox[1], params.bbox[2] : params.bbox[3]
            ].astype(image.dtype)
            bw.dtype = image_rotated.dtype
            bw.shape = image_rotated.shape
            bw[:] = image_rotated


        with BioWriter(out_dir_check.joinpath(file.name.replace(".tif", ".ome.tif"))) as bw:
            image_rotated = rotate(image, params.rotate, preserve_range=True)[
                params.bbox[0] : params.bbox[1], params.bbox[2] : params.bbox[3]
            ].astype(image.dtype)
            

            for x, y in itertools.product(params.X, params.Y):
                    draw_region(image_rotated,
                                x - params.bbox[2],
                                y - params.bbox[0],
                                radius=params.radii
                    )
            
            bw.dtype = image_rotated.dtype
            bw.shape = image_rotated.shape
            bw[:] = image_rotated
            

            # Visualize one of the plates with a well overlay
            # if index["index"] == 1:
            #     plt.figure()
            #     plt.imshow(
            #         image_rotated,
            #         cmap="gray",
            #         vmin=image_rotated.min(),
            #         vmax=image_rotated.max(),
            #     )

            #     plt.title(f"{index['index']}")
            #     ax = plt.gca()

            #     from itertools import product

                # for x, y in product(params.X, params.Y):
                #     circle = Circle(
                #         (x - params.bbox[2], y - params.bbox[0]),
                #         radius=params.radii,
                #         color="red",
                #     )
            #         ax.add_patch(circle)

# %%
