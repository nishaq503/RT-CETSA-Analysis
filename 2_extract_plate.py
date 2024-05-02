# %%

from pathlib import Path

from bfio import BioWriter
from filepattern import FilePattern
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from skimage.filters import threshold_otsu
from skimage.transform import rotate
from tifffile import imread

from lib import get_plate_params, get_wells

DATA_ROOT_DIR = "/Users/gerardinad/Documents/data/rt_cetsa_data"
IMAGE_PATH = Path(
    f"{DATA_ROOT_DIR}/20210318 LDHA compound plates/20210318 LDHA compound plate 1 6K cells/1.tif"
)

image = imread(IMAGE_PATH)

params = get_plate_params(image)

fp = FilePattern(str(IMAGE_PATH.parent), "{index:d+}.tif")

out_dir = Path(".data/output")
out_dir.mkdir(exist_ok=True)

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

            # Visualize one of the plates with a well overlay
            if index["index"] == 1:
                plt.figure()
                plt.imshow(
                    image_rotated,
                    cmap="gray",
                    vmin=image_rotated.min(),
                    vmax=image_rotated.max(),
                )

                plt.title(f"{index['index']}")
                ax = plt.gca()

                from itertools import product

                for x, y in product(params.X, params.Y):
                    circle = Circle(
                        (x - params.bbox[2], y - params.bbox[0]),
                        radius=params.radii,
                        color="red",
                    )
                    ax.add_patch(circle)
