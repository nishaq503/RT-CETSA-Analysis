# %%
from enum import Enum
from pathlib import Path

from bfio import BioWriter
from filepattern import FilePattern

import numpy as np
from tifffile import imread
from pydantic import BaseModel
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.transform import rotate

PLATE_SIZE = 384
PLATE_DIMS = (24, 16)
WELL_RADIUS = 10
WELL_RADIUS_TOLERANCE = 4


class PlateSize(Enum):
    SIZE_6 = 6
    SIZE_12 = 12
    SIZE_24 = 24
    SIZE_48 = 48
    SIZE_96 = 96
    SIZE_384 = 384
    SIZE_1536 = 1536


PLATE_DIMS = {
    PlateSize.SIZE_6: (2, 3),
    PlateSize.SIZE_12: (3, 4),
    PlateSize.SIZE_24: (4, 6),
    PlateSize.SIZE_48: (6, 8),
    PlateSize.SIZE_96: (16, 24),
    PlateSize.SIZE_48: (32, 48),
}

ROTATION = np.vstack(
    [
        -np.sin(np.arange(0, np.pi, np.pi / 180)),
        np.cos(np.arange(0, np.pi, np.pi / 180)),
    ]
)


class PlateParams(BaseModel):
    rotate: int
    """Counterclockwise rotation of image in degrees."""

    bbox: tuple[int, int, int, int]
    """Bounding box of plate after rotation, [ymin,ymax,xmin,xmax]."""

    size: PlateSize
    """The plate size, also determines layout."""

    radii: int
    """"""


def get_wells(image: np.ndarray) -> tuple[list[float], list[float], list[float], int]:
    """Get well locations and radii.

    Since RT-CETSA are generally high signal to noise, no need for anything fance
    to detect wells. Simple Otsu threshold to segment the well, image labeling,
    and estimation of radius based off of area (assuming the area is a circle).

    The input image is a binary image.
    """

    markers, n_objects = ndi.label(image)

    radii = []
    cx = []
    cy = []
    for s in ndi.find_objects(markers):
        cy.append((s[0].start + s[0].stop) / 2)
        cx.append((s[1].start + s[1].stop) / 2)
        radii.append(np.sqrt((markers[s] > 0).sum() / np.pi))

    return cx, cy, radii, n_objects


def get_plate_params(image: np.ndarray) -> PlateParams:

    # Calculate a simple threshold
    threshold = threshold_otsu(image)

    # Get initial well positions
    cx, cy, radii, n_objects = get_wells(image > threshold)

    # Calculate the counterclockwise rotations
    locations = np.vstack([cx, cy]).T
    transform = locations @ ROTATION

    # Find the rotation that aligns the long edge of the plate horizontally
    angle = np.argmin(transform.max(axis=0) - transform.min(axis=0))

    # Shortest rotation to alignment
    if angle > 90:
        angle -= 180

    # Rotate the plate and recalculate well positions
    image_rotated = rotate(image, angle, preserve_range=True)

    # Recalculate well positions
    cx, cy, radii, n_objects = get_wells(image_rotated > threshold)

    # Determine the plate layout
    n_wells = len(cx)
    plate_config = None
    for layout in PlateSize:
        error = abs(1 - n_wells / layout.value)
        if error < 0.05:
            plate_config = layout
            break
    if plate_config is None:
        raise ValueError("Could not determine plate layout")

    # Get the mean radius
    radii_mean = np.mean(radii)

    # Get the bounding box after rotation
    cx_min, cx_max = np.min(cx) - 2 * radii_mean, np.max(cx) + 2 * radii_mean
    cy_min, cy_max = np.min(cy) - 2 * radii_mean, np.max(cy) + 2 * radii_mean
    bbox = (int(cy_min), int(cy_max), int(cx_min), int(cx_max))

    return PlateParams(
        rotate=angle, size=plate_config, radii=int(radii_mean), bbox=bbox
    )


IMAGE_PATH = Path(
    "/home/jovyan/work/RT-CETSA-Analysis/.data/Data for Nick/20210318 LDHA compound plates/20210318 LDHA compound plate 1 6K cells/1.tif"
)
image = imread(IMAGE_PATH)

params = get_plate_params(image)

fp = FilePattern(str(IMAGE_PATH.parent), "{index:d+}.tif")

out_dir = Path(".data/output")
out_dir.mkdir(exist_ok=True)

from matplotlib import pyplot as plt

for index, files in fp(pydantic_output=False):
    for file in files:
        image = imread(file)
        with BioWriter(out_dir.joinpath(file.name.replace(".tif", ".ome.tif"))) as bw:
            image_rotated = rotate(image, params.rotate)[
                params.bbox[0] : params.bbox[1], params.bbox[2] : params.bbox[3]
            ].astype(image.dtype)
            bw.dtype = image_rotated.dtype
            bw.shape = image_rotated.shape
            bw[:] = image_rotated
            plt.figure()
            plt.imshow(
                rotate(image, params.rotate)[
                    params.bbox[0] : params.bbox[1], params.bbox[2] : params.bbox[3]
                ],
                cmap="gray",
            )
            plt.title(f"{index['index']}")

    # indexes.append(index["index"])

# %%
