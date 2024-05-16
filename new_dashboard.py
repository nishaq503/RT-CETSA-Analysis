import pathlib

import filepattern
import matplotlib.pyplot as plt
import pandas
import solara
import tifffile

from polus.images.segmentation.rt_cetsa_plate_extraction.__main__ import (
    main as extract_plates,
)
from polus.images.features.rt_cetsa_intensity_extraction.__main__ import (
    main as extract_intensities,
)
from polus.tabular.regression.rt_cetsa_moltprot.__main__ import (
    main as train_moltprot,
)


INP_DATA_DIR = "/home/nishaq/Documents/axle/data/Data for Nick/20210318 LDHA compound plates/20210318 LDHA compound plate 1 6K cells"

inp_dir_widget = solara.reactive(INP_DATA_DIR)
img_pattern_widget = solara.reactive("{id:d+}.tif")
inp_img_index_widget = solara.reactive(1)
ext_img_index_widget = solara.reactive(1)
row_idx_widget = solara.reactive(1)
col_idx_widget = solara.reactive(1)
continuous_update = solara.reactive(True)
show_ext_widget = solara.reactive(False)
show_int_widget = solara.reactive(False)
show_moltprot_widget = solara.reactive(False)
show_mask_widget = solara.reactive(False)
show_cell_widget = solara.reactive(False)


def show_plates(out_dir: pathlib.Path):
    out_pattern = "{id:d+}.ome.tiff"
    fp = filepattern.FilePattern(out_dir, out_pattern)
    img_paths = [f[1][0] for f in fp()]
    img_paths.sort(key=lambda x: int(x.name.split(".")[0]))
    solara.Markdown(f"**Number of images**: {len(img_paths)}")

    img_paths_dict = {
        i: img_path for i, img_path in enumerate(img_paths, start=1)
    }

    # Add a slider to navigate through the images
    solara.SliderInt(
        "Select image",
        value=ext_img_index_widget,
        min=1,
        max=len(img_paths),
    )
    # Display the selected image
    img_path = img_paths_dict[ext_img_index_widget.value]
    img = tifffile.imread(str(img_path))
    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.set_title(f"Extracted Plate: {img_path.name}")
    solara.display(fig)
    plt.close(fig)

    return


def show_intensities_df(out_dir: pathlib.Path):
    # load the dataframe
    plate_df = pandas.read_csv(out_dir / "plate.csv")
    solara.DataFrame(plate_df, scrollable=True, items_per_page=10)
    return


def show_moltprot_df(out_dir: pathlib.Path):
    # load the dataframe
    moltprot_df = pandas.read_csv(out_dir / "plate_moltprot.csv")
    solara.DataFrame(moltprot_df, scrollable=True, items_per_page=10)
    return


def show_mask(out_dir: pathlib.Path):
    mask_path = out_dir / "mask.ome.tiff"
    mask = tifffile.imread(str(mask_path))
    fig, ax = plt.subplots()
    ax.imshow(mask)
    ax.set_title("Mask")
    solara.display(fig)
    plt.close(fig)
    return


def show_cell_plot(plate_df: pandas.DataFrame):
    # Select a row and column
    solara.SliderInt(
        "Select column",
        value=col_idx_widget,
        min=1,
        max=24,
    )
    solara.SliderInt(
        "Select row",
        value=row_idx_widget,
        min=1,
        max=16,
    )
    col_idx = col_idx_widget.value
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    col = letters[col_idx - 1]
    row_idx = row_idx_widget.value
    cell_name = f"{col}{row_idx}"
    solara.Markdown(f"**Cell name**: {cell_name}")

    # Get column from plate.csv
    intensities = plate_df[cell_name].values
    temperatures = plate_df["Temperature"].values
    fig, ax = plt.subplots()
    ax.scatter(temperatures, intensities)
    ax.set_xlabel("Temperature")
    ax.set_ylabel("Intensity")
    ax.set_title(cell_name)
    solara.display(fig)
    plt.close(fig)

    return


@solara.component
def Page():

    solara.Checkbox(label="Continuous update", value=continuous_update)
    solara.InputText("Path to input directory ...", value=inp_dir_widget, continuous_update=continuous_update.value)
    with solara.Row():
        solara.Button("Reset", on_click=lambda: inp_dir_widget.set(INP_DATA_DIR))
    
    inp_dir = pathlib.Path(inp_dir_widget.value).resolve()
    if not inp_dir.exists():
        return solara.Markdown(f"**Path does not exist**: {inp_dir}")
    if not inp_dir.is_dir():
        return solara.Markdown(f"**Not a directory**: {inp_dir}")

    solara.Markdown(f"**Input Directory**: {inp_dir}")

    out_dir = inp_dir.parent / "out-plate-dashboard"
    out_dir.mkdir(exist_ok=True)

    solara.Markdown(f"**Output Directory**: {out_dir}")

    solara.InputText("Image pattern ...", value=img_pattern_widget, continuous_update=continuous_update.value)
    img_pattern = img_pattern_widget.value
    solara.Markdown(f"**Image Pattern**: {img_pattern}")

    inp_fp = filepattern.FilePattern(inp_dir, img_pattern)
    img_paths = [f[1][0] for f in inp_fp()]
    img_paths.sort(key=lambda x: int(x.name.split(".")[0]))

    solara.Markdown(f"**Number of images**: {len(img_paths)}")

    img_paths_dict = {
        i: img_path for i, img_path in enumerate(img_paths, start=1)
    }

    # Add a slider to navigate through the images
    solara.SliderInt(
        "Select image",
        value=inp_img_index_widget,
        min=1,
        max=len(img_paths),
    )
    # Display the selected image
    img_path = img_paths_dict[inp_img_index_widget.value]
    img = tifffile.imread(str(img_path))
    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.set_title(f"Raw Image: {img_path.name}")
    solara.display(fig)
    plt.close(fig)

    # Extract the plates
    solara.Button("Extract Plates", on_click=lambda: extract_plates(
        inp_dir=inp_dir,
        pattern=img_pattern,
        preview=False,
        out_dir=out_dir,
    ))
    solara.Checkbox(label="Show Extracted Plate", value=show_ext_widget)
    if show_ext_widget.value:
        show_plates(out_dir)

    # Extract the intensities
    solara.Button("Extract Intensities", on_click=lambda: extract_intensities(
        inp_dir=out_dir,
        pattern="{id:d+}.ome.tiff",
        preview=False,
        out_dir=out_dir,
    ))
    solara.Checkbox(label="Show Intensities DataFrame", value=show_int_widget)
    if show_int_widget.value:
        show_intensities_df(out_dir)

    # Run MoltProt regression
    solara.Button("Run MoltProt", on_click=lambda: train_moltprot(
        inp_dir=out_dir,
        pattern="plate.csv",
        preview=False,
        out_dir=out_dir,
    ))
    solara.Checkbox(label="Show MoltProt DataFrame", value=show_moltprot_widget)
    if show_moltprot_widget.value:
        show_moltprot_df(out_dir)
    
    # Show the mask
    solara.Checkbox(label="Show Mask", value=show_mask_widget)
    if show_mask_widget.value:
        show_mask(out_dir)

    # Show the cell plot
    solara.Checkbox(label="Show Cell Plot", value=show_cell_widget)
    if show_cell_widget.value:
        plate_df = pandas.read_csv(out_dir / "plate.csv")
        show_cell_plot(plate_df)
