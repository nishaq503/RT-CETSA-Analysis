import math
import pathlib

import filepattern
import matplotlib.pyplot as plt
import numpy
import pandas
import solara
import tifffile
from scipy.optimize import curve_fit

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


def model(
    x: float,
    kn: float,
    bn: float,
    ku: float,
    bu: float,
    dhm: float,
    tm: float,
) -> float:
    R = 8.31446261815324
    recip_t_diff = (1 / tm) - (1 / x)
    exp_term = numpy.exp((dhm / R) * recip_t_diff)
    # exp_inner = (x - tm * numpy.log(0.01/0.99)) / (t_onset - tm)
    # exp_term = numpy.exp(exp_inner)
    numerator = kn * x + bn + (ku * x + bu) * exp_term
    denominator = 1 + exp_term
    return numerator / denominator


def eq_two_state(
    temperatures: numpy.ndarray,
    intensities: numpy.ndarray,
    p0: numpy.ndarray,
) -> dict[str, float]:
    popt, _ = curve_fit(model, temperatures, intensities, method="trf")
    params = {
        "kn": float(popt[0]),
        "bn": float(popt[1]),
        "ku": float(popt[2]),
        "bu": float(popt[3]),
        "dhm": float(popt[4]),
        "tm": float(popt[5]),
    }
    return params


def show_cell_plot(plate_df: pandas.DataFrame, moltprot_df: pandas.DataFrame):
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

    moltprot_row = moltprot_df[moltprot_df["ID"] == cell_name]
    moltprot_row = moltprot_row.drop(columns=["ID"])

    param_names = [
        "kN_init", "bN_init", "kU_init", "bU_init", "dHm_init", "Tm_init",
        "kN_fit", "bN_fit", "kU_fit", "bU_fit", "dHm_fit", "Tm_fit",
        "S", "BS_factor", "T_onset", "dCp_component", "dG_std",
    ]
    params = {p: moltprot_row[p].values[0] for p in param_names}
    temperatures = plate_df["Temperature"].values + 273.15
    intensities = plate_df[cell_name].values
    # p0 = [
    #     params["kN_init"], params["bN_init"], params["kU_init"], params["bU_init"], params["dHm_init"], params["Tm_init"]
    # ]
    # params = eq_two_state(temperatures, intensities, numpy.asarray(p0))
    solara.Markdown(f"params: {params}")
    curve = [
        model(t, params["kN_fit"], params["bN_fit"], params["kU_fit"], params["bU_fit"], params["dHm_fit"], params["Tm_fit"])
        for t in temperatures
    ]
    # curve = [
    #     model(t, params["kn"], params["bn"], params["ku"], params["bu"], params["dhm"], params["tm"])
    #     for t in temperatures
    # ]

    # Get column from plate.csv
    intensities = plate_df[cell_name].values
    fig, ax = plt.subplots()
    ax.scatter(temperatures, intensities)
    ax.plot(temperatures, curve, color="red")
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
        baseline_fit=15,
        baseline_bounds=3,
        dCp=0,
        onset_threshold=0.01,
        savgol=10,
        trim_max=3,
        trim_min=3,
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
        moltprot_df = pandas.read_csv(out_dir / "plate_moltprot.csv")
        show_cell_plot(plate_df, moltprot_df)
