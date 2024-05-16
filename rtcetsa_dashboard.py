import warnings
from dataclasses import dataclass, replace
from functools import partial
from itertools import product
from io import BytesIO
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy
import solara
from matplotlib.patches import Circle
from skimage.transform import rotate
from solara.components.file_drop import FileInfo
from solara.lab import Ref
import pandas as pd
from tifffile import imread
from zipfile import ZipFile

from core import MoltenProtFit
from lib import extract_intensity
from lib import get_plate_params
from lib import index_to_battleship
from lib import PlateParams

# Silence pandas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

root_path = Path(__file__).parent


@dataclass(frozen=True)
class State:
    raw_images: dict[str, numpy.ndarray]
    plates: list[str, numpy.ndarray]
    plate_display: bytes | None
    plot_display: bytes | None
    plate_index: str | None
    params: PlateParams | None
    df: pd.DataFrame | None
    t_start: float = 37.0
    t_end: float = 95.0
    upload_progress: float = 0.0
    status: str = "Waiting for data"


def render_plate(
    plates: dict[str, numpy.ndarray], params: PlateParams, index: str
) -> bytes:

    plt.figure()
    plt.imshow(
        plates[index],
        cmap="gray",
        vmin=plates["1.tif"].min(),
        vmax=plates["1.tif"].max(),
    )

    plt.title(index)

    ax = plt.gca()

    for x, y in product(params.X, params.Y):
        circle = Circle(
            (x - params.bbox[2], y - params.bbox[0]),
            radius=params.radii,
            color="red",
            fill=False,
        )
        ax.add_patch(circle)

    bio = BytesIO()
    plt.savefig(bio, format="png")
    bio.seek(0)
    return bio.read()


def render_plot(x: int, y: int, state: State) -> bytes:

    params = state.value.params

    R = min(params.X[1] - params.X[0], params.Y[1] - params.Y[0]) // 2

    X = []
    Y = []
    for index, image in state.value.plates.items():
        X.append(int(index.split(".")[0]))
        Y.append(extract_intensity(image, params.X[x], params.Y[y], R))

    plt.figure()
    plt.scatter(x=X, y=Y)

    plt.title(index_to_battleship(x=x, y=y, size=params.size))

    bio = BytesIO()
    plt.savefig(bio, format="png")
    bio.seek(0)
    return bio.read()


def upload(file: FileInfo, state: solara.Reactive[State]):
    raw_images: dict[str, numpy.ndarray] = {}
    with ZipFile(file["file_obj"]) as zip_file:
        for f in zip_file.namelist():
            state.value = replace(state.value, status=f"Reading: {f}")
            if f.endswith(".tif") or f.endswith(".tiff"):
                with zip_file.open(f, "r") as fr:
                    raw_images[f] = imread(fr)

    raw_images = {
        r: raw_images[r] for r in sorted(raw_images, key=lambda x: int(x.split(".")[0]))
    }

    state.value = replace(state.value, status=f"Getting plate information...")
    params = get_plate_params(raw_images["1.tif"])

    plates: dict[str, numpy.ndarray] = {}
    for fname, img in raw_images.items():
        state.value = replace(state.value, status=f"Extracting plate: {f}")
        plates[fname] = rotate(img, params.rotate, preserve_range=True)[
            params.bbox[0] : params.bbox[1], params.bbox[2] : params.bbox[3]
        ].astype(img.dtype)

    plate_index = "1.tif"
    plate_display = render_plate(plates, params, index=plate_index)

    params = get_plate_params(plates["1.tif"])

    state.value = replace(
        state.value,
        raw_images=raw_images,
        plates=plates,
        status=f"Done!",
        params=params,
        plate_display=plate_display,
        plate_index=plate_index,
    )


def next_image(state: solara.Reactive[State]):
    next_val = False
    for i, im in state.value.plates.items():
        if next_val:
            plate_display = render_plate(state.value.plates, state.value.params, i)
            state.value = replace(
                state.value, plate_display=plate_display, plate_index=i
            )
            break
        if i == state.value.plate_index:
            next_val = True


def prev_image(state: solara.Reactive[State]):
    prev_val = False
    for i, im in reversed(list(state.value.plates.items())):
        if prev_val:
            plate_display = render_plate(state.value.plates, state.value.params, i)
            state.value = replace(
                state.value, plate_display=plate_display, plate_index=i
            )
            break
        if i == state.value.plate_index:
            prev_val = True


def select_well(x: int, y: int, state: solara.Reactive[State]):
    state.value = replace(state.value, plot_display=render_plot(x, y, state))


def analyze_data(state: solara.Reactive[State]):

    temperature_range = [state.value.t_start, state.value.t_end]
    params = state.value.params

    # Some wells seem overexposed, so we try to grab all light in each grid to compensate
    R = min(params.X[1] - params.X[0], params.Y[1] - params.Y[0]) // 2

    # Make an in memory csv
    csv = StringIO()

    # Write the headers
    headers = ["Temperature"]
    for x, y in product(range(len(params.X)), range(len(params.Y))):
        headers.append(index_to_battleship(x, y, params.size))
    csv.write(",".join(headers) + "\n")

    last_index = int(list(state.value.plates)[-1].split(".")[0])

    for name, image in state.value.plates.items():
        index = int(name.split(".")[0])
        temp = temperature_range[0] + (index - 1) / (last_index - 1) * (
            temperature_range[1] - temperature_range[0]
        )
        plate = [f"{temp:.1f}"]

        for x, y in product(range(len(params.X)), range(len(params.Y))):
            plate.append(extract_intensity(image, params.X[x], params.Y[y], R))

        csv.write(",".join([str(p) for p in plate]) + "\n")

    # Back to the beginning
    csv.seek(0)

    # Read in the dataframe
    df = pd.read_csv(csv)

    # Run rtcetsa analysis
    fit = MoltenProtFit(df, input_type="from_xlsx", parent_filename="raw_data.csv")

    fit.SetAnalysisOptions(
        model="santoro1988",
        baseline_fit=3,
        baseline_bounds=3,
        onset_threshold=0.01,
        savgol=10,
        blanks=[],
        exclude=[],
        invert=False,
        mfilt=None,
        shrink=None,
        trim_min=0,
        trim_max=0,
    )

    fit.PrepareData()
    fit.ProcessData()

    result = fit.plate_results.sort_values("BS_factor").reset_index()

    state.value = replace(state.value, df=result)


@solara.component
def Page():

    state = solara.use_reactive(
        State(
            raw_images=[],
            plates=[],
            plate_display=None,
            plate_index=None,
            params=None,
            plot_display=None,
            df=None,
        )
    )
    upload_progress = Ref(state.fields.upload_progress)
    status = Ref(state.fields.status)
    plate_display = Ref(state.fields.plate_display)
    t_start = Ref(state.fields.t_start)
    t_end = Ref(state.fields.t_end)

    with solara.Column(align="center", style={"width": "100%"}):
        with solara.Card(title="Step 1: Upload data"):
            solara.FileDrop(
                label=status.value,
                on_file=partial(upload, state=state),
                lazy=False,
                on_total_progress=upload_progress.set,
            )
            solara.ProgressLinear(upload_progress.value)
            solara.Markdown(status.value)

        with solara.Card(title="Step 2: Preview data"):
            if plate_display.value is not None:
                solara.Image(plate_display.value)

                with solara.Row():
                    solara.Button("Previous", on_click=lambda: prev_image(state))
                    solara.Button("Next", on_click=lambda: next_image(state))

        with solara.Card(title="Step 3: Plot Data"):
            if state.value.params is not None:
                for y in range(len(state.value.params.Y)):
                    with solara.Row(gap="2px", margin="2px", justify="center"):
                        for x in range(len(state.value.params.X)):
                            solara.Button(
                                index_to_battleship(
                                    x=x, y=y, size=state.value.params.size
                                ),
                                outlined=True,
                                style={
                                    "width": "10px",
                                    "margin": "1px",
                                    "padding": "1px",
                                },
                                on_click=lambda x=x, y=y: select_well(x, y, state),
                            )
            if state.value.plot_display is not None:
                solara.Image(state.value.plot_display)

        with solara.Card(title="Step 4: Analyze"):
            if state.value.params is not None:
                with solara.Row(justify="center"):
                    solara.InputFloat("Temperature Start", t_start)
                    solara.InputFloat("Temperature End", t_end)
                solara.Button("Analyze", on_click=lambda: analyze_data(state))
                if state.value.df is not None:
                    solara.DataFrame(state.value.df, scrollable=True)