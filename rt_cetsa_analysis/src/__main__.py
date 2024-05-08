"""CLI for rt-cetsa-moltprot-tool."""

import json
import logging
import os
import pathlib

import filepattern
import typer
from run_rscript import run_rscript

# get env
POLUS_LOG = os.environ.get("POLUS_LOG", logging.INFO)
POLUS_TAB_EXT = os.environ.get("POLUS_TAB_EXT", ".csv")

# Initialize the logger
logging.basicConfig(
    format="%(asctime)s - %(name)-8s - %(levelname)-8s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)
logger = logging.getLogger("rt_cetsa_analysis")
logger.setLevel(POLUS_LOG)

app = typer.Typer()

@app.command()
def main(
    inp_dir: pathlib.Path = typer.Option(
        ...,
        "--inpDir",
        help="Input directory containing the all data files.",
        exists=True,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    params_pattern: str = typer.Option(
        ".+",
        "--params",
        help="Match the molten fit params csv files in the input directory.",
    ),
    values_pattern: str = typer.Option(
        ".+",
        "--values",
        help="Match the baseline corrected values csv files in the input directory.",
    ),
    platemap_pattern: str = typer.Option(
        ".+",
        "--platemap",
        help="Match the platemap xlst files in the input directory.",
    ),
    preview: bool = typer.Option(
        False,
        "--preview",
        help="Preview the files that will be processed.",
    ),
    out_dir: pathlib.Path = typer.Option(
        ...,
        "--outDir",
        help="Output directory to save the results.",
        exists=True,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
) -> None:
    """CLI for rt-cetsa-moltprot-tool."""
    # TODO: Add to docs that input csv file should be sorted by `Temperature` column.
    logger.info("Starting the CLI for rt-cetsa-moltprot-tool.")

    logger.info(f"Input directory: {inp_dir}")
    logger.info(f"params_pattern: {params_pattern}")
    logger.info(f"values_pattern: {values_pattern}")
    logger.info(f"platemap_pattern: {platemap_pattern}")
    logger.info(f"Output directory: {out_dir}")

    fp_params = filepattern.FilePattern(inp_dir, params_pattern)
    fp_values = filepattern.FilePattern(inp_dir, values_pattern)
    fp_platemap = filepattern.FilePattern(inp_dir, platemap_pattern)
    
    params_files = [f for f in fp_params()]
    values_files = [f for f in fp_values()]
    platemap_files = [f for f in fp_platemap()]

    if preview:
        vals = fp_params.get_unique_values(fp_params.get_variables()[0])[fp_params.get_variables()[0]]
        out_json = {"files": [f"test_merged_molten_final_{exp}.csv" for exp in vals]}
        with (out_dir / "preview.json").open("w") as f:
            json.dump(out_json, f, indent=2)

    for params, values, platemap in zip(params_files, values_files, platemap_files):
        # TODO replace with exceptions
        if len(params[1]) != 1 or len(values[1]) != 1 or len(platemap[1]) != 1:
            raise Exception("patterns are incorrect. Each pattern is expected to match files uniquely.")

        run_rscript(params[1][0], values[1][0], platemap[1][0], out_dir)


if __name__ == "__main__":
    app()