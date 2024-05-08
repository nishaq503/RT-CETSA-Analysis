import subprocess
from pathlib import Path


def run_rscript(params_filepath : Path, values_filepath: Path, platemap_filepath: Path, out_dir: Path):
    print("run rscript with args: ", params_filepath, values_filepath, platemap_filepath, out_dir)

    cmd = [
        "Rscript",
         "./main.R",
         "--params",
         params_filepath.as_posix(),
         "--values",
         values_filepath.as_posix(),
         "--platemap",
         platemap_filepath.as_posix(),
         "--outdir",
         out_dir.as_posix()
    ]

    subprocess.run(
        args=cmd,
         cwd="src"
        ) 