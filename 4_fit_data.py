import warnings

import pandas as pd

from core import MoltenProtFit

# Silence pandas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

filename = ".data/output/plate.csv"

df = pd.read_csv(".data/output/plate.csv", index_col="Temperature", encoding="utf-8")

# IMPORTANT: data must be sorted by temperature in order to get correct readings
df = df.sort_index()

fit = MoltenProtFit(df, input_type="from_xlsx", parent_filename=filename)

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

print(fit.plate_results.sort_values("BS_factor"))
