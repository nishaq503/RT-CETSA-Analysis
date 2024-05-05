import warnings

import pandas as pd

from core import MoltenProtFit

# Silence pandas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

filename = ".data/output/plate.csv"
# filename = "data/cleaned_expt1.csv"
# filename="/Users/antoinegerardin/Documents/data/rt-cetsa/210318_MHR Analyzed/210318_Plate1/cleaned_expt1.csv"

# df = pd.read_csv(filename, index_col="Temperature", encoding="utf-8")

# IMPORTANT: data must be sorted by temperature in order to get correct readings
# df = df.sort_index()

# fit = MoltenProtFit(df, input_type="from_xlsx", parent_filename=filename, denaturant="C")
# fit = MoltenProtFit(df, input_type="csv", parent_filename=filename)
fit = MoltenProtFit(filename, input_type="csv", debug=True)


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

out = fit.plate_results
out['ind'] = out.index 
out['letter'] = out.apply(lambda row: row.ind[:1], axis=1)
out['number'] = out.apply(lambda row: row.ind[1:], axis=1).astype(int)
out.drop(columns=['ind'])
out = out.sort_values(['letter','number'])
out.to_csv('.data/output/molten_prot_out.csv')

from scipy.stats import gmean
g_mean= gmean(out['BS_factor'].values)
print("mean BS_factor: ", g_mean)