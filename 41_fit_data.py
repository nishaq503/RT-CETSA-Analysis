import warnings

import pandas as pd
from core import MoltenProtFit
from itertools import product
import concurrent.futures
from scipy.stats import gmean

# Silence pandas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

filename = ".data/output/plate.csv"

trim_min = range(0, 1, 3)
trim_max = range(0, 1, 3)
sav_gol = range(0, 10, 5)
baseline_fit = range(3 ,9, 3)
baseline_bounds = range(3, 9, 3)

settings_grid_search = product(baseline_bounds, baseline_fit, sav_gol, trim_max,trim_min)

def run_model(settings):
    baseline_bounds, baseline_fit, sav_gol, trim_max,trim_min = settings
    fit = MoltenProtFit(filename, input_type="csv", debug=True)
    fit.SetAnalysisOptions(
        model="santoro1988",
        baseline_fit=baseline_fit,
        baseline_bounds=baseline_bounds,
        onset_threshold=0.01,
        savgol=sav_gol,
        blanks=[],
        exclude=[],
        invert=False,
        mfilt=None,
        shrink=None,
        trim_min=trim_max,
        trim_max=trim_min,
    )
    fit.PrepareData()
    fit.ProcessData()
    out = fit.plate_results
    return out


def evaluate_model(settings):
    print(settings)
    out = run_model(settings)
    bs_factor= gmean(out['BS_factor'].values)
    return bs_factor

max_bs_factor = 0
best_settings = None
exceptions = []


if __name__ == '__main__':


    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
        # Start the load operations and mark each future with its URL
        future_to_settings = {executor.submit(evaluate_model, settings): settings for settings in settings_grid_search}
        for future in concurrent.futures.as_completed(future_to_settings):
            settings = future_to_settings[future]
            try:
                bs_factor = future.result()
                print(f"mean BS_factor for settings {settings}: {bs_factor}")
                if(bs_factor > max_bs_factor):
                    best_settings = settings
                    max_bs_factor = bs_factor
            except Exception as exc:
                exceptions.append((settings, exc))

    print(f"best settings: {best_settings}, bs_factor: {max_bs_factor}")

    out = run_model(best_settings)

    # sort output by row/col
    out['ind'] = out.index 
    out['letter'] = out.apply(lambda row: row.ind[:1], axis=1)
    out['number'] = out.apply(lambda row: row.ind[1:], axis=1).astype(int)
    out.drop(columns=['ind'])
    out = out.sort_values(['letter','number'])

    # save output
    output_path = '.data/output/molten_prot_out.csv'
    out.to_csv(output_path)
    
    print(f"results saved to {output_path}")