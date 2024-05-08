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
    return fit


def evaluate_model(settings):
    print(settings)
    out = run_model(settings)
    bs_factor= gmean(out['BS_factor'].values)
    return bs_factor

max_bs_factor = 0
best_settings = None
exceptions = []


if __name__ == '__main__':


    # with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
    #     # Start the load operations and mark each future with its URL
    #     future_to_settings = {executor.submit(evaluate_model, settings): settings for settings in settings_grid_search}
    #     for future in concurrent.futures.as_completed(future_to_settings):
    #         settings = future_to_settings[future]
    #         try:
    #             bs_factor = future.result()
    #             print(f"mean BS_factor for settings {settings}: {bs_factor}")
    #             if(bs_factor > max_bs_factor):
    #                 best_settings = settings
    #                 max_bs_factor = bs_factor
    #         except Exception as exc:
    #             exceptions.append((settings, exc))

    # print(f"best settings: {best_settings}, bs_factor: {max_bs_factor}")

    baseline_bounds, baseline_fit, sav_gol, trim_max,trim_min = 3, 3, 10, 0, 0 
    best_settings = baseline_bounds, baseline_fit, sav_gol, trim_max,trim_min
    fit = run_model(best_settings)
    out_fit_params = fit.plate_results


    # sort output by row/col
    out_fit_params['ind'] = out_fit_params.index 
    out_fit_params['letter'] = out_fit_params.apply(lambda row: row.ind[:1], axis=1)
    out_fit_params['number'] = out_fit_params.apply(lambda row: row.ind[1:], axis=1).astype(int)
    out_fit_params.drop(columns=['ind'])
    out_fit_params = out_fit_params.sort_values(['letter','number'])

    # save output
    output_path_fit_param = '.data/output/molten_prot_fit_params_out.csv'
    out_fit_params.to_csv(output_path_fit_param)

    output_path_baseline_corrected = '.data/output/molten_prot_baseline_corrected_out.csv'
    out_baseline_corrected = fit.plate_raw_corr 
    out_baseline_corrected.to_csv(output_path_baseline_corrected)

    print(f"fit params saved to {output_path_fit_param}")
    print(f"corrected baselines saved to {output_path_baseline_corrected}")