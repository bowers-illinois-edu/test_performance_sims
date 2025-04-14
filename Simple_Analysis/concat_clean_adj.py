#!/usr/bin/env python3

# A file to concatenate the many csv files produced by the simulation file.
# Currently must be hand edited, no command line arguments yet

import os
import glob
import pandas as pd
from pathlib import Path
# import concurrent.futures
from IPython.display import clear_output, Image, display

current_path = Path().absolute()
csv_path = current_path / 'Simple_Analysis' / 'CSVS_adj'

all_files = glob.glob(os.path.join(csv_path, "sim_[0-9]*.csv"))
## from https://github.com/softhints/python/blob/master/notebooks/python/Files/How_to_merge_multiple_CSV_files_with_Python.ipynb
all_df = []
for f in all_files:
    df = pd.read_csv(f, sep=',')
    df['file'] = f.split('/')[-1]
    all_df.append(df)

# with concurrent.futures.ProcessPoolExecutor() as executor:
#    all_df = executor.map(pd.read_csv, all_files)

merged_df = pd.concat(all_df, ignore_index=True)
# merged_df
merged_df.describe()

output_dir = current_path / 'Simple_Analysis'  # Combine path objects using /
output_file = output_dir / 'simple_sims_adj_merged_results.csv'

merged_df.to_csv(output_file)
