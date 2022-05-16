
import pandas as pd
import simulations.manifest as manifest
from simulations.helpers import load_coordinator_df


def load_sites():
    skipped_sites = list()

    coord_df = load_coordinator_df(characteristic=False, set_index=True)
    unfiltered_sites = coord_df.index.tolist()
    for site in unfiltered_sites:
        eir_df = pd.read_csv(manifest.input_files_path / coord_df.at[site, 'EIR_filepath'])
        if site not in eir_df.columns or "?" in site:
            skipped_sites.append(site)
    coord_df = coord_df[~coord_df.index.isin(skipped_sites)]

    sites = coord_df.index.tolist()
    nSims = coord_df['nSims'].tolist()
    return sites, nSims