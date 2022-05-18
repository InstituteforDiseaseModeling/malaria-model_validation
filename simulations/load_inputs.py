import pandas as pd
import simulations.manifest as manifest


def load_sites():
    skipped_sites = list()

    coord_df = pd.read_csv(manifest.simulation_coordinator_path)
    coord_df = coord_df.set_index('site')
    unfiltered_sites = coord_df.index.tolist()
    for site in unfiltered_sites:
        eir_df = pd.read_csv(manifest.input_files_path / coord_df.at[site, 'EIR_filepath'])
        if site not in eir_df.columns or "?" in site or not coord_df.at[site, 'include_site']:
            skipped_sites.append(site)
    coord_df = coord_df[~coord_df.index.isin(skipped_sites)]

    sites = coord_df.index.tolist()
    nSims = coord_df['nSims'].tolist()
    script_names = coord_df['run_script_name']
    return sites, nSims, script_names
