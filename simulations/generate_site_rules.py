from simulations.load_inputs import load_sites
from simulations.helpers import get_comps_id_filename

sites, nSims, script_names = load_sites()


def generate_rule(site, n, script_name="run_sims.py"):

    exp_id_file = get_comps_id_filename(site=site)
    analyzer_id_file = get_comps_id_filename(site=site, level=2)
    download_id_file = get_comps_id_filename(site=site, level=3)
    rule = f"""
rule {site}_run_sim:
    input: 
    output: '{exp_id_file}'
    priority: 50
    run:
        shell(get_command(script="{script_name}", site="{site}", n={n}))

rule {site}_analyzer:
    input: '{exp_id_file}'
    output: '{analyzer_id_file}'
    run:
        shell(get_command(script="run_analyzers.py", site="{site}"))

rule {site}_download:
    input: '{analyzer_id_file}'
    output: '{download_id_file}'
    run:
        shell(get_command(script="download_wi.py", site="{site}"))
    """
    return rule


def run(snakefile='snakefile_bak'):
    with open(snakefile, 'r') as file:
        snakefile_str = file.read()
    snakefile_str = delete_old_rules(snakefile_str)
    write_rules(snakefile_str, 'snakefile')


def delete_old_rules(snakefile_str):
    for site, n, script_name in zip(sites, nSims, script_names):
        if f'rule {site}' in snakefile_str:
            rule = generate_rule(site, n, script_name)
            snakefile_str = snakefile_str.replace(rule, '')
    return snakefile_str


def write_rules(snakefile_str, snakefile):
    new_rules = []
    for site, n, script_name in zip(sites, nSims, script_names):
        if f'rule {site}' not in snakefile_str:
            rule = generate_rule(site, n, script_name)
            new_rules.append(rule)
    with open(snakefile, 'w') as file:
        file.write(snakefile_str)
        for rule in new_rules:
            file.write(rule)


if __name__ == '__main__':
    run()
