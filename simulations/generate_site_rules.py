from simulations.load_inputs import load_sites

sites, nSims, script_names = load_sites()


def generate_rule(site, n, script_name="run_sims.py"):
    rule = f"""
rule {site}_run_sim:
    input: 
    output: "COMPS_ID/{site}_COMPS_ID_submit"
    priority: 50
    run:
        shell(get_command("{site}", script="{script_name}", n={n}))

rule {site}_analyzer:
    input: "COMPS_ID/{site}_COMPS_ID_submit"
    output: "COMPS_ID/{site}_COMPS_ID_done"
    run:
        shell(get_command("{site}", script="wait_for_experiment.py"))
    """
    return rule


def run(snakefile='snakefile_bak'):
    with open(snakefile, 'r') as file:
        snakefile_str = file.read()
    snakefile_str = delete_old_rules(sites, nSims, script_names, snakefile_str)
    write_rules(sites, nSims, script_names, snakefile_str, 'snakefile')


def delete_old_rules(sites, nSims, script_names, snakefile_str):
    for site, n, script_name in zip(sites, nSims, script_names):
        if f'rule {site}' in snakefile_str:
            rule = generate_rule(site, n, script_name)
            snakefile_str = snakefile_str.replace(rule, '')
    return snakefile_str


def write_rules(sites, nSims, script_names, snakefile_str, snakefile):
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
