from idmtools.core.platform_factory import Platform
from idmtools.core import ItemType

from COMPS.Data import Experiment
from COMPS.Data.WorkItem import RelationType

import simulations.manifest as manifest

import pandas as pd


def get_depends_on_ids(suite_ids: str) -> str:
    """
    Get IDs for all depends on entities in Comps for a given suite. Save them in a csv file organized by site names.
    Args:
        suite_ids (): Comps ID for a given suite.

    Returns: name of the csv file.

    """
    platform = Platform(manifest.platform_name)
    suite = platform.get_item(item_id=suite_ids, item_type=ItemType.SUITE)
    site_to_ids = {"site": [],
                   "exp_id": [],
                   "analyzer_wi_id": [],
                   "filter_wi_id": []}
    for exp in suite.experiments:
        exp_id = exp.id
        site_name = exp.name[len("validation_"):]
        site_to_ids["site"].append(site_name)
        site_to_ids["exp_id"].append(exp_id)

        e = Experiment.get(exp_id)
        analyzer_wi = e.get_parent_related_workitems(relation_type=RelationType.DependsOn)
        site_to_ids["analyzer_wi_id"].append(analyzer_wi[0].id)

        filter_wi = analyzer_wi[0].get_parent_related_workitems(relation_type=RelationType.DependsOn)
        site_to_ids["filter_wi_id"].append(filter_wi[0].id)

    df = pd.DataFrame(site_to_ids)
    df.to_csv(f"suite_{suite_ids}.csv")
    return f"suite_{suite_ids}.csv"


if __name__ == '__main__':
    with open(manifest.suite_id_file, 'r') as f:
        suite_ids = f.readline()
    csv_filename = get_depends_on_ids(suite_ids)
