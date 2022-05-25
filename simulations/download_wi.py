import argparse
import simulations.params as params
import simulations.manifest as manifest
from simulations.helpers import get_comps_id_filename, get_suite_id

from idmtools.core.platform_factory import Platform
from idmtools_platform_comps.utils.download.download import DownloadWorkItem
from idmtools.core import ItemType


def download_output(site: str, platform: Platform = None) -> bool:
    """
    Download output csv files to output folder from analyzer work item for given site.
    Args:
        site ():
        platform ():

    Returns: status of download work item

    """
    if not platform:
        platform = Platform(manifest.platform_name)
    analyzer_id_file = get_comps_id_filename(site, level=2)
    with open(analyzer_id_file, 'r') as id_file:
        wi_id = id_file.readline()
    wi_id = platform.get_item(wi_id, item_type=ItemType.WORKFLOW_ITEM)
    dl_wi = DownloadWorkItem(
        output_path="output",
        delete_after_download=False,
        extract_after_download=True,
        zip_name=f"{site}.zip",
        file_patterns=[f"{site}/**"]
    )
    dl_wi.related_work_items = [wi_id]
    suite_id = get_suite_id()
    dl_wi.tags['Suite'] = suite_id
    dl_wi.run(wait_until_done=True, platform=platform)

    # Check result
    if not dl_wi.succeeded:
        print(f"Download work item {dl_wi.uid} failed.\n")
    else:
        print(f"Download work item {dl_wi.uid} succeeded.")
        download_id_file = get_comps_id_filename(site, level=3)
        with open(download_id_file, 'w') as id_file:
            id_file.write(dl_wi.uid.hex)

    return dl_wi.succeeded


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process site name')
    parser.add_argument('--site', '-s', type=str, help='site name',
                        default=params.sites[0])  # not sure if we want to make this required argument
    args = parser.parse_args()
    download_output(args.site)
