#!/usr/bin/env python3
import argparse
import pathlib
import shutil

# emod_malaria
import emod_malaria.bootstrap as dtk

from simulations import manifest as manifest


def get_eradication(use_local_eradication, my_manifest=manifest):
    # Get Eradication and schema file
    msg = list()
    if not use_local_eradication:
        era_folder = pathlib.Path(my_manifest.eradication_path).parent
        shutil.copytree(era_folder, era_folder.parent / "Old_Eradication", dirs_exist_ok=True)
        msg.append(f"Copying old Eradication in local folder to : {era_folder.parent / 'Old_Eradication'}.\n")
        dtk.setup(era_folder)
        msg.append(f"Copying Eradication from emod_malaria to: {era_folder}.\n")

    else:
        msg.append(f"Using Eradication in local folder: {my_manifest.eradication_path}.\n")

    if my_manifest.eradication_path.is_file():
        with open(manifest.eradication_found, 'w') as file:
            for message in msg:
                file.write(message)
                print(message)
    else:
        raise FileNotFoundError(f'{my_manifest.eradication_path} is not found.')


if __name__ == "__main__":
    # TBD: user should be allowed to specify (override default) erad_path and input_path from command line 
    # plan = EradicationBambooBuilds.MALARIA_LINUX
    # print("Retrieving Eradication and schema.json from Bamboo...")
    # get_model_files( plan, manifest )
    # print("...done.")

    parser = argparse.ArgumentParser(description='Get Eradicaiton from emod_malaria')
    parser.add_argument('--use_local_eradication', '-u', action='store_true', help='using Eradication in local folder')

    args = parser.parse_args()
    get_eradication(use_local_eradication=args.use_local_eradication)
