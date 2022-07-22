import setuptools
import sys
import pathlib


def take_package_name(name):
    if not name.startswith("-"):
        return name.strip()


def load_requires_from_file(filepath):
    with open(filepath) as fp:
        return [take_package_name(pkg_name) for pkg_name in fp.readlines() if take_package_name(pkg_name)]


def load_arguments_from_file(filepath, arguments_list):
    with open(filepath) as fp:
        for pkg_name in fp.readlines():
            if pkg_name.startswith("-"):
                arguments_list.extend(pkg_name.strip().split(' '))


cur_path = pathlib.Path().parent.resolve()             
requirement_filepath = cur_path / "requirements.txt"          
requirements = load_requires_from_file(requirement_filepath)
arguments = []
load_arguments_from_file(requirement_filepath, arguments)

develop_install = 'develop' in sys.argv
if develop_install:
    sys.argv.extend(arguments)
    requirements.reverse()
    
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

authors = [
    ('Monique Ambrose', 'Monique.Ambrose@gatesfoundation.org'),
    ('Ye Chen', 'ye.chen@gatesfoundation.org')
]

setuptools.setup(
    name="malaria_model_validation",
    version="0.1.0",
    author=[author[0] for author in authors],
    author_email=[author[1] for author in authors],
    description="malaria-model_validation package",
    install_requires=requirements,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/InstituteforDiseaseModeling/malaria-model_validation",
    project_urls={
        "Bug Tracker": "https://github.com/InstituteforDiseaseModeling/malaria-model_validation/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"
    ],
    package_data={'simulation_inputs': ['*/*.csv', '*.csv'],
                  'simulations': ['snakefile_bak']},
    include_package_data=True,
    exclude_package_data={'': ['tests']},
    package_dir={},
    packages=setuptools.find_packages(exclude=['*tests*']),
    setup_requires=['wheel'],
    python_requires=">=3.9",
)