from setuptools import setup, find_packages
import codecs
import os.path


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="pydamage",
    version=get_version("pydamage/__init__.py"),
    description="Damage parameter estimation for ancient DNA",
    long_description=open("README.md").read(),
    url="https://github.com/maxibor/pydamage",
    long_description_content_type="text/markdown",
    license="GNU-GPLv3",
    python_requires=">=3.6",
    install_requires=[
        "click",
        "numpy",
        "pandas",
        "pysam",
        "scipy",
        "statsmodels",
        "matplotlib",
        "tqdm",
        "biopython"    
    ],
    packages=find_packages(include=["pydamage"]),
    entry_points={"console_scripts": ["pydamage = pydamage.cli:cli"]},
    include_package_data=True,
    package_data={"": ["models/accuracy_model_v2_python.pickle.gz"]},
)
