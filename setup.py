from setuptools import setup, find_packages

setup(
    name="sketchy",
    url="https://github.com/esteinig/sketchy",
    author="Eike J. Steinig",
    author_email="eikejoachim.steinig@my.jcu.edu.au",
    packages=find_packages(),
    include_package_data=True,
    version="0.4.0",
    license="MIT",
    description="Real-time lineage calling and genotyping of bacterial"
                "pathogens from uncorrected nanopore reads using MinHash",
)
