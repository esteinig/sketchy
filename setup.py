from setuptools import setup, find_packages

setup(
    name="sketchy",
    url="https://github.com/esteinig/sketchy",
    author="Eike J. Steinig",
    author_email="eikejoachim.steinig@my.jcu.edu.au",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "tqdm",
        "colorama",
        "pandas",
        "click",
        "pytest",
        "seaborn",
        "scipy",
        "biopython",
        "pypdf2",
        "deprecated"
    ],
    entry_points="""
        [console_scripts]
        sketchy=sketchy.terminal.client:terminal_client
    """,
    version="0.3",
    license="MIT",
    description="Real-time lineage hashing and genotyping of bacterial"
                "pathogens from uncorrected nanopore reads.",
)
