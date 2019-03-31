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
        "biopython",
        "google-cloud-storage",
        "click",
        "delegator.py",
        "pytest",
        "pre-commit",
        "scipy"
    ],
    entry_points="""
        [console_scripts]
        sketchy=sketchy.terminal.client:terminal_client
    """,
    version="0.1",
    license="MIT",
    description="Rapid lineage matching and antimicrobial susceptibility "
                "predictions using uncorrected nanopore reads. Uses MinHash "
                "sketches of species assembly databases genotyped with "
                "reproducible analysis pipelines from Pathfinder.",
)
