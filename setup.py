from setuptools import setup, find_packages

setup(
    name="sketchy",
    url="https://github.com/esteinig/sketchy",
    author="Eike J. Steinig",
    author_email="eikejoachim.steinig@my.jcu.edu.au",
    packages=find_packages(),
    include_package_data=True,
    version="0.4.4",
    license="MIT",
    entry_points="""
    [console_scripts]
    sketchy=sketchy.terminal.client:terminal_client
    """,
    install_requires=[
        'tqdm',
        'colorama',
        'pandas',
        'click',
        'pytest',
        'seaborn',
        'python-dateutil',
        'numpy',
        'pysam',
        'deprecation',
        'watchdog',
        'psutil',
        'pyfastx',
        'braceexpand',
        'networkx',
        'pyfastx',
        'watchdog',
        'psutil',
        'dendropy'
    ],
    description="Real-time lineage calling and genotyping of bacterial"
                "pathogens from uncorrected nanopore reads using "
                "genomic neighbor typing and MinHash",
)
