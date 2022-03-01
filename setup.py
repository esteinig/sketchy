from setuptools import setup, find_packages

setup(
    name="sketchy-utils",
    url="https://github.com/esteinig/sketchy",
    author="Eike Steinig",
    author_email="eike.steinig@unimelb.edu.au",
    packages=find_packages(),
    include_package_data=True,
    version="0.6.0",
    license="MIT",
    entry_points="""
    [console_scripts]
    sketchy-utils=sketchy_utils.terminal:app
    """,
    install_requires=[
        'typer-cli',
        'pyyaml',
        'pandas',
        'scikit-learn',
        'ijson',
        'mkdocs-markdownextradata-plugin',
        'mkdocs',
        'mkdocs-material',
        'mkdocs-material-extensions',
        'markdown-include'

    ],
    description="Lineage calling and genotyping of bacterial"
                "pathogens from uncorrected nanopore reads using "
                "genomic neighbor typing and MinHash",
)
