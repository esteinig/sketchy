from setuptools import setup, find_packages

setup(
    name="sketchy",
    url="https://github.com/esteinig/sketchy",
    author="Eike J. Steinig",
    author_email="eikejoachim.steinig@my.jcu.edu.au",
    packages=find_packages(),
    include_package_data=True,
    version="0.5.0",
    license="MIT",
    entry_points="""
    [console_scripts]
    sketchy-utils=sketchy.terminal.client:terminal_client
    """,
    install_requires=[
        'typer-cli',
        'pyyaml',
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
