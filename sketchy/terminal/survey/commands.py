import click
import json

from pathfinder.pipelines.data import SurveyData
from pathlib import Path


@click.command()
@click.option(
    '--directory', '-d', default=None, required=True, type=Path,
    help='Input directory containing summarized and filtered output'
         ' from pf-core/pf-survey for creating a sketch.'
)
@click.option(
    '--output', '-o', default='sketchy.survey.tsv', type=Path,
    help='Output CSV file with dataframe for Sketchy.'
)
@click.option(
    '--template', '-t', default=None, type=str,
    help='Use a preconfigured configuration template: mrsa, kleb'
)
@click.option(
    '--config', '-c', default=None, type=Path,
    help='Use a configuration file to extract specific columns for lineage, '
         'genotype and susceptibility from results of molecular typing '
         'with pf-core/pf-survey.'
)
def survey(directory, output, template, config):
    """ Create survey data frame for creating a sketch with MASH """

    s = SurveyData()
    s.read(results=directory)

    if template == "mrsa":
        # Default MRSA config
        cfg = dict(
            lineage=dict(
                mlst="sequence_type"
            ),
            genotype=dict(
                mykrobe_genotype=list()
            ),
            susceptibility=dict(
                mykrobe_phenotype=list()
            ),
        )
    elif template == "kleb":
        # Default Kleb config
        cfg = dict(
            lineage=dict(mlst='sequence_type'),
            genotype=dict(
                kleborate=[
                    'K_locus', 'O_locus',
                ]
            )
        )
    elif template == "tb":
        # Default Kleb config
        cfg = dict(
            lineage=dict(
                mykrobe_lineage=list()
            ),
            susceptibility=dict(
                mykrobe_phenotype=list()
            ),
            genotype=None,
        )
    else:
        # Default MLST config
        cfg = dict(
            lineage=dict(
                mlst="sequence_type"
            )
        )

    if config:
        if Path(config).exists():
            with Path(config).open('r') as json_cfg:
                cfg = json.load(json_cfg)


    df = s.sketchy(config=cfg)
    df.to_csv(output, sep='\t')