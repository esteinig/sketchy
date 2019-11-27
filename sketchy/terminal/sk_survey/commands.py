import click

from sketchy.sketchy import SketchySurvey
from pathlib import Path

DEFAULTS = {
    'kleb': dict(
        config=dict(
            kleborate=[
                'ST', 'virulence_score', 'resistance_score',
                'Yersiniabactin', 'K_locus', 'O_locus'
            ]
        ),
        binary=dict(
            Yersiniabactin='-'
        )
    ),
    'saureus': dict(
        config=dict(
            # sccion=[
            #     'mlst', 'meca', 'pvl', 'spa', 'scc',
            # ],
            mykrobe_phenotype=[
                'Clindamycin',
                'Rifampicin',
                'Ciprofloxacin',
                'Vancomycin',
                'Tetracycline',
                'Mupirocin',
                'Gentamicin',
                'Trimethoprim',
                'Penicillin',
                'Methicillin',
                'Erythromycin',
                'FusidicAcid'
            ]
        ),
        binary=dict(
        )
    )
}


@click.command()
@click.option(
    '--directory', '-d', default=None, required=True, type=Path,
    help='Input directory containing summarized and filtered output'
         ' from pf-core/pf-survey for creating a sketch.'
)
@click.option(
    '--output', '-o', default='survey.tsv', type=Path,
    help='Output file tab-delimited, genotype data for Sketchy.'
)
@click.option(
    '--template', '-t', default=None, type=str,
    help='Use a configuration template: saureus, kleb'
)
@click.option(
    '--config', '-c', default=None, type=Path,
    help='Use a configuration file to extract specific columns for lineage, '
         'genotype and susceptibility from results of molecular typing '
         'with pf-core/pf-survey.'
)
def sk_survey(directory, output, template, config):
    """ Create survey data frame for creating a sketch with MASH """

    survey = SketchySurvey(
        survey_directory=directory
    )

    data = survey.construct(
        config=DEFAULTS[template]['config'],
        binary=DEFAULTS[template]['binary']
    )

    # All columns should be lower-case
    data.columns = [c.lower() for c in data.columns]

    data.to_csv(output, sep='\t', index_label='uuid')

