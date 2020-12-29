import click

from sketchy.survey import SketchySurvey
from pathlib import Path

DEFAULTS = {
    'kpneumoniae': dict(
        merge=dict(
            kleborate=dict(
                rmp=['rmpA', 'rmpA2']
            )
        ),
        config=dict(
            kleborate=[
                'ST',
                'virulence_score',
                'resistance_score',
                'Yersiniabactin',
                'K_locus',
                'O_locus',
                'AGly',
                'Col',
                'Fcyn',
                'Flq',
                'Gly',
                'MLS',
                'Ntmdz',
                'Phe',
                'Rif',
                'Sul',
                'Tet',
                'Tmt',
                'Bla',
                'Bla_Carb',
                'Bla_ESBL',
                'Bla_broad'
            ]
        ),
        binary=dict(
            kleborate=[
                'Yersiniabactin',
                'rmp',
                'AGly',
                'Col',
                'Fcyn',
                'Flq',
                'Gly',
                'MLS',
                'Ntmdz',
                'Phe',
                'Rif',
                'Sul',
                'Tet',
                'Tmt',
                'Bla',
                'Bla_Carb',
                'Bla_ESBL',
                'Bla_broad'
            ]
        )
    ),
    'saureus': dict(
        config=dict(
            sccion=[
                'mlst',
                'meca',
                'pvl',
                'scc'
            ],
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
        merge=dict(),
        binary=dict()
    ),
    'mtuberculosis': dict(
        config=dict(
            mykrobe_phenotype=[
                'Ofloxacin',
                'Moxifloxacin',
                'Isoniazid',
                'Kanamycin',
                'Ethambutol',
                'Streptomycin',
                'Ciprofloxacin',
                'Pyrazinamide',
                'Rifampicin',
                'Amikacin',
                'Capreomycin'
            ],
            mykrobe_lineage=[
                'db_lineage'
            ]
        ),
        merge=dict(),
        binary=dict()
    ),
}


@click.command()
@click.option(
    '--directory', '-d', default=None, required=True, type=Path,
    help='Input directory with collected output from Pathfinder Survey'
)
@click.option(
    '--output', '-o', default='survey.tsv', type=Path,
    help='Tab-delimited genotype feature index for Sketchy'
)
@click.option(
    '--template', '-t', default=None, type=str,
    help='Use a configuration template: saureus, kpneumoniae, mtuberculosis'
)
@click.option(
    '--missing', '-m', default='-', type=str,
    help='Set a missing character [-]'
)
@click.option(
    '--intersect', '-i', is_flag=True,
    help='Take minimum intersection of all specified column values'
)
def construct(directory, output, template, intersect, missing):

    """ Construct genotype create data from Pathfinder Survey """

    survey = SketchySurvey(
        survey_directory=directory
    )

    survey.missing = missing

    if intersect:
        subset = list(
            DEFAULTS[template]['config'].keys()
        )
        survey.survey_data.isolate(*subset)
        survey.survey_data.intersection()

    data = survey.construct(
        config=DEFAULTS[template]['config'],
        binary=DEFAULTS[template]['binary'],
        merge=DEFAULTS[template]['db_merge']
    )

    # All columns should be lower-case
    data.columns = [c.lower() for c in data.columns]

    data.to_csv(output, sep='\t', index_label='uuid')

