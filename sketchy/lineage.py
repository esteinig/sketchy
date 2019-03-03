import delegator

from io import StringIO
from pathlib import Path

from pathfinder.results import SurveyResult


class LineageMatching:

    def __init__(
            self,
            survey: SurveyResult = None
    ):
        self.survey = survey

    def read_survey(
            self,
            survey_result: Path
    ) -> None:

        self.survey = SurveyResult(
            path=Path()
        )
        self.survey.data.read(
            survey_result
        )

    def select_sequence_types(
            self,
            outdir: Path = Path.home() / 'sketchy' / 'seqs',
            process: str = 'mlst',
            field: str = 'sequence_type',
            sample: int = None,
            values: list = None,
            atleast: int = None,
    ):

        """ Sequence type selection for database sketch with MinHash

        :param process:
        :param field:
        :param outdir:
        :param atleast:
        :param sample:
        :param values:
        :return:
        """

        data = self.survey.data.select(process, field, values=values,
                                       min_count=atleast, sample=sample)

        mlst = data.mlst.sequence_type.sort_index()
        pheno = data.mykrobe_phenotype.sort_index()

        sep = pandas.Series(
            ['_' for _ in data.iid.iid], index=mlst.index
        )

        index = mlst.index.to_series()

        index += sep + mlst.astype(str) + sep
        for column in pheno.columns:
            index += pheno[column]

        data.link_fasta(fdir=outdir, index=index.to_dict(), symlink=True)

    def dist(self, file, mashdb, ncpu=4, top=2):

        result = delegator.run(
            f'mash dist -p {ncpu} {mashdb} {file}'
        )

        df = pandas.read_csv(
            StringIO(result.out), sep='\t', header=None,
            names=[
                "id", 'file', 'dist', "p-value", "shared"
            ], index_col=False
        )

        shared = pandas.DataFrame(
            df.shared.str.split('/').tolist(), columns=['shared', 'total']
        )

        df.shared = shared.shared.astype(int)

        df = df.sort_values(by='shared', ascending=False)

        if top:
            df = df[:top]

        return df
