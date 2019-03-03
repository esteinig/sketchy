import delegator
import pandas

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

