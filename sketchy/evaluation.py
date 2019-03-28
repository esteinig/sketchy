from sketchy.sketchy import Sketchy
from pathlib import Path

import pandas

from tqdm import tqdm

import delegator


class Evaluator:
    """ Base access class of algorithm evaluation for the manuscript """

    def __init__(self, outdir: Path, sketch: str = None):

        self.sketchy = Sketchy()

        self.outdir = outdir
        self.tmpdir = outdir / '.tmp'

    def _identify_stages(self):
        """ First, first correct stable, isolate stable """
        pass

    def bootstrap(self, fastq, nb=1000, reads=1000) -> [Path]:

        print(f'Bootstrapping reads in: {fastq}')

        fname = fastq.stem

        self.sketchy.sort_fastq(
            fastq=fastq, shuffle=False, nbootstrap=nb, replacement=True,
            prefix=str(self.outdir / fname / 'boot_'), sample_size=reads, file=None
        )

        return list((self.outdir / fname).glob('*.fq'))

    def predict_bootstraps(self, bsfiles, sketch, data, cores=4, reads=1000) -> pandas.DataFrame:

        print(f'Predicting scores on bootstrap files.')

        score_data = []
        for fq in tqdm(bsfiles, desc='Bootstrap replicate'):
            replicate_name = fq.stem
            df = self.sketchy.predict_nanopore(
                fastq=fq,
                sketch=sketch,
                data=data,
                cores=cores,
                score=True,
                header=False,
                nreads=reads,
                top=1,
                out=self.outdir / f'{replicate_name}.csv',
                tmp=self.outdir / f'{replicate_name}_tmp',
                sort_by='shared'
            )
            # Bootstrap replicate IDs
            df['bootstrap'] = [replicate_name for _ in df.score]

        return pandas.concat(score_data)




