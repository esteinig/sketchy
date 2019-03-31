import click
import pandas
import shutil

from pathlib import Path


@click.command()
@click.option(
    '--dir', '-d', default=Path().cwd(), type=str, help='Directory containing Bootstrap replicate prediction results.'
)
@click.option(
    '--out', '-o', default='bootstrap_results', help='Bootstrap replicate prediction results.'
)
@click.option(
    '--recursive', '-r', is_flag=True, help='Summarize multiple directories by their names at depth of recursion = 1.'
)
@click.option(
    '--combine', '-c', is_flag=True, help='Summarize bootstrap replicate results from all sub-directories into a file.'
)
def psum(dir, out, recursive, combine):
    """ Predict on a bootstrap replicate (Nextflow) """

    if recursive:

        result_files = []
        pdir = Path(dir).glob('*')
        for d in pdir:
            bs_out = Path(out) / f'{Path(d).name}.tab'
            print(f'Producing bootstrap prediction summary: {bs_out}')
            if d.is_dir():
                results = Path(d).glob('*.tab')
                _cat_results(results, bs_out)
                result_files.append(
                    Path(bs_out)
                )
        if combine:
            all_results = []
            for result in result_files:
                df = pandas.read_csv(result, sep='\t', header=0, index_col=0)
                df['file'] = [result.stem for _ in df.bootstrap]
                all_results.append(df)

            pandas.concat(all_results).to_csv(
                Path(out) / 'combined_data.tab', sep='\t', header=True, index=True
            )
    else:
        results = Path(dir).glob('*.tab')
        _cat_results(results, Path(out) / f'{Path(dir).stem}.tab')


def _cat_results(results, out):

    return pandas.concat(
        [pandas.read_csv(result, sep='\t', header=0, index_col=0) for result in results]
    ).to_csv(out, sep='\t', header=True, index=True)

