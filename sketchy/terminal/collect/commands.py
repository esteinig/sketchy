import click
import pandas

from pathlib import Path


@click.command()
@click.option(
    '--directory', '-d', required=True, type=Path, help='Path to output directory from Nextflow'
)
@click.option(
    '--workflow', '-w', required=False, type=str, default="comparison", help='Nextflow subworkflow, one of: sketchy'
)
def collect(
    directory, workflow
):

    """ Collect predictions and summarize results from Nextflow """

    if workflow == "comparison":

        comparison_data = []
        for path in (directory / 'stream', directory / 'dist', directory / 'screen'):
            read_limit_paths = path.glob("*/")

            read_limit_data = []
            for read_limit_path in read_limit_paths:
                result_files = read_limit_path.glob("*.tsv")
                print(read_limit_path, result_files)
                result_data = []
                for file in result_files:
                    result_data.append(
                        pandas.read_csv(file, sep="\t")
                    )
                results = pandas.concat(result_data)
                results['read_limit'] = [read_limit_path.name for _ in results.iterrows()]
                read_limit_data.append(results)

            read_limits = pandas.concat(read_limit_data)
            read_limits['mode'] = [path.name for _ in read_limits.iterrows()]
            comparison_data.append(read_limits)

        comparison = pandas.concat(read_limits)

        print(comparison)
