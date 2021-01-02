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

        comparison_data = {}
        for path in (directory / 'stream', directory / 'dist', directory / 'screen'):
            database_paths = path.glob("*/")
            print(path.name)
            database_data = []
            for db_path in database_paths:
                db_header = pandas.read_csv(
                    db_path / 'header.txt', sep="\t", header=None, index_col=None
                ).iloc[0].tolist()

                read_limit_paths = db_path.glob("*/")

                read_limit_data = []
                for read_limit_path in read_limit_paths:
                    result_files = read_limit_path.glob("*.tsv")

                    result_data = []
                    for file in result_files:
                        try:
                            df = pandas.read_csv(file, sep="\t", header=None)
                            df.index = [file.name.strip(".tsv") for _ in df.iterrows()]

                            if path.name == "stream":
                                df.columns = ["read"] + db_header
                            elif path.name == "dist":
                                df.columns = ["rank", "distance", "shared_hashes"] + db_header + ["id"]
                            elif path.name == "screen":
                                df.columns = ["rank", "identity", "shared_hashes"] + db_header + ["id"]
                            else:
                                raise ValueError("Something went seriously wrong, dude! Get your shit together.")

                        except pandas.errors.EmptyDataError:
                            # This can happen to 'screen' if very few reads are used
                            print(f"Could not read results from: {file} - skipping ...")
                            continue

                        result_data.append(df)

                    results = pandas.concat(result_data)
                    results['read_limit'] = [read_limit_path.name for _ in results.iterrows()]
                    read_limit_data.append(results)

                read_limits = pandas.concat(read_limit_data)
                read_limits['db'] = [db_path.name for _ in read_limits.iterrows()]
                database_data.append(read_limits)

            dbs = pandas.concat(database_data)
            dbs['mode'] = [path.name for _ in dbs.iterrows()]
            comparison_data[path.name] = dbs

        for k, v in comparison_data.items():
            print(k)
            print(v)

