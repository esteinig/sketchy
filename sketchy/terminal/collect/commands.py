import click
import pandas

from pathlib import Path


@click.command()
@click.option(
    '--directory', '-d', required=True, type=Path, help='Path to output directory from Nextflow'
)
@click.option(
    '--workflow', '-w', required=False, type=str, default="comparison", help='Nextflow subworkflow, one of: comparison'
)
@click.option(
    '--outdir', '-o', required=False, type=Path, default="nxf-results", help='Nextflow summary output directory'
)
@click.option(
    '--id_last', '-i', required=False, is_flag=True, help='Add identifier column last into the database header (saureus db)'
)
def collect(
    directory, workflow, outdir, id_last
):

    """ Collect predictions and summarize results from Nextflow """

    if workflow == "comparison":

        outdir.mkdir(parents=True, exist_ok=True)

        comparison_data = {}
        for path in (directory / 'dist', directory / 'screen', directory / 'screen_winner', directory / 'stream'):
            database_paths = [p for p in path.glob("*") if p.is_dir()]

            database_data = []
            for db_path in database_paths:
                db_header = pandas.read_csv(
                    db_path / 'header.txt', sep="\t", header=None, index_col=None
                ).iloc[0].tolist()

                if id_last:
                    db_header = db_header + ['id', 'replicate']  # use for saureus
                else:
                    db_header = ['id'] + db_header + ['replicate']  # quick fix

                read_limit_paths = [p for p in db_path.glob("*") if p.is_dir()]

                read_limit_data = []
                for read_limit_path in read_limit_paths:
                    result_files = read_limit_path.glob("*.tsv")

                    result_data = []

                    for file in result_files:

                        try:
                            df = pandas.read_csv(file, sep="\t", header=None)
                            name = file.name.strip(".tsv").split("_")
                            df.index = ["_".join(name[:-1]) for _ in df.iterrows()]
                            df['replicate'] = [name[-1] for _ in df.iterrows()]

                            print(df)

                            if path.name == "stream":
                                df.columns = ["read"] + db_header
                            elif path.name == "dist":
                                df.columns = ["rank", "distance", "shared_hashes"] + db_header
                            elif path.name == "screen" or path.name == "screen_winner":
                                df.columns = ["rank", "identity", "shared_hashes"] + db_header
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
            v.to_csv(f"{outdir / k}.tsv", sep="\t", index=True, header=True)

