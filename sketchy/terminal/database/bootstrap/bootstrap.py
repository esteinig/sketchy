import click

from pathlib import Path
from sketchy.sketchy import SketchyDatabase


@click.command()
@click.option(
    '--fasta_directory', '-f', type=Path, required=True,
    help='Path to directory containing reference sketch assembly files (.fasta)'
)
@click.option(
    '--reference_database', '-r', type=Path, required=True,
    help='Reference database sketch directory containing the original database to be bootstrapped '
)
@click.option(
    '--bootstrap_samples', '-b', type=int, required=False, default=1000,
    help='Number of samples to draw from the reference assemblies with replacement (size of database) [1000]'
)
@click.option(
    '--outdir', '-o', type=Path, required=False, default=Path('bootstrap_fasta'),
    help='Output path to symlink assembly files for sketch construction in Mash [bootstrap_fasta]'
)
@click.option(
    '--outdb', '-o', type=Path, required=False, default=Path.cwd(),
    help='Output path to database directory with genotype file sub-set to bootstrap sample [bootstrap_db]'
)
def bootstrap(fasta_directory, reference_database, bootstrap_samples, outdir, outdb):

    """ Bootstrap sample a reference database for Sketchy """

    db = SketchyDatabase(db_path=reference_database)

    db.bootstrap_sample(fasta_dir=fasta_directory, samples=bootstrap_samples, outdir=outdir, outdb=outdb)





