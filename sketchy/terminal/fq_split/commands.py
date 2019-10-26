import click
import pysam

from pathlib import Path


@click.command()
@click.option(
    "--fpath",
    "-f",
    type=Path,
    help="Path to Fast{a,q} file to split.",
    required=True,
)
@click.option(
    "--prefix",
    "-p",
    type=str,
    help="Prefix for split output in format: {prefix}.{split}.fq",
    required=True,
)
@click.option(
    "--size",
    "-s",
    default=10,
    help="Maximum size to split reads into files in GB [10]",
    type=int,
)
def fq_split(fpath: Path, prefix: str, size: int):
    """ Split Fast{a,q} into separate files with approx. maximum size """

    chunk = 1
    giga = 1 << 30
    with pysam.FastxFile(fpath) as fastq_in:
        output = Path(f'{prefix}.{chunk}.fq')

        fh = output.open('w')
        for i, read in enumerate(fastq_in):
            strout = str(read)
            if not strout.endswith('\n'):
                strout += '\n'
            fh.write(strout)
            # Check every 100 reads if output size exceeds maximum size
            if i % 10000 == 0:
                fsize = output.stat().st_size / giga
                if fsize > size:
                    print(f'Reads ({fsize}GB) written to: {output}')
                    fh.close()
                    chunk += 1
                    output = Path(f'{prefix}.{chunk}.fq')
                    fh = output.open('w')
                else:
                    pass
            if i % 10000 == 0:
                print(i)
    fh.close()
