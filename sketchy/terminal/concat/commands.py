import click

from pathlib import Path
from PyPDF2 import PdfFileMerger

@click.command()
@click.option(
    '--path', '-p', default=Path.cwd(), type=Path,
    help='Path to PDFs, default: PWD '
)
@click.option(
    '--glob', '-g', default='*.pdf',
    help='Glob to concatenate, default: *.pdf '
)
@click.option(
    '--output', '-o', default=Path.cwd() / 'concatenated.pdf', type=Path,
    help='Output concatenated PDF'
)
def concat(path, glob, output):
    """ Concatenate PDFs """

    pdfs = [str(f) for f in path.glob(glob)]

    merger = PdfFileMerger()

    for pdf in pdfs:
        merger.append(pdf)

    merger.write(
        str(output)
    )
    merger.close()

