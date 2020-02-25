import tarfile
import logging

from pathlib import Path
from sketchy.utils import PoreLogger, run_cmd
from colorama import Fore

from google.cloud import storage

C = Fore.CYAN
G = Fore.GREEN
Y = Fore.YELLOW
R = Fore.RED
B = Fore.LIGHTBLUE_EX
RE = Fore.RESET

SPECIES = dict(
    saureus='Staphylococcus aureus',
    kpneumoniae='Klebsiella pneumoniae',
    mtuberculosis='Mycobacterium tuberculosis',
    ecoli='Escherichia coli'
)


class GoogleCloudSketch:

    def __init__(
        self, sketch_path=Path.home() / '.sketchy', verbose: bool = True
    ):

        ########################################
        # Public Google Cloud Storage Settings #
        ########################################

        self.bucket_name = 'sketchy-sketch'

        self.sketches = ['kpneumoniae', 'saureus', 'mtuberculosis']

        self.pl = PoreLogger(logging.INFO if verbose else logging.ERROR)
        self.sketch_path = sketch_path

    def pull(self):

        """ Download all sketch archives """

        for name in self.sketches:
            self.download_sketch_archive(archive_name=name)

        self.pl.logger.info(
            f'Set SKETCHY_PATH={self.sketch_path} or --home {self.sketch_path} '
            f'for access to databases in task: sketchy run'
        )

    def list_sketches(self):

        """ List cached or remote reference sketch collection paths """

        print(f'{"collection":<20}{"k-mer":<10}{"size":<10}{"alias":<30}{"run":<30}')
        for name in self.sketches:
            dir_path = self.sketch_path / f'{name}'

            if not dir_path.exists():
                dir_path = Path(f'gs://{self.bucket_name}/{name}.tar.gz')

            for f in Path(dir_path).glob("*.msh"):
                name, kmer_size, sketch_size = f.name.strip(".msh").split("_")
                alias = f.name.strip('.msh')
                run = f"sketchy run -s {alias}"
                print(
                    f"{G}{name:<20}{RE}{C}{kmer_size:<10}"
                    f"{sketch_size:<10}{RE}{Y}{alias:<30}{RE}{B}{run:<30}{RE}"
                )

    def download_sketch_archive(self, archive_name: str):

        self.pl.logger.info(
            f'Download collection to: {self.sketch_path / archive_name}'
        )

        (self.sketch_path / archive_name).mkdir(parents=True, exist_ok=True)

        archive_file = archive_name+'.tar.gz'
        archive_file = self.sketch_path / archive_name / archive_file

        try:
            self.download_blob(
                self.bucket_name, archive_file, archive_name
            )

            with tarfile.open(self.sketch_path / archive_name / archive_file) as tar:
                tar.extractall(path=self.sketch_path / archive_name)

            archive_file.unlink()
        except tarfile.ReadError:
            pass

        self.pl.logger.info(
            f'Downloads complete, use `sketchy list` to list local sketches'
        )

    def download_blob(self, bucket_name, archive_file, archive_name):

        run_cmd(
            f'wget https://storage.googleapis.com/{bucket_name}/{archive_file.name} '
            f'-O {self.sketch_path / archive_name / archive_file}'
        )

    @staticmethod
    def download_blob_(bucket_name, source_blob_name, destination_file_name):
        """Downloads a blob from the bucket."""

        # bucket_name = "your-bucket-name"
        # source_blob_name = "storage-object-name"
        # destination_file_name = "local/path/to/file"

        storage_client = storage.Client()

        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(source_blob_name)
        blob.download_to_filename(destination_file_name)

        print(
            "Blob {} downloaded to {}.".format(
                source_blob_name, destination_file_name
            )
        )
