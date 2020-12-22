import tarfile
import logging

from pathlib import Path
from sketchy.utils import PoreLogger, run_cmd
from colorama import Fore

C = Fore.CYAN
G = Fore.GREEN
Y = Fore.YELLOW
R = Fore.RED
B = Fore.LIGHTBLUE_EX
RE = Fore.RESET


class GoogleCloudSketch:

    def __init__(
        self,
        sketch_path=Path.home() / '.sketchy',
        full: bool = False,
        verbose: bool = True
    ):

        ########################################
        # Public Google Cloud Storage Settings #
        ########################################

        self.bucket_name = 'sketchy-sketch'

        self.sketches = ['kpneumoniae', 'saureus']
        self.full = full

        self.pl = PoreLogger(logging.INFO if verbose else logging.ERROR)
        self.sketch_path = sketch_path

    def pull(self):

        """ Download all sketch archives """

        for name in self.sketches:
            self.download_sketch_archive(archive_name=name)

    def list_sketches(self):

        """ List cached or remote reference sketch collection paths """

        print(f'{"collection":<20}')
        for f in self.sketch_path.glob("*"):
            if f.is_dir():
                file_name = f.name
                print(
                    f"{G}{file_name:<20}{RE}"
                )

    def download_sketch_archive(self, archive_name: str):

        self.pl.logger.info(
            f'Download archive {archive_name} to: {self.sketch_path}'
        )

        self.sketch_path.mkdir(parents=True, exist_ok=True)

        ext = '.tar.gz' if self.full else '.min.tar.gz'
        archive_file = archive_name+ext

        archive_file = self.sketch_path / archive_file

        try:
            self.download_blob(self.bucket_name, archive_file)

            with tarfile.open(self.sketch_path / archive_file) as tar:
                tar.extractall(path=self.sketch_path)

            archive_file.unlink()
        except tarfile.ReadError:
            pass

    def download_blob(self, bucket_name, archive_file):

        run_cmd(
            f'wget -q https://storage.googleapis.com/{bucket_name}/{archive_file.name} '
            f'-O {self.sketch_path / archive_file}'
        )