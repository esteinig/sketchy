from pathlib import Path
from sketchy.utils import PoreLogger, run_cmd

from colorama import Fore

C = Fore.CYAN
G = Fore.GREEN
Y = Fore.YELLOW
R = Fore.RED
RE = Fore.RESET


class GoogleCloudSketch:

    def __init__(self, sketch_path=Path.home() / '.sketchy' / 'db'):

        ########################################
        # Public Google Cloud Storage Settings #
        ########################################

        self.base_url = 'https://storage.googleapis.com/np-core-sketchy/'

        self.sketch_files = {
            'kleb': 'kleb.default.msh',
            'mrsa': 'mrsa.default.msh',
            'tb': 'tb.default.msh'
        }

        self.pl = PoreLogger()
        self.sketch_path = sketch_path

    def list_dbs(self):

        for sketch in self.sketch_files.keys():
            file = self.sketch_files[sketch]
            file_path = self.sketch_path / file

            if not file_path.exists():
                file_path = self.base_url + file

            if sketch == 'mrsa':
                species = 'Staphylococcus aureus'
            elif sketch == 'kleb':
                species = 'Klebsiella pneumoniae'
            elif sketch == 'tb':
                species = 'Mycobacterium tuberculosis'
            else:
                species = '-'

            color = Y if isinstance(file_path, str) else G
            print(
                f'{C}{species:<35}{RE}',
                f'{Y}{sketch:<10}{RE}',
                f'{color}{str(file_path):<45}{RE}'
            )

    def download(self):

        self.pl.logger.info(f'Initiating  download to: {self.sketch_path}')

        self.sketch_path.mkdir(parents=True, exist_ok=True)

        db_paths = []
        for sketch in self.sketch_files.keys():
            file = self.sketch_files[sketch]
            file_path = self.sketch_path / file

            cmd = f'wget {self.base_url + file} -O {file_path}'
            if not file_path.exists():
                self.pl.logger.info(f'Downloading {file} from {self.base_url}')
                run_cmd(cmd)
            else:
                self.pl.logger.info(f'File exists: {file_path}')

            db_paths.append(
                self.sketch_path / file
            )

        self.pl.logger.info(
            f'Downloads complete, use `sketchy db-list` to view local sketches.'
        )

        return db_paths
