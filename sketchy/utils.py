import logging


class PoreLogger:

    def __init__(self):

        self.logger = logging.getLogger(self.__class__.__name__)
        logging.basicConfig(filename=f"{self.__class__.__name__}.log")