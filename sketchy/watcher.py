"""
Code adapted from Pomoxis: https://github.com/nanoporetech/pomoxis

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
(c) 2016 Oxford Nanopore Technologies Ltd.
"""

import os

from pathlib import Path
from sketchy.utils import PoreLogger
from watchdog.observers import Observer
from watchdog.events import RegexMatchingEventHandler

EVENT_TYPE_MOVED = "moved"
EVENT_TYPE_DELETED = "deleted"
EVENT_TYPE_CREATED = "created"
EVENT_TYPE_MODIFIED = "modified"


class Watcher(PoreLogger):
    def __init__(self):

        PoreLogger.__init__(self)
        self.watcher = None  # active watcher

    def watch_path(
        self,
        path: Path,
        callback: callable,
        recursive: bool = False,
        wait: bool = True,
        regexes: list = (r".*\.fastq$",),
        **kwargs,
    ):
        """Watch a filepath indefinitely for new files, applying callback to files.

        :param path: path to watch.
        :param callback: callback to apply to newly found files.
        :param recursive: recursively watch path?
        :param wait: wait for file completion
        :param regexes: filter files by applying regexes
        :param kwargs: additional arguments to pass to callback
        """

        self.logger.info(
            f"Initiate event handler for regexes: {' '.join(regexes)}"
        )
        handler = StandardRegexMatchingEventHandler(
            callback=callback, regexes=regexes, wait=wait, **kwargs
        )

        self.watcher = MonitorLizard(
            str(path), event_handler=handler, recursive=recursive
        )

        self.logger.info(f"Directory to watch: {path}")
        self.logger.info(f"Night gathers, and now my watch begins ...")

        self.watcher.start()


# Helper functions

def wait_for_file(fname):
    """Block until a filesize remains constant."""
    size = None
    while True:
        try:
            newsize = os.path.getsize(fname)
        except:
            newsize = None
        else:
            if newsize is not None and size == newsize:
                break

        size = newsize


# EventHandler classes

class StandardRegexMatchingEventHandler(RegexMatchingEventHandler):
    def __init__(
        self,
        callback: callable,
        regexes: list,
        wait: bool = True,
        **kwargs
    ):
        RegexMatchingEventHandler.__init__(self, regexes=regexes)

        self.callback_arguments = kwargs
        self.callback = callback
        self.wait = wait

    def _process_file(self, event):
        """Process an event when a file is created (or moved).
        :param event: watchdog event.
        :returns: result of applying `callback` to watched file.
        """
        if event.event_type == EVENT_TYPE_CREATED:
            fpath = Path(event.src_path)
        else:
            fpath = Path(event.dest_path)

        # need to wait for file to be closed
        if self.wait:
            wait_for_file(fpath)

        return self.callback(fpath, **self.callback_arguments)



class MonitorLizard(object):
    def __init__(self, path, event_handler, recursive=False):
        """Wrapper around common watchdog idiom.
        :param path: path to watch for new files.
        :param event_handler: subclass of watchdog.events.FileSystemEventHandler.
        :param recursive: watch path recursively?
        """
        self.observer = Observer()
        self.observer.schedule(event_handler, path, recursive)

    def start(self):
        """Start observing path."""

        self.observer.start()

    def stop(self):
        """Stop observing path."""
        self.observer.stop()
        self.observer.join()

