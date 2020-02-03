import sys
import logging
import datetime
import collections
from swamp import __version__
from swamp.logger.colorformat import ColorFormatter


class SwampLogger(logging.Logger):
    """Class that extends logging.Logger and assist with swamp logging messages.

    :param name: name to identify the logger
    :type name: str
    :param silent
    :type silent: bool
    :param args: arguments passed to logging.Logger
    :type args: list, tuple
    :param kwargs: arguments passed to logging.Logger
    :type kwargs: dict
    """

    def __init__(self, name, silent=False, *args, **kwargs):
        self._silent = silent
        self._msg_record = []
        super(SwampLogger, self).__init__(name, *args, **kwargs)

    def __repr__(self):
        return 'SwampLogger(name=%s, silent=%s)' % (self.name, self.silent)

    # ------------------ Properties ------------------

    @property
    def greeting_msg(self):
        """First message to appear when the logger initiates"""

        return """\n#########################################################################
#########################################################################
#########################################################################
# SWAMP - Solving structures With Alpha helical Membrane Pairs          #
#########################################################################\n\n
SWAMP version %s
Current time: %s\n
""" % (__version__, self.now)

    @property
    def error_header(self):
        """Header for error messages"""

        return """
\033[91m\033[1m**********************************************************************
*******************           SWAMP ERROR          *******************
**********************************************************************
\033[0m
"""

    @property
    def logging_levels(self):
        """A dictionary with the different logger levels that can be set"""

        return {
            "notset": logging.NOTSET,
            "info": logging.INFO,
            "debug": logging.DEBUG,
            "warning": logging.WARNING,
            "error": logging.ERROR,
            "critical": logging.CRITICAL,
        }

    @property
    def silent(self):
        """If set to True messages are not logged"""
        return self._silent

    @silent.setter
    def silent(self, value):
        self._silent = value

    @property
    def msg_record(self):
        """List containing the recorded messages"""
        return self._msg_record

    @msg_record.setter
    def msg_record(self, value):
        self._msg_record = value

    @property
    def record(self):
        """Named tuple to be used as a template for records"""
        return collections.namedtuple('record', ['msg', 'level'])

    @property
    def msg(self):
        """Format a given message before logging it"""
        return '%s {}: {}' % self.now

    @property
    def now(self):
        """Current time"""
        now = datetime.datetime.now()
        return now.strftime("%b-%d-%Y %H:%M:%S")

    # ------------------ Methods ------------------

    def error(self, msg, *args, **kwargs):
        """Extend logging.logger method to include error message header and check if the logger is silenced"""

        msg = self.msg.format('ERROR', msg)

        if not self.silent:
            super(SwampLogger, self).error(self.error_header)
            super(SwampLogger, self).error(msg, *args, **kwargs)

        self.msg_record.append(self.record(msg=self.error_header, level=logging.ERROR))
        self.msg_record.append(self.record(msg=msg, level=logging.ERROR))

    def warning(self, msg, *args, **kwargs):
        """Extend logging.logger method to check if the logger is silenced"""

        msg = self.msg.format('WARNING', msg)

        if not self.silent:
            super(SwampLogger, self).warning(msg, *args, **kwargs)

        self.msg_record.append(self.record(msg=msg, level=logging.WARNING))

    def info(self, msg, *args, **kwargs):
        """Extend logging.logger method to check if the logger is silenced"""

        if isinstance(msg, collections.Iterable) and '****' not in msg and '####' not in msg:
            msg = self.msg.format('INFO', msg)

        if not self.silent:
            super(SwampLogger, self).info(msg, *args, **kwargs)

        self.msg_record.append(self.record(msg=msg, level=logging.INFO))

    def debug(self, msg, *args, **kwargs):
        """Extend logging.logger method to check if the logger is silenced"""

        msg = self.msg.format('DEBUG', msg)

        if not self.silent:
            super(SwampLogger, self).debug(msg, *args, **kwargs)

        self.msg_record.append(self.record(msg=msg, level=logging.DEBUG))

    def add_console_handler(self, level='info'):
        """Add the console handler to the logger"""

        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(self.logging_levels.get(level, logging.INFO))
        ch.setFormatter(ColorFormatter())
        self.addHandler(ch)

    def add_file_handler(self, fname, level='debug'):
        """Add the file handler to the logger"""

        fh = logging.FileHandler(fname)
        fh.setLevel(self.logging_levels.get(level, logging.INFO))
        fh.setFormatter(logging.Formatter())
        self.addHandler(fh)

    def init(self, console_level="info", logfile_level="debug", logfile=None, use_console=True):
        """Method to initiate the logger

        :param console_level: indicate the console logging level (default: 'info')
        :type console_level: str
        :param logfile_level: indicate the logfile logging level (default: 'debug')
        :type logfile_level: str
        :param logfile: file name to create a log file (default: None)
        :type logfile: str, None
        :param use_console: indicate whether or not to use the console while logging (default: True)
        :type use_console: bool
        """

        if use_console:
            self.add_console_handler(console_level)

        if logfile is not None:
            self.add_file_handler(fname=logfile, level=logfile_level)

        self.info(self.greeting_msg)
        self.debug("File logger level: %s", self.logging_levels.get(logfile_level, logging.INFO))
        self.debug("Console logger level: %s", self.logging_levels.get(console_level, logging.INFO))

    def dump_records(self):
        """Log all the records in the message record list"""

        current_list = list(self.msg_record)
        for record in current_list:
            self.log(record.level, record.msg)
