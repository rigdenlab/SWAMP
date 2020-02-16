import sys
import logging
import datetime
import collections
from swamp import __version__
from swamp.logger.colorformat import ColorFormatter


class SwampLogger(logging.Logger):
    """Class that implements methods to assist with SWAMP logging messages.

    :param str name: name to identify the logger
    :param bool silent: if True, messages are not logged
    :param list args: arguments passed to :py:obj:`~logging.Logger`
    :param dict kwargs: arguments passed to :py:obj:`~logging.Logger`
    :ivar list msg_record: a list with all the messaged recorded into the log (even if `silent` is set to True)
    """

    def __init__(self, name, silent=False, *args, **kwargs):
        self.silent = silent
        self.msg_record = []
        super(SwampLogger, self).__init__(name, *args, **kwargs)

    def __repr__(self):
        return 'SwampLogger(name=%s, silent=%s)' % (self.name, self.silent)

    # ------------------ Properties ------------------

    @property
    def greeting_msg(self):
        """A greeting message to appear when the logger initiates, informs of the current time and SWAMP version"""

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
        """Header to highlight error messages"""

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
    def record(self):
        """Named tuple to be used as a template for the records stored at \
        :py:attr:`~swamp.logger.swamplogger.SwampLogger.msg_record`"""
        return collections.namedtuple('record', ['msg', 'level'])

    @property
    def msg(self):
        """Property to be used as the format for logged messages"""
        return '%s {}: {}' % self.now

    @property
    def now(self):
        """Current time"""
        now = datetime.datetime.now()
        return now.strftime("%b-%d-%Y %H:%M:%S")

    # ------------------ Methods ------------------

    def error(self, msg, *args, **kwargs):
        """Extend :py:func:`logging.logger.error` to include the \
        :py:attr:`~swamp.logger.swamplogger.SwampLogger.error_header`, check if \
        :py:attr:`~swamp.logger.swamplogger.SwampLogger.silent` and format the message"""

        msg = self.msg.format('ERROR', msg)

        if not self.silent:
            super(SwampLogger, self).error(self.error_header)
            super(SwampLogger, self).error(msg, *args, **kwargs)

        self.msg_record.append(self.record(msg=self.error_header, level=logging.ERROR))
        self.msg_record.append(self.record(msg=msg, level=logging.ERROR))

    def warning(self, msg, *args, **kwargs):
        """Extend :py:func:`logging.logger.warning` to check if \
        :py:attr:`~swamp.logger.swamplogger.SwampLogger.silent` and format the message"""
        msg = self.msg.format('WARNING', msg)

        if not self.silent:
            super(SwampLogger, self).warning(msg, *args, **kwargs)

        self.msg_record.append(self.record(msg=msg, level=logging.WARNING))

    def info(self, msg, *args, **kwargs):
        """Extend :py:func:`logging.logger.info` to check if \
        :py:attr:`~swamp.logger.swamplogger.SwampLogger.silent` and  and format the message"""

        if isinstance(msg, collections.Iterable) and '****' not in msg and '####' not in msg:
            msg = self.msg.format('INFO', msg)

        if not self.silent:
            super(SwampLogger, self).info(msg, *args, **kwargs)

        self.msg_record.append(self.record(msg=msg, level=logging.INFO))

    def debug(self, msg, *args, **kwargs):
        """Extend :py:func:`logging.logger.debug` to check if \
        :py:attr:`~swamp.logger.swamplogger.SwampLogger.silent` and  and format the message"""

        msg = self.msg.format('DEBUG', msg)

        if not self.silent:
            super(SwampLogger, self).debug(msg, *args, **kwargs)

        self.msg_record.append(self.record(msg=msg, level=logging.DEBUG))

    def add_console_handler(self, level='info'):
        """Add the console handler to the :py:obj:`~logging.Logger`"""

        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(self.logging_levels.get(level, logging.INFO))
        ch.setFormatter(ColorFormatter())
        self.addHandler(ch)

    def add_file_handler(self, fname, level='debug'):
        """Add the file handler to the :py:obj:`~logging.Logger`"""

        fh = logging.FileHandler(fname)
        fh.setLevel(self.logging_levels.get(level, logging.INFO))
        fh.setFormatter(logging.Formatter())
        self.addHandler(fh)

    def init(self, console_level="info", logfile_level="debug", logfile=None, use_console=True):
        """Method to initiate the :py:obj:`~swamp.logger.swamplogger.SwampLogger`

        :param str console_level: indicate the console logging level (default: 'info')
        :param str logfile_level: indicate the logfile logging level (default: 'debug')
        :param str logfile: file name to create a log file (default: None)
        :param bool use_console: indicate whether or not to use the console while logging (default: True)
        """

        if use_console:
            self.add_console_handler(console_level)

        if logfile is not None:
            self.add_file_handler(fname=logfile, level=logfile_level)

        self.info(self.greeting_msg)
        self.debug("File logger level: %s", self.logging_levels.get(logfile_level, logging.INFO))
        self.debug("Console logger level: %s", self.logging_levels.get(console_level, logging.INFO))

    def dump_records(self):
        """Log all the records in the :py:attr:`~swamp.logger.swamplogger.SwampLogger.msg_record`"""

        current_list = list(self.msg_record)
        for record in current_list:
            self.log(record.level, record.msg)
