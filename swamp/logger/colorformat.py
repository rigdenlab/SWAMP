import logging


class ColorFormatter(logging.Formatter):
    """Class to format the logging messages with the appropiate color according to the logging level"""

    @property
    def colors(self):
        """Dictionary with the colors corresponding with each logging level"""
        return {
            logging.DEBUG: 34,  # blue
            logging.WARNING: 33,  # yellow
            logging.ERROR: 31,  # red
            logging.CRITICAL: 31,  # red
        }

    def format(self, record):
        """Format a given record with the appropiate color"""
        if record.levelno in self.colors:
            prefix = "\033[1;{}m".format(self.colors[record.levelno])
            postfix = "\033[0m"
            record.msg = "\n".join([prefix + l + postfix for l in str(record.msg).splitlines()])
        return logging.Formatter.format(self, record)
