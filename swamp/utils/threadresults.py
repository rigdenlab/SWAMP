import threading
import logging


class ThreadResults(object):
    """Class to hold the results from a multi-threaded process

       Implements `threading.semaphore` to regulate thread I/O into a result list

       :param `~swamp.logger.swamplogger.SwampLogger` logger: logger instance to record log messages
       :ivar `threading.lock` lock: lock to control I/O to the result instance
       :ivar `pandas.DataFrame` value: dataframe with the results of the grid search
       """

    def __init__(self, logger=None):
        self.lock = threading.Lock()
        self.value = []
        if logger is not None:
            self.logger = logger
        else:
            self.logger = logging.getLogger(__name__)

    def register(self, new_results):
        """Register a set of new results into :py:attr:`swamp.utils.threadresults.ThreadResults.value`

        :param list new_results: the new results to be registered
        """

        self.logger.debug('Waiting for lock')
        self.lock.acquire()
        self.logger.debug('Acquired lock')
        self.value.append(new_results)
        self.lock.release()
        self.logger.debug('Released lock')
