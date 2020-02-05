import threading
import logging


class ThreadResults(object):
    """Class to hold the results from a multi-threaded process

       Implements semaphore methods to regulate thread I/O into the result list

       :param :obj:`swamp.logger.swamplogger` logger: logger instance to record log messages
       :ivar :obj:`threading.lock` lock: lock to control I/O to the result instance
       :ivar :obj:`pandas.DataFrame` value: dataframe with the results of the grid search
       """

    def __init__(self, logger=None):
        self.lock = threading.Lock()
        self.value = []
        if logger is not None:
            self.logger = logger
        else:
            self.logger = logging.getLogger(__name__)

    def register(self, new_results):
        self.logger.debug('Waiting for lock')
        self.lock.acquire()
        self.logger.debug('Acquired lock')
        self.value.append(new_results)
        self.lock.release()
        self.logger.debug('Released lock')
