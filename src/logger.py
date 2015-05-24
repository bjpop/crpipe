'''
Initialisation and use of concurrency-friendly logging facility.
'''

import ruffus.cmdline as cmdline


class Logger(object):
    '''Concurrency friendly logging facility'''
    def __init__(self, prog_name, log_file, verbosity):
        proxy, mutex = cmdline.setup_logging(__name__, log_file, verbosity)
        self.proxy = proxy
        self.mutex = mutex


    def info(self, message):
        '''Display an informational message to the log file'''
        with self.mutex:
            self.proxy.info(message)
