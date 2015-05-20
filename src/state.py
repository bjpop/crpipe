'''
Global state of the pipeline:
    - options: the command line arguments of the pipeline program
    - config: the parsed contents of the pipeline configuration file
    - logger: the concurrency friendly logging facility
    - drmaa_session: the DRMAA session for running jobs on the cluster
'''

from collections import namedtuple

State = namedtuple("State", ["options", "config", "logger", "drmaa_session"])
