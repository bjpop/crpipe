'''
A Bioinformatics pipeline based on Ruffus.
Author: Bernie Pope (bjpope@unimelb.edu.au).
'''

# Remember to do export DRMAA_LIBRARY_PATH=/usr/local/slurm_drmaa/1.0.7-gcc/lib/libdrmaa.so

from ruffus import *
import ruffus.cmdline as cmdline
import drmaa
from version import version
import sys
from config import Config
from state import State
from logger import Logger
from pipeline import make_pipeline

# default place to save cluster job scripts
# (mostly useful for post-mortem debugging)
DEFAULT_JOBSCRIPT_DIR = 'jobscripts'
# default name of the pipeline configuration file
DEFAULT_CONFIG_FILE = 'pipeline.config'


def parse_command_line():
    '''Parse the command line arguments of the pipeline'''
    parser = cmdline.get_argparse(description='Colorectal cancer pipeline',
        ignored_args = ["version"] )
    parser.add_argument('--config', type=str, default=DEFAULT_CONFIG_FILE,
        help='Pipeline configuration file in YAML format, defaults to {}'.format(DEFAULT_CONFIG_FILE))
    parser.add_argument('--jobscripts', type=str, default=DEFAULT_JOBSCRIPT_DIR,
        help='Directory to store cluster job scripts created by the pipeline, defaults to {}'.format(DEFAULT_JOBSCRIPT_DIR))
    parser.add_argument('--version', action='version',
        version='%(prog)s ' + version)
    return parser.parse_args()


def main():
    '''Initialise the pipeline, then run it'''
    # Parse command line arguments
    options = parse_command_line()
    # Initialise the logger
    logger = Logger(__name__, options.log_file, options.verbose)
    # Log the command line used to run the pipeline
    logger.info(' '.join(sys.argv))
    # Set up the DRMAA session for running cluster jobs
    drmaa_session = drmaa.Session()
    drmaa_session.initialize()
    # Parse the configuration file, and initialise global state
    config = Config(options.config)
    state = State(options=options, config=config, logger=logger,
                  drmaa_session=drmaa_session)
    # Build the pipeline workflow
    pipeline = make_pipeline(config, state)
    # Run (or print) the pipeline
    cmdline.run(options)
    # Shut down the DRMAA session
    drmaa_session.exit()


if __name__ == '__main__':
    main()
