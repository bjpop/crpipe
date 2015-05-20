'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import *
from stages import make_stage, fastqc


def make_pipeline(config, runner):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name="crpipe")
    # Get a list of paths to all the FASTQ files
    fastq_paths = config.get_option('fastqs')

    # Run fastQC on the FASTQ files
    pipeline.transform(task_func=make_stage(runner, fastqc),
        input=fastq_paths,
        filter=suffix('.fastq.gz'),
        output='_fastqc')

    return pipeline
