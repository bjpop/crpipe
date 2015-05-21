'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix
from stages import make_stage, fastqc, index_reference_bwa, index_reference_samtools


def make_pipeline(config, runner):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name="crpipe")
    # Get a list of paths to all the FASTQ files
    fastq_files = config.get_option('fastqs')

    # Run fastQC on the FASTQ files
    pipeline.transform(task_func=make_stage(runner, fastqc),
        input=fastq_files,
        filter=suffix('.fastq.gz'),
        output='_fastqc')

    reference_file = config.get_option('reference')

    # Index the reference using BWA 
    pipeline.transform(task_func=make_stage(runner, index_reference_bwa),
        input=reference_file,
        filter=suffix('.fa'),
        output=['.fa.amb', '.fa.ann', '.fa.pac', '.fa.sa', '.fa.bwt'])
    
    # Index the reference using samtools 
    pipeline.transform(task_func=make_stage(runner, index_reference_samtools),
        input=reference_file,
        filter=suffix('.fa'),
        output='.fa.fai')

    return pipeline
