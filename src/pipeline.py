'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name="crpipe")
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')

    stages = Stages(state)

    # Run fastQC on the FASTQ files
    pipeline.transform(
        task_func=stages.fastqc,
        name="fastqc",
        input=fastq_files,
        filter=suffix('.fastq.gz'),
        output='_fastqc')

    reference_file = state.config.get_option('reference')

    # Index the reference using BWA 
    pipeline.transform(
        task_func=stages.index_reference_bwa,
        name="index_reference_bwa",
        input=reference_file,
        filter=suffix('.fa'),
        output=['.fa.amb', '.fa.ann', '.fa.pac', '.fa.sa', '.fa.bwt'])
    
    # Index the reference using samtools 
    pipeline.transform(
        task_func=stages.index_reference_samtools,
        name="index_reference_samtools",
        input=reference_file,
        filter=suffix('.fa'),
        output='.fa.fai')

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name="align_bwa",
        input=fastq_files,
        # Match the R1 (read 1) FASTQ file and grab the path and sample name. 
        # This will be the first input to the stage.
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+)_R1.fastq.gz'),
        # Add two more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        #    2. The reference genome file
        add_inputs=add_inputs(['{path[0]}/{sample[0]}_R2.fastq.gz', reference_file]),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='{path[0]}/{sample[0]}.bam') \
        .follows("index_reference_bwa") \
        .follows("index_reference_samtools")

    return pipeline
