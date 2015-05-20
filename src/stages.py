'''
Individual stages of the pipeline implemented as functions from
input files to output files.
'''

from utils import safe_make_dir

def make_stage(runner, function):
    '''Convenience wrapper to build an arity 2 function from a stage function.
    This just makes a closure which is closed over the runner paramater.
    Ruffus stages should be functions with 2 arguments, but our functions have
    3 arguments due to the need to pass the runner as a parameter.
    This avoids making the runner a global variable.
    The runner knows everything about submitting jobs and has full access to
    the state of the pipeline, such as config, options, DRMAA and the logger.
    '''
    closure = lambda input, output: function(runner, input, output)
    # Hack to make the closure have the same name as the input function.
    # Ruffus uses the func_name property to identify stages.
    closure.func_name = function.func_name
    return closure 

def fastqc(runner, fastq_in, dir_out):
    '''Quality check fastq file using fastqc'''
    safe_make_dir(dir_out)
    command = "fastqc --quiet -o {dir} {fastq}".format(dir=dir_out, fastq=fastq_in)
    runner.run_stage('fastqc', command)

def index_reference_bwa(runner, reference_in, index_file_out):
    '''Index the reference genome using BWA'''
    command = "bwa index -a bwtsw {ref}".format(ref=reference_in)
    runner.run_stage('index_reference_bwa', command)

def index_reference_samtools(runner, reference_in, index_file_out):
    '''Index the reference genome using samtools'''
    command = "samtools faidx {ref}".format(ref=reference_in)
    runner.run_stage('index_reference_samtools', command)
