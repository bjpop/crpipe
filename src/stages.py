'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to
the state of the pipeline, such as config, options, DRMAA and the logger.
'''

from utils import safe_make_dir
from runner import run_stage

def make_stage(state, function):
    '''Convenience wrapper to build an arity 2 function from a stage function.
    This just makes a closure which is closed over the state paramater.
    Ruffus stages should be functions with 2 arguments, but our functions have
    3 arguments due to the need to pass the state as a parameter.
    This avoids making the state a global variable.
    '''
    closure = lambda input, output: function(state, input, output)
    # Hack to make the closure have the same name as the input function.
    # Ruffus uses the func_name property to identify stages.
    closure.func_name = function.func_name
    return closure 

def fastqc(state, fastq_in, dir_out):
    '''Quality check fastq file using fastqc'''
    safe_make_dir(dir_out)
    command = "fastqc --quiet -o {dir} {fastq}".format(dir=dir_out, fastq=fastq_in)
    run_stage(state, 'fastqc', command)

def index_reference_bwa(state, reference_in, index_file_out):
    '''Index the reference genome using BWA'''
    command = "bwa index -a bwtsw {ref}".format(ref=reference_in)
    run_stage(state, 'index_reference_bwa', command)

def index_reference_samtools(state, reference_in, index_file_out):
    '''Index the reference genome using samtools'''
    command = "samtools faidx {ref}".format(ref=reference_in)
    run_stage(state, 'index_reference_samtools', command)

'''
Could write like so by passing state in extra parameters
def index_reference_samtools(reference_in, index_file_out, state):
    command = "samtools faidx {ref}".format(ref=reference_in)
    runner.run_stage('index_reference_samtools', command)
'''
