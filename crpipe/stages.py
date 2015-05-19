'''
Individual stages of the pipeline implented as functions from
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
    return lambda input, output: function(runner, input, output)

def fastqc(runner, fastq, outdir):
    '''Quality check fastq file using fastqc'''
    safe_make_dir(outdir)
    command = "fastqc --quiet -o {outdir} {fastq}".format(outdir=outdir, fastq=fastq)
    runner.run_stage('fastqc', command)
