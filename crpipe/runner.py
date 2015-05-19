'''
Run a pipeline state either locally or on the cluster,
taking into account the configuration settings,
and command line options of the pipeline.
'''

# Remember to do export DRMAA_LIBRARY_PATH=/usr/local/slurm_drmaa/1.0.7-gcc/lib/libdrmaa.so

from ruffus.drmaa_wrapper import run_job, error_drmaa_job

# slurm memory is requested in MB, but the config file specifies in GB
MEGABYTES_IN_GIGABYTE = 1024

'''
SLURM options:

--account=name          Charge job to specified accounts
--exclusive             Allocate nodenumber of tasks to invoke on each nodes 
                        in exclusive mode when cpu consumable resource is enabled
--mem=MB                Minimum amount of real memory
--mem-per-cpu=MB        Maximum amount of real memory per allocated cpu required by a job
--mincpus=n             Minimum number of logical processors (threads) per node
--nodes=N               Number of nodes on which to run (N = min[-max])
--ntasks-per-node=n     Number of tasks to invoke on each node
--partition=partition   Partition requested
--reservation=name      Allocate resources from named reservation
--time=hours:minutes    Set a maximum job wallclock time
--ntasks=n              Number of tasks
--mail-type             Notify user by email when certain event types occur. Valid 
                        type values are BEGIN, END, FAIL, REQUEUE, and ALL 
                        (any state change)
'''

class Runner(object):
    def __init__(self, state):
        # State contains:
        #    - command line options
        #    - configuration options
        #    - logger
        #    - DRMAA session
        self.state = state

    def run_stage(self, stage, command):

        # Grab the configuration options for this stage
        config = self.state.config
        modules = config.get_stage_option(stage, 'modules')
        mem = config.get_stage_option(stage, 'mem') * MEGABYTES_IN_GIGABYTE 
        account = config.get_stage_option(stage, 'account')
        queue = config.get_stage_option(stage, 'queue')
        walltime = config.get_stage_option(stage, 'walltime')
        run_local = config.get_stage_option(stage, 'local')
    
        # Generate a "module load" command for each required module
        module_loads = '\n'.join(['module load ' + module for module in modules])
        cluster_command = '\n'.join([module_loads, command])
    
        # Specify job-specific options for SLURM
        job_options = '--time={} --mem={} --partition={} --account={}'.format(walltime, mem, queue, account)

        # Log a message about the job we are about to run
        log_messages = ['Running stage: {}'.format(stage),
                        'Command: {}'.format(command)]
        if not run_local:
            log_messages.append('Job options: {}'.format(job_options))
        self.state.logger.info('\n'.join(log_messages))
    
        # Run the job, capturing stdout and stderr
        stdout_res, stderr_res = None, None
        try:
            stdout_res, stderr_res = \
                run_job(cmd_str=cluster_command,
                    job_name = stage,
                    logger = self.state.logger.proxy,
                    drmaa_session = self.state.drmaa_session,
                    # Determines whether to run the command on the local
                    # machine or run it on the cluster
                    run_locally = run_local,
                    # Keep a copy of the job script for diagnostic purposes
                    retain_job_scripts = True,
                    job_script_directory = self.state.options.jobscripts, 
                    job_other_options = job_options)
        except error_drmaa_job as err:
            raise Exception("\n".join(map(str, ["Failed to run:", command, err, stdout_res, stderr_res])))
