# A bioinformatics pipeline based on [Ruffus](http://www.ruffus.org.uk/)

Author: Bernie Pope (bjpope@unimelb.edu.au)

Crpipe is based on the [Ruffus](http://www.ruffus.org.uk/) library for writing bioinformatics pipelines. Its features include:

 * Job submission on a cluster using DRMAA (currently only tested with SLURM).
 * Job dependency calculation and checkpointing.
 * Pipeline can be displayed as a flowchart.
 * Re-running a pipeline will start from the most up-to-date stage. It will not redo previously completed tasks.

## License

3 Clause BSD License. See LICENSE.txt in source repository.

## Installation

I recommend using a virtual environment:

```
cd /place/to/install
virtualenv crpipe
source crpipe/bin/activate
pip install -U https://github.com/bjpop/crpipe
```

If you don't want to use a virtual environment then you can just install with pip:

```
pip install -U https://github.com/bjpop/crpipe
```

## Worked example

The `example` directory in the source distribution contains a small dataset to illustrate the use of the pipeline.

#### Get a copy of the source distribution:

```
cd /path/to/test/directory
git clone https://github.com/bjpop/crpipe
```

#### Install `crpipe` as described above.

#### Get a reference genome.

```
cd crpipe/example
mkdir reference
# copy your reference into this directory, or make a symbolic link
# call it reference/genome.fa
```

####  Run `crpipe` and ask it what it will do next.

```
crpipe -n --verbose 3
```

#### Run the pipeline.

```
crpipe --use_threads --log_file pipeline.log --jobs 2 --verbose 3
```

## Usage

You can get a summary of the command line arguments like so:

```
crpipe -h
usage: crpipe [-h] [--verbose [VERBOSE]] [-L FILE] [-T JOBNAME] [-j N]
              [--use_threads] [-n] [--touch_files_only] [--recreate_database]
              [--checksum_file_name FILE] [--flowchart FILE]
              [--key_legend_in_graph] [--draw_graph_horizontally]
              [--flowchart_format FORMAT] [--forced_tasks JOBNAME]
              [--config CONFIG] [--jobscripts JOBSCRIPTS] [--version]

Colorectal cancer pipeline

optional arguments:
  -h, --help            show this help message and exit
  --config CONFIG       Pipeline configuration file in YAML format, defaults
                        to pipeline.config
  --jobscripts JOBSCRIPTS
                        Directory to store cluster job scripts created by the
                        pipeline, defaults to jobscripts
  --version             show program's version number and exit

Common options:
  --verbose [VERBOSE], -v [VERBOSE]
                        Print more verbose messages for each additional
                        verbose level.
  -L FILE, --log_file FILE
                        Name and path of log file

pipeline arguments:
  -T JOBNAME, --target_tasks JOBNAME
                        Target task(s) of pipeline.
  -j N, --jobs N        Allow N jobs (commands) to run simultaneously.
  --use_threads         Use multiple threads rather than processes. Needs
                        --jobs N with N > 1
  -n, --just_print      Don't actually run any commands; just print the
                        pipeline.
  --touch_files_only    Don't actually run the pipeline; just 'touch' the
                        output for each task to make them appear up to date.
  --recreate_database   Don't actually run the pipeline; just recreate the
                        checksum database.
  --checksum_file_name FILE
                        Path of the checksum file.
  --flowchart FILE      Don't run any commands; just print pipeline as a
                        flowchart.
  --key_legend_in_graph
                        Print out legend and key for dependency graph.
  --draw_graph_horizontally
                        Draw horizontal dependency graph.
  --flowchart_format FORMAT
                        format of dependency graph file. Can be 'pdf', 'svg',
                        'svgz' (Structured Vector Graphics), 'pdf', 'png'
                        'jpg' (bitmap graphics) etc
  --forced_tasks JOBNAME
                        Task(s) which will be included even if they are up to
                        date.
```

## Configuration file

You must supply a configuration file for the pipeline in YAML format.

Here is an example:

```
# Default settings for the pipeline stages.
# These can be overridden in the stage settings below.

defaults:
    # Number of CPU cores to use for the task
    cores: 1
    # Maximum memory in gigabytes for a cluster job
    mem: 4
    # VLSCI account for quota
    account: VR0002
    queue: main
    # Maximum allowed running time on the cluster in Hours:Minutes
    walltime: '1:00'
    # Load modules for running a command on the cluster.
    modules: 
    # Run on the local machine (where the pipeline is run)
    # instead of on the cluster. False means run on the cluster.
    local: False

# Stage-specific settings. These override the defaults above.
# Each stage must have a unique name. This name will be used in
# the pipeine to find the settings for the stage.

stages:
    # Run quality checks on the FASTQ files using fastQC
    fastqc:
        walltime: '10:00'
        mem: 8
        modules:
            - 'fastqc/0.10.1'

    # Index the hg19 human genome reference with BWA
    index_reference_bwa:
        walltime: '10:00'
        mem: 8
        modules:
            - 'bwa-intel/0.7.12' 
    
    # Index the hg19 human genome reference with samtools
    index_reference_samtools:
        walltime: '10:00'
        mem: 8
        modules:
            - 'samtools-intel/1.1'
    
    # Align paired end FASTQ files to the reference
    align_bwa:
        cores: 8
        walltime: '48:00'
        mem: 32
        modules:
            - 'bwa-intel/0.7.12'
            - 'samtools-intel/1.1'

# The Human Genome in FASTA format

reference: /path/to/reference/genome.fa 

# The input FASTQ files.

fastqs:
   - /path/to/fastqs/sample1_R1.fastq.gz
   - /path/to/fastqs/sample1_R2.fastq.gz
   - /path/to/fastqs/sample2_R1.fastq.gz
   - /path/to/fastqs/sample2_R2.fastq.gz

read_groups:
   'sample1': '@RG\tID:id1\tPU:pu1\tSM:sample1\tPL:ILLUMINA\tLB:lib_sample1'
   'sample2': '@RG\tID:id2\tPU:pu2\tSM:sample2\tPL:ILLUMINA\tLB:lib_sample2'
