# A bioinformatics pipeline based on Ruffus 

Author: Bernie Pope (bjpope@unimelb.edu.au)

Crpipe is based on the Ruffus library for writing bioinformatics pipelines. Its features include:

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

## Usage

You can get a summary of the command line arguments like so:

```
crpipe -h
```

## Configuration file

You must supply a configuration file for the pipeline in YAML format.

Here is an example:

```
```
