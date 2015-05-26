'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such 
as config, options, DRMAA and the logger.
'''

from utils import safe_make_dir
from runner import run_stage
import os

class Stages(object):
    def __init__(self, state):
        self.state = state


    def fastqc(self, fastq_in, dir_out):
        '''Quality check fastq file using fastqc'''
        safe_make_dir(dir_out)
        command = "fastqc --quiet -o {dir} {fastq}".format(dir=dir_out, fastq=fastq_in)
        run_stage(self.state, 'fastqc', command)
    

    def index_reference_bwa(self, reference_in, index_file_out):
        '''Index the reference genome using BWA'''
        command = "bwa index -a bwtsw {ref}".format(ref=reference_in)
        run_stage(self.state, 'index_reference_bwa', command)
    

    def index_reference_samtools(self, reference_in, index_file_out):
        '''Index the reference genome using samtools'''
        command = "samtools faidx {ref}".format(ref=reference_in)
        run_stage(self.state, 'index_reference_samtools', command)
    

    def align_bwa(self, inputs, bam_out, sample):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, [fastq_read2_in, reference_in] = inputs
        # Get the read group information for this sample from the configuration file
        read_group = self.state.config.get_read_group(sample)
        # Get the number of cores to request for the job, this translates into the
        # number of threads to give to bwa's -t option
        cores = self.state.config.get_stage_option('align_bwa', 'cores')
        # Run bwa and pipe the output through samtools view to generate a BAM file
        command = 'bwa mem -t {cores} -R "{read_group}" {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -S -b - > {bam}' \
                  .format(cores=cores,
                      read_group=read_group,
                      fastq_read1=fastq_read1_in,
                      fastq_read2=fastq_read2_in,
                      reference=reference_in,
                      bam=bam_out)
        run_stage(self.state, 'align_bwa', command)
 

    def bamtools_stats(self, bam_in, stats_out):
        '''Generate alignment stats with bamtools'''
        command = 'bamtools stats -in {bam} > {stats}' \
                  .format(bam=bam_in, stats=stats_out)
        run_stage(self.state, 'bamtools_stats', command)


    def extract_discordant_alignments(self, bam_in, discordants_bam_out):
        '''Extract the discordant paired-end alignments using samtools'''
        command = 'samtools view -b -F 1294 {input_bam} > {output_bam}' \
                  .format(input_bam=bam_in, output_bam=discordants_bam_out)
        run_stage(self.state, 'extract_discordant_alignments', command)


    def extract_split_read_alignments(self, bam_in, splitters_bam_out):
        '''Extract the split-read alignments using samtools'''
        command = ('samtools view -h {input_bam} | ' \
                   'extractSplitReads_BwaMem -i stdin | ' \
                   'samtools view -Sb - > {output_bam}' 
                   .format(input_bam=bam_in, output_bam=splitters_bam_out))
        run_stage(self.state, 'extract_split_read_alignments', command)

    # Samtools annoyingly takes the prefix of the output bam name as its argument.
    # So we pass this as an extra argument. However Ruffus needs to know the full name
    # of the output bam file, so we pass that as the normal output parameter.
    def sort_bam(self, bam_in, sorted_bam_out, sorted_bam_prefix):
        '''Sort the reads in a bam file using samtools'''
        command = 'samtools sort {input_bam} {output_bam_prefix}' \
                  .format(input_bam=bam_in, output_bam_prefix=sorted_bam_prefix)
        run_stage(self.state, 'sort_bam', command)


    def structural_variants_lumpy(self, inputs, vcf_out):
        '''Call structural variants with lumpy'''
        sample_bam, [splitters_bam, discordants_bam] = inputs
        command = 'lumpyexpress -B {sample_bam} -S {splitters_bam} ' \
                  '-D {discordants_bam} -o {vcf}' \
                  .format(sample_bam=sample_bam, splitters_bam=splitters_bam,
                          discordants_bam=discordants_bam, vcf=vcf_out)
        run_stage(self.state, 'structural_variants_lumpy', command)


    def genotype_svtyper(self, inputs, vcf_out):
        '''Call genotypes on lumpy output using SVTyper'''
        vcf_in, [sample_bam, splitters_bam] = inputs
        command = 'svtyper -B {sample_bam} -S {splitters_bam} ' \
                  '-i {vcf_in} -o {vcf_out}' \
                  .format(sample_bam=sample_bam, splitters_bam=splitters_bam,
                          vcf_in=vcf_in, vcf_out=vcf_out)
        run_stage(self.state, 'genotype_svtyper', command)


    def index_bam(self, bam_in, index_out):
        '''Index a bam file with samtools'''
        command = 'samtools index {bam}'.format(bam=bam_in)
        run_stage(self.state, 'index_bam', command)
