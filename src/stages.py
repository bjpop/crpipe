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


    def original_reference(self, output):
        '''Original reference file'''
        pass


    def original_fastqs(self, output):
        '''Original fastq files'''
        pass


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
        command = 'samtools faidx {ref}'.format(ref=reference_in)
        run_stage(self.state, 'index_reference_samtools', command)


    def index_reference_bowtie2(self, reference_in, index_file_out, output_prefix):
        '''Index the reference genome using bowtie2'''
        command = 'bowtie2-build {ref} {output_prefix}' \
                      .format(ref=reference_in, output_prefix=output_prefix)
        run_stage(self.state, 'index_reference_bowtie2', command)


    def reference_dictionary_picard(self, reference_in, dict_file_out):
        '''Create a FASTA sequence dictionary for the reference using picard'''
        command = 'java -jar $PICARD_HOME/lib/CreateSequenceDictionary.jar ' \
                  'R={ref} O={dict_file}'.format(ref=reference_in, dict_file=dict_file_out)
        run_stage(self.state, 'reference_dictionary_picard', command)


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


    def extract_genes_bedtools(self, bam_in, bam_out):
        '''Extract MMR genes from the sorted BAM file'''
        bed_file = self.state.config.get_stage_option('extract_genes_bedtools', 'bed') 
        command = 'bedtools intersect -abam {bam_in} -b {bed_file} > {bam_out}' \
                  .format(bam_in=bam_in, bed_file=bed_file, bam_out=bam_out)
        run_stage(self.state, 'extract_genes_bedtools', command)


    def extract_chromosomes_samtools(self, bam_in, bam_out):
        '''Extract selected chomosomes from the bam files'''
        command = 'samtools view -h -b {bam_in} chr2 chr3 chr7 > {bam_out}' \
                  .format(bam_in=bam_in, bam_out=bam_out)
        run_stage(self.state, 'extract_chromosomes_samtools', command)


    #def alignment_coverage_gatk(self, inputs, summary_out, output_prefix):
    #    '''Compute depth of coverage of the alignment with GATK DepthOfCoverage'''
    #    bam_in, [reference_in] = inputs
    #    # Give the Java runtime 1GB less than requested, for the max heap size
    #    mem = int(self.state.config.get_stage_option('alignment_coverage_gatk', 'mem')) - 1 
    #    command = 'java -Xmx{mem}g -jar $GATK_HOME/GenomeAnalysisTK.jar -T DepthOfCoverage -R {ref} -o {output_base} -I {bam}' \
    #              .format(mem=mem, ref=reference_in, output_base=output_prefix, bam=bam_in)
    #    run_stage(self.state, 'alignment_coverage_gatk', command)


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

    def sort_bam_sambamba(self, bam_in, sorted_bam_out):
        '''Sort the reads in a bam file using sambamba'''
        cores = self.state.config.get_stage_option('sort_bam_sambamba', 'cores')
        # Get the tmp directory
        tmp = self.state.config.get_option('tmp') 
        # Get the amount of memory requested for the job
        mem = int(self.state.config.get_stage_option('sort_bam_sambamba', 'mem'))
        mem_limit = max(mem - 4, 1)
        command = 'sambamba sort --nthreads={cores} --memory-limit={mem}GB --tmpdir={tmp} --out={output_bam} {input_bam}' \
                  .format(cores=cores, mem=mem_limit, tmp=tmp, input_bam=bam_in, output_bam=sorted_bam_out)
        run_stage(self.state, 'sort_bam_sambamba', command)


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


    def structural_variants_socrates(self, bam_in, variants_out, sample_dir):
        '''Call structural variants with Socrates'''
        threads = self.state.config.get_stage_option('structural_variants_socrates', 'cores') 
        # jvm_mem is in gb
        jvm_mem = self.state.config.get_stage_option('structural_variants_socrates', 'jvm_mem') 
        bowtie2_ref_dir = self.state.config.get_stage_option('structural_variants_socrates', 'bowtie2_ref_dir') 
        output_dir = os.path.join(sample_dir, 'socrates')
        safe_make_dir(output_dir)
        command = \
        '''
cd {output_dir}
export _JAVA_OPTIONS='-Djava.io.tmpdir={output_dir}'
Socrates all -t {threads} --bowtie2_threads {threads} --bowtie2_db {bowtie2_ref_dir} --jvm_memory {jvm_mem}g {bam}
        '''.format(output_dir=output_dir, threads=threads, bowtie2_ref_dir=bowtie2_ref_dir, jvm_mem=jvm_mem, bam=bam_in)
        run_stage(self.state, 'structural_variants_socrates', command)
