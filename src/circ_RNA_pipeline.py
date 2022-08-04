#!/usr/bin/env python3
import argparse
import logging
import sys


sys.path.append("/src")
from fastq_to_ubam import *
from run_STAR_js import index_samtools, sort_samtools 

################################################################################
# Setup 
################################################################################



# TODO argparse to get all the input
parser = argparse.ArgumentParser(description='Run circular RNAseq pipeline')
# to do maybe make this one argument that sometimes needs two values
parser.add_argument('-r1', '--raw_input', help='Path to raw data (Read 1 for PE data)', required = True)
parser.add_argument('-r2', '--input_read_2', help = 'Path to raw data read 2 file for PE data')
parser.add_argument('--sample', help = 'Sample name', required = True)
# TODO double check that the three options here makes sense -- maybe add a crammed option? 
parser.add_argument('--file_type', help = 'Input file type', choices = ["ubam", "fastq", "bam"], required = True) 
parser.add_argument('--tmp_dir', help ='Path to directory for temporary files', required = True)
parser.add_argument('--read_type', help = 'Paired end or Single end reads', choices = ['PE', 'SE'], required = True)
parser.add_argument('--stranded', help = 'Include for strand specific libraries', default = "NONE", choices = ["NONE", "FIRST_READ_TRANSCRIPTION_STRAND", "SECOND_READ_TRANSCRIPTION_STRAND"])
parser.add_argument('--cohort', help = 'Full RNAseq cohort name', required = True )
parser.add_argument('--tissue', help = 'Sample tissue type', choices = ['brain', 'blood_pax', 'ipsc', 'plasma', 'csf'])
parser.add_argument('--STAR_index', help = "Path to the directory with the STAR genome index refernce", required = True)
parser.add_argument('--out_dir', help = 'directory for output', required = True)
parser.add_argument('--out_struct', help = "Directory structure of output", default = 'none', choices = ['none', 'hydra'])

args = parser.parse_args()


if args.read_type == 'PE' and args.file_type == 'fastq' and args.input_read_2 is None:
    parser.error("--input_read_2 required if input is fastq and read type is PE")

if args.out_struct == 'hydra' and args.tissue is None: 
    parser.error("--tissue required if output structure is hydra")


# Maybe make this an error later 
if args.raw_input == args.input_read_2: 
    parser.error('--input_read_2 and --raw_input cannot be the same file. Please check input.')

# TODO step to move from our server to storage1 (perhaps make this flexible?)


# TODO step to move from storage1 to scratch? 


# set up logs
log_path = args.cohort + '_circular_pipeline.log'
# Create logger 
circ_logs = logging.getLogger('circular_RNA_pipeline')
circ_logs.setLevel(logging.DEBUG)

# Handler 1: file for all info+ logs
log_file = logging.FileHandler(log_path)
file_format = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")
log_file.setLevel(logging.INFO)
log_file.setFormatter(file_format)

# Handler 2: stream for all warning +
stream = logging.StreamHandler()
streamformat = logging.Formatter("%(levelname)s:%(module)s:%(message)s")
stream.setLevel(logging.WARNING)
stream.setFormatter(streamformat)

# Add handlers to logs 
circ_logs.addHandler(log_file)
circ_logs.addHandler(stream)

# Pipeline Steps 
################################################################################

def create_out_dir(dir_to_create):
    already_exists = os.path.exists(dir_to_create)
    if not already_exists: 
        os.makedirs(dir_to_create, exist_ok = True) #note this will make any dirs missing on the path

def setup_output_dirs(output_struct, out_dir, cohort_name, tissue, sample_id):
     if output_struct == "hydra": 
          if tissue == "brain":
               tissue_folder = "01-Brain"
          elif tissue == "blood_pax": 
               tissue_folder = "02-Blood_PAXgene"
          elif tissue == "ipsc":
               tissue_folder = "03-iPSC"
          elif tissue == "plasma":
               tissue_folder = "04-Plasma"
          else: 
               tissue_folder = "05-CSF"
          circ_logs.info(f'Setting up output structure to comply with Hydras dir structure and placing in {out_dir}')
          # 02-Processed/01-Brain/03-Circular/${SAMPLEID} -.bam, .bai, SJ.out.tab, chimeric.out.juction
          circ_proccessed_dir = os.path.join(out_dir, '02-Processed/02-GRCh38', tissue_folder, '03-Circular', cohort_name, sample_id)
          print(f'Circ_processed dir: {circ_proccessed_dir}')
          create_out_dir(circ_proccessed_dir)
          #  03-AnalysisReady/01-Brain/03-Circular/${SAMPLEID} -- DCC ouoput -- CircRNACount, CircCoordinates, Linear_Count, Circ Skip junctions
          circ_analysis_dir = os.path.join(out_dir, '03-AnalysisReady/01-Brain/03-Circular', cohort_name)
          print(f'Circ analaysis dir: {circ_analysis_dir}')
          create_out_dir(circ_analysis_dir)
          output_dirs = { "circ_bams": circ_proccessed_dir,  "DCC_out": circ_analysis_dir}
     else: 
          print(f'Placing all output in dir:{out_dir}')
          output_dirs = {"circ_bams": out_dir, "DCC_out": out_dir }
     return output_dirs

def align_STAR_chimeric(input_1, input_2, bam_out_dir, sample_name, file_type, read_type, star_genome, tmp):
     if file_type == 'fastq': 
          if read_type == 'PE':
               star_chimeric_fastq_PE(input_1, input_2, sample_name, bam_out_dir, star_genome)
          else: 
               star_chimeric_fastq_SE(input_1, sample_name, bam_out_dir, star_genome)
     elif file_type == "bam":
          if read_type == 'PE':
               #TODO star_chimeric_bam_PE
               revert_bam(input_1, bam_out_dir) #TODO not sure if this needs to be different to separate out reads 
               star_chimeric_bam_PE(input_1, bam_out_dir)
          else:
               #TODO test star_chimeric_bam_SE
               #get ubam
               ubam_out = os.path.join(bam_out_dir, f"{sample_name}_unmapped.bam")
               bam_to_ubam(input_1, ubam_out, tmp , by_readgroup = 'false')
               # align with star
               star_chimeric_bam_SE(ubam_out, bam_out_dir, sample_name, star_genome)
     else: 
          if read_type == 'PE':
               #TODO star_chimeric_bam_PE
               star_chimeric_bam_PE(input_1, bam_out_dir)
          else:
               #TODO test star_chimeric_bam_SE
               star_chimeric_bam_SE(input_1, bam_out_dir, sample_name, star_genome)

def star_chimeric_fastq_PE(input_1, input_2, sample_name, bam_out_dir, genome_dir, thread_N = 12):
     # run for both mates toghter
     output_prefix = os.path.join(bam_out_dir, f'{sample_name}_Unified.')
     cmd_both_mates = (
          f'/opt/STAR-2.7.8a/bin/Linux_x86_64/STAR --runThreadN {thread_N}'
          f' --genomeDir {genome_dir}'
          f' --outSAMtype BAM SortedByCoordinate'
          f' --readFilesIn {input_1} {input_2}'
          f' --readFilesCommand zcat' # not sure if this works if they aren't zipped
          f' --outFileNamePrefix {output_prefix}'
          f' --outReadsUnmapped Fastx'
          f' --outSJfilterOverhangMin 15 15 15 15'
          f' --alignSJoverhangMin 15'
          f' --alignSJDBoverhangMin 15'
          f' --outFilterMultimapNmax 20'
          f' --outFilterScoreMin 1'
          f' --outFilterMatchNmin 1'
          f' --outFilterMismatchNmax 2'
          f' --chimSegmentMin 15'
          f' --chimScoreMin 15'
          f' --chimScoreSeparation 10'
          f' --chimJunctionOverhangMin 15'
     )
     print ('Command:' + cmd_both_mates) 
     logger.info(cmd_both_mates)
     cmd_to_call_both = cmd_both_mates.split()
     subprocess.check_call(cmd_to_call_both)

     output_prefix_R1 = os.path.join(bam_out_dir, f'{sample_name}_R1.')
     cmd_mate1 = (
          f'/opt/STAR-2.7.8a/bin/Linux_x86_64/STAR --runThreadN {thread_N}'
          f' --genomeDir {genome_dir}'
          f' --outSAMtype None'
          f' --readFilesIn {input_1}'
          f' --readFilesCommand zcat' #TODO not sure if this works if they aren't zipped
          f' --outFileNamePrefix {output_prefix_R1}'
          f' --outReadsUnmapped Fastx'
          f' --outSJfilterOverhangMin 15 15 15 15'
          f' --alignSJoverhangMin 15'
          f' --alignSJDBoverhangMin 15'
          f' --seedSearchStartLmax 30'
          f' --outFilterMultimapNmax 20'
          f' --outFilterScoreMin 1'
          f' --outFilterMatchNmin 1'
          f' --outFilterMismatchNmax 2'
          f' --chimSegmentMin 15'
          f' --chimScoreMin 15'
          f' --chimScoreSeparation 10'
          f' --chimJunctionOverhangMin 15'
     )
     print ('Command:' + cmd_mate1) 
     logger.info(cmd_mate1)
     cmd_to_call_mate1 = cmd_mate1.split()
     subprocess.check_call(cmd_to_call_mate1)

     output_prefix_R2 = os.path.join(bam_out_dir, f'{sample_name}_R2.')
     cmd_mate2 = (
          f'/opt/STAR-2.7.8a/bin/Linux_x86_64/STAR --runThreadN {thread_N}'
          f' --genomeDir {genome_dir}'
          f' --outSAMtype None'
          f' --readFilesIn {input_2}'
          f' --readFilesCommand zcat' #TODO not sure if this works if they aren't zipped
          f' --outFileNamePrefix {output_prefix_R2}'
          f' --outReadsUnmapped Fastx'
          f' --outSJfilterOverhangMin 15 15 15 15'
          f' --alignSJoverhangMin 15'
          f' --alignSJDBoverhangMin 15'
          f' --seedSearchStartLmax 30'
          f' --outFilterMultimapNmax 20'
          f' --outFilterScoreMin 1'
          f' --outFilterMatchNmin 1'
          f' --outFilterMismatchNmax 2'
          f' --chimSegmentMin 15'
          f' --chimScoreMin 15'
          f' --chimScoreSeparation 10'
          f' --chimJunctionOverhangMin 15'
     )
     print ('Command:' + cmd_mate2) 
     logger.info(cmd_mate2)
     cmd_to_call_mate2 = cmd_mate2.split()
     subprocess.check_call(cmd_to_call_mate2)

def star_chimeric_fastq_SE(input_1, sample_name, bam_out_dir, genome_dir, thread_N = 12):
     output_prefix = os.path.join(bam_out_dir, f'{sample_name}_R1.')
     cmd = (
          f'/opt/STAR-2.7.8a/bin/Linux_x86_64/STAR --runThreadN {thread_N}'
          f' --genomeDir {genome_dir}'
          f' --outSAMtype BAM SortedByCoordinate'
          f' --readFilesIn {input_1}'
          f' --readFilesCommand zcat' #TODO not sure if this works if they aren't zipped
          f' --outFileNamePrefix {output_prefix}'
          f' --outReadsUnmapped Fastx'
          f' --outSJfilterOverhangMin 15 15 15 15'
          f' --alignSJoverhangMin 15'
          f' --alignSJDBoverhangMin 15'
          f' --outFilterMultimapNmax 20'
          f' --outFilterScoreMin 1'
          f' --outFilterMatchNmin 1'
          f' --outFilterMismatchNmax 2'
          f' --chimSegmentMin 15'
          f' --chimScoreMin 15'
          f' --chimScoreSeparation 10'
          f' --chimJunctionOverhangMin 15'
     )
     print ('Command:' + cmd) 
     logger.info(cmd)
     cmd_to_call = cmd.split()
     subprocess.check_call(cmd_to_call)

def index(sample_name, out_dir):
     out_prefix = os.path.join(out_dir, sample_name)
     if args.read_type == 'PE': 
          sorted_bam = f"{out_prefix}_Unified.Aligned.sortedByCoord.out.bam"
     else: 
          sorted_bam = f"{out_prefix}.Aligned.sortedByCoord.out.bam"
     index_samtools(sorted_bam)

def revert_bam(input_aligned_bam, out_dir):
     #TODO write this 
     print("Still need to write this")

def star_chimeric_bam_PE(input_ubam, out_dir):
     #TODO I need to write this 
     print('Still need to write this')

def star_chimeric_bam_SE(input_ubam, out_dir, sample_name, genome_dir, \
     read_cmd = 'samtools view -h', thread_N = 12):
     output_prefix = os.path.join(out_dir, f'{sample_name}.')
     STAR_file_input = "SAM SE"
     cmd = (
          f'/opt/STAR-2.7.8a/bin/Linux_x86_64/STAR --runThreadN {thread_N}'
          f' --genomeDir {genome_dir}'
          f' --outSAMtype BAM SortedByCoordinate'
          f' --readFilesIn {input_ubam}'
          f' --readFilesType {STAR_file_input}'
          f' --readFilesCommand {read_cmd}'
          f' --outFileNamePrefix {output_prefix}'
          f' --outReadsUnmapped Fastx'
          f' --outSJfilterOverhangMin 15 15 15 15'
          f' --alignSJoverhangMin 15'
          f' --alignSJDBoverhangMin 15'
          f' --outFilterMultimapNmax 20'
          f' --outFilterScoreMin 1'
          f' --outFilterMatchNmin 1'
          f' --outFilterMismatchNmax 2'
          f' --chimSegmentMin 15'
          f' --chimScoreMin 15'
          f' --chimScoreSeparation 10'
          f' --chimJunctionOverhangMin 15'
     )
     print ('Command:' + cmd) 
     logger.info(cmd)
     cmd_to_call = cmd.split()
     subprocess.check_call(cmd_to_call)

# create folders : 
out_dirs = setup_output_dirs(args.out_struct, args.out_dir, args.cohort, args.tissue, args.sample)

# Align chimerically with star (Sort & index with samtools & convert to Ubam if needed)
align_STAR_chimeric(args.raw_input, args.input_read_2, out_dirs['circ_bams'], args.sample, args.file_type, args.read_type, args.STAR_index, args.tmp_dir)

# currently sorting with star -- but might take up less memory to sort with samtools 
# index bam with samtools 
index(args.sample, out_dirs['circ_bams'])

# TODO Post Alignment Picard QC ( Collect RNAseq metrics, Collect Alignemt summary metrics, Mark dups) -- can this be done here? 

# TODO run mulitqc??

# TODO quantify with DCC -- perhaps leave this for a seperate image since it needs all samples