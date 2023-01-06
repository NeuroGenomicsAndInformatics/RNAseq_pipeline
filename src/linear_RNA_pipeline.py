#!/usr/bin/env python3
import argparse
import logging
import sys
import shutil
import re

sys.path.append("/src")
from fastq_to_ubam import *
from run_STAR_js import align_STAR_normal, cram_samtools, index_samtools, sort_samtools, cram_samtools 
from run_Picard_QC import picard_collect_RNA_metrics, picard_collect_alignment_metrics, picard_mark_dups
from run_Salmon import salmon_quant
from run_tin import calc_tin
from run_fastqc import fastqc_fastq_PE, fastqc_SE
################################################################################
# Setup 
################################################################################



# Argparse to get all the input
parser = argparse.ArgumentParser(description='Run linear RNAseq pipeline')
# to do maybe make this one argument that sometimes needs two values
parser.add_argument('-r1', '--raw_input', help='Path to raw data (Read 1 for PE data)', nargs = '+', required= True)
parser.add_argument('-r2', '--input_read_2', help = 'Path to raw data read 2 file for PE data', nargs = '*')
parser.add_argument('--sample', help = 'Sample name', required = True)
# TODO double check that the three options here makes sense -- maybe add a crammed option? 
parser.add_argument('--file_type', help = 'Input file type for -r1 and -r2', choices = ["ubam", "fastq", "bam"], required = True) 
parser.add_argument('--tmp_dir', help ='Path to directory for temporary files', required = True)
parser.add_argument('--read_type', help = 'Paired end or Single end reads', choices = ['PE', 'SE'], required = True)
parser.add_argument('--stranded', help = 'Include for strand specific libraries', default = "NONE", choices = ["NONE", "FIRST_READ_TRANSCRIPTION_STRAND", "SECOND_READ_TRANSCRIPTION_STRAND"])
parser.add_argument('--cohort', help = 'Full RNAseq cohort name', required = True )
parser.add_argument('--tissue', help = 'Sample tissue type', choices = ['brain', 'blood_pax', 'ipsc', 'plasma', 'csf'])
parser.add_argument('--STAR_index', help = "Path to the directory with the STAR genome index refernce", required = True)
parser.add_argument('--ref_flat', help = "Path to reference annotation in flat format", required = True)
parser.add_argument('--annotation', help = "Path to annotaiton file for Salmon", required = True)
parser.add_argument('--annote_bed', help = "Path to annotation in bed format for TIN", required = True)
parser.add_argument('--rib_int', help = "Ribosomal (rRNA) interval list", required = True)
parser.add_argument('--transcripts', help = "Path to Salmon transcripts index", required = True)
parser.add_argument('--adapter_seq', default = 'null', help ="Adapter sequence for picard")
parser.add_argument('--out_dir', help = 'Directory for output', required = True)
parser.add_argument('--out_struct', help = "Directory structure of output", default = 'none', choices = ['none', 'hydra'])
parser.add_argument( '-d', '--delete', action ='store_true',  help = 'Include to delete the following files: .Aligned.out.bam, Aligned.sortedByCoord.out.bam, Aligned.sortedByCoord.out.bam.bai, ._STARgenome, ._STARpass1, _unmapped.bam and Aligned.sortedByCoord.out.md.bam.bai')
parser.add_argument('-c', '--cram', action = 'store_true', help = 'Cram the following bams: Aligned.sortedByCoord.out.md.bam')
parser.add_argument('--ref', help = 'Path to reference genome fasta (for cramming)', required= True) # TODO change this to be required only when cramming
parser.add_argument('--input_to_merge', help = 'Path to secondary input file to concatenate with primary input (ie MSBB which is separated into aligned and unaligend reads)')
parser.add_argument('--merge_file_type', help = 'Input file type for --input_to_merge', choices = ["ubam", "fastq", "bam"], default = "fastq")
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
log_file_name = args.cohort + '_linear_pipeline.log'
log_path = os.path.join(args.out_dir, log_file_name)
# Create logger 
linear_logs = logging.getLogger('linear_RNA_pipeline')
linear_logs.setLevel(logging.DEBUG)

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
linear_logs.addHandler(log_file)
linear_logs.addHandler(stream)


################################################################################
# Pipeline Steps 
################################################################################
def check_refs_exist(star_ref, flat_ref, annot_ref, annot_bed_ref, rib_int_ref, transcript_ref):
    star_exists = os.path.exists(star_ref)
    flat_exists = os.path.exists(flat_ref)
    annot_exists = os.path.exists(annot_ref)
    annot_bed_exists = os.path.exists(annot_bed_ref)
    rib_int_exists = os.path.exists(rib_int_ref)
    transcript_exists = os.path.exists(transcript_ref)
    if not star_exists: 
        sys.exit(f'STAR reference not found: {star_ref}, exiting')
    if not flat_exists: 
        sys.exit(f'flat reference not found: {flat_ref}, exiting')
    if not annot_exists: 
        sys.exit(f'Annotation not found: {annot_ref}, exiting')
    if not annot_bed_exists: 
        sys.exit(f'Bed format annotation not found: {annot_bed_ref}, exiting')
    if not rib_int_exists: 
        sys.exit(f'Ribosomal interval list not found: {rib_int_ref}, exiting')
    if not transcript_exists: 
        sys.exit(f'Salmon transcripts index not found: {transcript_ref}, exiting')

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
        linear_logs.info(f'Setting up output structure to comply with Hydras dir structure and placing in {out_dir}')
        #If outstructure = hydra create paths for :
        # 02-Processed/01-Brain/01-FastQC/${SAMPLEID} -- fastqc.htlm, fastqc.zip, picard qc
        fastqc_out_dir = os.path.join(out_dir, '02-Processed/02-GRCh38/', tissue_folder, '01-FastQC', cohort_name, sample_id)
        create_out_dir(fastqc_out_dir)
        # 02-Processed/01-Brain/02-Linear_TIN/${SAMPLEID} -- .bam, .bai, .tin.csv
        linear_tin_processed_dir = os.path.join(out_dir, '02-Processed/02-GRCh38/', tissue_folder, '02-Linear_TIN', cohort_name, sample_id)
        create_out_dir(linear_tin_processed_dir)
        # 03-AnalysisReady/01-Brain/01-Linear/${SAMPLEID} -- quant
        quant_dir = os.path.join(out_dir, '03-AnalysisReady/02-GRCh38', tissue_folder, '01-Linear', cohort_name, sample_id)
        create_out_dir(quant_dir)
        # 03-AnalysisReady/01-Brain/02-TIN/${SAMPLEID} - tin summary
        tin_analysis_dir = os.path.join(out_dir, '03-AnalysisReady/02-GRCh38', tissue_folder, '02-TIN', cohort_name, sample_id)
        create_out_dir(tin_analysis_dir)
        # 03-AnalysisReady/01-Brain/04-MultiQC/${POOLID}-- multiquc
        multiqc_dir = os.path.join(out_dir, '03-AnalysisReady/02-GRCh38', tissue_folder, '04-MultiQC', cohort_name)
        create_out_dir(multiqc_dir)
        output_dirs = {"fastqc": fastqc_out_dir, "linear_tin": linear_tin_processed_dir, "quant": quant_dir, "tin_summary": tin_analysis_dir, "multiqc": multiqc_dir}
    else: 
        print(f'Placing all output in dir: {out_dir}')
        output_dirs = {"fastqc": out_dir, "linear_tin": out_dir, "quant": out_dir, "tin_summary": out_dir, "multiqc":out_dir }
    return output_dirs

def convert_to_ubam(out_dir, sample_name, file_type, read_type, raw_input, input_read_2, tmp_dir, rg_name  = 'A'):
    ubam_out = os.path.join(out_dir, f"{sample_name}_unmapped.bam")
    if file_type == 'fastq': 
        print('Not converting to ubam, just using fastqs for STAR')
        #print('Converting fastq to ubam')
        #if read_type == 'PE':
        #    print('Using ' + raw_input + ' and ' + input_read_2 +' as input for ummaped bam')
        #    fastq_to_ubam_PE(raw_input, input_read_2, ubam_out, sample_name, rg_name)
        #    star_input = ubam_out
        #else: 
        #    print('Converting single fq to unmapped bam')
        #    fastq_to_ubam_SE(raw_input, ubam_out, sample_name, rg_name)
        #star_input = ubam_out

    elif args.file_type == 'bam':  
        bam_to_ubam(raw_input, ubam_out, tmp_dir, by_readgroup = 'false'),
        star_input = ubam_out
    else:
        print('Input data already ubam, proceeding to alignment')
        star_input = raw_input
    return star_input

def align_with_star(star_input, sample_name, read_type, STAR_index, out_dir):
    out_prefix = os.path.join(out_dir, f'{sample_name}.')
   # if read_type == 'bam':
   #     print('Aligning ubams with star')
    # TODO Align with star (Sort & index with samtools) -- for reverted bams since the code here is slighlty different 
   # else:   
    print('Aligning ubam with STAR')
    STAR_file_input = "SAM "+ read_type
    align_STAR_normal( star_input, out_prefix, STAR_file_input, STAR_index)
    aligned_bam_out =  f"{out_prefix}Aligned.out.bam"
    aligned_transcript_bam = f"{out_prefix}Aligned.toTranscriptome.out.bam"
    return (aligned_bam_out, aligned_transcript_bam)

def sort(out_dir, sample_name, aligned_bam_in):
    print('Sorting with samtools')
    out_prefix = os.path.join(out_dir, sample_name)
    sort_samtools(aligned_bam_in, out_prefix)
    sorted_bam_out = f"{out_prefix}.Aligned.sortedByCoord.out.bam"
    return sorted_bam_out

def index(sorted_bam_out):
    print('Indexing with samools')
    index_samtools(sorted_bam_out)
    indexed_bam_out = f"{sorted_bam_out}.bai"
    return indexed_bam_out

def fastqc(out_dir, sample_name, input_file_type, input_read_type, raw_input, raw_input2):
    if input_file_type == 'fastq' and input_read_type == 'PE': 
        fastqc_fastq_PE(raw_input, raw_input2, sample_name, out_dir)
    else: 
        fastqc_SE(raw_input, sample_name, out_dir)

def delete_extras(delete, sample_name, aligned_bam, sorted_bam, indexed_bam, STAR_dir, indexed_md_bam, merged_ubam, input_2):
    if delete is True: 
        print("Deleting extra files")
        os.remove(aligned_bam)
        os.remove(sorted_bam)
        os.remove(indexed_bam)
        os.remove(indexed_md_bam)
        unmapped_bam = os.path.join(STAR_dir, f"{sample_name}_unmapped.bam")
        os.remove(unmapped_bam)
        # remove the input to merge if it exists 
        if input_2 is not None:
            os.remove(merged_ubam)
            unmapped_bam_input2 = os.path.join(STAR_dir, f"{sample_name}_input2_unmapped.bam")
            os.remove(unmapped_bam_input2)
        STAR_pass1 = os.path.join(STAR_dir, f"{sample_name}._STARpass1")
        STAR_genome = os.path.join(STAR_dir,f"{sample_name}._STARgenome" )
        shutil.rmtree(STAR_pass1)
        shutil.rmtree(STAR_genome)
        linear_logs.info(f'Deleting files: {aligned_bam}, {sorted_bam}, {indexed_bam}, {indexed_md_bam}, {unmapped_bam}, {STAR_pass1}, and {STAR_genome}')

def cram_bams(cram, mark_dups_bam, ref, sample_name, out_dir ):
    if cram is True: 
        print(f"Cramming {mark_dups_bam_out}")
        out_prefix = os.path.join(out_dir, sample_name)
        md_cram_out = f"{out_prefix}.Aligned.sortedByCoord.out.md.cram"
        cram_samtools(mark_dups_bam, md_cram_out, ref )
        linear_logs.info(f"Cramming: {mark_dups_bam}")
        #then delete the orginal bam 
        linear_logs.info(f"Deleting {mark_dups_bam}")
        os.remove(mark_dups_bam)

def get_MSBB_read_group(sample_name):
    #MSBB read group ID is just the sample name UNLESS it has a third _ in name -- then it is everything before this underscore
    regex = re.compile('^[^_]+_[^_]+_[^_]+')
    read_group_id = regex.findall(sample_name)[0]
    return read_group_id



#  Note: currently only setup to work with SE reads -- for MSBB 
def merge_files(out_dir, sample_name, input_1_ubam, input_2, input_2_file_type, read_type, read_2, tmp_dir):
    if input_2 is not None:
        # fix sample names so they are different 
        sample_name_input2 = f"{sample_name}_input2"
        # convert second file to ubam (first file should already be)
        #get the read group for MSBB -- TODO probably should make this check if this is MSBB
        msbb_read_group = get_MSBB_read_group(sample_name)
        input_2_ubam = convert_to_ubam(out_dir, sample_name_input2, input_2_file_type, read_type, input_2, read_2, tmp_dir, msbb_read_group)
        # concatenate with samtools 
        unmapped_cat_bam = os.path.join(out_dir, f"{sample_name}_input1_input2_cat.bam")
        concat_ubams(input_1_ubam, input_2_ubam, unmapped_cat_bam)
        return unmapped_cat_bam
    else: 
        return input_1_ubam




# TODO decide if we want to log the size of the input file(s)

# 0.1 check all references exist (so that don't have to quit half way thorugh): 
check_refs_exist(args.STAR_index, args.ref_flat, args.annotation, args.annote_bed, args.rib_int, args.transcripts)

# 0.2 Set up output paths -- have an option to have this automatically structure like hydra
out_dirs = setup_output_dirs(args.out_struct, args.out_dir, args.cohort, args.tissue, args.sample)
print(f'''Output locations: \n fastqc: {out_dirs["fastqc"]} \n linear and tin processed: {out_dirs["linear_tin"]} 
salmon quant: {out_dirs["quant"]} \n multiqc: {out_dirs["multiqc"]} \n tin_summary: {out_dirs["tin_summary"]}''')

# 1. Run fastqc (if we don't need to merge files first)
if args.input_to_merge is None: 
    fastqc(out_dirs["fastqc"], args.sample, args.file_type, args.read_type, args.raw_input, args.input_read_2) 

# 2. convert to Ubam (from fastq or aligned bam)
ubam_out = convert_to_ubam(out_dirs["linear_tin"], args.sample, args.file_type, args.read_type, args.raw_input, args.input_read_2, args.tmp_dir)

# 2.1 concatenate if needed 
merged_ubam_out = merge_files(out_dirs["linear_tin"], args.sample, ubam_out, args.input_to_merge, args.merge_file_type, args.read_type, args.input_read_2, args.tmp_dir)

# 2.1.2 run fastqc on concatenated bams 
if args.input_to_merge is not None: 
    fastqc(out_dirs["fastqc"], args.sample, "ubam", args.read_type, merged_ubam_out, args.input_read_2)

# 3. Align with STAR
aligned_bam_out, aligned_transcript_bam_out = align_with_star(merged_ubam_out, args.sample, args.read_type, args.STAR_index, out_dirs["linear_tin"])

# 4. samtools sort
sorted_bam_out = sort(out_dirs["linear_tin"], args.sample, aligned_bam_out)

# 5. samtools index
indexed_bam_out = index(sorted_bam_out)

# 6. Post Alignment Picard QC ( Collect RNAseq metrics, Collect Alignemt summary metrics, Mark dups) 
RNAseq_metrics_out = os.path.join(out_dirs["fastqc"], f"{args.sample}.RNA_Metrics.txt")
picard_collect_RNA_metrics(sorted_bam_out, RNAseq_metrics_out, args.ref_flat, args.rib_int, args.tmp_dir, args.stranded)
Alignment_metrics_out = os.path.join(out_dirs["fastqc"], f"{args.sample}.Summary_metrics.txt")
picard_collect_alignment_metrics(sorted_bam_out, Alignment_metrics_out, args.tmp_dir, args.adapter_seq)
mark_dups_bam_out = os.path.join(out_dirs["linear_tin"], f"{args.sample}.Aligned.sortedByCoord.out.md.bam")
mark_dups_txt_out = os.path.join(out_dirs["fastqc"],f"{args.sample}.marked_dup_metrics.txt")
picard_mark_dups(sorted_bam_out, mark_dups_bam_out, mark_dups_txt_out, args.tmp_dir )

# 7. quantify with salmon 
salmon_out = os.path.join(out_dirs["quant"],f'{args.sample}.salmon')
salmon_quant(aligned_transcript_bam_out, salmon_out, args.annotation, args.transcripts)

# 8. TIN
# first index md bam
indexed_md_bam_out = index(mark_dups_bam_out)
#then run tin.py
calc_tin(mark_dups_bam_out, args.annote_bed)
# move summary tin ouput
tin_summary_current = f"{args.sample}.Aligned.sortedByCoord.out.md.summary.txt"
tin_summary_new_loc = os.path.join(out_dirs["tin_summary"], f"{args.sample}.Aligned.sortedByCoord.out.md.summary.txt")
tin_xls_current = f"{args.sample}.Aligned.sortedByCoord.out.md.tin.xls"
tin_xls_new_loc = os.path.join(out_dirs["linear_tin"], f"{args.sample}.Aligned.sortedByCoord.out.md.tin.xls")
shutil.move(tin_summary_current, tin_summary_new_loc)
shutil.move(tin_xls_current, tin_xls_new_loc)

# TODO run mulitqc -- this is normally run on all the samples together, so perhaps not?? -- maybe have last all complete function? 

## 9. Clean up
#  Remove raw data #TODO maybe only do this when delete flag is turned on -- might be misleading 
os.remove(args.raw_input)
if args.input_read_2 is not None:
    os.remove(args.input_read_2)

# delete the intermediate files we don't normally keep
delete_extras(args.delete, args.sample, aligned_bam_out, sorted_bam_out, indexed_bam_out, out_dirs["linear_tin"], indexed_md_bam_out, merged_ubam_out, args.input_to_merge)

# cram bams 
cram_bams(args.cram, mark_dups_bam_out, args.ref, args.sample, out_dirs["linear_tin"])

# TODO clean this up: copy everthing back to storage1
# TODO change this from hard coding to accepting an argument for storage locaion
# TODO look into using glob here otherwise make this it's own fuction
#storage1_path = os.path.join('/storage1/fs1/cruchagac/Active/jksanford/', args.cohort , '02.-ProcessedData')
#storage1_qc = os.path.join(storage1_path, '01.-QC')
#storage1_bams = os.path.join(storage1_path, '02.-bams')
#storage1_star = os.path.join(storage1_path, '03.-STAR/Hg38' )
#storage1_salmon = os.path.join(storage1_path, '04.-Salmon/Hg38')
#storage1_tin = os.path.join(storage1_path, '05.-PostAlign/TIN')
#storage1_circ = os.path.join(storage1_path, '06.-circRNA/Hg38')
# everything (with matching sample ID) in 01-QC
#shutil.move(os.path.join(args.storage_qc,args.sample_id+'.marked_dup_metrics.txt'), os.path.join(storage1_qc,args.sample_id+'.marked_dup_metrics.txt'))
#shutil.move(os.path.join(args.storage_qc,args.sample_id+'.RNA_Metrics.txt'), os.path.join(storage1_qc,args.sample_id+'.RNA_Metrics.txt'))
#shutil.move(os.path.join(args.storage_qc,args.sample_id+'.Summary_Metrics.txt'), os.path.join(storage1_qc,args.sample_id+'.Summary_Metrics.txt'))

# everything in 02-bams 
#shutil.move(os.path.join(args.storage_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.bam'), os.path.join(storage1_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.bam'))
#shutil.move(os.path.join(args.storage_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.bam.bai'), os.path.join(storage1_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.bam.bai'))

# everything in 03.-STAR
#shutil.move(os.path.join(args.storage_aligned,args.sample_id), os.path.join(storage1_star,args.sample_id))

# everything in 04.-Salmon 
#shutil.move(os.path.join(args.storage_quant,args.sample_id), os.path.join(storage1_salmon,args.sample_id))

# everything in 05.-PostAlign
#shutil.move(os.path.join(args.storage_tin,args.sample_id+'.Aligned.sortedByCoord.out.md.tin.xls'), os.path.join(storage1_tin,args.sample_id+'.Aligned.sortedByCoord.out.md.tin.xls'))
#shutil.move(os.path.join(args.storage_tin,args.sample_id+'.Aligned.sortedByCoord.out.md.summary.txt'), os.path.join(storage1_tin,args.sample_id+'.Aligned.sortedByCoord.out.md.summary.txt'))

# everything in  06.-circRNA/Hg38
#shutil.move(os.path.join(args.storage_circ,args.sample_id), os.path.join(storage1_circ,args.sample_id))
