#!/usr/bin/env python3
import logging 
import subprocess
import os

logger = logging.getLogger('linear_RNA_pipeline')

def align_STAR_normal(input, file_out_prefix,  file_read_type,  genome_dir, n_thrd = '12', two_pass_mode = 'Basic', \
    mismatch_max = '999',  multimap_max =  '20', \
    SJ_overhang_min = '8', SJDB_overhang_min = '1', \
    mismatch_n_over_l_max = '0.1', align_intron_min = '20', \
    align_intron_max = '1000000', align_mates_gap_max = '1000000', \
    filter_type = 'BySJout', score_min_over_l_read = '0.33', \
    match_n_min_over_l_read = '0.33', limit_sjbd_insert_nsj = '1200000', \
    read_cmd = 'samtools view -h', SAM_strand_field = 'intronMotif', \
    filter_intron_motifs = 'None', soft_clip_at_ref_ends = 'Yes' , \
    quant_mode = ['TranscriptomeSAM', 'GeneCounts'], \
    out_sam_type = ['BAM', 'Unsorted'], \
    out_sam_mapped = 'Within', genome_load = 'NoSharedMemory', \
    limit_genome_gen_ram = '31000000000', genome_chr_bin_Nbits = '14', \
    chim_seg_min = '15', chim_junct_overhang_min = '15', \
    chim_out_type = ['WithinBAM', 'SoftClip'], chim_main_seg_mult_nmax = '1', 
    out_SAM_attributes = ['NH', 'HI', 'AS', 'nM', 'NM', 'ch']):

    logger.info('Starting STAR alingment')
    input_exists = os.path.exists(input)
    if (input_exists ):
        #TODO note: file_read_type needs to be either SAM SE or SAM PE I think
        cmd = (f'/opt/STAR-2.7.8a/bin/Linux_x86_64/STAR --runMode alignReads'
            f' --runThreadN {n_thrd}'
            f' --genomeDir {genome_dir}'
            f' --twopassMode {two_pass_mode}'
            f' --outFilterMismatchNmax {mismatch_max}'
            f' --outFilterMultimapNmax {multimap_max}'
            f' --alignSJoverhangMin {SJ_overhang_min}'
            f' --alignSJDBoverhangMin {SJDB_overhang_min}'
            f' --outFilterMismatchNoverLmax {mismatch_n_over_l_max}'
            f' --alignIntronMin {align_intron_min}'
            f' --alignIntronMax {align_intron_max}'
            f' --alignMatesGapMax {align_mates_gap_max}'
            f' --outFilterType {filter_type}'
            f' --outFilterScoreMinOverLread {score_min_over_l_read}'
            f' --outFilterMatchNminOverLread {match_n_min_over_l_read}'
            f' --limitSjdbInsertNsj {limit_sjbd_insert_nsj}'
            f' --readFilesIn {input}'
            f' --readFilesType {file_read_type}'
            f' --readFilesCommand {read_cmd}'
            f' --outFileNamePrefix {file_out_prefix}'
            f' --outSAMstrandField {SAM_strand_field}'
            f' --outFilterIntronMotifs {filter_intron_motifs}'
            f' --alignSoftClipAtReferenceEnds {soft_clip_at_ref_ends}'
            f' --quantMode {" ".join(quant_mode)}'
            f' --outSAMtype {" ".join(out_sam_type)}'
            f' --outSAMunmapped {out_sam_mapped}'
            f' --genomeLoad {genome_load}'
            f' --limitGenomeGenerateRAM {limit_genome_gen_ram}'
            f' --genomeChrBinNbits {genome_chr_bin_Nbits}'
            f' --chimSegmentMin {chim_seg_min}'
            f' --chimJunctionOverhangMin {chim_junct_overhang_min}'
            f' --chimOutType {" ".join(chim_out_type)}'
            f' --chimMainSegmentMultNmax {chim_main_seg_mult_nmax}'
            f' --outSAMattributes {" ".join(out_SAM_attributes)}'
        )
        print ('Command:' + cmd) 
        logger.info(cmd)
        cmd_to_call = cmd.split()
        subprocess.check_call(cmd_to_call)
    else: 
        logger.error('uBam file (STAR input) not found')
        #TODO add in
        #quit()

def align_STAR_fastq(input_read1, input_read2, file_out_prefix, genome_dir,  \
    file_read_type = 'zcat', n_thrd = '12', two_pass_mode = 'Basic', \
    mismatch_max = '999',  multimap_max =  '20', \
    SJ_overhang_min = '8', SJDB_overhang_min = '1', \
    mismatch_n_over_l_max = '0.1', align_intron_min = '20', \
    align_intron_max = '1000000', align_mates_gap_max = '1000000', \
    filter_type = 'BySJout', score_min_over_l_read = '0.33', \
    match_n_min_over_l_read = '0.33', limit_sjbd_insert_nsj = '1200000', \
    read_cmd = 'samtools view -h', SAM_strand_field = 'intronMotif', \
    filter_intron_motifs = 'None', soft_clip_at_ref_ends = 'Yes' , \
    quant_mode = ['TranscriptomeSAM', 'GeneCounts'], \
    out_sam_type = ['BAM', 'Unsorted'], \
    out_sam_mapped = 'Within', genome_load = 'NoSharedMemory', \
    limit_genome_gen_ram = '31000000000', genome_chr_bin_Nbits = '14', \
    chim_seg_min = '15', chim_junct_overhang_min = '15', \
    chim_out_type = ['WithinBAM', 'SoftClip'], chim_main_seg_mult_nmax = '1', 
    out_SAM_attributes = ['NH', 'HI', 'AS', 'nM', 'NM', 'ch']):

    logger.info('Starting STAR alingment')
    input1_exists =[]
    for R1_fastq in input_read1:
        input1_exists = os.path.exists(R1_fastq)
    input2_exists = []
    for R2_fastq in input_read2:
        input2_exists = os.path.exists(R2_fastq)
    if (all(input1_exists) and all(input2_exists)):
        cmd = (f'/opt/STAR-2.7.8a/bin/Linux_x86_64/STAR --runMode alignReads'
            f' --runThreadN {n_thrd}'
            f' --genomeDir {genome_dir}'
            f' --twopassMode {two_pass_mode}'
            f' --outFilterMismatchNmax {mismatch_max}'
            f' --outFilterMultimapNmax {multimap_max}'
            f' --alignSJoverhangMin {SJ_overhang_min}'
            f' --alignSJDBoverhangMin {SJDB_overhang_min}'
            f' --outFilterMismatchNoverLmax {mismatch_n_over_l_max}'
            f' --alignIntronMin {align_intron_min}'
            f' --alignIntronMax {align_intron_max}'
            f' --alignMatesGapMax {align_mates_gap_max}'
            f' --outFilterType {filter_type}'
            f' --outFilterScoreMinOverLread {score_min_over_l_read}'
            f' --outFilterMatchNminOverLread {match_n_min_over_l_read}'
            f' --limitSjdbInsertNsj {limit_sjbd_insert_nsj}'
            f' --readFilesIn {",".join(input_read1)} {",".join(input_read2)}'
            f' --readFilesType {file_read_type}'
            f' --readFilesCommand {read_cmd}'
            f' --outFileNamePrefix {file_out_prefix}'
            f' --outSAMstrandField {SAM_strand_field}'
            f' --outFilterIntronMotifs {filter_intron_motifs}'
            f' --alignSoftClipAtReferenceEnds {soft_clip_at_ref_ends}'
            f' --quantMode {" ".join(quant_mode)}'
            f' --outSAMtype {" ".join(out_sam_type)}'
            f' --outSAMunmapped {out_sam_mapped}'
            f' --genomeLoad {genome_load}'
            f' --limitGenomeGenerateRAM {limit_genome_gen_ram}'
            f' --genomeChrBinNbits {genome_chr_bin_Nbits}'
            f' --chimSegmentMin {chim_seg_min}'
            f' --chimJunctionOverhangMin {chim_junct_overhang_min}'
            f' --chimOutType {" ".join(chim_out_type)}'
            f' --chimMainSegmentMultNmax {chim_main_seg_mult_nmax}'
            f' --outSAMattributes {" ".join(out_SAM_attributes)}'
        )
        print ('Command:' + cmd) 
        logger.info(cmd)
        cmd_to_call = cmd.split()
        subprocess.check_call(cmd_to_call)
    else: 
        logger.error('uBam file (STAR input) not found')
        #TODO add in
        #quit()


def sort_samtools(input_align_bam, sample_name, threads = 12 ):
    logger.info('Starting samtools sort')
    input_exists = os.path.exists(input_align_bam)
    if(input_exists):
        cmd = (f'samtools sort --threads {threads}'
            f' -o {sample_name}.Aligned.sortedByCoord.out.bam'
            f' {input_align_bam}'
        )
        print ('Command:' + cmd)
        logger.info(f'Command: {cmd}')
        cmd_to_call = cmd.split()
        subprocess.check_call(cmd_to_call)
    else: 
        
        logger.error(f'Aligned bam {input_align_bam} (samtools input) not found')
        #TODO add in
        #quit()

def index_samtools(input_sorted_bam):
    input_exits = os.path.exists(input_sorted_bam)
    if(input_exits):
        cmd = f'samtools index {input_sorted_bam}'
        print('Command:' + cmd)
        logger.info(f'Command: {cmd}')
        cmd_to_call = cmd.split()
        subprocess.check_call(cmd_to_call)
    else: 
        logger.error(f'Sorted bam {input_sorted_bam} (samtools index input) not found')
        #TODO add in 
        #quit()

def cram_samtools(bam_to_cram, output_cram, ref): 
    cmd = (f'samtools view -C'
    f' -T {ref}'
    f' -o {output_cram}' 
    f' {bam_to_cram}')
    logger.info(f'Command: {cmd}')
    cmd_to_call = cmd.split()
    subprocess.check_call(cmd_to_call)