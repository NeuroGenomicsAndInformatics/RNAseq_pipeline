#!/usr/bin/env python3
import logging 
import subprocess
import os
import sys



logger = logging.getLogger('linear_RNA_pipeline')

# define function to convert fastq to ubam
def fastq_to_ubam_SE(input_fq, output, sample_name, memory = 8):
    input_exists = os.path.exists(input_fq)
    logger.info('Running fastq to ubam')
    if (input_exists): 
        cmd = (f"java -Dpicard.useLegacyParser=false -jar -Xmx{memory}g /opt/picard-tools/picard.jar FastqToSam"
            f" -FASTQ {input_fq}"
            f" -O {output}"
            f" -SM {sample_name}"
            )
        print('Command:'+ cmd)
        logger.info('Command: {cmd}')
        cmd_to_call = cmd.split()
        subprocess.check_call(cmd_to_call)
    else: 
        logger.error('No fastq file found. Exiting.')
        sys.exit("No fastq file found. Exiting.")
        

def fastq_to_ubam_PE(input_fq_1, input_fq_2, output, sample_name, memory = 8):
    logger.info('Running fastq to ubam')
    fq1_exists = os.path.exists(input_fq_1)
    fq2_exists = os.path.exists(input_fq_2)
    if (fq1_exists and fq2_exists):
        cmd = (f"java -jar -Xmx{memory}g /opt/picard-tools/picard.jar FastqToSam"
            f' F1={input_fq_1} F2={input_fq_2}'
            f' O={output}'
            f' SM={sample_name}'
        )
        print ('Command:' + cmd) 
        logger.info(cmd)
        cmd_to_call = cmd.split()
        subprocess.check_call(cmd_to_call)

    else: 
        if fq1_exists:
            logger.error(f'Fastq file {input_fq_2} not found, exiting')
            sys.exit(f'Fastq file {input_fq_2} not found, exiting')
        elif fq2_exists:
            logger.error(f'Fastq file {input_fq_1} not found, exiting')
            sys.exit(f'Fastq file {input_fq_2} not found, exiting')
        else: 
            logger.error (f'Fastq files {input_fq_1} and {input_fq_2} not found, exiting')
            sys.exit(f'Fastq files {input_fq_1} and {input_fq_2} not found, exiting')

 
def bam_to_ubam(input_bam, output, tmp, \
    by_readgroup = 'true', sanitize = 'true', \
    max_discard = '0.005', max_ram_rec = '150000', \
    val_stringency = 'LENIENT', \
    clear_attr = ['NH', 'HI', "nM", "NM", 'ch'], memory = 8): # fix clear attributes here
    logger.info('Reverting bam')
    bam_exists = os.path.exists(input_bam)
    print('Reverting bam to ubam')
    if(bam_exists):
        cmd = (f"java -Dpicard.useLegacyParser=false -jar -Xmx{memory}g /opt/picard-tools/picard.jar RevertSam"
            f' -I {input_bam}'
            f' -OUTPUT_BY_READGROUP {by_readgroup}' 
            f' -SANITIZE {sanitize}'
            f' -MAX_DISCARD_FRACTION {max_discard}'
            f' -MAX_RECORDS_IN_RAM {max_ram_rec}'
            f' -VALIDATION_STRINGENCY {val_stringency}'
            f' -ATTRIBUTE_TO_CLEAR {" -ATTRIBUTE_TO_CLEAR ".join(clear_attr)}'
            f' -TMP_DIR {tmp}'
            f' -O {output}'
        ) 


        print('Command:'+ cmd)
        cmd_to_call = cmd.split()
        logger.info(cmd)
        subprocess.check_call(cmd_to_call)
    else: 
        logger.error('No input bam file found')
        sys.exit('No input bam file found')
def concat_ubams(input_bam_1, input_bam_2, output, threads = 12):
    logger.info('Starting samtools concatenate')
    cmd = (f"samtools cat -b {input_bam_1} {input_bam_2} -o {output} -@ {threads}" )
    cmd_to_call  = cmd.split()
    logger.info(cmd)
    subprocess.check_call(cmd_to_call)