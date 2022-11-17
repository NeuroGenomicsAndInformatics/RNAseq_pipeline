#!/usr/bin/env python3
import logging 
import subprocess

logger = logging.getLogger('linear_RNA_pipeline')

def fastqc_fastq_PE(R1_fastq, R2_fastq, sample_id, fastqc_out_dir):
    print('Running fastqc')
    cmd1 = (f"/opt/FastQC-0.11.9/fastqc"
        f" {R1_fastq}" 
        f" -o {fastqc_out_dir}"
    )
    cmd2 = (f"/opt/FastQC-0.11.9/fastqc"
        f" {R2_fastq}" 
        f" -o {fastqc_out_dir}"
    )

    logger.info('Running fastqc')
    logger.info(f"Command: {cmd1}")
    logger.info(f"Command: {cmd2}")
    cmd_to_call1 = cmd1.split()
    cmd_to_call2 = cmd2.split()
    subprocess.check_call(cmd_to_call1)
    subprocess.check_call(cmd_to_call2)

def fastqc_SE(input, sample_id, fastqc_out_dir):
    print('Running fastqc')
    # TODO -- check if input exits here since this step is first 
    cmd = (f"/opt/FastQC-0.11.9/fastqc"
        f" {input}" 
        f" -o {fastqc_out_dir}"
    )


    logger.info('Running fastqc')
    logger.info(f"Command: {cmd}")
    cmd_to_call = cmd.split()
    subprocess.check_call(cmd_to_call)