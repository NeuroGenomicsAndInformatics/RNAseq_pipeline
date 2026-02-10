#!/usr/bin/env python3
import logging 
import subprocess

logger = logging.getLogger('linear_RNA_pipeline')
def calc_tin(input_md_bam, reference):
    cmd = (f'tin.py -i {input_md_bam}'
    f' -r {reference}')
    print ('Command:' + cmd) 
    logger.info('Starting TIN')
    logger.info(f'Command:{cmd}')
    cmd_to_call = cmd.split()
    subprocess.check_call(cmd_to_call)

