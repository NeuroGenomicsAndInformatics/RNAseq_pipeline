#!/usr/bin/env python3
import logging 
import subprocess

logger = logging.getLogger('linear_RNA_pipeline')

def salmon_quant(aligned_transciptome_bam, quant_output, annotation, transcripts_path, library_type = 'A', threads = 2):
    cmd = (f'/opt/salmon-1.7.0_linux_x86_64/bin/salmon quant'
    f' --threads {threads}'
    f' -t {transcripts_path}'
    f' -a {aligned_transciptome_bam}' 
    f' --libType {library_type}' # note: this infers the read type on the first few thousand reads
    f' --seqBias --gcBias'
    f' --geneMap {annotation}'
    f' -o {quant_output}')
    logger.info('Starting Salmon quantification')
    logger.info(f'Command:{cmd}')
    cmd_to_call = cmd.split()
    subprocess.check_call(cmd_to_call)