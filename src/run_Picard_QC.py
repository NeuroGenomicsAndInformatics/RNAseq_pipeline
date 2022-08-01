#!/usr/bin/env python3
import logging 
import subprocess

logger = logging.getLogger('linear_RNA_pipeline')
#TODO figure out if strand specficity has a default here
def picard_collect_RNA_metrics(input, output, ref_flat, ribosomal_int, tmp, strand_specficity, memory = 8 ):
    #TODO if else for input not found
    print('Running Picard QC: RNA metrics')

    cmd = (f'java -Dpicard.useLegacyParser=false -jar -Xmx{memory}g /opt/picard-tools/picard.jar CollectRnaSeqMetrics'
        f" -I {input}" 
        f" -O  {output}"
        f" -REF_FLAT {ref_flat}"
        f" -RIBOSOMAL_INTERVALS {ribosomal_int}"
        f" -STRAND_SPECIFICITY {strand_specficity}"
        f" -TMP_DIR {tmp}"
    )
    logger.info('Running Picard QC: RNA metrics')
    logger.info(f"Command: {cmd}")
    cmd_to_call = cmd.split()
    subprocess.check_call(cmd_to_call)

def picard_collect_alignment_metrics(input, output, tmp, adapter_seq, memory = 8):
    print('Running Picard QC: Collect Alignment metrics')
    logger.info('Running Picard QC: Collect Alignment metrics')
    cmd = (
        f"java -Dpicard.useLegacyParcer=false -jar -Xmx{memory}g /opt/picard-tools/picard.jar CollectAlignmentSummaryMetrics"
        f" -I {input}"
        f" -ADAPTER_SEQUENCE {adapter_seq}"
        f" -O {output}"
        f" -TMP_DIR {tmp}"
    )
    logger.info(f"Command {cmd}")
    cmd_to_call = cmd.split()
    subprocess.check_call(cmd_to_call)

def picard_mark_dups(input, output, metrics_out, tmp, prgm_id = 'null', \
     assume_sort_order = 'coordinate', optical_duplicate_pixel_distance = 100,  max_in_ram = 150000,  memory = 8 ):
    print('Running Mark Duplicates')
    logger.info('Running Mark Duplicates')
    cmd = (
        f"java -Dpicard.useLegacyParcer=false -jar -Xmx{memory}g  /opt/picard-tools/picard.jar MarkDuplicates"
        f" -I {input}"
        f" -O {output}"
        f" -PROGRAM_RECORD_ID  {prgm_id}"
        f" -M {metrics_out}"
        f" -TMP_DIR {tmp}"
        f" -ASSUME_SORT_ORDER  {assume_sort_order}"
        f" -MAX_RECORDS_IN_RAM {max_in_ram}"
        f" -OPTICAL_DUPLICATE_PIXEL_DISTANCE {optical_duplicate_pixel_distance}"
    )
    logger.info(f"Command {cmd}")
    cmd_to_call = cmd.split()
    subprocess.check_call(cmd_to_call)