"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

MxFinder.py: main script. Check github for more details (https://github.com/JLTrincado/MxFinder)
"""

import sys
import time
import re

from argparse import ArgumentParser, RawTextHelpFormatter

import logging, sys, os
import subprocess
from lib.format_gtf_refseq import *
from lib.create_paths_v2 import *
from lib.get_distance_to_ss import *
from lib.test_NMDs import *

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# create console handler and set level to info
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


def main():
    try:

        logger.info("Starting execution")

        gtf_path = sys.argv[1]
        gene_of_interest = sys.argv[2]
        python2 = sys.argv[3]
        mosea = sys.argv[4]
        fast_genome = sys.argv[5]
        bedtools = sys.argv[6]
        output_path = sys.argv[7]

        # gtf_path = "/home/shinoda/Desktop/Florida/annotation/refseq_mm10_full.formatted.gtf"
        # gene_of_interest = "Mbnl1"
        # python2 = "/home/shinoda/miniconda3/envs/python2.7/bin/python2.7"
        # mosea = "/home/shinoda/Software/MoSEA-master/mosea.py"
        # fast_genome = "/home/shinoda/Software/MoSEA-master/test_files/genome/mm10.fa"
        # output_path = "/home/shinoda/Desktop/Florida/results/mm10"

        # 1. Obtain the exon coordinates from the GTF. Format it in a bed format
        output_path_aux =  output_path+"/"+gene_of_interest+"_exons.bed"
        format_gtf_refseq(gtf_path, gene_of_interest, output_path_aux)

        # 2. Run MoSEA (python2)
        subprocess.call([python2, mosea, 'getfasta','--bedfile', output_path_aux,'--genome',fast_genome,
                         '--bedtoolspath', bedtools, '--output', output_path_aux+'.fa'])

        # 3. Get all the possible paths with the info given. It returns the fasta sequence associated to each combination
        output_path_aux2 = output_path+"/"+gene_of_interest+"_possible_transcripts.fa"
        create_paths(output_path_aux, output_path_aux+'.fa', output_path_aux2)

        # 4. Run extract_orfs (python2)
        output_path_aux3 = output_path+"/"+gene_of_interest+"_possible_ORFs.fa"
        orfs_scripts_path = os.path.dirname(os.path.abspath(__file__))+"/extract_orfs.py"
        with open(output_path_aux3, "w") as outfile:
            subprocess.call([python2,orfs_scripts_path,output_path_aux2,'75'], stdout=outfile)
        logger.info("Saved "+output_path_aux3)

        # 5. Given the ORFs per transcript (take the longest per transcript), the sequences of each transcript and the
        # position of the ss, get the relative distance to this ss
        output_path_aux4 = output_path+"/"+gene_of_interest+"_evaluated.paths"
        get_distance_to_ss(output_path_aux3, output_path_aux2+".paths", output_path_aux4)

        # 6. The following script do a Fisher test for each pair of exons, in order to test independence between all the exons respect to the NMD condition
        output_path_aux5 = output_path+"/"+gene_of_interest+"_evaluated_p_value.paths"
        test_NMDs(output_path_aux4, output_path_aux5)

        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()