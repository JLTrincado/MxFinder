"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

format_gtf_refseq.py: given a gtf and a gene, extract the info assocated to this gene
and transform it to bed file. Just using the refseq annotation
"""

import sys
import time
import re

from argparse import ArgumentParser, RawTextHelpFormatter

import logging, sys, os

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


def format_gtf_refseq(gtf_path, gene_of_interest, output_path):
    try:

        # gtf_path = sys.argv[1]
        # gene_of_interest = sys.argv[2]
        # output_path = sys.argv[3]

        # gtf_path = "/home/shinoda/projects_rg/SUPPA2/annotation/refseq_hg19.formatted.gtf"
        # gene_of_interest = "MBNL1"
        # output_path = "/home/shinoda/Desktop/Florida/annotation/MBNL1_exons_hg19_refseq.bed"

        #1. Get the information of the gene and transform it ot bed
        outFile = open(output_path, 'w')
        with open(gtf_path) as f:
            logger.info("Processing "+output_path+"...")
            for line in f:
                #Skip the comments and all the lines not related with the gene of interest
                if(not(re.search("#", line)) and re.search("\""+gene_of_interest+"\"", line) and re.search("exon", line)):
                    tokens = line.rstrip().split("\t")
                    #Extract just the number of the exon
                    exon = "0"
                    transcript_id = str(tokens[8].split("transcript_id")[1])[2:-2]
                    gene_id = str(tokens[8].split(";")[0].split("gene_id")[1])[2:-1]
                    outFile.write(tokens[0]+"\t"+str(int(tokens[3])-1)+"\t"+tokens[4]+"\t"+gene_id+":"+transcript_id+
                        ":"+"exon_"+exon+":"+str(int(tokens[3])-1)+"_"+tokens[4]+":"+tokens[6]+"\t"+str(0)+"\t"+tokens[6]+"\n")

        outFile.close()
        logger.info("Saved " + output_path)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
