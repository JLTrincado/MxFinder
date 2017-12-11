"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

check_results.py: check the results
"""

import logging, sys, os, re
from copy import deepcopy
from argparse import ArgumentParser, RawTextHelpFormatter
import glob

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

        path = sys.argv[1]
        output_path = sys.argv[2]
        #
        # path = "/home/shinoda/Desktop/Florida/annotation/tables_introns"
        # output_path = "/home/shinoda/Desktop/Florida/annotation/tables_introns/results.tab"

        files = glob.glob(path+"/*.paths.sorted")
        outFile = open(output_path, 'w')
        outFile.write("Gene\tPair_exons\todds_ratio\tp_value\tcontingency_table\n")
        for file in files:
            gene = file.split("/")[-1].split("_")[0]
            with open(file) as f:
                logger.info("Processing file "+file)
                next(f)
                for line in f:
                    #Check if there is any signficant line. If that's the case, output it
                    tokens = line.rstrip().split("\t")
                    if(tokens[2]!="nan"):
                        p_value = float(tokens[2])
                        if(p_value<=0.05):
                            outFile.write(gene+"\t"+"\t".join(tokens)+"\n")
                        elif(p_value==1.0):
                            break

        outFile.close()
        logger.info("Generated "+output_path)
        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()