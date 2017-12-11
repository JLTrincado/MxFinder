"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

test_NMDs.py: The following script do a Fisher test for each pair of exons, in order to test
independence between all the pairs of exons respect to the NMD condition.
Output: p_value of the Fisher_test and odds_ratio. The greater the odds_ratio the greater the association
of coocurrence with nonNMD. The smaller the odds_ratio the greater the association of coocurrence with nonNMD
(therefore, there is mutual exclusion)
"""

import logging, sys, os, re
import scipy.stats
from statsmodels.sandbox.stats.multicomp import multipletests
from itertools import combinations
from copy import deepcopy
from argparse import ArgumentParser, RawTextHelpFormatter

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

        paths_file = sys.argv[1]
        output_path = sys.argv[2]

        # paths_file = "/home/shinoda/Desktop/Florida/annotation/MBNL1_TEST/MBNL1_evaluated_refseq.paths"
        # output_path = "/home/shinoda/Desktop/Florida/annotation/MBNL1_TEST/MBNL1_evaluated_p_value_refseq.paths"

        # 1. Get the whole list of exons in the file
        exons_list = [] #list of exons
        with open(paths_file) as f:
            logger.info("Getting the whole list of exons in the file...")
            for line in f:
                tokens = line.rstrip().split("\t")
                if(tokens[1]!="No ORF"):
                    aux_list = re.split('\W+', tokens[0])[2:-1]
                    #Transform all elements to integer
                    aux_list_int = list(map(int, aux_list))
                    for x in aux_list_int:
                        #Put the exon in the exons_lis if it's not already
                        if(x not in exons_list):
                            exons_list.append(x)

        # 2. Read again the file, count by pairs how many times appears depending of the NMD
        exon_counts = {} # counts of NMD or not NMD per exon_pair, and coocurrence
        # (e.g: [2,3] -> [coocurrence:[NMD,noNMD], nocoocurrence:[NMD,noNMD]]
        with open(paths_file) as f:
            logger.info("Counting ocurrences of exon pairs...")
            for line in f:
                tokens = line.rstrip().split("\t")
                if(tokens[1]!="No ORF"):
                    if(int(re.split('\W+', tokens[3])[1])<50):
                        NMD_flag = False
                    else:
                        NMD_flag = True
                    aux_list = re.split('\W+', tokens[0])[2:-1]
                    # Transform all elements to integer
                    aux_list_int = list(map(int, aux_list))
                    # Get all the combinations for exons_list
                    combinations_list = list(combinations(exons_list, 2))
                    # Check if each pair show up in the transcript
                    for x in combinations_list:
                        #coocurrence
                        if((x[0] in aux_list_int and x[1] in aux_list_int) or
                               (x[0] not in aux_list_int and x[1] not in aux_list_int)):
                            if(x not in exon_counts):
                                #NMD
                                if(NMD_flag):
                                    exon_counts[x] = [[0,0],[1,0]]
                                #no NMD
                                else:
                                    exon_counts[x] = [[0,0],[0,1]]
                            else:
                                #NMD
                                if(NMD_flag):
                                    exon_counts[x][1][0] += 1
                                #no NMD
                                else:
                                    exon_counts[x][1][1] += 1
                        # no coocurrence
                        else:
                            if(x not in exon_counts):
                                #NMD
                                if(NMD_flag):
                                    exon_counts[x] = [[1,0],[0,0]]
                                #no NMD
                                else:
                                    exon_counts[x] = [[0,1],[0,0]]
                            else:
                                #NMD
                                if(NMD_flag):
                                    exon_counts[x][0][0] += 1
                                #no NMD
                                else:
                                    exon_counts[x][0][1] += 1


        # 3. Get all the combinations between pairs of exons, and do a fisher test.
        # Save the p-values for applying an FDR correction
        outFile = open(output_path+".temp", 'w')
        outFile.write("Pair_exons\todds_ratio\tp_value\tcontingency_table\n")
        raw_pvals = []
        for key, values in exon_counts.items():
            #Do a fisher test with the associated values
            oddsratio, pvalue = scipy.stats.fisher_exact([values[0],values[1]])
            raw_pvals.append(pvalue)
            #Output the results
            outFile.write(str(key)+"\t"+str(oddsratio)+"\t"+str(pvalue)+"\t["+
                          str(values[0])+","+str(values[1])+"]\n")
        outFile.close()

        # 4. Apply the FDR correction and output it in the file
        _, pvals_corrected, _, _ = multipletests(raw_pvals, method='fdr_bh', alpha=0.05)
        pvals_corrected_list = pvals_corrected.tolist()
        outFile = open(output_path, 'w')
        outFile.write("Pair_exons\todds_ratio\tp_value\tFDR\tcontingency_table\n")
        cont = 0
        with open(output_path+".temp") as f:
            next(f)
            logger.info("Applying multiple testing correction...")
            for line in f:
                tokens = line.rstrip().split("\t")
                outFile.write(tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+str(pvals_corrected_list[cont])+"\t"+tokens[3]+"\n")
                cont += 1
        outFile.close()

        # 5. Sort the output file
        outFile = open(output_path+".sorted", 'w')
        lines = open(output_path, 'r').readlines()[1:]
        outFile.write("Pair_exons\todds_ratio\tp_value\tFDR\tcontingency_table\n")
        for line in sorted(lines, key=lambda line: float(line.split("\t")[2])):
            outFile.write(line)
        outFile.close()

        # Remove the auxiliary file
        os.remove(output_path)
        os.remove(output_path+".temp")


        logger.info("Saved " + output_path+".sorted")
        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()