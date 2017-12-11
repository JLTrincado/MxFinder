"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

check_introns_3UTR.py: check which genes have introns in the 3' UTR region
"""

import logging, sys, os, re
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

        gtf_path = sys.argv[1]
        output_path = sys.argv[2]
        output_path = sys.argv[3]
        #

        # gtf_path = "/home/shinoda/Desktop/Florida/annotation/refseq_hg19_full.formatted.gtf"
        # output_path = "/home/shinoda/Desktop/Florida/annotation/genes_introns_3UTR.tab"
        # output_path2 = "/home/shinoda/Desktop/Florida/annotation/genes_no_introns_3UTR.tab"


        # 1. Get all the information for each gene
        transcript_gene, transcripts_exons, transcripts_CDS, transcripts_strand, transcripts_stop_codon = {}, {}, {}, {}, {}
        with open(gtf_path) as f:
            logger.info("Processing gtf file...")
            for line in f:
                tokens = line.rstrip().split("\t")
                #We are gonna add the chr to the transcript_id, cause there are some transcripts repeated in different
                # chromosomes (e.g: chrX and Y)
                transcript_id = tokens[0]+"|"+str(tokens[8].split("transcript_id")[1])[2:-2]
                gene_id = str(tokens[8].split(";")[0].split("gene_id")[1])[2:-1]
                #Save the association transcript - gene
                transcript_gene[transcript_id] = gene_id
                # Save the strand
                if(transcript_id in transcripts_strand):
                    if(transcripts_strand[transcript_id]!=tokens[6]):
                        raise Exception("Opposite strands associated to the same transcript: "+transcript_id)
                    else:
                        pass
                else:
                    transcripts_strand[transcript_id] = tokens[6]
                # Save the exon cordinates
                if(tokens[2]=="exon"):
                    if (transcript_id in transcripts_exons):
                        transcripts_exons[transcript_id].append([tokens[3],tokens[4]])
                    else:
                        transcripts_exons[transcript_id] = [[tokens[3],tokens[4]]]
                # Save the stop codon cordinates
                if(tokens[2]=="stop_codon"):
                    if (transcript_id in transcripts_stop_codon):
                        if(tokens[6]=="+" and tokens[4]>transcripts_stop_codon[transcript_id][1]):
                            transcripts_stop_codon[transcript_id] = transcripts_stop_codon[transcript_id][0]+"-"+tokens[4]
                        elif(tokens[6]=="-" and tokens[3]<transcripts_stop_codon[transcript_id][0]):
                            transcripts_stop_codon[transcript_id] = tokens[3]+"-"+transcripts_stop_codon[transcript_id][1]
                    else:
                        transcripts_stop_codon[transcript_id] = [tokens[3],tokens[4]]

        # 2. Sort the exons cordinates, depending on the strand
        logger.info("Sorting exons...")
        for key,values in transcripts_exons.items():
            if(transcripts_strand[key] == "+"):
                transcripts_exons[key] = sorted(values, key=lambda x: (x[1], x[0]))
            else:
                transcripts_exons[key] = sorted(values, key=lambda x: (x[1], x[0]), reverse=True)

        # 3. Check for each transcript_id if the after the stop_codon there are more exons
        logger.info("Checking the stop codon localizations...")
        list_genes_introns, list_genes_no_introns = [], []
        for key, values in transcripts_exons.items():
            flag_stop_codon = False
            gene_id = transcript_gene[key]
            if(key in transcripts_stop_codon):
                stop_codon = transcripts_stop_codon[key]
                for element in values:
                    #Check if the stop codon is in any of the exons
                    if(transcripts_strand[key] == "+"):
                        if((not flag_stop_codon) and (element[0]<=transcripts_stop_codon[key][1] and
                        element[1]>=transcripts_stop_codon[key][1])):
                            flag_stop_codon = True
                            #If we have reached the last exon, put it in the list of genes with no introns in 3UTR
                            if(element[0]==values[-1][0] and element[1]==values[-1][1]):
                                if (gene_id not in list_genes_no_introns):
                                    list_genes_no_introns.append(gene_id)
                                break
                        elif(flag_stop_codon):
                            # There are introns in the 3UTR region. Output this gene
                            if(gene_id not in list_genes_introns):
                                list_genes_introns.append(gene_id)
                            break
                    else:
                        if ((not flag_stop_codon) and (element[0] <= transcripts_stop_codon[key][0] and
                                                               element[1] >= transcripts_stop_codon[key][0])):
                            flag_stop_codon = True
                            #If we have reached the last exon, put it in the list of genes with no introns in 3UTR
                            if(element[0]==values[-1][0] and element[1]==values[-1][1]):
                                if (gene_id not in list_genes_no_introns):
                                    list_genes_no_introns.append(gene_id)
                                break
                        elif (flag_stop_codon):
                            # There are introns in the 3UTR region. Output this gene
                            if(gene_id not in list_genes_introns):
                                list_genes_introns.append(gene_id)
                            break

        # 4. Output the list of genes sorted with introns
        list_genes_introns_sorted = sorted(list_genes_introns)
        outFile = open(output_path, 'w')
        for x in list_genes_introns_sorted:
            outFile.write("\""+x+"\""+"\n")
        outFile.close()
        logger.info("Generated "+output_path)

        # 5. Output the list of genes sorted without introns
        list_genes_no_introns_sorted = sorted(list_genes_no_introns)
        outFile = open(output_path2, 'w')
        for x in list_genes_no_introns_sorted:
            outFile.write("\""+x+"\""+"\n")
        outFile.close()
        logger.info("Generated "+output_path2)

        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()