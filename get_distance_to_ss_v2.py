"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_distance_to_ss_v2.py: Given the ORFs per transcript, the sequences of each transcript and
the position of the ss, get the relative distance to this ss. Then, take the longest that starts
with the given start codon
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
        #
        orfs_path = sys.argv[1]
        start_codon_user = sys.argv[2]
        transcripts_paths = sys.argv[3]
        output_path = sys.argv[4]
        #
        # orfs_path = "/home/shinoda/Desktop/Florida/annotation/MBNL1_possible_ORFs_refseq.fa"
        # start_codon_user = "152017982"
        # transcripts_paths = "/home/shinoda/Desktop/Florida/annotation/MBNL1_possible_transcripts_refseq.fa.paths"
        # output_path = "/home/shinoda/Desktop/Florida/annotation/MBNL1_evaluated_refseq.paths"

        # # 1. Extract just the first ORF (the longest) from each transcript
        # outFile = open(orfs_path+".unique", 'w')
        # with open(orfs_path) as f:
        #     logger.info("Processing ORFs file...")
        #     transcript_id, transcript_id_old = "Transcript_0", ""
        #     for line in f:
        #         # If its header
        #         if(re.search(">", line)):
        #             #If the transcript id changes, output the next 2 lines
        #             if(not(re.search(transcript_id, line))):
        #                 tokens = line.rstrip().split(":")
        #                 transcript_id = tokens[0][1:]
        #                 pass_flag = False
        #                 outFile.write(line)
        #             else:
        #                 pass_flag = True
        #         # If its sequence
        #         else:
        #             if(pass_flag):
        #                 pass
        #             else:
        #                 outFile.write(line)
        # outFile.close()

        # 1. Using the info in the header of the ORFs file, we can obtain the information of
        # whether our stop codon is falling
        results = {}
        with open(orfs_path + ".unique") as f:
            logger.info("Processing ORFs unique file...")
            for line in f:
                # If its header
                if (re.search(">", line)):
                    # Two things: obtain the differences for each exon associated (with the cordinates)
                    # and obtain the length of the ORF
                    tokens = line.rstrip().split(":")
                    strand = tokens[1]
                    transcript_id = tokens[0][1:]
                    exons = tokens[2:-1]
                    ORF = tokens[-1].split("-")
                    # Get the difference for each exon. Go substracting the sequences according to the distance
                    # to the ORF end
                    first_pos = int(ORF[0])
                    second_pos = int(ORF[1])
                    pos = 0
                    flag_1, flag_2, flag_3 = True, False, False
                    i = 0
                    for x in exons:
                        coords = x.split("_")
                        length_exon = int(coords[1]) - int(coords[0])
                        if (flag_1):
                            # coords = x.split("_")
                            # length_exon = int(coords[1]) - int(coords[0])
                            if ((pos + length_exon) < first_pos):
                                pass
                                # pos += length_exon
                            # Reached start codon
                            else:
                                offset = first_pos - pos
                                if(strand=="+"):
                                    start_codon = int(coords[0]) + offset
                                else:
                                    start_codon = int(coords[1]) - offset
                                flag_1 = False
                                flag_2 = True
                                if ((pos + length_exon) < second_pos):
                                    pass
                                    # pos += length_exon
                                # Reached stop codon
                                else:
                                    # offset = second_pos - (pos + length_exon)
                                    offset = second_pos - pos
                                    if (strand == "+"):
                                        stop_codon = int(coords[0]) + offset
                                    else:
                                        stop_codon = int(coords[1]) - offset
                                    flag_2 = False
                                    flag_3 = True
                                    # pos += length_exon
                                    # if (i == len(exons) - 1):
                                    #     distance = pos - second_pos
                                    # else:
                                    #     pass
                        if (flag_2):
                            # coords = x.split("_")
                            # length_exon = int(coords[1]) - int(coords[0])
                            if ((pos + length_exon) < second_pos):
                                pass
                                # pos += length_exon
                            # Reached stop codon
                            else:
                                offset = second_pos - pos
                                if (strand == "+"):
                                    stop_codon = int(coords[0]) + offset
                                else:
                                    stop_codon = int(coords[1]) - offset
                                flag_2 = False
                                flag_3 = True
                                # pos += length_exon
                                # If reached last exon, the stop codon is in the last exon. Therefore, there won't be NMD
                                # if (i==len(exons)-1):
                                #     # distance = (pos + length_exon) - second_pos
                                #     distance = 0
                                #     break
                                # else:
                                #     pass
                        if (flag_3):
                            # If not reached last exon, get the distance to the 3'ss
                            if (i != len(exons) - 1):
                                # coords = x.split("_")
                                # length_exon = int(coords[1]) - int(coords[0])
                                # pos += length_exon
                                distance = (pos + length_exon) - second_pos
                                break
                            # If reached last exon, the stop codon is in the last exon. Therefore, there won't be NMD
                            else:
                                distance = 0
                                break
                        pos += length_exon
                        i += 1
                    results[transcript_id] = [start_codon, stop_codon, distance]
                else:
                    pass

        # 2. Get the longer transcript that starts with the given start codon
        results_filtered = {}
        for key, values in results.items():
            transcript_id = "_".join(key.split("_")[0:2])
            #Save it if the start_codon coincides with the one that the user specified
            if(values[0]==int(start_codon_user)):
                if(transcript_id not in results_filtered):
                    results_filtered[transcript_id] = values
                #If there is already a value, get the one with greater distance
                else:
                    if(results_filtered[transcript_id][2]<values[2]):
                        results_filtered[transcript_id] = values

        # 3. Associate the info to the transcripts_paths file
        outFile = open(output_path, 'w')
        with open(transcripts_paths) as f:
            logger.info("Processing transcripts_paths file...")
            for line in f:
                tokens = line.rstrip().split(":")
                transcript_id = tokens[0][1:]
                #Get the information from the previous results
                if(transcript_id in results_filtered):
                    info = results_filtered[transcript_id]
                    outFile.write(":".join(tokens)+"\tStart_codon: "+str(info[0])+
                                  "\tStop_codon: "+str(info[1])+
                                  "\tDistance: "+str(info[2])+"\n")
                else:
                    outFile.write(":".join(tokens)+"\tNo ORF\n")
        outFile.close()
        logger.info("Saved " + output_path)
        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()