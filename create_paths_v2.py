"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

create_paths.py: create all the possible transcripts from the list of possible exons.
It returns the fasta sequence associated to each combination
V2: this version is without restricting the paths to a unique stop codon
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

# def get_unique_tuples(elements):
#     '''
#     Given a list of tuples (of 2 cordinates), return a unique list of elements
#     '''
#     aux_dict = {}
#     for tuple in elements:
#         if(tuple[0] not in aux_dict):
#             aux_dict[tuple[0]] = [tuple[1]]
#         else:
#             #Save the element if is not included in the associated list
#             if(tuple[1] not in aux_dict[tuple[0]]):
#                 aux_dict[tuple[0]].append(tuple[1])
#
#     #Return the unique list of elements
#     list_out = []
#     for key, value in aux_dict.items():
#         for x in value:
#             list_out.append([key,x])
#
#     return list_out

def get_pos(list, exon):
    '''
    Returns the position of the exon in the given list
    '''
    i = 0
    for x in list:
        if(x[0]==exon[0] and x[1]==exon[1]):
            return i
        i += 1

def check_all_possible_paths(matrix, local_path, pos):
    '''
    Using recursion, it generates all possible paths according to the transcipts-exon map
    This function call to an auxiliary one (check_all_possible_paths_aux)
    '''
    paths_list = []
    for j in range(0,len(matrix[pos])):
        if (matrix[pos][j] == 1):
            check_all_possible_paths_aux(matrix, paths_list, deepcopy(local_path), j+1)
    return paths_list

def check_all_possible_paths_aux(matrix, paths_list, local_path, pos):
    '''
    Auxiliary function of check_all_possible_paths
    '''
    # Put the value in local_path and keep on looking for more exons, from pos
    local_path.append(pos-1)
    # If it have reached the the exon_stop, save local_path in paths_dict if it doesn't exist and return
    if(pos-1==len(matrix[pos])-1):
        if (local_path not in paths_list):
            paths_list.append(local_path)
            return None
    else:
        for j in range(0,len(matrix[pos])):
            if (matrix[pos][j] == 1):
                if(j!=len(matrix[pos])-1):
                    check_all_possible_paths_aux(matrix, paths_list, deepcopy(local_path), j+1)
                # If it have reached the end, save local_path in paths_dict if it doesn't exist and return
                else:
                    if (local_path not in paths_list):
                        paths_list.append(local_path)
                        return None
        return None

def main():
    try:

        bed_path = sys.argv[1]
        fasta_path = sys.argv[2]
        output_path = sys.argv[3]

        # bed_path = "/home/shinoda/Desktop/Florida/annotation/A1BG_exons.bed"
        # fasta_path = "/home/shinoda/Desktop/Florida/annotation/A1BG_exons.bed.fa"
        # output_path = "/home/shinoda/Desktop/Florida/annotation/A1BG_possible_transcripts_refseq.fa"

        # exon_stop_coords = exon_stop.split("-")

        # 1. Sort the exon information. Also, save the information of the transcripts and their associated exons
        with open(bed_path) as f:
            logger.info("Processing bed file...")
            exons_list = []
            transcripts_dict = {}
            for line in f:
                tokens = line.rstrip().split("\t")
                chr = tokens[0]
                # Save the cordinates in the list
                if([tokens[1],tokens[2]] not in exons_list):
                    exons_list.append([tokens[1],tokens[2]])
                transcript = str(tokens[3].split(":")[1])
                if(transcript not in transcripts_dict):
                    transcripts_dict[transcript] = [[tokens[1],tokens[2]]]
                else:
                    transcripts_dict[transcript].append([tokens[1], tokens[2]])

        # If there is only 1 transcript, stop the execution.
        # Raise an ERROR exception (we are interested in cases with various transcripts)
        if(len(transcripts_dict)!=1):
            raise Exception("Only 1 transcript associated to this gene. Stop execution.")

        # Sort the list of exons
        exons_list_sorted = sorted(exons_list, key=lambda x: (x[1], x[0]))
        # # Get the position of the exon_stop in exons_list_sorted
        # exon_stop_coords_pos = exons_list_sorted.index(exon_stop_coords)

        # Sort the exons associated to each transcript
        for key, values in transcripts_dict.items():
            transcripts_dict[key] = sorted(values, key=lambda x: (x[1], x[0]))
        # # Remove all the transcripts that doesn't include the exon_end
        # # Copy the dictionary for iterating on it
        # transcripts_dict_copy = deepcopy(transcripts_dict)
        # for key, values in transcripts_dict_copy.items():
        #     if(exon_stop_coords not in values):
        #         del transcripts_dict[key]

        # 2. Associate the sequence to each exon
        exons_sequence = {}
        i = 0
        header_id = ""
        with open(fasta_path) as f:
            logger.info("Processing fasta file...")
            for line in f:
                if(re.search(">", line)):
                    prev_header_id = header_id
                    header_id = str(line.rstrip().split(":")[3])
                    header_flag = True
                # If we have reach a new header, save the stored sequence in
                if(i!=0):
                    if(header_flag):
                        exons_sequence[prev_header_id] = sequence
                        sequence = ""
                        header_flag = False
                    #Save the sequence lines until the next header
                    else:
                        sequence = sequence + line.rstrip()
                else:
                    sequence = ""
                    header_flag = False
                i += 1
        exons_sequence[header_id] = sequence

        # 3. Build the connection matrix
        width = len(exons_list_sorted)
        height = len(transcripts_dict)
        # Leave an extra row and column for the start and end positions
        matrix = [[0 for x in range(width+1)] for y in range(width+1)]
        # For each transcript, check the associated exons and mark them in the matrix
        i = 0
        for key, exon_list in transcripts_dict.items():
            j = 0
            pos_previous = 0
            for exon in exon_list:
                #Check the associated position in exons_list_sorted
                pos = get_pos(exons_list_sorted,exon)
                # Mark the first exon of this transcript as connected with start
                if(j==0):
                    matrix[0][pos] = 1
                    pos_previous = pos
                # Mark the last exon of this transcript as connected with end
                else:
                    matrix[pos_previous+1][pos] = 1
                    pos_previous = pos
                if (j == len(exon_list) - 1):
                    matrix[pos + 1][len(exons_list_sorted)] = 1
                j += 1
            i += 1

        # 4. For each row, get all the possible pathways from left to right
        paths_list = check_all_possible_paths(matrix, [], 0)

        # 5. Obtain the sequence for each pathway
        i = 1
        outFile = open(output_path, 'w')
        # Save all the combinations of the exons in a separate file
        outFile2 = open(output_path+".paths", 'w')
        for element in paths_list:
            header_id = ">Transcript_"+str(i)
            outFile2.write(header_id+": "+str(element) + "\n")
            full_seq = ""
            for x in element:
                # Check the coordinates associated to the number of the exon
                coords = exons_list_sorted[x]
                id_coords = coords[0]+"_"+coords[1]
                # Add this cordinates to the header_id
                header_id = header_id + ":" + id_coords
                # Get the associated sequence
                seq = exons_sequence[id_coords].rstrip()
                full_seq = full_seq + seq
            #Output the results
            outFile.write(header_id+"\n"+full_seq+"\n")
            i += 1
        outFile.close()
        outFile2.close()
        logger.info("Number of possible transcripts: "+str(len(paths_list)))
        logger.info("Saved " + output_path)
        logger.info("Saved " + output_path+".paths")


        # 6. Generate a bed track with the annotation of all the possible exons
        # Obtain all the unique possible exons from paths_list
        unique_exon_list = []
        for element in paths_list:
            for x in element:
                if(x not in unique_exon_list):
                    unique_exon_list.append(x)
        # Sort the list
        unique_exon_list_sorted = sorted(unique_exon_list)
        # Get the cordinates for each exon and output the results in an external bed tracks file
        outFile = open(output_path+".bedtrack", 'w')
        outFile.write("track name=exon_number color=0,0,0\n")
        for x in unique_exon_list_sorted:
            coords = exons_list_sorted[x]
            outFile.write(chr+"\t"+coords[0]+"\t"+coords[1]+"\t"+str(x)+"\n")
        outFile.close()
        logger.info("Generated also bed track with the info of all the possible exons: " + output_path+".bedtrack")

        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()