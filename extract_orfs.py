	#! /usr/bin/python
__author__="jruiz"
__date__ ="$Sep 08, 2015 12:24:43 PM$"
'''Extract all ORFs in a transcript FASTA
'''

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

fasta = sys.argv[1]
try:
	threshold = sys.argv[2] #nucleotides
except:
	threshold = 75

#OBJECTS
class orf_object:
	def __init__(self, sequence, start, end):
		self.sequence = sequence
		self.start = start
		self.end = end

#FUNCTIONS
def find_all(sequence, subsequence):
	''' Returns a list of indexes within sequence that are the start of subsequence'''
	start = 0
	idxs = []
	next_idx = sequence.find(subsequence, start)
 
	while next_idx != -1:
		idxs.append(next_idx)
		start = next_idx + 1# Move past this on the next time around
		next_idx = sequence.find(subsequence, start)
 
	return idxs
		
def find_orfs(sequence, threshold):
	""" Finds all valid open reading frames in the string 'sequence', and
		returns them as a list"""
 
	starts = find_all(sequence, 'ATG')
	stop_amber = find_all(sequence, 'TAG')
	stop_ochre = find_all(sequence, 'TAA')
	stop_umber = find_all(sequence, 'TGA')
	stops = stop_amber + stop_ochre + stop_umber
	stops.sort()

	orfs = []

	for start in starts:
		for stop in stops:
			if start < stop \
				and (start - stop) % 3 == 0:  # Stop is in-frame
					if len(sequence[start:stop+3]) >= int(threshold):
						orf_obj = orf_object(sequence[start:stop+3], start, stop+3)
						orfs.append(orf_obj)
					break

	orfs.sort(key=lambda x: len(x.sequence), reverse=True)
	return orfs

cdna = SeqIO.index(fasta, "fasta")
for sequence in cdna:
	ends = []
	orfs = find_orfs(cdna[sequence].seq, threshold)
	n = 1
	for orf in orfs:
		if orf.end in ends:
			continue
		ends.append(orf.end)
		print(">" + sequence + "_" + str(n) + ":" + str(orf.start) + "-" + str(orf.end) + "\n" + str(orf.sequence))
		n += 1

exit(0)