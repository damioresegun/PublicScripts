#!/usr/bin/env python3
import sys
'''
This script is a small function that works specifically with the bmreads list output from BMTAGGER. 
The function takes in the list and the fastq files containing the human contamination and filters the identified reads
from the reads. This script is implemented in full in another package (AuOF). 
USAGE: ./Filter_BMTaggerList.py bmreads.list reads_1.fastq > reads_1.clean.fastq
'''
# Load in the human reads:
humanRead={}
for line in open(sys.argv[1]):
        humanRead[line.strip()] = None
# Print out the fastq file line by line unless the read is human:
skip = False
for i, line in enumerate(open(sys.argv[2])):
    if i%4 == 0:
        if line[1:].split("/")[0].split()[0] in humanRead: 
            skip = True
    else: 
        skip = False

    if skip == False:
        print(line.rstrip() + "\n")