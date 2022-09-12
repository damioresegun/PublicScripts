#!/usr/bin/env python3
'''
THIS SCRIPT IS INCOMPLETE!. IT WORKS TO A CERTAIN EXTENT BUT NEEDS EXTENSIVE EXPANSION


Script to carry out QIIME2 analysis 
for sequences derived from soil microbiomes
- Requirements:
    - qiime2 installed via conda
        - Would recommend following the instructions here:
        - https://docs.qiime2.org/2022.2/install/native/#install-qiime-2-within-a-conda-environment

IMPORTANT NOTE: The format of the input files is extremely important. Due to the nature of 
QIIME2, I have had to make this script extremely limited and specific for a particular format
type. This has been tested on the Casava 1.8 demultiplexed read format and this means that
the naming of the file has to conform to this format as well. This is just how QIIME2 reads
in the file. So the proper naming convention for the Casava format is 
P35791_lib581607_7976_R1_001.fastq.gz i.e. sample_barcode_lane_R1(or R2)_001.fastq.gz
If you choose to use other format types, please make sure they conform to the expected format and 
naming conventions!

Usage: 
Basic: python3 Workflow.py -i Demultiplexed/ -o Project_Out/ -t paired
'''
# Author: Dami Oresegun (2022)
# load in the modules
import os
import sys
import argparse
import subprocess
from pathlib import Path
##########################################################################
# set the needed arguments
def get_args():
    parser = argparse.ArgumentParser(description="Workflow to carry out" +
                                    "QIIME2 analysis on soil-derived " +
                                    "whole-genome metagenomic " + 
                                    "Illumina sequence data.")
    ##########################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-i", "--input", 
                               dest="Input_DIR",
                               action="store",
                               type=str,
                               help="The path to the directory " +
                               "holding the demultiplexed FASTQ " +
                               "files. FASTQ files can be gzipped " +
                               "or left uncompressed. "+ 
                               "Note: The fastq files have to be " +
                               "named with _1 and _2",
                               required=True)
    required_args.add_argument("-o", "--output",
                               dest="Output_DIR",
                               action="store",
                               type=str,
                               help="Path to the directory to save " +
                               "analysis outputs",
                               required=True)
    required_args.add_argument("-t", "--type",
                               dest="sequence_type",
                               action="store",
                               choices=['single','paired','other'],
                               default="paired",
                               type=str,
                               help="The type of read data you have. " +
                               "Can be single or paired end. " +
                               "Default is paired",
                               required=True)
    ##########################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-f", "--format",
                                dest="Input_FORMAT",
                                action="store",
                                type=str,
                                default="CasavaOneEightLanelessPerSampleDirFmt",
                                help="The format used to carry out the " +
                                "demultiplexing and quality encoding of " +
                                "of your samples. Only two choices shown " +
                                "however others can be added at your " +
                                "discretion, although this script has only "+
                                "been tested with the paired end format. " +
                                "Other formats can be seen in the " +
                                "import_formats.txt file of this package")
    optional_args.add_argument("-p", "--processes",
                               dest='Threads',
                               action="store",
                               type=int,
                               default="24",
                               help="Number of threads. Default is 24")
    optional_args.add_argument("-ot", "--other_type",
                                dest="Other_type",
                                action="store",
                                type=str,
                                help="Only to be used if using the [other] " +
                                "sequence type option. If chosen, the user " +
                                "will need to provide the chosen sequence " +
                                "type. Can be found in QIIME documentation " +
                                "or in the import_types.txt file in this " +
                                "package")
    ##########################################################################
    args = parser.parse_args()
    return args
##############################################################################
# set global variables
args = get_args()
INPDIR = args.Input_DIR
OUTDIR = args.Output_DIR
STYPE = args.sequence_type
SFORMAT = args.Input_FORMAT
THREADS = args.Threads
OTYPE = args.Other_type
##############################################################################
# check if the input exists
if os.path.exists(INPDIR):
    pass
else:
    print('Are you sure your input directory exists?')
    print('Input directory not found. Please check again')
    sys.exit(1)
# check if output exists
if os.path.exists(OUTDIR):
    print('Output folder already exists, outputs will be save here')
    pass
else:
    print('Output folder does not exist, creating it now')
    os.makedirs(OUTDIR)
# check sequence types and formats
if STYPE == "single":
    print('You have given single end reads')
    STYPEI = "SampleData[SequencesWithQuality]"
    pass
elif STYPE == "paired":
    print('You have given paired-end reads')
    STYPEI = "SampleData[PairedEndSequencesWithQuality]"
    pass
elif STYPE == "single" and SFORMAT == "PairedEndFastqManifestPhred33V2":
    print('You have given single end reads but gave a paired end format. ' +
    'These are incompatible. Please try again')
    sys.exit(1)
elif STYPE == "paired" and SFORMAT == "SingleEndFastqManifestPhred33V2":
    print('You have given paired end reads but gave a single end format. ' +
    'These are incompatible. Please try again')
    sys.exit(1)
elif STYPE == "other":
    print('You have provided a type or format that is not present in the ' +
    'default options. We hope this is compatible with the sequence ' +
    'type you have selected. However, be aware that there may be some ' +
    'failure due to this. If this occurs, please check on the QIIME2 ' +
    'documentation website')
    STYPEI = OTYPE
##############################################################################
'''Import the data into QIIME2'''
def importData(inpt,outpt,typer,formatr):
    # check and make imported folder
    impO = os.path.join(outpt, "ImportedData")
    if os.path.exists(impO):
        pass
    else:
        os.makedirs(impO)
    for i in Path(inpt).glob('*'):
        fname = os.path.basename(i).split("_")[0]
        # build the import command
        impOf = os.path.join(impO, fname)
        qimp = ("qiime tools import --input-path", str(i), "--type", typer,
                "--input-format", formatr, "--output-path", impOf)
        runQimp = ' '.join(qimp)
        print(runQimp)
        #subprocess.call(runQimp, shell=True)
        print('Import complete')
        # summarise the outputs
        qimv = ("qiime demux summarize --i-data", impOf+".qza", 
                "--o-visualization", impOf+".qzv")
        runQimv = ' '.join(qimv)
        print(runQimv)
        #subprocess.call(runQimv, shell=True)
        print('Demultiplexed reads have been summarised as qzv files')
        print('The best way to view them is to download the files and ' +
                'view them on www.view.qiime2.org')
        return impO
##############################################################################
def fastqc(inpt,outpt,threads):
    for i in Path(inpt).glob('*'):
        fname = os.path.basename(i).split("_")[0]
        fastOut = os.path.join(outpt, "FastQC", fname)
        if os.path.exists(fastOut):
            pass
        else:
            os.makedirs(fastOut)
        fasfi1 = os.path.join(i, "*_R1*")
        print(fasfi1)
        fasfi2 = os.path.join(i, "*_R2*")
        fatq1 = ("fastqc", "-t", str(threads), "-o", fastOut, fasfi1)
        runFatq1 = ' '.join(fatq1)
        print(runFatq1)
        fatq2 = ("fastqc", "-t", str(threads), "-o", fastOut, fasfi2)
        runFatq2 = ' '.join(fatq2)
        print(runFatq2)
        # run the fastqc command
        subprocess.call(runFatq1, shell=True)
        subprocess.call(runFatq2, shell=True)
##############################################################################
# first check quality of data
fastqc(INPDIR,OUTDIR,THREADS)
# call the import function
#impData = importData(INPDIR,OUTDIR,STYPEI,SFORMAT)
#print(impData)
