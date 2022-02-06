#!/bin/bash
#$ -cwd

# Script to take in raw nanopore files and take them through basecalling. Set up for Albacore, Guppy barcoding, guppy demultiplexer and qcat demultiplexer. Uncomment for the ones where necessary. Ensure that this is run in environments where the necessary tools are. 

# This currently set for qcat demultiplexing
set -e
#Set Paths
Raw_IN="$HOME/All_RawFiles/June2019_Data/LWBPBK_4_06_19_CLINICALS/LWBPBK_4_06_19_CLINICALS/20190604_2030_MN17366_FAJ02943_3b1eee27/fast5/"
Base_OUT="$HOME/All_GuppyBasecalled/June2019_Data1"   # output folder for the basecalling
Dem_OUT="$HOME/All_QcatDemultiplex/June2019_Data1"    # output folder for the demultiplexing
THREADS="32"   # number of threads
CALLERS="4"  # number of callers
PERCAL="8"   # number of threads per caller
FLOWCELL="FLO-MIN106"
KIT="SQK-RBK004"  # sequencing kit used. this version is for guppy
KIT2="RBK004"    # sequencing kit used. this version is for qcat

mkdir -p $Base_OUT
mkdir -p $Dem_OUT

#Albacore command
#sf5_read_fast5_basecaller.py -i ${input} -t 14 -s ${output} -f FLO-MIN106 -k SQK-RBK004 --barcoding -o fast5,fastq -q 0

#guppy basecaller
guppy_basecaller --flowcell $FLOWCELL --kit $KIT --fast5_out -i ${Raw_IN} -s ${Base_OUT} -v -r -q 0 --qscore_filtering --num_callers $CALLERS --cpu_threads_per_caller $PERCAL

#guppy barcoder
#ginput=$HOME/ReBaseCall/Basecalls/April2019/pass
#goutput=$HOME/ReBaseCall/Demultiplex/April2019
#guppy_barcoder -i ${ginput} -s ${goutput} --barcode_kits "SQK-RBK004" -t 10 -q 0

#qcat
cat $Base_OUT/pass/* | qcat -b $Dem_OUT --detect-middle -t $THREADS --trim -k $KIT2 --dual

echo "All done"
