#!/bin/bash
## Automation script nanopore sequencing analysis. 
## Version 0.1.0 
## Works to take the data through albacore, poretools, porechop, assembly stats, nanoQC, nanoFilt and NanoPlot
## Requires (all future implementations) directory for data, directory to generate outputs
## Requires all tools to be available in the bin folder. If not, enter in configuration file (to be implemented if necessary)
## Future implements: flags?
## NOTE: Will need to change how this is set up. Rather than carrying out each sub step and checking for barcodes. Check for barcodes if yes take through the entire process. If not then take through the entire process. Should make it faster than checking for barcode each and every step
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#configFile=""
#WPATH="" # path to where you want to be the working directory
DPATH="$HOME/5xBarcodingJCS1734/RawData/reads/20171219_1833_barcodex5JCS34_17/downloads/pass/" # path to MinION data.
#DPATH="$HOME/All"
#DPATH="$HOME/LaptopData_preBarcode/RawData/20171016_1550_SKS201a/downloads/pass/" # path to MinION data
SPATH="$HOME/Scripts/Testing/"  # path to save destination
PTOC="SQK-RAD002"    # State the protocol used e.g SQK-RAD002
FLOCL="FLO-MIN106" # State the flow cell used
THRD="6" # State the number of threads you wish to use
PORTOL="poretools" # Path to poretools if not in $PATH
#PLTFM="y" # Is this run on Marvin?
FINAME="All_Test" # Give a name for this run; all files will be saved with this prefix
BACODS="y" # Was this a barcoding run?
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                            SET TOOLS, PATHS, DOWNLOAD DATA
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################################################################################################################
## Get the path to save all results in a created folder
echo $SPATH
cd $SPATH
## Make a folder to hold all outputs from this script
mkdir -p "$SPATH/PreProcessing"
currPATH="$SPATH/PreProcessing"
echo
## Change into new working directory
cd $currPATH
echo
echo "You are in the working directory: $currPATH"
echo
## Determine if user is working on Marvin in order to load necessary modules. Will be needed later for NanoFilt
echo "Is this run on marvin? Enter Y or N:"
read PLTFM
#................................................................................................................................................................
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                            RAW MINKNOW STATS
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################################################################################################################
# Check if the data is barcoded or not.
if [ $BACODS = "Y" ] || [ $BACODS = "y" ] || [ $BACODS = "yes" ] || [ $BACODS = "Yes" ]
then
# If data is barcoded; make directory to hold the individual barcode stats of each stage
mkdir -p $currPATH/Statistics/RawMinKNOWStats
echo
echo "Running pipeline for barcoded MinION data..."
# for each barcode in the provided path:
for i in $DPATH/*; do
echo
#echo $i
echo
# Check if the platform being used is marvin, if yes:
if [ $PLTFM = "Y" ] || [ $PLTFM = "y" ] || [ $PLTFM = "yes" ] || [ $PLTFM = "Yes" ]
then
echo "Running Poretools on Marvin..."
echo
namr=$(basename $i)
# Carry out stats analysis using poretools and save in the created directory with the barcode as the filename
echo "poretools stats $i > $currPATH/Statistics/${namr}_RawMinKNOWStats.txt"
else
# if platform is not marvin then use the provided path to poretools and save each barcode in the created directory
echo "${PORTOL} stats $i > $currPATH/Statistics/${i}_RawMinKNOWStats.txt"
fi   # if marvin statement is completed here
done   # for statement is completed here
# Else if the data is not barcoded, check if the platform is marvin or not
elif [ $PLTFM = "Y" ] || [ $PLTFM = "y" ] || [ $PLTFM = "yes" ] || [ $PLTFM = "Yes" ]
then
echo "Running Poretools on Marvin..."
echo
# Run poretools on Marvin on sequencing run MinION data saving the results under the provided filename
echo "poretools stats $DPATH > $currPATH/Statistics/${FINAME}_RawMinKNOWStats.txt"
else
echo "Running Poretools on your local system..."
# Else use the path to poretools to run poretools and save the results under the provided user name
${PORTOL} stats $DPATH > $currPATH/Statistics/${FINAME}_RawMinKNOWStats.txt
echo "Poretools complete. The statistics are saved in you directory."
echo "Moving to next step..."
fi   # if barcoded loop completed here
exit
#................................................................................................................................................................
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                            ALBACORE BASECALLING
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################################################################################################################
#................................................................................................................................................................
echo "Starting Albacore Basecalling..."
## Make directory for albacore output
mkdir -p $currPATH/Albacore_Output
## Based on answer of previous question; carry out albacore
if [ $PLTFM = "Y" ] || [ $PLTFM = "y" ] || [ $PLTFM = "yes" ] || [ $PLTFM = "Yes" ]
then
echo "Running Albacore on Marvin..."
## Load module required
module load python/3.4
## Display the module loaded
python -V
#ytp='python -V'
echo "Your python module is loaded"
## Run albacore on marvin
#read_fast5_basecaller.py -f $FLOCL -k $PTOC -o fast5,fastq -i $DPATH -r -s $currPATH/Albacore_Output/ -t $THRD -q 0

#Checking to see if the command fails or succeeds
#if [ $? -eq 0 ]; then
#    echo OK
#else
#    echo FAIL
#fi
echo
else
## If running on a local PC; it is assumed that the user has got the files in $PATH. If not; paths will be in config file (TO-DO: Make config file)
echo "Running Albacore on your local system..."
echo
## Run albacore on local system
#read_fast5_basecaller.py -f $FLOCL -k $PTOC -o fast5,fastq -i $DPATH -r -s $currPATH/Albacore_Output/ -t $THRD -q 0
fi
echo
echo "Albacore finished. Moving to the next step... "
#................................................................................................................................................................
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                            ALBACORE STATS
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################################################################################################################
#................................................................................................................................................................
#
#................................................................................................................................................................
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                             PORETOOLS CONVERSION AND PORECHOP DEMULTIPLEXING
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################################################################################################################
#................................................................................................................................................................
echo
echo "Beginning Poretools analysis of your basecalled data..."
pwd
# make directory for output
mkdir -p $currPATH/Poretools_Output
cd Poretools_Output
BACODS="y"
# check if data is barcoded
if [ $BACODS = "Y" ] || [ $BACODS = "y" ] || [ $BACODS = "yes" ] || [ $BACODS = "Yes" ]
then
echo "Running Poretools demultiplexing to fastq on Marvin..."
#### CHANGE SOURCE PATH TO ALBACORE OUTPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DPATH="$HOME/5xBarcodingJCS1734/RawData/reads/20171219_1833_barcodex5JCS34_17/downloads/pass/"
echo $DPATH
# get the names of the subdirectories present in the path
coda=($DPATH/*/)
#iterate through the path folder to get folder names
for coda in "${coda[@]}"
do
    echo "$coda"
    codar=$(basename "$coda")    
    echo $codar
    mkdir -p $codar
    #cd $codar
    poretools fastq ${coda}/*/ > $codar/${codar}.fastq
    gzip $codar/${codar}.fastq
done
else 
echo "Your data is not barcoded..."
echo
echo "Running Poretools" 
fi
##Note: Combining into clusters begins right before porechop

