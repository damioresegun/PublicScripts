#!/bin/bash
## Automation script nanopore sequencing analysis. 
## Version 0.1.1 
## Works to take the data through albacore, poretools, porechop, assembly stats, nanoQC, nanoFilt and NanoPlot
## Requires (all future implementations) directory for data, directory to generate outputs
## Requires all tools to be available in the bin folder. If not, enter in configuration file (to be implemented if necessary)
## Future implements: flags?
## NOTE: Will need to change how this is set up. Rather than carrying out each sub step and checking for barcodes. Check for barcodes if yes take through the entire process. If not then take through the entire process. Should make it faster than checking for barcode each and every step
## If barcode pipeline works then might just make two separate files that are sourced into this one
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#configFile=""
#WPATH="" # path to where you want to be the working directory
DPATH="$HOME/5xBarcodingJCS1734/RawData/reads/20171219_1833_barcodex5JCS34_17/downloads/pass/" # path to MinION data.
#DPATH="$HOME/All"
#DPATH="$HOME/LaptopData_preBarcode/RawData/20171016_1550_SKS201a/downloads/pass/" # path to MinION data
SPATH="$HOME/Scripts/Testr/"  # path to save destination
ALBAC="read_fast5_basecaller.py"
PTOC="SQK-RAD002"   # State the protocol used e.g SQK-RAD002
FLOCL="FLO-MIN106" # State the flow cell used
THRD="6" # State the number of threads you wish to use
PORTOL="poretools" # Path to poretools if not in $PATH
PLTFM="y" # Is this run on Marvin?
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
#read PLTFM
echo "Did your protocol involve barcoding your reads? Enter Y or N:"
#read BACODS
#................................................................................................................................................................
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                            BARCODED PIPELINE
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
###################################################################################################################################################################
#                                                                 Albacore Basecalling and Demultiplexing                
###################################################################################################################################################################
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
read_fast5_basecaller.py -f $FLOCL -k $PTOC -o fast5,fastq -i $DPATH -r -s $currPATH/Albacore_Output/NewOutr -t $THRD -q 0
echo
else
## If running on a local PC; it is assumed that the user has got the files in $PATH. If not; paths will be in config file (TO-DO: Make config file)
echo "Running Albacore on your local system..."
echo
## Run albacore on local system
${ALBAC} -f $FLOCL -k $PTOC -o fast5,fastq -i $DPATH -r -s $currPATH/Albacore_Output/ -t $THRD -q 0
#this ran albacore worked on the local
fi
echo
echo "Albacore finished. Moving to the next step... "
###################################################################################################################################################################
#                                                                              Raw MinKNOW Statistics
###################################################################################################################################################################
if [ $BACODS = "Y" ] || [ $BACODS = "y" ] || [ $BACODS = "yes" ] || [ $BACODS = "Yes" ]
then
# If data is barcoded; make directory to hold the individual barcode stats of each stage
mkdir -p $currPATH/Statistics/RawMinKNOWStats
echo
echo "Running pipeline for barcoded MinION data..."
# get the names of the subdirectories present in the path
coda=($DPATH/*/)
#iterate through the path folder to get folder names
for coda in "${coda[@]}"
do
    echo "$coda"
    codar=$(basename "$coda")    
    echo $codar
    if [ $PLTFM = "Y" ] || [ $PLTFM = "y" ] || [ $PLTFM = "yes" ] || [ $PLTFM = "Yes" ]
    then
    echo "Assessing the run statistics for ${codar} on Marvin..."
    echo
    # Carry out stats analysis using poretools and save in the created directory with the barcode as the filename
    poretools stats -q ${coda}/*/ > $currPATH/Statistics/RawMinKNOWStats/${codar}_RawMinKNOWStats.txt
    else
    # if platform is not marvin then use the provided path to poretools and save each barcode in the created directory
    echo "Assessing the run statistics for ${codar} on your local system..."
    ${PORTOL} stats -q ${coda}/*/ > $currPATH/Statistics/RawMinKNOWStats/${codar}_RawMinKNOWStats_local.txt
    fi   # if marvin statement is completed here
done
fi
echo "Poretools complete. The statistics text files are saved in you directory."
echo "Moving to next step..."
###################################################################################################################################################################
#                                                                              Albacore Output Statistics
###################################################################################################################################################################
echo "Beginning Poretools analysis of your basecalled data..."
pwd
# make directory for output
mkdir -p $currPATH/Poretools_Output
cd Poretools_Output
# check if data is barcoded
if [ $BACODS = "Y" ] || [ $BACODS = "y" ] || [ $BACODS = "yes" ] || [ $BACODS = "Yes" ]
then
echo "Running Poretools demultiplexing to fastq on Marvin..."
#### CHANGE SOURCE PATH TO ALBACORE OUTPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#DPATH="$HOME/5xBarcodingJCS1734/RawData/reads/20171219_1833_barcodex5JCS34_17/downloads/pass/"
echo $DPATH
echo $coda
# get the names of the subdirectories present in the path
coda=($currPATH/Albacore_Output/*/)
#iterate through the path folder to get folder names
for barz in "${barz[@]}"
do
    echo "$barz"
    barza=$(basename "$barz")    
    echo $barza
    mkdir -p $barza
    #cd $codar
    poretools fastq ${barz}/*/ > $barza/${barza}.fastq"
    gzip $barz/${barza}.fastq
done
else 
echo "Your data is not barcoded..."
echo
cd $currPATH
pwd
echo "Running Poretools" 

fi
