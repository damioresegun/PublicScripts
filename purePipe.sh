#!/bin/bash 
## Automation script nanopore sequencing analysis. 
## Version 0.1.1 
## Works to take the data through the pipeline from albacore, through poretools, porechop, assembly stats, the nanopack suite, minimap2 and then further downstream filtering
## (Future implementations) Requires directory for data, directory to generate outputs
## Requires all tools to be available in the bin folder. If not, enter in configuration file (to be implemented if necessary)
## Future implements: flags?
#configFile=""
#WPATH="" # path to where you want to be the working directory
#DPATH="$HOME/5xBarcodingJCS1734/RawData/reads/20171219_1833_barcodex5JCS34_17/downloads/pass/" # path to MinION data. Path here is to the barcoded
DPATH="$HOME/WBC_Depleted_Clinical_MinION/NewRaw/reads/20181115_1455_SKS/fast5/pass/"
#DPATH="$HOME/LaptopData_preBarcode/RawData/20171016_1550_SKS201a/downloads/pass/" # path to MinION data. Here is the 201a data
SPATH="$HOME/OldClinicals/NewAlbacore"  # path to save destination
#SPATH="$HOME/Scripts/Testr/"  # path to save destination
ALBAC="read_fast5_basecaller.py"  # Which albacore function to use
#PTOC="SQK-RAD002"   # State the protocol used e.g SQK-RAD002
PTOC="SQK-RBK004"
#PTOC="SQK-RBK001"
FLOCL="FLO-MIN106" # State the flow cell used
THRD="16" # State the number of threads you wish to use
PORTOL="poretools" # Path to poretools if not in $PATH
#PLTFM="y" # Is this run on Marvin?
FINAME="All_Test" # Give a name for this run; all files will be saved with this prefix
BACODS="y" # Was this a barcoding run?
set -e
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                            SET TOOLS, PATHS, DOWNLOAD DATA
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################################################################################################################
## Get the path to save all results in a created folder
if [ -d "$SPATH" ]; then
	echo "${SPATH} exists. Subsequent folders will be made here"
else
echo "${SPATH} does not exist. Folder will be made to hold subsequent folders and data"
mkdir -p $SPATH
echo "Made ${SPATH}"
fi
currPATH="$SPATH"
cd $currPATH
pwd
#................................................................................................................................................................
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                            PIPELINE
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
###################################################################################################################################################################
#                                                                 Albacore Basecalling and Demultiplexing                
###################################################################################################################################################################
echo "Starting Albacore Basecalling..."
## Make directory for albacore output
mkdir -p $currPATH/Albacore_Output/
## Load module required
module load python/3.6.4
## Display the module loaded
python3 -V
echo "Your python module is loaded"
## Based on answer of previous question; carry out albacore
if [ $BACODS = "Y" ] || [ $BACODS = "y" ] || [ $BACODS = "yes" ] || [ $BACODS = "Yes" ]
then
echo "Running Albacore on barcoded samples ..."
echo
## Run albacore on barcoded samples with the barcoding flag
echo "Are you sure you want to do Albacore? Y or N..."
read $albaGO
echo "Beginning.."
read_fast5_basecaller.py -f $FLOCL -k $PTOC --barcoding -o fast5,fastq -i $DPATH -r -s $currPATH/Albacore_Output/ -t $THRD -q 0
echo
echo "Albacore demultiplexing complete. Proceeding to next step..."
echo
else
echo "Running Albacore on non-barcoded samples ..."
## Run albacore on non-barcoded samples with the batching disabled
echo "Albacore will be run on non barcoded pipeline. Proceed?: "
read $nonBarcAlbaGo
read_fast5_basecaller.py -f $FLOCL -k $PTOC -o fast5,fastq -i $DPATH -r -s $currPATH/Albacore_Output/ -t $THRD -q 0 -n 0
echo
echo "Albacore basecalling complete. Proceeding to next step..."
fi
echo
module unload python/3.6.4
exit
###################################################################################################################################################################
#                                                                              Albacore Output Statistics
###################################################################################################################################################################
pwd
echo "Beginning Poretools statistical analysis of your basecalled data..."
module load python/2.7
mkdir -p $currPATH/Statistics/Albacore_Stats
if [ $BACODS = "Y" ] || [ $BACODS = "y" ] || [ $BACODS = "yes" ] || [ $BACODS = "Yes" ]
then
# If data is barcoded; make directory to hold the individual barcode stats of each stage
echo "Running pipeline for barcoded MinION data..."
# get the names of the subdirectories present in the path
coda=($DPATH/*/)
echo "${coda}"
#iterate through the path folder to get folder names
for coda in "${coda[@]}"
do
    echo "$coda"
    codar=$(basename "$coda")    
    echo $codar
    #if [ $PLTFM = "Y" ] || [ $PLTFM = "y" ] || [ $PLTFM = "yes" ] || [ $PLTFM = "Yes" ]
    #then
    echo "Assessing the run statistics for ${codar} on Marvin..."
    echo
    # Carry out stats analysis using poretools and save in the created directory with the barcode as the filename
    echo "poretools stats -q ${coda}/*/ > $currPATH/Statistics/RawMinKNOWStats/${codar}_RawMinKNOWStats.txt"
done
else
    # if data is not barcoded, then use the name of the sample to make folder and carry on poretools
    echo "Assessing the run statistics for ${codar} on your local system..."
    poretools stats $DPATH > $currPATH/Statistics/RawMinKNOWStats/${FINAME}_RawMinKNOWStats.txt
    #cat "${PORTOL} stats -q ${coda}/*/ > $currPATH/Statistics/RawMinKNOWStats/${codar}_RawMinKNOWStats_local.txt"
    #fi   # if statement is completed here
fi
echo "Poretools complete. The statistics text files are saved in you directory."
echo "Moving to next step..."

#echo "Beginning Poretools analysis of your basecalled data..."
#pwd
# make directory for output
#mkdir -p $currPATH/Poretools_Output
#cd Poretools_Output
# check if data is barcoded
#if [ $BACODS = "Y" ] || [ $BACODS = "y" ] || [ $BACODS = "yes" ] || [ $BACODS = "Yes" ]
#then
#echo "Running Poretools demultiplexing to fastq on Marvin..."
#### CHANGE SOURCE PATH TO ALBACORE OUTPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#DPATH="$HOME/5xBarcodingJCS1734/RawData/reads/20171219_1833_barcodex5JCS34_17/downloads/pass/"
#echo $DPATH
#echo $coda
# get the names of the subdirectories present in the path
#coda=($currPATH/Albacore_Output/*/)
#iterate through the path folder to get folder names
#for barz in "${barz[@]}"
#do
#    echo "$barz"
  #  barza=$(basename "$barz")    
 #   echo $barza
   # mkdir -p $barza
    #cd $codar
    #poretools fastq ${barz}/*/ > $barza/${barza}.fastq"
    #gzip $barz/${barza}.fastq
#done
#else 
#echo "Your data is not barcoded..."
#echo
#cd $currPATH
#pwd
#echo "Running Poretools" 

#fi

exit
