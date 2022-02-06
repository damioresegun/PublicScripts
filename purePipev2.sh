#!/bin/bash
## Automation
## Automation script for nanopore sequencing data analysis
## Version 0.1.2
## Takes input data through the albacore (and sf5 albacore for multi-read fast5 inputs), porechop, asseembly-stats, the nanopack suite of tools, minimap2 and then further downstream analysis of the alignment file
## Future implementations -- Requires directory for data, directory to generate outputs, configuration file, Flags?
## Requirements: All tools need to be in available in their bin folder and/or on the PATH. If not, enter absolute path in the configuration file
###############################################################################################
##                                     Requirements
#
###############################################################################################
# The path to the config file if the config file flag is indicated. Config file will be in the package folder
#configFile" == "
# Path to the working directory
WPATH=""
# Path to the MinION data. This can be raw MinKNOW output or already MinKNOW basecalled
DPATH=""
# Path to save outputs of the script
SPATH="Testing/amazo"
# Which albacore function to use. 1= Albacore 2= sf5 albacore. Both need python 3.6.4
whichAlba="1" 
ALBAC="read_fast5_basecaller.py"
SF5ALBAC="sf5_read_fast5_basecaller.py"
# Which protocol was used for the sequencing run
PTOC="SQK-RBK004"
# Number of threads needed
THRD="14"
# Porechop command to remove adapters. Has to be called in python 3.6.4
PRCHP="porechop"
# Was this a barcoding run?
BACD="y"
# Stop execution of script if an error is detected
set -e 

# Prior checks for save destination
if [ -d "$SPATH" ]; then
	echo "${SPATH} exists. Subsequent folders will be made here"
else
	echo "{SPATH} does not exist. The folder and parent folders in the path will be made to hold subsequent folders and data"
mkdir -p $SPATH
echo "The save folder and parent folders [${SPATH}] have been created"
fi
echo
echo
##################################################################################################
#
#                    PIPELINE - ALBACORE BASECALLING AND DEMULTIPLEXING
#
##################################################################################################
# Change into the save folder and show that it is the current working directory
cd $SPATH
currPATH=`pwd`
echo "Your working directory is now ${currPATH}"
echo
echo "Starting Albacore basecalling..."
# Make directory for albacore outputs
mkdir -p $currPATH/Albacore_Output/
# Load the module needed and display the module is loaded
module load python/3.6.4
python3 -V
echo
echo "The python module is loaded. Albacore will begin now ..."
echo
# Check if the sequencing run was barcoded
if [ "$BACD" == "Y" ] || [ "$BACD" == "y" ] || [ "$BACD" == "YES" ] || [ "$BACD" == "yes" ] || [ "$BACD" == "Yes" ]
then
# Confirm that albacore should be run. Necessary due to the high memory usage required
	read -p "Albacore will run on barcoded sequences... Continue? [y or n] `echo $'\n> '`" albaGo
	echo
	#read $albaGO
	if [ "$albaGO" == "Yes" ] || [ "$albaGo" == "y" ]
	then
		echo "Beginning..."
		echo $whichAlba
		if [ "$whichAlba" == "1" ]
		then
			echo "read_fast5_basecaller.py -f $FLOCL -k $PTOC --barcoding -o fast5,fastq -i $DPATH -r -s $currPATH/Albacore_Output/ -t $THRD -q 0"
			echo
			echo "Albacore demultiplexing complete. Proceeding to next step..."
			echo
		else
			if [ "$whichAlba" == "2" ]
			then
				echo "sf5_read_fast5_basecaller.py -f $FLOCL -k $PTOC --barcoding -o fast5,fastq -i $DPATH -r -s $currPATH/Albacore_Output/ -t $THRD -q 0"
				echo
				echo "Sf5 Albacore demultiplexing complete. Proceeding to the next step..."	
			else 
				echo "No Albacore protocol was chosen. Exiting..."
				exit
		#fi
	#else
		if [ "$albaGO" == "n" ] || [ "$albaGO" == "No" ] || [ "$albaGO" == "no" ]
		then	
			read -p "Albacore will not proceed. Continue to next step [y]? Or exit [n]... `echo $'\n> '`" skipAlba
	#read $skipAlba
				if [ "$skipAlba" == "n" ] || [ "$skipAlba" == "no" ] || [ "$skipAlba" == "No" ]
				then
					echo "Script ended, Have a good day!"
				else
					echo "Albacore has been skipped. Proceeding to the next step..."
				fi
		else
			echo "No choice given. Exiting..."
		fi

	fi
fi
fi
fi
echo "done"
