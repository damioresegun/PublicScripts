#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -r y
#$ -j y
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -N VsPk_FlyePilonBusco
#$ -M dro@st-andrews.ac.uk
##Script made to take input data through busco. This is specifically made to receive data from the output of pilon. PLEASE BE CARE ABOUT INPUTS AND HOW THE BUSCO COMMAND IS CALLED. DOUBLE CHECK!!!!
## In order to run this, it is necessary to have the busco docker image installed. Version 3.0.0 for busco and version 3.2.2 for Augustus. It is also important that the AUGUSTUS_CONFIG_PATH is set here and you are able to write to this folder. So it is advised that the folder should be in a location you have control over rather than a shared folder


set -e
#set the variables needed
ref1="$HOME/Index/PainRefV2.fasta"
ref2="$HOME/Index/PkLappRefSeq.fna"
input="$HOME/QuickPilon/VsPk/Flye"
output="$HOME/All_AssemblyQC/Pilon_BUSCO/VsPk/Flye"
augConfig="$HOME/Tools/augustusConfig/config/"
busDock="$HOME/Tools/buscoDocker/busco-docker_3.0.0.with_r.sif"
lineage="$HOME/Tools/augustusConfig/protists_ensembl/"

THREADS="48"
SPEC="pfalciparum"
#set the location of the config file
export AUGUSTUS_CONFIG_PATH=$augConfig
echo $AUGUSTUS_CONFIG_PATH
#make the output where necessary
mkdir -p $output

echo "Moving on to the isolates"
for f in $input/* ; 
	do
	echo "You are in ${f}" 
	#cd $f
	fname=$(basename "$f")    
    echo $fname
	echo "power dirc is:"
	pwd
		#run busco command
		mkdir -p ${output}/${fname}
		
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"
		
		cd ${output}/${fname}/
		singularity run ${busDock} -i ${input}/${fname} -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		echo "${PREFIX} done successfully"

done

