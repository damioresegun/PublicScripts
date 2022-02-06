#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -r y
#$ -j y
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -N VsPk_FlyeBusco_Medak
#$ -M dro@st-andrews.ac.uk
##Script made to take data that has been polished by medaka through busco.
## In order to run this, it is necessary to have the busco docker image installed. Version 3.0.0 for busco and version 3.2.2 for Augustus. It is also important that the AUGUSTUS_CONFIG_PATH is set here and you are able to write to this folder. So it is advised that the folder should be in a location you have control over rather than a shared folder
## layout of the script: for every isolate in the input folder, run busco, move the run output folder to the output folder, delete the tmp folder and move to the next isolate

set -e
#set the variables needed

input="$HOME/All_Medaka/VsPk/Canu"
output="$HOME/All_AssemblyQC/Medaka_BUSCO/VsPk/Canu"
augConfig="$HOME/Tools/augustusConfig/config/"
busDock="$HOME/Tools/buscoDocker/busco-docker_3.0.0.with_r.sif"
lineage="$HOME/Tools/augustusConfig/protists_ensembl/"

THREADS="16"
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
		# busco command
		cd ${output}/${fname}
		singularity run ${busDock} -i ${input}/${fname}/${fname}_consensus.fasta -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		echo "${PREFIX} done successfully"

done
#run_busco -i ~/All_De_Novo/sks339/sks339.contigs.fasta -o ./ -l ~/Index/busco_eukaryota_odb9/dataset.cfg -m geno -c 12

