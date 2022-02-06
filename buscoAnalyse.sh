#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -r y
#$ -j y
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -N VsPk_FlyeGenomeBusco
#$ -M dro@st-andrews.ac.uk
##Script made to take input data through busco. This is specifically made to receive data from the output of the wtdbg2 parameter search. It is likely that it will need to be adapted to other use-cases
## In order to run this, it is necessary to have the busco docker image installed. Version 3.0.0 for busco and version 3.2.2 for Augustus. It is also important that the AUGUSTUS_CONFIG_PATH is set here and you are able to write to this folder. So it is advised that the folder should be in a location you have control over rather than a shared folder
## layout of the script: for every isolate in the wtdbg2 output, run busco, move the run output folder to the output folder, delete the tmp folder and move to the next isolate

set -e
#set the variables needed
ref1="$HOME/Index/PainRefV2.fasta"
ref2="$HOME/Index/PkLappRefSeq.fna"
input="$HOME/All_CompanionGenomes/VsPk/Flye"
output="$HOME/All_AssemblyQC/Genomes_BUSCO/VsPk/Flye"
augConfig="$HOME/Tools/augustusConfig/config/"
busDock="$HOME/Tools/buscoDocker/busco-docker_3.0.0.with_r.sif"
lineage="$HOME/Tools/augustusConfig/protists_ensembl/"

THREADS="12"
SPEC="pfalciparum"
#set the location of the config file
export AUGUSTUS_CONFIG_PATH=$augConfig
echo $AUGUSTUS_CONFIG_PATH
#make the output where necessary
mkdir -p $output

#run busco for the reference files to see how our data match up
#mkdir -p $HOME/All_AssemblyQC/BUSCO/References
#cd $HOME/All_AssemblyQC/BUSCO/References
#singularity run ${busDock} -i $ref1 -l ${lineage} -o PkPain -f -m geno -c $THREADS -sp ${SPEC} --long
#echo "Reference 1 done!"
#singularity run ${busDock} -i $ref2 -l ${lineage} -o PkLapp -f -m geno -c $THREADS -sp ${SPEC} --long
#echo
#echo "Reference 2 done!"
echo "Moving on to the isolates"
for f in $input/* ; 
	do
	echo "You are in ${f}" 
	#cd $f
	fna=$(basename "$f")
	fname=${fna%*.*}
    echo $fname
	echo "power dirc is:"
	pwd
		#run busco command
		mkdir -p ${output}/${fname}
		#mkdir -p ${output}/${fname}/Clean
		#mkdir -p ${output}/${fname}/Raw
		
		PREFIX="$fname"
		echo $PREFIX
		echo "current working directory $(pwd)"
		
		#canu -- normal
		#cd ${output}/${fname}
		#singularity run ${busDock} -i ${input}/${fname}/${fname}.contigs.fasta -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		
		
		#canu raw
		#cd ${output}/${fname}/Raw
		#singularity run ${busDock} -i ${input}/${fname}/${fname}.contigs.fasta -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		
		
		#canu clean
		#cd ${output}/${fname}/Clean
		#singularity run ${busDock} -i ${input}/${fname}/${fname}_clean.contigs.fasta -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		
		#flye -- normal
		#cd ${output}/${fname}/
		#singularity run ${busDock} -i ${input}/${fname}/assembly.fasta -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		
		#flye raw
		#cd ${output}/${fname}/Raw
		#singularity run ${busDock} -i ${input}/${fname}/assembly.fasta -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		
		
		# flye clean
		#cd ${output}/${fname}/Clean
		#singularity run ${busDock} -i ${input}/${fname}/clean_assembly.fasta -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		
		
		#####GENOME 
		cd ${output}/${fname}/
		singularity run ${busDock} -i ${input}/${fname}.fasta -l ${lineage} -o $PREFIX -f -m geno -c $THREADS -sp ${SPEC} --long -z
		
		echo "${PREFIX} done successfully"
		
done
#run_busco -i ~/All_De_Novo/sks339/sks339.contigs.fasta -o ./ -l ~/Index/busco_eukaryota_odb9/dataset.cfg -m geno -c 12
