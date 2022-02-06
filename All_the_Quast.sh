#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -r y
#$ -j y
#$ -N VsPk_FlyeQuastGenomes
#$ -q all.q
#$ -l hostname="node8"
#$ -M dro@st-andrews.ac.uk
#Has to be run in the QUAST CONDA env
set -e
#exec 2>&1 | tee VsPk_Flye_Quastlog.txt
Assem="$HOME/All_CompanionGenomes/VsPk/Flye"
output="$HOME/AssemblyQc/Genomes_Quast/VsPk/Flye"
lappRef="$HOME/Index/PkLappRefSeq.fna"
lapGf="$HOME/Index/PkLappRefSeq.gff"
THREADS="12"
cd $Assem

#cd /storage/home/users/dro/Assembly_QC/sks339/newAssem/quast/VsErnest/

for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	fname="${fna%*.*}"
	#fname="${fnamer%*_*.*}"
	
    echo $fname
    #cd $folder
    #mkdir -p $HOME/Assembly_QC/Quast/Guppy/${fname}
	#mkdir -p $output/${fname}/Raw
	#mkdir -p $output/${fname}/Clean
    pwd
    #cd ${Assem}/${fname}
    #echo "Running quast Vs Ernest Reference..."
    #/shelf/apps/dro/conda/envs/quast/lib/python3.6/site-packages/quast-5.0.2-py3.6.egg-info/scripts/quast.py -t 12 /storage/home/users/dro/All_De_Novo/newAssem/sks339/sks339.contigs.fasta -o /storage/home/users/dro/Assembly_QC/newAssemblyRun/${fname}/quast/VsErnest/ -r /storage/home/users/dro/Index/PkErnestRefSeq.fasta -g /storage/home/users/dro/Index/PkErnestRefSeq.gff3 -e -f
    #echo "Quast vs Ernest Reference done..."
    echo
    echo "Running Quast on Lapp reference..."
    echo 
    #/shelf/apps/dro/conda/envs/quast/lib/python3.6/site-packages/quast-5.0.2-py3.6.egg-info/scripts/quast.py -t 12 /storage/home/users/dro/All_De_Novo/newAssem/sks339/sks339.contigs.fasta -o /storage/home/users/dro/Assembly_QC/newAssemblyRun/${fname}/quast/VsLapp/ -r /storage/home/users/dro/Index/PkLappRefSeq.fna -g /storage/home/users/dro/Index/PkErnestRefSeq.gff3 -e -f
    
	##real quast command
	##generic
	quast -t $THREADS ${Assem}/${fname}.fasta -o ${output}/${fname}/ -r ${lappRef} -g ${lapGf} -e -f --circos
		
	#canu
	###raw assem
	#quast -t $THREADS ${Assem}/${fname}/${fname}.contigs.fasta -o ${output}/${fname}/Raw/ -r ${lappRef} -g ${lapGf} -e -f
	###clean
	#quast -t $THREADS ${Assem}/${fname}/${fname}_clean.contigs.fasta -o ${output}/${fname}/Clean/ -r ${lappRef} -g ${lapGf} -e -f
	
	
	
	#flye
	###raw
	#quast -t $THREADS ${Assem}/${fname}/assembly.fasta -o ${output}/${fname}/Raw/ -r ${lappRef} -g ${lapGf} -e -f
	###clean
	#quast -t $THREADS ${Assem}/${fname}/clean_assembly.fasta -o ${output}/${fname}/Clean/ -r ${lappRef} -g ${lapGf} -e -f
	
	
    cd $Assem
	
	### quast command for companion outputs
	#quast -t 16 $Assem/${fname}*.fasta -o $output/${fname} -r $lappRef -g $lapGf -e -f 
done