#!/bin/bash
#$ -m e
#$ -V
#$ -r y
#$ -j y
#$ -N VsPk_Flye_Decontam
## add $ after first # to activate these
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -M dro@st-andrews.ac.uk

#Script to take data through part of decontamination after blobtools. Should go through the blobDB text file to keep only entries that have Apicomplexa, no-hit, undef and # thus theoretically keeping only the phylum for Apicomplexa and other undefined sequences. After this, the script will get the names of the remaining sequences and put them in a list to then use to filter the contigs file. The resulting contig file will then be taken through blobtools again for confirmation of decontamination before being taken forward for polishing
#inputFol="$HOME/All_BlobTools/Canu"
inputFol="$HOME/All_BlobTools/VsPk/Flye/"
#path to the assembly fassta file
fastaPath="$HOME/All_De_Novo/VsPk/Flye/"

module load python/3.6.4

for i in $inputFol/* ;
do
ls $i
cd $i
pwd

# take out the lines in the table that are not the header, Apicomplexa, no-hits and undefined
grep -E 'Apicomplexa|#|no-hit|undef' blobDB.table.txt > clean.blobDB.table.txt

#make a list of the contigs we want to keep in our fasta
cut -f 1,1 clean.blobDB.table.txt > nodes.txt

# take off the header in the nodes list
grep -v '^#' nodes.txt > node_names.txt

# remove the first nodes files
rm nodes.txt
cd
done
echo "Nodes cleaned. Moving on to cleaning the FASTA..."

for fas in $inputFol/* ;
do
ls $fas
fname=$(basename "$fas")
echo $fname
# call the fasta tool for the cleaning of the fasta
#canu
#python3 $HOME/Scripts/FastaTool.py ${fastaPath}/${fname}/${fname}.contigs.fasta ${inputFol}/${fname}/node_names.txt > ${fastaPath}/${fname}/${fname}_clean.contigs.fasta
#flye default
python3 $HOME/Scripts/FastaTool.py ${fastaPath}/${fname}/assembly.fasta ${inputFol}/${fname}/node_names.txt > ${fastaPath}/${fname}/clean_assembly.fasta
done
echo "All cleaning done. Felicitations!"

