#$ -cwd

# Script for SSPACE long reads. The tools uses BLASR for alignment which has been shown to be consistently lower precision, quality and running time for nanopore reads. So this script is done to test this on our data. Most likely will not be used
#
## Script is set up to take in the clean contigs from All_Denovo folder, the raw reads and place it in the output folder. The parameters stated here are necessary to be kept as is. They will be hard coded in. Only changeable variables are the inputs and output folder. 

contigsPath="$HOME/All_De_Novo/Qcat"
rawsPath="$HOME/All_Isolate_HumanUnMapped/Qcat"
outputPath="$HOME/All_SSPACE/Qcat"
spacePath="$HOME/Tools/SSPACE-LongRead_v1-1"
perlPath="/usr/bin/perl"

for i in $contigsPath/* ;
do
fname=$(basename "$i")
echo "You are working on the ${fname} folder"
mkdir -p $outputPath/$fname
#cd $outputPath/$fname
pwd
echo "perl"
perl -v
echo
echo
echo
echo "perlpath"
$perlPath -v
run=`$perlPath ${spacePath}/SSPACE-LongRead.pl -c ${contigsPath}/${fname}/${fname}_clean.contigs.fasta -p ${rawsPath}/${fname}.fastq -b $outputPath/$fname -t 12 -k 1 -i 70 -a 1500 -g -5000 -l 10`
echo "${fname} is complete"
eval $run
done
echo "All done"
