#$ -cwd
# script to take you through finisher. You need to provide the path to all your contigs folders (i.e. the folder holding all your isolates) and raw reads folder. 
# The script will make a new directory based on your input and then copy the contigs and reads files in there, The files will be renamed and then taken through for finisher

#### RUN IN THE PYTHON2.7 (PYTHON27) CONDA ENVIRONMENT!!!!!!!!!!!!!!!!!!!!!!1

contigsPath="$HOME/All_De_Novo/Qcat"
rawsPath="$HOME/All_Isolate_HumanUnMapped/Qcat"
outputPath="$HOME/All_Finisher/Qcat"
finshPath="$HOME/ToolsAgain/finishingTool"

for i in $contigsPath/* ;
do
fname=$(basename "$i")
echo "You are working on the ${fname} folder"
mkdir -p $outputPath/$fname
cd $outputPath/$fname
pwd
# copy the contigs into the output file
cp ${contigsPath}/${fname}/${fname}_clean.contigs.fasta contigs.fasta
# convert the raw reads into a fasta in the output directory
sed -n '1~4s/^@/>/p;2~4p' ${rawsPath}/${fname}.fastq > raw_reads.fasta
# change the formatting of the raw reads and contigs file to make them easy to run on finisher
perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' raw_reads.fasta > newRaw_reads.fasta
cp newRaw_reads.fasta raw_reads.fasta
perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' contigs.fasta > newContigs.fasta
cp newContigs.fasta contigs.fasta
pwd
python -V
wait
# run the finisher tool
echo "python2 ${finshPath}/finisherSC.py -par 12 ${outputPath}/${fname} /shelf/apps/dro/conda/envs/python27/bin/"

python2 ${finshPath}/finisherSC.py -par 12 ${outputPath}/${fname} $HOME/Tools/MUMmer3.23/
cd
done

