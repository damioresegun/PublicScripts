#!/bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -N PkAligning
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -M dro@st-andrews.ac.uk
#for aligning and sorting SAM files to BAM files for all sequences within the input folder and to the output folder. Alignment done against the Pk reference. 
# intermediate files will be deleted i.e. SAM files and unsorted BAM files. Sorted BAM files will be renamed to just .bam rather than _sorted.bam
##NOTE: MIGHT HAVE TO RUN THIS IN A SAMTOOLS CONDA ENVIRONMENT
set -e
#set input and output folders
inputFol="$HOME/All_AdapterRemoved"
outSave="$HOME/All_Aligned/VsPk"
fastqSave="$HOME/All_Isolate_PkMapped"
THREADS="32"
statsFol="$HOME/All_Stats/"
fastqForm="$HOME/All_Isolate_PkMapped_Formatted"
# load necessary modules
#module load samtools/1.6
mkdir -p $fastqSave
mkdir -p $statsFol/postPk_Align_Extraction
#gupsta=$(dirname $inputFol)
#cd $gupsta
#state the current working directory
hre=$(pwd)
echo $hre
#load necessary modules
module load assembly-stats/gitv0_003a372
module load python/3.6.4
#make an unmapped folder
mkdir -p $outSave/Unmapped
mkdir -p $outSave/Mapped
unmpd=$outSave/Unmapped
mpd=$outSave/Mapped
#mkdir -p $unmpd/VsLapp
#lapr=$unmpd/VsLapp
for f in ${inputFol}/*
do
#get the filename
file=$(basename $f)
echo $file
#remove the extension
fname="${file%.*}"
echo $fname
echo
#set the alignment name to the path and the file
alname=${outSave}/${fname}
echo $alname
#align to Pk reference
minimap2 -ax map-ont $HOME/Index/PkLappRefSeq.fna.mmi $inputFol/${fname}.fastq > ${alname}VsPk.sam -t $THREADS
# now convert to bam, sort and delete intermediate files
echo "alignment done"
samtools view --threads $THREADS -bS ${alname}VsPk.sam > ${alname}VsPk.bam
echo
echo "bam conversion done"
samtools sort --threads $THREADS -o ${alname}VsPk_sorted.bam ${alname}VsPk.bam 
rm ${alname}VsPk.sam
rm ${alname}VsPk.bam
mv ${alname}VsPk_sorted.bam ${alname}VsPk.bam
echo ${fname} >> All_FlagstatMappedVsPk_stats.txt
samtools flagstat --threads $THREADS ${alname}VsPk.bam >> ${statsFol}/postAlign_Extraction/All_FlagstatMappedVsPk_stats.txt
samtools view --threads $THREADS -F 0x04 -b ${alname}VsPk.bam > ${mpd}/${fname}VsPk_mapped.bam
bedtools bamtofastq -i ${mpd}/${fname}VsPk_mapped.bam -fq ${mpd}/${fname}VsPk_mapped.fastq
echo "Converted to fastq"
# extract unmapped reads
cd $unmpd
pwd
samtools view --threads $THREADS -f 4 -b ${alname}VsPk.bam > ${unmpd}/${fname}VsPk_unmapped.bam
echo "unmapped extraction done"
#convert unmapped reads to fastq
bedtools bamtofastq -i ${unmpd}/${fname}VsPk_unmapped.bam -fq ${unmpd}/${fname}VsPk_unmapped.fastq
echo
#check the stats
assembly-stats ${mpd}/${fname}VsPk_mapped.fastq > ${statsFol}/postPk_Align_Extraction/${fname}_VsPkMappedStats.txt

assembly-stats ${unmpd}/${fname}VsPk_unmapped.fastq > ${statsFol}/postPk_Align_Extraction/${fname}_VsPkUnmappedStats.txt
#move fastqs to their own folder
mv ${mpd}/${fname}VsPk_mapped.fastq $fastqSave
echo "all done"
#cd $hre
done
# format the fastq
~/Scripts/fastqformatter.sh ${mpd} ${fastqForm}


cat ${statsFol}/postPk_Align_Extraction/*_VsPkMappedStats.txt >> ${statsFol}/postPk_Align_Extraction/All_VsPkMappedStats.txt
cat ${statsFol}/postPk_Align_Extraction/*_VsPkUnmappedStats.txt >> ${statsFol}/postPk_Align_Extraction/All_VsPkUnmappedStats.txt
cat ${statsFol}/postPk_Align_Extraction/All_VsPkMappedStats.txt | tr -s '[:blank:]' ',' > ${statsFol}/postPk_Align_Extraction/All_VsPkMappedStats.csv
cat ${statsFol}/postPk_Align_Extraction/All_VsPkUnmappedStats.txt | tr -s '[:blank:]' ',' > ${statsFol}/postPk_Align_Extraction/All_VsPkUnmappedStats.csv
rm ${statsFol}/postPk_Align_Extraction/All_VsPkUnmappedStats.txt
exit

