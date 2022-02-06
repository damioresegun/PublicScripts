#!/bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -r y
#$ -j y
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -N Chromosomes_Align
#$ -M dro@st-andrews.ac.uk
#for aligning and sorting SAM files to BAM files for all sequences within the input folder and to the output folder
# intermediate files will be deleted i.e. SAM files and unsorted BAM files. Sorted BAM files will be renamed to just .bam rather than _sorted.bam
##NOTE: MIGHT HAVE TO RUN THIS IN A SAMTOOLS CONDA ENVIRONMENT
set -e
#set input and output folders
inputFol="$HOME/FinalGenomes"
outSave="$HOME/Genome_Chromosomes_Aligned"
THREADS="48"
statsFol="$HOME/All_AssemblyStats"
# load necessary modules
mkdir -p $statsFol/Chromosome_Aligned
#state the current working directory
hre=$(pwd)
echo $hre
#load necessary modules
module load assembly-stats/gitv0_003a372
#make an unmapped folder
mkdir -p $outSave/Unmapped
mkdir -p $outSave/Mapped
unmpd=$outSave/Unmapped
mpd=$outSave/Mapped
for f in ${inputFol}/*
do
#get the filename
echo $f
filer=$(basename $f)
echo "this is $filer"
#remove the extension
file="${filer%.*}"
echo $file
echo $f
echo
#set the alignment name to the path and the file
alname=${outSave}/${file}
echo $alname
for iN in ${f}/*
do
for u in {01..14};
do
echo "in here ${iN}"
newn=$(basename "$iN")
echo $newn
noona="${newn%.*}"
echo $noona

#align to Pain reference
echo "minimap2 -ax map-ont $HOME/Index/PkPainRefSeq.fna.mmi ${iN} > ${alname}/${noona}VsRef.sam -t $THREADS"
done
done
exit
exit
# now convert to bam, sort and delete intermediate files
echo "alignment done"
samtools view --threads $THREADS -bS ${alname}/${noona}VsRef.sam > ${alname}/${noona}VsRef.bam
echo
echo "bam conversion done"
samtools sort --threads $THREADS -o ${alname}/${noona}VsRef_sorted.bam ${alname}/${noona}VsRef.bam 
#rm ${alname}VsRef.sam
rm ${alname}/${noona}VsRef.bam
mv ${alname}/${noona}VsRef_sorted.bam ${alname}/${noona}VsRef.bam
# index the bam file
samtools index ${alname}/${noona}VsRef.bam
# flagstat
echo ${noona} >> ${statsFol}/Chromosome_Aligned/Chromosomes_FlagstatMappedVsRef_stats.txt
samtools flagstat --threads $THREADS ${alname}/${noona}VsRef.bam >> ${statsFol}/Chromosome_Aligned/Chromosomes_FlagstatMappedVsRef_stats.txt
# extract mapped reads
samtools view --threads $THREADS -F 0x04 -b ${alname}/${noona}VsRef.bam > ${mpd}/${noona}VsRef_mapped.bam
echo "mapped extraction done"
bedtools bamtofastq -i ${mpd}/${noona}VsRef_mapped.bam -fq ${mpd}/${noona}VsRef_mapped.fastq
# extract unmapped reads
#cd $unmpd
#pwd
samtools view --threads $THREADS -f 4 -b ${alname}/${noona}VsRef.bam > ${unmpd}/${noona}VsRef_unmapped.bam
echo "unmapped extraction done"
#convert unmapped reads to fastq
bedtools bamtofastq -i ${unmpd}/${noona}VsRef_unmapped.bam -fq ${unmpd}/${noona}VsRef_unmapped.fastq
echo
#check the stats
#assembly-stats ${unmpd}/${fname}VsRef_unmapped.fastq > ${statsFol}/postAlign_Extraction/${fname}_VsRefUnmappedStats.txt
# get coverage information
cvg=$outSave/Coverage
mkdir -p $outSave/Coverage
bedtools genomecov -ibam ${alname}/${noona}VsRef.bam > $outSave/Coverage/${noona}VsRef_wholeCov.txt
bedtools genomecov -ibam ${alname}/${noona}VsRef.bam -bga > $outSave/Coverage/${noona}VsRef_bgaCov.txt
#bedtools genomecov -ibam ${alname}VsRef.bam -d > $outSave/Coverage/${file}VsRef_perBaseCov.txt
echo "all done"
#cd $hre
done
done
exit

