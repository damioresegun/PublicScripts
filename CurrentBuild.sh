#!/bin/bash
## This script will take you through the mapping of your data to determine the visualise it against the reference genome data.
## Requires: the config file containing the different paths required within this script
copycatpro=$1
source $copycatpro
echo "$RPATH"
echo $startPos
echo $copyCat
echo $configFile
exit
if(($csvFile == "blank"))
then
echo "csvFile is blank"
fi
#if (($downloadSAMTOOLS = "Yes"))
#then
#echo "This is crap $downloadSAMTOOLS"
#wget $samLINK -P $toolPATH
#eval $great
#fi
#if (($downloadSAMTOOLS = "No"))
#then
#echo "Or else moron $downloadSAMTOOLS"
#fi
#echo $downloadSAMTOOLS
#echo $samLINK
 
exit
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                            SET TOOLS, PATHS, DOWNLOAD DATA
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################################################################################################################
#................................................................................................................................................................
## Get the path to save all results in a created folder
echo -n "Please enter the path for you wish to save your results. Please note that a folder name 'Mapping' will be made to hold your files. Press [ENTER] when finished:"
    read SPATH
		echo "The path you have entered is $SPATH"
    echo
    ## Make the folder to hold all results
    mkdir -p "./Mapping"

## Set the working directory to the created folder
currPATH="$SPATH/Mapping"
echo
#mkdir -p $currPATH/test
cd $currPATH
echo
echo "Working directory is:"
pwd
echo "You are now in the Mapping folder $currPATH"
#................................................................................................................................................................
## Get the full path to the reference file
echo
echo "Please enter the path to your reference genome file:"
#read RPATH
echo
echo "The path you have entered for your reference genome is $RPATH"
## Get the filename of the reference file
fullRefName=$(basename "$RPATH")
refName="${fullRefName%.*}"
echo
echo "The name of the reference file is $refName"
#................................................................................................................................................................
## Get the path to the data being tested
echo
echo "Please enter the path to your dataset csv file"
#read csvFile
#................................................................................................................................................................
echo "Please enter the path to the folder where your dataset are saved:"
#read DPATH
echo
#................................................................................................................................................................
### Extract the first column from the csv file
cut -f 1 -d ',' $csvFile > accession_code.txt
### Skip the first line of the column
#tail -n +2 accession_code.txt > codes.txt
### download the data from the website
#while IFS=, read -r line || [[ -n "$line" ]]; do
sed 1d accession_code.txt | tr -d '\015' | while read line; do
    ##get each line from the text file
    echo "Text read from file: $line"
    ## get the first six letters of the 
    #echo "First six is: ${line:0:6}"
    preName=${line:0:6}
    url1="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$preName/$line/${line}_1.fastq.gz"
    url2="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$preName/$line/${line}_2.fastq.gz"
    urlcmd1="wget $url1 -P $DPATH"
    #eval $urlcmd1
    urlcmd2="wget $url2 -P $DPATH"
    #eval $urlcmd2  
done
echo "Your datasets are saved in $DPATH"

#...............................................................................................................................................................
echo "Enter your start position"
#read startPos
echo "Enter your stop position"
#read stopPos
echo "Enter the chromosome of interest e.g. Pf3D7_12_v3"
#echo chrDes
#................................................................................................................................................................
echo
echo "Please enter the path to your samtools package executable:"
#read SAMTOOLS
#................................................................................................................................................................
echo
#................................................................................................................................................................
echo
echo "Please enter the path to your bedtools package executable:"
#read BEDPATH
#................................................................................................................................................................
echo
#................................................................................................................................................................
echo
echo "Please enter the path to your delly package executable:"
#read dPATH
#................................................................................................................................................................
echo
# Set FreeC path
echo "Enter your FreeC executable path"
#read FPATH
## Parent directory of freeC
FsPATH=$(dirname $FPATH)
echo $FsPATH
#................................................................................................................................................................
## Parent directory of 
# Get the paths to the bcftools and vcfutils from the provided samtools path
samDir=$(dirname $SAMTOOLS)
echo $samDir
bcfTools=$samDir/bcftools/bcftools
echo $bcfTools
vcfScript=$samDir/bcftools/vcfutils.pl
echo $vcfScript
echo
#................................................................................................................................................................
#Make an index directory to hold index files
mkdir -p "./Index"
indexPATH="./Index"
cd $indexPATH
#................................................................................................................................................................
# copy reference file into current directory
copText="cp $RPATH ./"
#echo $copText > copy.txt
eval $copText
rPATH=`readlink -f ./$fullRefName`
#echo $rPATH > new.txt
echo $rPATH
echo "The new path to the reference genome is: $rPATH"
echo
#................................................................................................................................................................
# Carry out bwa and samtools indexing
indexer="bwa index $rPATH" 
samIndexer="$SAMTOOLS faidx $rPATH"
echo $indexer > index.txt
echo $samIndexer > samIndex.txt
eval $indexer
eval $samIndexer
#................................................................................................................................................................
#### Split the reference into its chromosomes #################
splitr=`perl ~/split.pl $rPATH`
echo $splitr
mkdir -p "./chromosomes"
mv ./*.fa ./chromosomes
cPATH=`readlink -f ./chromosomes`
echo $cPATH
#................................................................................................................................................................
# Get the length of each chromosome and save in a file
for g in $(find $cPATH -name *.fa); do
  echo $g
  countr=`cat $g | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' >> Pf.fasta.len`
  eval $countr
done
lenFILE=`readlink -f ./Pf.fasta.len`
echo "The length file is $lenFILE"
#................................................................................................................................................................
# Check the start and Stop positions entered by user
# Start Position
if (($startPos-10000 <= 0))
then
echo "Starting position entered is less than length of chromosome. Starting position set to 0"
startPos="0"
else
echo "Starting position remains $startPos"
fi
#................................................................................................................................................................
# Stop position
# extract length of chromosome from file
grep $chrDes $lenFILE > chrStuff.txt
chr12="chrStuff.txt"
fa6t=`cut -f 2 $chr12`
echo $fa6t
rm $chr12
if (($stopPos+10000 >= $fat6t))
then
echo "Stop position exceeds length of chromosome. Stop position set to $fa6t"
stopPos=$fa6t
echo $stopPos
else
echo "Stop position remains $stopPos"
fi
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
cd ..
echo "Working is:"
nPATH=`pwd`
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                             MAPPING
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
## For each downloaded data file:
for fastq1 in $(find $DPATH -name *1.fastq.gz); do
	echo $fastq1
 ## Extract the data filename to name the folder.
	folderPATH="${fastq1%_*.*}"
	echo $folderPATH
        fileName=$(basename "$folderPATH")
        echo $fileName
        
        ## make a folder for each datafile
	      mkdir -p "./${fileName}_genome"
        aPATH=`readlink -f ./${fileName}_genome`
       ## Change into each recurring folder to carry out mapping
             cd ./${fileName}_genome
             
       ##second file
             fastq2=$(echo ${fastq1} | sed s/_1\.fastq\.gz/_2\.fastq\.gz/g)
             echo $fastq2
#................................................................................................................................................................  
       ## Carry out mapping using bwa
             #testingCall="bwa mem $rPATH ${folderPATH}_1.fastq.gz ${folderPATH}_2.fastq.gz > ${fileName}.sam";
             testingCall="bwa mem -t 70 $rPATH $fastq1 $fastq2 > ${fileName}.sam";
             #echo $testingCall > command_prep.txt
             eval $testingCall
#................................................................................................................................................................  
       ## Call samtools
             samCall="$SAMTOOLS import ${rPATH}.fasta.fai ${fileName}.sam ${fileName}.bam ; $SAMTOOLS sort ${fileName}.bam ${fileName}.sorted ; $SAMTOOLS index ${fileName}.sorted.bam ; $SAMTOOLS mpileup -q 10 -Q 23 -d 2000 -C 50 -ugf $rPATH ${fileName}.sorted.bam | $bcfTools view -bcvg - > ${fileName}.raw.tmp.bcf ; $bcfTools view ${fileName}.raw.tmp.bcf | $vcfScript varFilter -d 10 -D 2000 > ${fileName}.var.filt.vcf ; gzip -9 ${fileName}.var.filt.vcf"
             #echo $samCall > sam_prep.txt
             eval $samCall
#................................................................................................................................................................
             ##Run FreeC for copy variant analysis
         ##Set up the config file
             cp $configFile ./${fileName}_configFile.txt
             bamr=`readlink -f ${fileName}.sorted.bam`
             echo $bamr
             cmd=`sed -e s+length+$lenFILE+g -e s+CHROMOSOMES+$cPATH+g -e s+SAMTOOLS+$SAMTOOLS+g -e s+matefile+$bamr+g ./${fileName}_configFile.txt > ./${fileName}_config.txt`
             eval $cmd
             rm ./${fileName}_configFile.txt
             freeC="$FPATH -conf ./${fileName}_config.txt"
             #echo $freeC > free.txt
             eval $freeC
             echo "$fileName is done"
#................................................................................................................................................................
       ## Create a BED format file for the results
             BEDMAKE="${FsPATH}/scripts/freec2bed.pl -f ${fileName}.sorted.bam_ratio.txt > ${fileName}.bed"
             echo $BEDMAKE
             #eval $BEDMAKE
       ## Create a window of 1000 bp and a slider of 100bp in that window
             makeCOVER="$BEDPATH coverage -abam ${fileName}.sorted.bam -b ${fileName}.bed > ${fileName}_cov.bed"
             echo $makeCOVER
             #eval $makeCOVER
#................................................................................................................................................................
       ## Use delly to determine coverage
             makeCOVER="${dPATH}/cov -s 500 -o 250 -f ${fileName}.cov.gz -g $rPATH ./${fileName}.sorted.bam ; zcat ${fileName}.cov.gz > ./${fileName}.cov"
             echo $makeCOVER
             eval $makeCOVER
#................................................................................................................................................................
       ## Plot the graphs 
             pRun=`Rscript ~/coverage_grapher.R $aPATH ${fileName} $chrDes $startPos $stopPos`
               eval $pRun
               #mv ~/*pdf ~/Graphs
             rm ./${fileName}.sam
       ## Extract the CNVs from the CNV file of each sample to find duplicate in range between start and stop provided
             cnv="${fileName}.sorted.bam_CNVs"
             eRun="Rscript ~/cnvextract.r ${fileName} ${cnv} ${chrDes} ${startPos} ${stopPos}"
             eval $eRun
             
       ## Call estMOI             
             #estCall="$estPATH ${fileName}.sorted.bam ${fileName}.var.filt.vcf.gz $rPATH --out=${fileName}"
             #echo $estCall > estCaller.txt
             #eval $estCall
             cd ..
done