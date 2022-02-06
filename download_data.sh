#!/bin/bash
### Extract the first column from the csv file
cut -f 1 -d ',' $1 > accession_code.txt
### Skip the first line of the column
tail -n +2 accession_code.txt > codes.txt
### download the data from the website
while IFS='' read -r line || [[ -n "$line" ]]; do
    ##get each line from the text file
    echo "Text read from file: $line"
    ## get the first six letters of the 
    echo "First six is: ${line:0:6}"
    preName=${line:0:6}
    url1="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$preName/$line/${line}_1.fastq.gz"
    url2="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$preName/$line/${line}_2.fastq.gz"
    cmd1=`wget $url1`
    cmd2=`wget $url2`
    
done < codes.txt
