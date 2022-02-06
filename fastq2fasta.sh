#!/bin/bash
path="LaptopData_preBarcode/RawData/20171016_1550_SKS201a/fastq/pass/*"
for f in $path; do
echo $f
drec=$(dirname "$f")
echo $drec
fullfilename=$(basename -- "$f")
echo $fullfilename
filename="${fullfilename%.*}"
echo $filename
sed -n '1~4s/^@/>/p; 2~4p' $f > "${drec}/${filename}.fasta"
done 
