#!/bin/bash
startPos=10000
stopPos=2500000
lenFILE="/home/dami/Pf.fasta.len"
chr="Pf3D7_12_v3"
echo $startPos

if (($startPos-5000 <= 0))
then
startPos="0"
echo $startPos
else
echo "Starting position remains $startPos"
fi
grep $chr $lenFILE > chrStuff.txt
chr12="chrStuff.txt"
fat=`cut -f 2 $chr12`
echo $fat
rm $chr12
if (($stopPos+10000 >= $fat))
then
echo "Stop position exceeds length of chromosome. Stop position set to $fa6t"
stopPos=$fat
echo $stopPos
else 
echo "Starting position remains $stopPos" 
fi
#length=`cut -f 2 $chr12`
#echo $length