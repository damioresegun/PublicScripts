#!/bin/bash

#for batch copying demultiplexed output from qcat from individual fastq outputs from guppy
guppyOut="$HOME/ReBaseCall/QcatDemultiplex/5xBarc/qcat_0/"
batchSave="$HOME/ReBaseCall/QcatDemultiplex/5xBarc/Consensus/"
cp ${guppyOut}/* $batchSave
gupsta=$(dirname $guppyOut)
cd $gupsta
pwd
for f in ${gupsta}/qcat_1/*
do
file=$(basename $f)
fname="${file%.*}"
echo $fname
cat $f >> ${batchSave}/${fname}.fastq
done
cd ${gupsta}
echo
echo "qcat1 done"
echo
for f in ${gupsta}/qcat_2/*
do
file=$(basename $f)
fname="${file%.*}"
echo $fname
cat $f >> ${batchSave}/${fname}.fastq
done
echo
echo "qcat2 done"
echo
for f in ${gupsta}/qcat_3/*
do
file=$(basename $f)
fname="${file%.*}"
echo $fname
cat $f >> ${batchSave}/${fname}.fastq
done
echo
echo "qcat3 done"
echo
for f in ${gupsta}/qcat_4/*
do
file=$(basename $f)
fname="${file%.*}"
echo $fname
cat $f >> ${batchSave}/${fname}.fastq
done
echo
echo "qcat4 done"