#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -j y
#$ -N MaxGuppywtdbg2
#$ -pe multi 32
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -M dro@st-andrews.ac.uk

# Script to take data through parameter searching for wtdbg2. It will take in the error corrected reads from canu and then take the data through the different parameter
set -e
reads="$HOME/All_De_Novo/MaxGuppy/quick"
output="$HOME/All_Wtdbg2/MaxGuppy"
tool="$HOME/Tools/wtdbg2-master"
mkdir -p $output
THREADS="32"
#parameters to be used are: 
	## -p : Kmer psize, 0 <= p <= 25, [21]        k + p <= 25, seed is <k-mer>+<p-homopolymer-compressed>	
	## -t : number of threads
	## -L : drop reads shorter than this integer [5000 recommended for PacBio]
	## -k : kmer fsize, 0 <= k <= 25, [0]
	## -AS : Keep contained reads during alignment; Subsampling kmers, 1/(<-S>) kmers are indexed, [4.00]
            ## -S is very useful in saving memeory and speeding up
            ## please note that subsampling kmers will have less matched length
	## -K : Filter high frequency kmers, maybe repetitive, [1000.05]
            ## >= 1000 and indexing >= (1 - 0.05) * total_kmers_count
	## -s : in similarity, calculated by kmer matched length / aligned length, [0.05]
	## -fo : Force to overwrite output files; -o=Prefix of output files (REQUIRED), []
	## -e : Min read depth of a valid edge, [3]
	## -g : approximate genome size
	## -x : presets. for ont: -x ont
module load oraclejava/jdk1.8.0_74

for f in $reads/* 
do
	echo $f
	fname=$(basename "$f")
	echo $fname
	mkdir -p $output/$fname
	cd $output/$fname
	pwd
	# for each file, go through the following parameters
		## wtdbg2 preset for ont and corrected with l 1000
		$tool/wtdbg2 -t $THREADS -L 1000 -p 21 -k 0 -AS 4 -K 0.05 -s 0.5 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN1000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN1000.ctg.lay.gz -fo ${fname}_MIN1000.ctg.preset_L1000.lay.fa
	## wtdbg2 preset for ont and corrected with l 2000
		$tool/wtdbg2 -t $THREADS -L 2000 -p 21 -k 0 -AS 4 -K 0.05 -s 0.5 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.preset_L2000.lay.fa
	## changing the kmer size parameter
		#default kmer size
		$tool/wtdbg2 -t $THREADS -L 2000 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.pdefault_L2000.lay.fa
		## kmer size of 19 and min length 2000
		$tool/wtdbg2 -t $THREADS -L 2000 -p 19 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p19_L2000.lay.fa
		## kmer size of 20 and min length 2000
		$tool/wtdbg2 -t $THREADS -L 2000 -p 20 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p20_L2000.lay.fa
		## kmer size of 21 and min length 2000
		$tool/wtdbg2 -t $THREADS -L 2000 -p 21 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p21_L2000.lay.fa
		## kmer size of 22 and min length 2000
		$tool/wtdbg2 -t $THREADS -L 2000 -p 22 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p22_L2000.lay.fa
		## kmer size of 23 and min length 2000
		$tool/wtdbg2 -t $THREADS -L 2000 -p 23 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p23_L2000.lay.fa
		## kmer size of 24 and min length 2000
		#$tool/wtdbg2 -t $THREADS -L 2000 -p 24 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		#$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p24_L2000.lay.fa
		## kmer size of 25 and min length 2000
		#$tool/wtdbg2 -t $THREADS -L 2000 -p 25 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		#$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p25_L2000.lay.fa		
	## changing the read depth
		#default kmer size and changed read depth
		$tool/wtdbg2 -t $THREADS -L 2000 -e2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.pdefault_e2_L2000.lay.fa
		## kmer size of 19 and min length 2000 and changed read depth
		$tool/wtdbg2 -t $THREADS -L 2000 -p 19 -e2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p19_e2_L2000.lay.fa
		## kmer size of 20 and min length 2000 and changed read depth
		$tool/wtdbg2 -t $THREADS -L 2000 -p 20 -e2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p20_e2_L2000.lay.fa
		## kmer size of 21 and min length 2000 and changed read depth
		$tool/wtdbg2 -t $THREADS -L 2000 -p 21 -e2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p21_e2_L2000.lay.fa
		## kmer size of 22 and min length 2000 and changed read depth
		$tool/wtdbg2 -t $THREADS -L 2000 -p 22 -e2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p22_e2_L2000.lay.fa
		## kmer size of 23 and min length 2000 and changed read depth
		$tool/wtdbg2 -t $THREADS -L 2000 -p 23 -e2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p23_e2_L2000.lay.fa
		## kmer size of 24 and min length 2000 and changed read depth
		#$tool/wtdbg2 -t $THREADS -L 2000 -p 24 -e2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		#$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p24_e2_L2000.lay.fa
		## kmer size of 25 and min length 2000 and changed read depth
		#$tool/wtdbg2 -t $THREADS -L 2000 -p 25 -e2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		#$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p25_e2_L2000.lay.fa

    ## changing the subsampling
		#default kmer size and changed read depth and subsampling
		$tool/wtdbg2 -t $THREADS -L 2000 -e2 -AS2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.pdefault_e2_AS2_L2000.lay.fa
		## kmer size of 19 and min length 2000 and changed read depth and subsampling
		$tool/wtdbg2 -t $THREADS -L 2000 -p 19 -e2 -AS2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p19_e2_AS2_L2000.lay.fa
		## kmer size of 20 and min length 2000 and changed read depth and subsampling
		$tool/wtdbg2 -t $THREADS -L 2000 -p 20 -e2 -AS2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p20_e2_AS2_L2000.lay.fa
		## kmer size of 21 and min length 2000 and changed read depth and subsampling
		$tool/wtdbg2 -t $THREADS -L 2000 -p 21 -e2 -AS2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p21_e2_AS2_L2000.lay.fa
		## kmer size of 22 and min length 2000 and changed read depth and subsampling
		$tool/wtdbg2 -t $THREADS -L 2000 -p 22 -e2 -AS2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p22_e2_AS2_L2000.lay.fa
		## kmer size of 23 and min length 2000 and changed read depth and subsampling
		$tool/wtdbg2 -t $THREADS -L 2000 -p 23 -e2 -AS2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p23_e2_AS2_L2000.lay.fa
		## kmer size of 24 and min length 2000 and changed read depth and subsampling
		#$tool/wtdbg2 -t $THREADS -L 2000 -p 24 -e2 -AS2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		#$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p24_e2_AS2_L2000.lay.fa
		## kmer size of 25 and min length 2000 and changed read depth and subsampling
		#$tool/wtdbg2 -t $THREADS -L 2000 -p 25 -e2 -AS2 -i ${f}/${fname}.correctedReads.fasta.gz -fo ${fname}_MIN2000
		#$tool/wtpoa-cns -t $THREADS -i ${fname}_MIN2000.ctg.lay.gz -fo ${fname}_MIN2000.ctg.p25_e2_AS2_L2000.lay.fa
	cd $output
done



