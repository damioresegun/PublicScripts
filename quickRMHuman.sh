#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -r y
#$ -j y
#$ -N RepeatMaskHumOct19
#$ -q all.q
#$ -l hostname="node1"
#$ -M dro@st-andrews.ac.uk

set -e

# for i in $HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/ReWrittenAssemblies/*;
# do foi=$(basename "$i");
# foiu=${foi%.*};
# echo $foiu; 
# yup="/storage/home/users/dro/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/RepeatModeling/${foiu}";
# mkdir -p $yup;
# cd $yup;
# $HOME/PipelineScripts/RM/run_pipeline_Human.py $foiu $i;
# cd;
# done

# for i in $HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/ReWrittenAssemblies/*; 
# do foi=$(basename "$i"); 
# foiu=${foi%.*}; 
# echo $foiu; 
# mkdir -p $HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/RepeatMasker/${foiu}; 
# $HOME/PipelineScripts/RM/run_pipeline_Human.py $foiu $i $HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/CDHIT/Pk_RepeatLib.fa $HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/RepeatMasker/${foiu}; 
# done

for i in $HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/ReWrittenAssemblies/*; 
do foit=$(basename "$i");
foitr=${foit%.*};
echo $foitr;
outr="/storage/home/users/dro/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/GFFs";
mkdir -p ${outr}/${foitr};
cd ${outr}/${foitr};
~/PipelineScripts/RM/run_pipeline_Human.py $foitr $i \
$HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/RepeatModeling/CombinedIsolates_consensi.fa.censor \
$HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/RepeatMasker/${foitr} \
$HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/LTRHarvest/${foitr}/${foitr}_ltr.out \
$HOME/PipelineAnalyses/HumanAnalysis/Testing/RepeatMasking/TransposonPSI/${foitr}/${foitr}.fasta.TPSI.allHits.chains.bestPerLocus
cd;
done
