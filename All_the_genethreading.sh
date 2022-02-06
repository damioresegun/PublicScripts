#$ -cwd

Assem="$HOME/All_De_Novo/newCanu"
cd $Assem
pred="/storage/home/users/dro/All_GenePredict/newCanu"
#cd $pred
#echo `pwd`
#exit

#/storage/home/users/dro/StudentWork/gth-1.7.1-Linux_x86_64-64bit/bin/gth -genomic /storage/home/users/dro/StudentWork/sks339.contigs.fasta.gz -cdna /storage/home/users/dro/StudentWork/Plasmodium_knowlesi_strain_h_gca_900004885.PKNA1-C.2.cdna.all.fa.gz -gff3out -skipalignmentout -o /storage/home/users/dro/StudentWork/Genethreading/sks339_myGenePredict.gff3.gz

#"gth -genomic /storage/home/users/dro/StudentWork/sks339.contigs.fasta -protein /storage/home/users/dro/StudentWork/Plasmodium_knowlesi_strain_h_gca_900004885.PKNA1-C.2.pep.all.fa -skipalignmentout -gff3out -o /storage/home/users/dro/StudentWork/NewRun/MyGenomeThreader_sks339.aln -force -paralogs -prseedlength 20 -prhdist 2 -gcmincoverage 80 -prminmatchlen 20"


#echo "/storage/home/users/dro/StudentWork/gm_et_linux_64/gmes_petap/gmes_petap.pl --verbose --sequence=/storage/home/users/dro/StudentWork/sks339.contigs.fasta --EP=/storage/home/users/dro/StudentWork/NewRun/MyGenomeThreader_sks339.aln --cores=8"

for folder in ${Assem}/*
do
    echo "The folder you are working on is ${folder}"
    fname=$(basename "$folder")    
    echo $fname
    #cd $folder
    mkdir -p ${pred}/${fname}
    pwd
    cd ${Assem}/${fname}
    pwd
    echo "Running quast with the Lapp Reference..."
    
/shelf/apps/dro/conda/envs/quast/lib/python3.6/site-packages/quast-5.0.2-py3.6.egg-info/scripts/quast.py -t 14 ${Assem}/${fname}/${fname}.contigs.fasta -o /storage/home/users/dro/Assembly_QC/newCanuRun/${fname}/quast/ -r /storage/home/users/dro/Index/PkLappRefSeq.fna -g /storage/home/users/dro/Index/PkLappRefSeq.gff -e -f
   echo 
       echo "Running genethreading with the Lapp Reference..."
       echo 
gth -genomic ${Assem}/${fname}/${fname}.contigs.fasta -protein $HOME/GenethreadingIndex/PkLappRefSeq_AA.faa -skipalignmentout -gff3out -o ${pred}/${fname}/${fname}.aln -force -paralogs -prseedlength 20 -prhdist 2 -gcmincoverage 80 -prminmatchlen 20
    
    echo "Genethreading done..."

 #   echo "/shelf/apps/dro/conda/envs/quast/lib/python3.6/site-packages/quast-5.0.2-py3.6.egg-info/scripts/quast.py -t 12 /storage/home/users/dro/All_De_Novo/newAssem/sks339/sks339.contigs.fasta -o /storage/home/users/dro/Assembly_QC/newAssemblyRun/${fname}/quast/VsLapp/ -r /storage/home/users/dro/Index/PkLappRefSeq.fna -g /storage/home/users/dro/Index/PkErnestRefSeq.gff3 -e -f"
    cd $Assem
done