#$ -cwd

cd /storage/home/users/dro/StudentWork/NewRun

#echo `pwd`
#exit

#/storage/home/users/dro/StudentWork/gth-1.7.1-Linux_x86_64-64bit/bin/gth -genomic /storage/home/users/dro/StudentWork/sks339.contigs.fasta.gz -cdna /storage/home/users/dro/StudentWork/Plasmodium_knowlesi_strain_h_gca_900004885.PKNA1-C.2.cdna.all.fa.gz -gff3out -skipalignmentout -o /storage/home/users/dro/StudentWork/Genethreading/sks339_myGenePredict.gff3.gz

/storage/home/users/dro/StudentWork/gth-1.7.1-Linux_x86_64-64bit/bin/gth -genomic /storage/home/users/dro/StudentWork/sks339.contigs.fasta -protein /storage/home/users/dro/StudentWork/Plasmodium_knowlesi_strain_h_gca_900004885.PKNA1-C.2.pep.all.fa -skipalignmentout -gff3out -o /storage/home/users/dro/StudentWork/NewRun/MyGenomeThreader_sks339.aln -force -paralogs -prseedlength 20 -prhdist 2 -gcmincoverage 80 -prminmatchlen 20


/storage/home/users/dro/StudentWork/gm_et_linux_64/gmes_petap/gmes_petap.pl --verbose --sequence=/storage/home/users/dro/StudentWork/sks339.contigs.fasta --EP=/storage/home/users/dro/StudentWork/NewRun/MyGenomeThreader_sks339.aln --cores=8

