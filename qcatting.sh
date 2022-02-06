#$ -cwd

echo "Beginning master..."

cat ~/All_Basecalled_Data/MarchGuppy/pass/* | qcat -b ~/Demultiplex/qcat/ --detect-middle -t 16 --trim -k RBK004 --dual

echo "We are done out here"

