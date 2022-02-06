 #!/bin/bash

#Method for checking if the software is installed on the machine (0 = installed, 1 = non-installed)
test_install()
{
	SOFTWARE_NAME="$1"

	STATUS=$(dpkg-query -W -f='${Status}' "$SOFTWARE_NAME")

	if [ "$STATUS" = "install ok installed" ]; then
		echo 0 #True
	else
		echo 1 #False
	fi
}



#Method is asking user if he wants to install a software (name as parameter) and give back the user's response (0 = yes, 1 = no)
yes_no()
{
	SOFTWARE_NAME="$1"
	
	continue="true"

	while [ $continue = "true" ];
	do
	    echo "Do you wish to install "$SOFTWARE_NAME" ? [y/n]"
	    read yn

	    if [ $yn = "y" ]; then
		return 0 #yes
		continue="false"
	    elif [ $yn = "n" ]; then
		return 1 #no
		continue="false"
	    else
		echo "Please answer with y or n"
	    fi

	done
}

#Method for creating directory based on the user inpit and unzip fastleopard inside it
FastLeopard_dir_chooser()
{
	continue="true"

	while [ $continue = "true" ];
	do
		echo -n "Please enter the folder in which FastLeopard will be installed > "
		read FastLeopard_PATH
		#Try to cd to the folder
		mkdir -p "$FastLeopard_PATH"
		#If the command returns 0 -> success
		ret=$?
		#If the command returns 0 -> success, we stop the while
		if [ $ret = 0 ]; then
			continue="false"
		else
			echo "Please enter an existing directory"
		fi
	done

#tar xvfz "FastLeopard.tar.gz" --directory "$FastLeopard_PATH";
}



# ----------------------------------------
# 	Dependencies for FastLeopard
# ----------------------------------------

#Ask for folder
FastLeopard_dir_chooser

# Set directory
FULL_PATH=$PWD/$FastLeopard_PATH

#Perl installation
if [ $(test_install "perl") = 0 ]; then
	echo "Perl already installed"
else
	echo "\n############################"
	echo "Perl installation"
	echo "############################\n"
	sudo apt-get install perl
	echo "\n############################"
	echo "Perl installation finished"
	echo "############################\n"
fi

sudo apt-get install build-essential
sudo apt-get install liblist-moreutils-perl
sudo apt-get install zlib1g-dev

#R installation
if [ $(test_install "r-base") = 0 ]; then
	echo "R already installed"
else
	
	echo "\n####################################################"
	echo "R installation and required libraries and packages"
	echo "#####################################################\n"
	sudo apt-get install r-base
	
	sudo su - -c "R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('Biobase')\""
	sudo su - -c "R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('edgeR')\""
	sudo su - -c "R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('DESeq2')\""
	sudo su - -c "R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('limma')\""
	sudo su - -c "R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('qvalue')\""
	sudo su - -c "R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('cluster')\""
	sudo su - -c "R -e \"source('http://bioconductor.org/biocLite.R'); biocLite('EBSeq')\""
        
	echo "\n############################"
	echo "R installation finished"
	echo "############################\n"

fi

#Blast installation
if [ $(test_install "ncbi-blast+") = 0 ]; then
	echo "Blast+ already installed"
else
	echo "Installation of blast+"
	sudo apt-get install ncbi-blast+
fi


#Blast databases
	DB_PATH="$FastLeopard_PATH/blast/db"
	cd "$FastLeopard_PATH/blast/"
	mkdir "./db/"

		#Swissprot
		yes_no "swissprot"
		ret=$?
		if [ $ret = 0 ]; then
			echo "Swissprot DB Downloading... (in $DB_PATH)"
			wget -N -P "$DB_PATH" "ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz"
		        tar xvfz "$DB_PATH/swissprot.tar.gz" --directory "$DB_PATH";
		fi

		#nt
		yes_no "nt"
		ret=$?
		if [ $ret = 0 ]; then
			echo "nt Downloading... (in $DB_PATH)"
			cd $DB_PATH
			NT=which update_blastdb --passive nt
			for f in nt.*
			do
				tar zxvf "$f" "$DB_PATH"
			done
		fi

		#nr
		yes_no "nr"
		ret=$?
		if [ $ret = 0 ]; then
			echo "nr Downloading... (in $DB_PATH)"
			cd $DB_PATH
			NR=which update_blastdb --passive nr
		fi

	#Update all the paths created
	echo "updating paths..."
	sudo updatedb


#Bowite installation
if [ $(test_install "bowtie") = 0 ]; then
	echo "Bowtie already installed"
else
	echo "\n############################"
	echo "Bowtie installation"
	echo "############################\n"

	sudo apt-get install bowtie

	echo "\n############################"
	echo "R installation finished"
	echo "############################\n"
fi

#Bowtie2 installation
if [ $(test_install "bowtie2") = 0 ]; then
	echo "Bowtie2 already installed"
else
	echo "\n############################"
	echo "Bowtie2 installation"
	echo "############################\n"

	sudo apt-get install bowtie2

	echo "\n############################"
	echo "Bowtie2 installation finished"
	echo "############################\n"
fi

#Samtools
if [ $(test_install "samtools") = 0 ]; then
	echo "Samtools already installed"
else
	echo "\n############################"
	echo "Samtools installation"
	echo "############################\n"

	sudo apt-get install samtools

	echo "\n############################"
	echo "Samtools installation finished"
	echo "############################\n"
fi

#Trinity
	cd $FULL_PATH
	
	Trinity_PATH="trinity"	
	if [ ! -d $Trinity_PATH ]; then
		mkdir -p "trinity"		
	fi	

	cd $Trinity_PATH

	yes_no "trinity"
	ret=$?
	if [ $ret = 0 ]; then
		url="https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.3.2.tar.gz"
		base=$(basename ${url})
		filename="trinity.tar.gz"				
			echo "\n############################"
			echo "Trinity is downloading... (in $Trinity_PATH)"
			echo "############################\n"
			wget -P "$Trinity_PATH" -O "$filename" "$url" 
									
			tar xvfz "trinity.tar.gz";
			cd "./trinityrnaseq-Trinity-v2.3.2/"
			
			sudo make
			sudo make plugins
				
	fi

	
## RSEM

	cd $FULL_PATH
	
	RSEM_PATH="rsem"
		
	if [ ! -d $RSEMm_PATH ]; then
		mkdir -p "rsem"
	fi
	
	cd $RSEM_PATH

		yes_no "RSEM"
		ret=$?
		if [ $ret = 0 ]; then
			url="https://github.com/deweylab/RSEM/archive/master.zip"
			base=$(basename ${url})
			filename="rsem.zip"
			
				echo "\n############################"
				echo "RSEM is downloading... (in $RSEM_PATH)"
				echo "############################\n"

				wget -P "$RSEM_PATH" -O "$filename" "$url" 
						
		        unzip $filename
			
			cd "./RSEM-master/"
			
			sudo make
			
		fi

## QUAST

	cd $FULL_PATH
	
	QUAST_PATH="quast"
		
	if [ ! -d $QUAST_PATH ]; then
		mkdir -p "quast"
	fi
	
	cd $QUAST_PATH

		yes_no "quast"
		ret=$?
		if [ $ret = 0 ]; then
			url="https://downloads.sourceforge.net/project/quast/quast-4.5.tar.gz"
			base=$(basename ${url})
			filename="quast.tar.gz"
				echo "\n############################"
				echo "QUAST is downloading... (in $QUAST_PATH)"
				echo "############################\n"
				sudo apt-get install -y pkg-config libfreetype6-dev libpng-dev python-matplotlib
				wget -P "$QUAST_PATH" -O "$filename" "$url"
									
		        	
				tar xvfz $filename;
		fi


## TRIMMOMATIC

cd $FULL_PATH
	
	TRIM_PATH="trimmomatic"
		
	if [ ! -d $TRIM_PATH ]; then		
		mkdir -p "trimmomatic"
	fi
	
	cd $TRIM_PATH

		yes_no "trimmomatic"
		ret=$?
		if [ $ret = 0 ]; then			
			url="https://github.com/timflutre/trimmomatic/archive/master.zip"
			base=$(basename ${url})
			filename="trimmomatic.zip"
			
				echo "\n############################"
				echo "Trimmomatic is downloading... (in $TRIM_PATH)"
				echo "############################\n"

				wget -P "$TRIM_PATH" -O "$filename" "$url" 
			
						
		        unzip $filename
			
			cd "./trimmomatic-master/"
			
			sudo make
		fi


## Printing all of the installed tools
echo "\n\n\n\n\n"
echo $FULL_PATH/$Trinity_PATH"/trinityrnaseq-Trinity-v2.3.2/Trinity"
echo $FULL_PATH/$RSEM_PATH
echo $FULL_PATH/$TRIM_PATH"/classes/trimmomatic.jar"
echo $FULL_PATH/$QUAST_PATH"/quast-master/quast.py"

