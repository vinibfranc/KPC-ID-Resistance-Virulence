#!/bin/bash
set -e
set -u
set -o pipefail

echo "Configurando dependências para rodar o pipeline..."

# Para baixar os arquivos FASTQ do Sequence Read Archive (SRA)
mkdir -p $1
cd $1
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6-1/sratoolkit.2.9.6-1-ubuntu64.tar.gz
tar xzvf sratoolkit.2.9.6-1-ubuntu64.tar.gz
cd sratoolkit.2.9.6-1-ubuntu64/
echo 'PATH=$PATH:'$(pwd)/bin/ >> ~/.bashrc
cd ..

# Instalação do FASTQC
sudo apt-get install -y fastqc

# Instalação do cutadapt
sudo apt-get -y install python3
sudo apt install -y python3-pip
pip3 install --user cutadapt

# Instalação do TrimGalore
sudo apt-get -y install curl
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.0.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
cd TrimGalore-0.6.0/
echo 'PATH=$PATH:'$(pwd) >> ~/.bashrc
cd ..
sudo apt-get install -y pigz

# Instalação do Bowtie2
sudo apt-get install -y unzip
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download
unzip download
cd bowtie2-2.3.3.1-linux-x86_64/
echo 'PATH=$PATH:'$(pwd) >> ~/.bashrc
cd ..

# Dependências do SAMTOOLS e BEDTOOLS
sudo apt-get -y install gcc make libncurses5-dev libncursesw5-dev libbz2-1.0 libz-dev libbz2-dev libbz2-ocaml libbz2-ocaml-dev liblzma-dev

# Instalação do SAMTOOLS
sudo apt-get -y install samtools

# Instalação do  BEDTOOLS
sudo apt-get -y install bedtools

# Instalação do megahit
wget https://github.com/voutcn/megahit/releases/download/v1.2.8/MEGAHIT-1.2.8-Linux-x86_64-static.tar.gz
tar zvxf MEGAHIT-1.2.8-Linux-x86_64-static.tar.gz
cd MEGAHIT-1.2.8-Linux-x86_64-static/bin/
echo 'PATH=$PATH:'$(pwd) >> ~/.bashrc
cd ../..

# Dependências do Kraken2
sudo apt-get -y install build-essential
sudo apt-get -y install ncbi-blast+

# Instalação do Kraken2
wget https://github.com/DerrickWood/kraken2/archive/master.zip
unzip master.zip
cd kraken2-master
mkdir bin
./install_kraken2.sh bin
echo 'PATH=$PATH:'$(pwd)/bin/ >> ~/.bashrc
cd ..

# Instalação do KRONA
wget https://github.com/marbl/Krona/archive/master.zip
unzip master.zip.2
cd Krona-master/
KronaTools/install.pl --prefix $1/Krona-master/KronaTools/
echo 'PATH=$PATH:'$(pwd)/KronaTools/scripts/ >> ~/.bashrc