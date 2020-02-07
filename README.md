# KPC-ID-Resistance-Virulence

Neste tutorial faremos uma análise metagenômica para detecção de patógenos e identificação de genes bacterianos (KPC) de resistência e virulência em duas amostras de líquor de pacientes com neuroinfecções.

## Visão geral do pipeline

Nosso tutorial irá conter as seguintes etapas:

1. Configuração das ferramentas de bioinformática
2. Download das amostras
3. Controle de qualidade
4. Remoção de reads humanos
5. Classificação taxonômica, filtragem e visualização de abundância microbiana
6. Montagem de metagenoma
7. Identificação de genes de resistência
8. Identificação de genes de virulência

----------------------------------------------------------------

## Pré-requisitos

Para realizar esse tutorial, você precisa ter um computador com o [Ubuntu 18.04](https://ubuntu.com/) instalado. Caso não tenha, você pode seguir os passos abaixo para usar uma máquina virtual:

- Baixar e instalar a [VirtualBox](https://www.virtualbox.org/).
- Baixar e carregar a imagem (ISO) do [Ubuntu 18.04](http://releases.ubuntu.com/18.04/) no VirtualBox. Para isso, você pode seguir este [tutorial](https://www.youtube.com/watch?v=zsqJhle7CXE).

## Pipeline

Para que possamos ter acesso aos arquivos deste tutorial, iremos baixar os dados dele. Para isso, vamos clonar o repositório do Github para ter os mesmos arquivos na nossa máquina.

Primeiro instalamos o git:

```
sudo apt-get install git
```

Depois, clonamos o repositório:

```
git clone https://github.com/vinibfranc/KPC-ID-Resistance-Virulence
```

Pronto! Agora já podemos iniciar nosso tutorial!

### 1. Configuração das ferramentas de bioinformática

Primeiramente iremos configurar as ferramentas necessárias para a execução do pipeline para não nos preocuparmos com isso mais tarde. 

O primeiro passo é dar permissão para esse script ```configuracao.sh``` ser executado:

```
$ chmod +x configuracao.sh
```

Para rodar o script, basta criar uma variável com o caminho para baixar e instalar as ferramentas e em seguida rodar o script:

```
$ MEU_CAMINHO="/home/vinibfranc/TCC/AR_Bac/KPC-ID-Resistance-Virulence/ferramentas"
$ ./configuracao.sh $MEU_CAMINHO
```

OBS.: O caminho irá variar de acordo com o computador, por isso você precisa escolher um caminho válido no seu computador.

Depois, precisamos salvar as alterações no arquivo ```~/.bashrc```, fazendo:

```
$ source ~/.bashrc
```


### 2. Download das amostras

Inicialmente, iremos criar uma pasta para colocar nossas amostras:

```
$ mkdir -p amostras
$ cd amostras
```

As amostras no formato FASTQ podem ser baixadas com a ferramenta [sra-tools](https://github.com/ncbi/sra-tools) utilizando o seguinte comando:

```
$ fasterq-dump SRR8580960
$ fasterq-dump SRR8580963
```

Para visualizar alguns reads da amostra com código de acesso SRR8580963, podemos executar o comando: 

```
$ head SRR8580963.fastq
```

Para visualizar todos, rodamos:

```
$ cat SRR8580963.fastq
```

### 3. Controle de qualidade

Para visualizar como estão nossos reads após o sequenciamento utilizamos o [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/):

```
$ fastqc *.fastq -t 4
```

Para visualizar os relatórios de qualidade (QC) das amostras, podemos abrir e analisar os arquivos ```SRR8580960_fastqc.html``` e ```SRR8580963_fastqc.html``` gerados na pasta ```amostras/```.

Para remover reads de baixa qualidade e adaptadores das amostras, utilizamos as ferramentas [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) e [cutadapt](https://cutadapt.readthedocs.io/en/stable/):


```
$ trim_galore --quality 30 --phred33 --fastqc_args "-t 4" SRR8580960.fastq
$ trim_galore --quality 30 --phred33 --fastqc_args "-t 4" SRR8580963.fastq
```

Agora iremos visualizar os relatórios de qualidade (QC) das amostras após o controle de qualidade, acessíveis nos arquivos ```SRR8580960_trimmed_fastqc.html``` e ```SRR8580963_trimmed_fastqc.html``` gerados na pasta ```amostras/```.

Os arquivos com os reads que passaram no controle de qualidade e podem ser usados nas etapas posteriores são: ```SRR8580960_trimmed.fq``` e ```SRR8580963_trimmed.fq```.

### 4. Remoção de reads humanos

Essa etapa irá requerer o download do genoma humano (GRCh38), a construção de uma hash table para o alinhamento e o alinhamento propriamente dito utilizando [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Por ser bastante demorado, irei disponibilizar os arquivos resultantes no [link](https://mega.nz/#F!8QYnWC7D!ZX6EkNGuJ5wDN838oWC45w).

Baixar os arquivos ```SRR8580960_trimmed.sam``` e ```SRR8580963_trimmed.sam```, criar pasta em ```results/bowtie2/sam``` e copiar os arquivos para ela.

<b>-----> Você não precisa executar os comandos abaixo, pois os arquivos já foram baixados! <-----</b>

De qualquer forma, os passos para replicação dessa etapa são:

#### 4.1. Download do genoma humano já indexado

```
$ cd ..
$ mkdir -p "ref_dbs/human_db"
$ cd "ref_dbs/human_db"
$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
$ tar xvzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
```

#### 4.2. Alinhamento ao genoma humano

```
$ cd ../..
$ mkdir -p results/bowtie2/sam
$ bowtie2 --threads 4 -x ref_dbs/human_db/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -U amostras/SRR8580960_trimmed.fq -S results/bowtie2/sam/SRR8580960_trimmed.sam
$ bowtie2 --threads 4 -x ref_dbs/human_db/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -U amostras/SRR8580963_trimmed.fq -S results/bowtie2/sam/SRR8580963_trimmed.sam
```

Será gerado o arquivo SAM, o qual armazena as sequências alinhadas à sequência de referência, bem como suas coordenadas genômicas.

OBS.: É esperado que nada (ou praticamente nada) alinhe ao genoma humano, pois os reads humanos já haviam sido retirados antes da submissão ao SRA.

<b>-----> Agora você pode voltar a executar os comandos normalmente! <-----</b>

#### 4.3. Remoção de reads humanos

Inicialmente, convertemos o arquivo SAM para um arquivo binário (BAM):

```
$ mkdir -p results/bowtie2/bam
$ samtools view -bS results/bowtie2/sam/SRR8580960_trimmed.sam > results/bowtie2/bam/SRR8580960_trimmed.bam
$ samtools view -bS results/bowtie2/sam/SRR8580963_trimmed.sam > results/bowtie2/bam/SRR8580963_trimmed.bam
```

Depois, pegamos os reads não mapeados ao genoma humano: 

```
$ mkdir -p results/bowtie2/unmapped
$ samtools view -u -f 4 results/bowtie2/bam/SRR8580960_trimmed.bam > results/bowtie2/unmapped/SRR8580960_trimmed.bam
$ samtools view -u -f 4 results/bowtie2/bam/SRR8580963_trimmed.bam > results/bowtie2/unmapped/SRR8580963_trimmed.bam
```

Finalmente, transformamos o arquivo BAM para FASTQ para fazer a classificação taxonômica posteriormente: 

```
$ mkdir -p results/bowtie2/fastq
$ samtools fastq results/bowtie2/unmapped/SRR8580960_trimmed.bam > results/bowtie2/fastq/SRR8580960_trimmed.fastq
$ samtools fastq results/bowtie2/unmapped/SRR8580963_trimmed.bam > results/bowtie2/fastq/SRR8580963_trimmed.fastq
```

### 5. Classificação taxonômica, filtragem e visualização de abundância microbiana

Nessa etapa, os reads serão comparados contra um abrangente banco de dados utilizando o [Kraken2](https://ccb.jhu.edu/software/kraken2/), a fim de identificar os micro-organismos presentes neste metagenoma. 

Essa etapa irá requerer o download de genomas de vírus, bactérias, fungos e parasitas do NCBI, a construção de uma hash table para o alinhamento e o alinhamento propriamente dito. Por ser bastante demorado, irei disponibilizar os arquivos resultantes no [link](https://mega.nz/#F!8QYnWC7D!ZX6EkNGuJ5wDN838oWC45w).

Baixar os arquivos da pasta ```kraken2``` no link disponibilizado e inserir dentro da pasta ```results``` do seu computador.

<b>-----> Você não precisa executar os comandos abaixo, pois os arquivos já foram baixados! <-----</b>

De qualquer forma, os passos para replicação dessa etapa são:

#### 5.1. Download dos genomas de referência de micro-organismos

Inicialmente, fazemos o download dos genomas presentes no Kraken2:

```
$ kraken2-build --download-taxonomy --threads 4 --db ref_dbs/pathogens_db
$ kraken2-build --download-library viral --threads 4 --db ref_dbs/pathogens_db
$ kraken2-build --download-library archaea --threads 4 --db ref_dbs/pathogens_db
$ kraken2-build --download-library fungi --threads 4 --db ref_dbs/pathogens_db
$ kraken2-build --download-library protozoa --threads 4 --db ref_dbs/pathogens_db
$ kraken2-build --download-library bacteria --threads 4 --db ref_dbs/pathogens_db
```

Depois, fazemos o download do genoma de outros micro-organismos que já foram identificados como causadores de neuroinfecções, mas não estavam no banco de dados anterior. Para isso rodamos o script ```download_extra.sh```:

```
$ chmod +x download_extra.sh
$ ./download_extra.sh
```

#### 5.2. Construção do index

```
$ kraken2-build --build --threads 4 --max-db-size 13000000000 --db ref_dbs/pathogens_db
```

#### 5.3. Classificação taxonômica

```
$ DB_PATH="ref_dbs/pathogens_db"
$ KRAKEN2="results/kraken2"
$ mkdir -p $KRAKEN2
$ mkdir -p $KRAKEN2/classified
$ mkdir -p $KRAKEN2/unclassified
$ mkdir -p $KRAKEN2/tabular
$ mkdir -p $KRAKEN2/report
$ kraken2 --db $DB_PATH --threads 4 \
            --report "$KRAKEN2/report/SRR8580960_trimmed.kreport" \
            --classified-out "$KRAKEN2/classified/SRR8580960_trimmed.fastq" \
            --unclassified-out "$KRAKEN2/unclassified/SRR8580960_trimmed.fastq" \
            --output "$KRAKEN2/tabular/SRR8580960_trimmed.txt" results/bowtie2/fastq/SRR8580960_trimmed.fastq
$ kraken2 --db $DB_PATH --threads 4 \
            --report "$KRAKEN2/report/SRR8580963_trimmed.kreport" \
            --classified-out "$KRAKEN2/classified/SRR8580963_trimmed.fastq" \
            --unclassified-out "$KRAKEN2/unclassified/SRR8580963_trimmed.fastq" \
            --output "$KRAKEN2/tabular/SRR8580963_trimmed.txt" results/bowtie2/fastq/SRR8580963_trimmed.fastq
```

<b>-----> Agora você pode voltar a executar os comandos normalmente! <-----</b>

Analise os resultados gerados em ```results/kraken2/classified```, ```results/kraken2/unclassified```, ```results/kraken2/report``` e ```results/kraken2/tabular```.

#### 5.4. Filtragem de resultados para incluir somente patógenos de neuroinfecções

Para reduzir nosso escopo de análise, podemos utilizar um script em Python para filtrar os resultados para considerar somente patógenos de neuroinfecções:

```
$ python3 filtrar_patogenos.py 
```

Analise os resultados gerados em ```results/kraken_genus_species_strains``` e ```kraken_pathogens```.

#### 5.5. Estimativa de abundância (KRONA plot)

Para saber a abundância total de micro-organismos nas amostras, podemos utilizar a ferramenta [KRONA](https://github.com/marbl/Krona):

```
$ mkdir -p ferramentas/Krona-master/KronaTools/taxonomy
$ ferramentas/Krona-master/KronaTools/updateTaxonomy.sh
$ mkdir -p results/krona
$ ImportTaxonomy.pl -o results/krona/SRR8580960_trimmed_krona.html -t 3 -s 4 results/kraken2/tabular/SRR8580960_trimmed.txt
$ ImportTaxonomy.pl -o results/krona/SRR8580963_trimmed_krona.html -t 3 -s 4 results/kraken2/tabular/SRR8580963_trimmed.txt
```

Agora podemos visualizar o gráfico gerado, que mostra a composição microbiana de cada amostra, nos arquivos ```SRR8580960_trimmed_krona.html``` e ```SRR8580963_trimmed_krona.html``` que estão na pasta ```results/krona```.


### 6. Montagem de metagenoma

Para montar os reads em fragmentos maiores com o objetivo posterior de identificar genes de resistência e virulência, utilizamos a ferramenta [megahit](https://github.com/voutcn/megahit):

```
$ mkdir -p results/megahit
$ megahit -r results/bowtie2/fastq/SRR8580960_trimmed.fastq -o results/megahit/SRR8580960_trimmed.out
$ megahit -r results/bowtie2/fastq/SRR8580963_trimmed.fastq -o results/megahit/SRR8580963_trimmed.out
```

O resultado final das montagens pode ser encontrado nos arquivos ```final.contings.fa```, dentro da pasta ```results/megahit```.

### 7. Identificação de genes de resistência

Para identificar genes de resistência, iremos utilizar o banco de dados "The Comprehensive Antibiotic Resistance Database" ([CARD](https://card.mcmaster.ca/analyze)).

Inicialmente vamos acessar https://card.mcmaster.ca/download e procurar por ```Download CARD Data```. Em seguida, vamos fazer o download do banco de dados.

Em seguida vamos extrair o arquivo e mover para a pasta ```ref_dbs```.

Após ler o arquivo ```CARD-Download-README.txt```, podemos perceber que o arquivo ```nucleotide_fasta_protein_homolog_model.fasta``` possui as sequências de referência que desejamos. Então iremos construir uma banco de dados de referência para depois fazer um alinhamento utilizando o [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) em linha de comando.

Para construir o banco de dados, digite:

```
$ makeblastdb -in ref_dbs/card-data/nucleotide_fasta_protein_homolog_model.fasta -title Resistance_db -dbtype nucl -out ref_dbs/card-data/Resistance_db
```

Depois, iremos usar o BLASTN para alinhar nossas sequências montadas ao banco de dados CARD:

```
$ mkdir -p results/blastn_card_montados
$ blastn -query results/megahit/SRR8580960_trimmed.out/final.contigs.fa -db ref_dbs/card-data/Resistance_db -out results/blastn_card_montados/SRR8580960_trimmed.tab -evalue 1e-5 -outfmt 6 -max_target_seqs 1
$ blastn -query results/megahit/SRR8580963_trimmed.out/final.contigs.fa -db ref_dbs/card-data/Resistance_db -out results/blastn_card_montados/SRR8580963_trimmed.tab -evalue 1e-5 -outfmt 6 -max_target_seqs 1
```

Para analisar os dados apresentados em ```results/blastn_card_montados```, consideremos a seguinte tabela, que apresenta o significado das colunas em ordem:

Field | Description
| --- | --- |
| qseqid | query (e.g., gene) sequence id
| sseqid | subject (e.g., reference genome) sequence id
| pident | percentage of identical matches
| length | alignment length
| mismatch | number of mismatches
| gapopen | number of gap openings
| qstart | start of alignment in query
| qend | end of alignment in query
| sstart | start of alignment in subject
| send | end of alignment in subject
| evalue | expect value
| bitscore | bit score

Para analisarmos mais especificamente nossos resultados, podemos acessar a interface Web do CARD (https://card.mcmaster.ca/) e fazer as buscas. O que precisamos para isso é ir em nossos arquivos tabulados presentes em ```results/blastn_card_montados``` e procurar pelos acession numbers na segunda coluna do arquivo. Por exemplo, ```ARO:3004122```. Depois de encontrá-los, podemos ir na barra de busca no canto superior direito da página e colar esse código. Finalmente, vamos ter acesso a uma página como esta: https://card.mcmaster.ca/ontology/41247, com dados relevantes desta sequência associada a resistência que foi encontrada na nossa amostra. 

Abaixo, uma lista com os códigos encontrados nesta análise, para facilitar nossas buscas:

Para a amostra SRR8580960:

```
ARO:3004122 (Klebsiella)
ARO:3001328 (Escherichia)
ARO:3000216
ARO:3000518
ARO:3000793
ARO:3002985
ARO:3000491
ARO:3000216
ARO:3003056
ARO:3004580 (Klebsiella)
ARO:3003952
ARO:3003923
ARO:3000830
ARO:3004588 (Klebsiella)
ARO:3003209
ARO:3000676
ARO:3000796
ARO:3001084
ARO:3004122 (Klebsiella)
ARO:3000237
```

Para a amostra SRR8580963:

```
ARO:3002683
ARO:3002683
ARO:3000167
```

### 8. Identificação de genes de virulência

Para identificar genes de virulência, iremos utilizar o banco de dados "Virulence Factors of Bacterial Pathogens" ([VFDB](http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi)).

Inicialmente vamos acessar http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi e procurar por ```Download``` no canto inferior esquerdo. Em seguida, vamos fazer o download do banco de dados clicando em ```Download sequences of full dataset```.

Em seguida vamos extrair o arquivo, criar uma pasta em ```ref_dbs/vfdb``` e mover o arquivo para a pasta ```ref_dbs/vfdb```.

Para construir o banco de dados, digite:

```
$ makeblastdb -in ref_dbs/vfdb/VFDB_setB_nt.fas -title Virulence_db -dbtype nucl -out ref_dbs/vfdb/Virulence_db
```

Depois, iremos usar o BLASTN para alinhar nossas sequências montadas ao banco de dados VFDB:

```
$ mkdir -p results/blastn_vfdb_montados
$ blastn -query results/megahit/SRR8580960_trimmed.out/final.contigs.fa -db ref_dbs/vfdb/Virulence_db -out results/blastn_vfdb_montados/SRR8580960_trimmed.tab -evalue 1e-5 -outfmt 6 -max_target_seqs 1
$ blastn -query results/megahit/SRR8580963_trimmed.out/final.contigs.fa -db ref_dbs/vfdb/Virulence_db -out results/blastn_vfdb_montados/SRR8580963_trimmed.tab -evalue 1e-5 -outfmt 6 -max_target_seqs 1
```

Para analisar os dados apresentados em ```results/blastn_vfdb_montados```, consideremos a mesma tabela apresentada logo acima.

Para analisarmos mais especificamente nossos resultados...