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

## Pipeline

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

Depois, precisamos salvar as alterações no arquivo ```~/.bashrc```, fazendo:

```
source ~/.bashrc
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

Essa etapa irá requerer o download do genoma humano (GRCh38), a construção de uma hash table para o alinhamento e o alinhamento propriamente dito utilizando [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Por ser bastante demorado, irei disponibilizar os arquivos resultantes no [link](https://mega.nz/#F!8QYnWC7D).

Baixar os arquivos ```SRR8580960_trimmed.sam``` e ```SRR8580963_trimmed.sam```, criar pasta em ```results/bowtie2/sam``` e copiar os arquivos para ela.

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

Essa etapa irá requerer o download de genomas de vírus, bactérias, fungos e parasitas do NCBI, a construção de uma hash table para o alinhamento e o alinhamento propriamente dito. Por ser bastante demorado, irei disponibilizar os arquivos resultantes no [link](https://mega.nz/#F!8QYnWC7D).

Baixar os arquivos da pasta ```kraken2```.

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

### 6. Montagem de metagenoma

### 7. Identificação de genes de resistência

### 8. Identificação de genes de virulência

