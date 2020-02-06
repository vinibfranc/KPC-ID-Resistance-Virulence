#!/bin/bash
set -e
set -u
set -o pipefail


DB_PATH="ref_dbs/pathogens_db"
mkdir -p "$DB_PATH/extra_taxids/"
cd "$DB_PATH/extra_taxids/"

mkdir -p parasites
cd parasites

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/870/725/GCA_001870725.1_MEX_genome_complete.1-6-13/GCA_001870725.1_MEX_genome_complete.1-6-13_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/524/195/GCF_000524195.1_ASM52419v1/GCF_000524195.1_ASM52419v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/469/725/GCA_000469725.3_EMULTI002/GCA_000469725.3_EMULTI002_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/368/765/GCA_006368765.1_ASM636876v1/GCA_006368765.1_ASM636876v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/237/925/GCF_000237925.1_ASM23792v2/GCF_000237925.1_ASM23792v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/699/445/GCF_000699445.1_SchHae_1.0/GCF_000699445.1_SchHae_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/803/305/GCA_000803305.1_Toxocara_canis_adult_r1.0/GCA_000803305.1_Toxocara_canis_adult_r1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/884/285/GCA_001884285.1_ASM188428v1/GCA_001884285.1_ASM188428v1_genomic.fna.gz
gunzip *.gz
find . -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db "../../../../$DB_PATH"

cd ..
mkdir -p virus
cd virus

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/861/205/GCF_000861205.1_ViralProj15297/GCF_000861205.1_ViralProj15297_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14518/GCF_000845085.1_ViralProj14518_genomic.fna.gz
gunzip *.gz
find . -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db "../../../../$DB_PATH"

cd ..
mkdir -p bacteria
cd bacteria

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/308/635/GCF_000308635.1_ASM30863v1/GCF_000308635.1_ASM30863v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/160/075/GCF_000160075.2_ASM16007v2/GCF_000160075.2_ASM16007v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/160/675/GCF_000160675.1_ASM16067v1/GCF_000160675.1_ASM16067v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/296/505/GCF_000296505.1_Acti_turi_ACS-279-V-COL4_V1/GCF_000296505.1_Acti_turi_ACS-279-V-COL4_V1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/307/545/GCF_001307545.1_ASM130754v1/GCF_001307545.1_ASM130754v1_genomic.fna.gz
gunzip *.gz
find . -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db "../../../../$DB_PATH"

cd ..
mkdir -p protozoa
cd protozoa

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/185/145/GCA_001185145.1_ASM118514v1/GCA_001185145.1_ASM118514v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/499/105/GCA_000499105.1_Naegleria_fowleri_1.0/GCA_000499105.1_Naegleria_fowleri_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/208/925/GCF_000208925.1_JCVI_ESG2_1.0/GCF_000208925.1_JCVI_ESG2_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/065/GCF_000209065.1_ASM20906v1/GCF_000209065.1_ASM20906v1_genomic.fna.gz
gunzip *.gz
find . -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db "../../../../$DB_PATH"

cd ..
mkdir -p fungi
cd fungi

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/585/GCF_000149585.1_ASM14958v1/GCF_000149585.1_ASM14958v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/835/475/GCF_000835475.1_Clad_bant_CBS_173_52_V1/GCF_000835475.1_Clad_bant_CBS_173_52_V1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/835/455/GCF_000835455.1_Fons_pedr_CBS_271_37_V1/GCF_000835455.1_Fons_pedr_CBS_271_37_V1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/525/GCA_000003525.2_BD_ER3_V1/GCA_000003525.2_BD_ER3_V1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/335/GCF_000149335.2_ASM14933v2/GCF_000149335.2_ASM14933v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/335/GCF_000151335.2_JCVI-cpa1-1.0/GCF_000151335.2_JCVI-cpa1-1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/477/535/GCF_001477535.1_Pneu_jiro_RU7_V2/GCF_001477535.1_Pneu_jiro_RU7_V2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/985/GCF_000001985.1_JCVI-PMFA1-2.0/GCF_000001985.1_JCVI-PMFA1-2.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/735/GCF_000150735.1_Paracocci_br_Pb18_V2/GCF_000150735.1_Paracocci_br_Pb18_V2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/275/GCF_000006275.2_JCVI-afl1-v2.0/GCF_000006275.2_JCVI-afl1-v2.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/276/285/GCA_002276285.1_Lprolificans_pilon/GCA_002276285.1_Lprolificans_pilon_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/230/625/GCF_000230625.1_Exop_derm_V1/GCF_000230625.1_Exop_derm_V1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/835/555/GCF_000835555.1_Rhin_mack_CBS_650_93_V1/GCF_000835555.1_Rhin_mack_CBS_650_93_V1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/710/275/GCA_000710275.1_ASM71027v1/GCA_000710275.1_ASM71027v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/901/145/GCA_002901145.1_ASM290114v1/GCA_002901145.1_ASM290114v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/615/GCF_000149615.1_ASM14961v1/GCF_000149615.1_ASM14961v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/642/055/GCF_001642055.1_Altal1/GCF_001642055.1_Altal1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/836/295/GCF_000836295.1_O_gall_CBS43764/GCF_000836295.1_O_gall_CBS43764_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/022/145/GCF_004022145.1_Paevar1/GCF_004022145.1_Paevar1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/221/725/GCA_002221725.1_ScBoyd1.0/GCA_002221725.1_ScBoyd1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/696/995/GCA_000696995.1_ApoEleB7760-1.0/GCA_000696995.1_ApoEleB7760-1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/738/825/GCA_000738825.1_ASM73882v1/GCA_000738825.1_ASM73882v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/599/735/GCA_001599735.1_JCM_2334_assembly_v001/GCA_001599735.1_JCM_2334_assembly_v001_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/335/GCF_000006335.3_ASM633v3/GCF_000006335.3_ASM633v3_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/961/545/GCF_000961545.1_S_schenckii_v1/GCF_000961545.1_S_schenckii_v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/855/GCF_000002855.3_ASM285v2/GCF_000002855.3_ASM285v2_genomic.fna.gz
gunzip *.gz
find . -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db "../../../../$DB_PATH"