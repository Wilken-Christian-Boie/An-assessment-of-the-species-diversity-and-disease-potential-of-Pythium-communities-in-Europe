#!/bin/bash

# Run the script: "./pipeline.sh | tee ../output/pipeline_log.txt"
###############################################################################
# Procedure of microbiome analysis (Pythium Survey Europe)
#[>] 1. Download databases (Pythium|Fungi|Fusarium|Rhizoctonia)
#[>] 2. Prepare NGS reads
#[>] 3. ASV analyis using 'dada2' R package
#[>] 4. Create multiple sequence alignment (MSA) using 'mafft' software
#[>] 5. Create phylogenetic tree using 'FastTree' software
#[>] 6. Specific microbiome analysis using 'phyloseq' R package
###############################################################################

######### Prerequisites ##########
# set color for output
CYAN='\033[0;36m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# create log file of pipeline
rm ../logs/pipeline_status.txt
touch ../logs/pipeline_status.txt

# set length to truncate forward and reverse reads
pythium_length_f=200
pythium_length_r=200

######## Start analysis ##########
echo -e "[${GREEN}+${NC}] Start analysis: $(date)"
echo -e "[${GREEN}+${NC}] Download database" 
#[>] 1. Download databases
python3 download_bold.py
echo

########## Prepare reads #########
echo -e "[${GREEN}+${NC}] Prepare reads (pick and cut sequences with cutapadt out of raw sequences)"
#[>] 2. Prepare NGS reads
CLEAN=1 ./prepare_reads.sh
echo

## Get read length distribution ##
echo -e "[${GREEN}+${NC}] Get read length distribution of clean reads"
echo -e "[>] ${CYAN}Pythium${NC}" && Rscript read_length_distribution.r -f ../logs/Pythium_R1_length_distribution.txt -r ../logs/Pythium_R2_length_distribution.txt --trimmf $pythium_length_f --trimmr $pythium_length_r -o ../output/Pythium_read_length.pdf

###### Start dada2 analysis ######
echo -e "[${GREEN}+${NC}] Run dada2 analysis (main analysis)"
#[>] 3. ASV analyis using 'dada2' R package
echo -e "[>] ${CYAN}Pythium${NC}" && Rscript run_dada2.r --input ../2_clean_data --output ../output --trimmf $pythium_length_f --trimmr $pythium_length_r --project "Pythium" --database ../database/Pythium_Globisporangium_Sparolegnia_BOLD.fa

#### Assign taxonomy via blast ###
echo -e "[>] ${CYAN}Assign taxonomy${NC}" && CLEAN=1 ./assign_taxonomy.sh
echo

## Multiple sequence alignment ###
#[>] 4. Create multiple sequence alignment (MSA) using 'mafft' software
echo -e "[>] ${CYAN}Pythium${NC}" && mafft --maxiterate 1000 --localpair --thread 15 ../output/ASVs_Pythium.fa > ../output/ASVs_Pythium.aln

#### Create phylogenetic tree ####
#[>] 5. Create phylogeneic tree using 'FastTree' software
echo -e "[>] ${CYAN}Pythium${NC}" && FastTree -nt ../output/ASVs_Pythium.aln > ../output/ASVs_Pythium.nwk

#### Create phyloseq object ####
#[>] 6. Create phyloseq object
echo -e "[>] ${CYAN}Pythium${NC}" && ./create_phyloseq_object.r
echo -e "[${GREEN}+${NC}] ${GREEN}Analysis finished:${NC} $(date)"

