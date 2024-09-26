#!/bin/bash

# set color for output
CYAN='\033[0;36m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

if [[ $CLEAN ]]; then
	rm ../2_clean_data/*.fastq
	rm ../logs/*summary.txt
	rm ../logs/*.html
	rm ../logs/*.zip
	rm -rf ../logs/multiqc*
	rm ../logs/*length_distribution.txt
fi

echo "[+] Start preparing raw reads: $(date)"
# store primer sequence for cutadapt
pythium_f="TGCGGAAGGATCATTACCACAC" # Pythium (taxid:4797) -> 304 bp fragment
pythium_r="GCGTTCAAAATTTCGATGACTC"

# filter reads function 
filter_reads(){
  COUNT=1
  FASTA_COUNT=$(ls ../1_data/*.fastq.gz | wc -l)
  for prefix in $(ls $1 | sed -r 's/_R[12][.]fastq.gz//' | uniq); do
    echo -e "[${GREEN}+${NC}] ${GREEN}Sample: $5 ${NC} ($COUNT/$(expr $FASTA_COUNT / 2)) $(basename $prefix)"
    cutadapt --cores 15 -g $3 -G $4 -o $2$(basename $prefix)_adaptertrimmed_$5_R1.fastq -p $2$(basename $prefix)_adaptertrimmed_$5_R2.fastq "$prefix"_R1.fastq.gz "$prefix"_R2.fastq.gz --discard-untrimmed > ../logs/$(basename $prefix)_adaptertrimmed_$5_summary.txt
    (( COUNT++ ))
    echo
  done
}

# filter all reads
echo "[+] Filter all reads based on primer sequence"
echo "[+] Info: Untrimmed reads were discarded"
filter_reads "../1_data/*.fastq.gz" "../2_clean_data/" $pythium_f $pythium_r Pythium

echo -e "[${GREEN}+${NC}] Run FastQC"
fastqc ../2_clean_data/*.fastq -q -t 15 -o ../logs/

echo -e "[${GREEN}+${NC}] Get read length distribution of all filtered files"
read_length(){
  awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}}' $1 > $2
}
read_length "../2_clean_data/*R1*.fastq" "../logs/Pythium_R1_length_distribution.txt"
read_length "../2_clean_data/*R2*.fastq" "../logs/Pythium_R2_length_distribution.txt"

echo -e "[${GREEN}+${NC}] Generate multiqc report: $(date)"
multiqc  ../logs/ --outdir ../logs/

echo
echo "[+] Preparing reads done: $(date)"
echo -e "prepare_reads.sh \n Preparing reads done successfully: $(date)" >> ../logs/pipeline_status.txt
