#!/bin/bash
if [[ $CLEAN ]]; then
	rm ../output/ASVs_Pythium_taxonomy*.txt 
fi

assign() {
	##### Create blast database ######
	echo -e "\n[+] Create blastdb $1 $(date)"
	makeblastdb -in $1 -dbtype nucl -out ../database/blastdb/pythium_database

	######### Run localblast #########
	echo -e "\n[+] Run localblast"
	blastn -task megablast -query ../output/ASVs_Pythium.fa -db ../database/blastdb/pythium_database -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp" -num_threads 15 -out ../output/ASVs_Pythium_taxonomy_blast_raw.txt

	# filter blast output by highest bitscore and similarity >= 95 % otherwise -> "NA"
	echo -e "\n[+] Format blast output"
	python3 parse_blast.py ../output/ASVs_Pythium_taxonomy_blast_raw.txt ../output/ASVs_Pythium.fa

	# cleanup
	rm ../output/ASVs_Pythium_taxonomy_blast_raw.txt
}

assign ../database/Pythium_robideau_BOLD.fa

echo -e "assign_taxonomy.sh\nAssignment of taxonomy successfully: $(date)\n" >> ../logs/pipeline_status.txt
