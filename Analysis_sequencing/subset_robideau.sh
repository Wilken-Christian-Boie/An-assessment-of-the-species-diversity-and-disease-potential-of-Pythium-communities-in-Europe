#!/bin/bash

echo "[+] Get BOLD IDs of Robideau et al. $(date)"
cat ../database/database_robideau.txt | cut -f6 | grep "\S" > ../database/robideau_bold_ids

echo "[+] Get Genbank IDs of Robideau et al. $(date)"
cat ../database/database_robideau.txt | cut -f7,8,9 | sed s/\,\\t/\\n/g > ../database/robideau_genbank_ids
