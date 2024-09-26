#!/bin/bash

if [[ $CLEAN ]]; then
	rm ../3_cooccurrence/*.tsv
	rm ../3_cooccurrence/bootstrap_correlation/*.tsv
	rm ../3_cooccurrence/bootstrap_counts/*.tsv
fi

BOOTSTRAPS=10000
ITERATIONS=100

# create summarised ASV counts
# outputs: ASVs_Pythium_counts_summarized.txt
echo -e "[+] Create summarised ASV counts $(date)"
Rscript summarise_before_fastspar.r

echo -e "[+] Run fastspar $(date)"
echo "NOTE: Added 'OTU_ID' to header as requirement to fastspar"
fastspar --yes --iterations $ITERATIONS --threads 15 --otu_table ../3_cooccurrence/ASVs_Pythium_counts_summarized.tsv --correlation ../3_cooccurrence/correlation_fastspar.tsv --covariance ../3_cooccurrence/covariance_fastspar.tsv

echo -e "[+] Calculate exact pvalues $(date)"
echo -e "[+] First, perform bootstrapping $(date)"
fastspar_bootstrap --threads 15 --otu_table ../3_cooccurrence/ASVs_Pythium_counts_summarized.tsv --number $BOOTSTRAPS --prefix ../3_cooccurrence/bootstrap_counts/bootstraps

echo -e "[+] Second, calculate correlation values $(date)"
parallel fastspar --yes --otu_table {} --correlation ../3_cooccurrence/bootstrap_correlation/cor_{/} --covariance ../3_cooccurrence/bootstrap_correlation/cov_{/} --iterations $ITERATIONS ::: ../3_cooccurrence/bootstrap_counts/*

echo -e "[+] Finally, calculate pvalues $(date)"
fastspar_pvalues --threads 15 --otu_table ../3_cooccurrence/ASVs_Pythium_counts_summarized.tsv --correlation ../3_cooccurrence/correlation_fastspar.tsv --prefix ../3_cooccurrence/bootstrap_correlation/cor_ --permutations $BOOTSTRAPS --outfile ../3_cooccurrence/pvalues_fastspar.tsv
