#!/usr/bin/env python3

import os
import sys
import pickle
from collections import defaultdict

BLAST = sys.argv[1]
FASTA = sys.argv[2]

class OUTPUT:
	summary = "../output/ASVs_Pythium_taxonomy_decision.txt"
	taxonomy = "../output/ASVs_Pythium_taxonomy.txt"

class THRESHOLD:
	""" define blast filter criteria """
	show = 50
	similarity = 95

def highest_abundant_species(hits):
	""" select highest abundant species """
	if(len(hits) == 0): return None
	species = [x["tax"] for x in hits]
	occurrences = dict(sorted({x:species.count(x) for x in species}.items(), key = lambda x:x[1], reverse = True))
	max_occ = list(occurrences.values()).count(max(occurrences.values()))
	if max_occ > 1:
		return(list(sorted(dict(filter(lambda x:x[1] == max(occurrences.values()), occurrences.items())).keys())))
	return(list(occurrences.keys()))

def load_blast(filepath):
	""" read and parse blast output """
	blast_dict = defaultdict(list)
	with open(filepath, "r") as file:
		for line in file:
			elements = line.strip().split("\t")
			blast_dict[elements[0]].append({"tax": elements[1], "similarity": float(elements[2]), "eval": float(elements[10]), "bitscore": float(elements[11])}) 
	return(blast_dict)

def filter_blast(filepath):
	""" load blast output and select most occurrent hit based on highest bitscore, similarity and quantity if multiple hits """
	blast = load_blast(filepath)
	hits = {}		
	for k,v in blast.items():
		max_ = max([x["bitscore"] for x in v])
		best_hits = [x for x in v if x["bitscore"] == max_ and x["similarity"] >= THRESHOLD.similarity]
		if len(best_hits) != 0:
			best_hits_tax_list = [x["tax"] for x in best_hits]
			best_hits_quant = {x:best_hits_tax_list.count(x) for x in best_hits_tax_list}
			top_hit = highest_abundant_species(best_hits)[0].replace(";"," ")
			hits[k] = top_hit
			save_summary(k,top_hit,best_hits,best_hits_quant,v)
	return(hits)

def save_summary(asv,top,best,best_quant,all_):
	with open(OUTPUT.summary, "a") as summary:
		summary.write(f"\n{'#'*10} {asv} {'#'*10}\n")
		summary.write("\n--- TOP HIT --- \n")
		summary.write(f"{top}\n")
		summary.write(f"\n--- GROUP OF BEST HITS --- ({len(best)})\n")
		summary.write(f"Quantity:\n")
		for k,v in best_quant.items(): summary.write(f"{k}\t found {v} time(s)\n")
		summary.write(f"\nBlast hits:\n")
		for x in best: summary.write(f"{str(x)}\n")
		summary.write(f"\n--- ALL HITS --- ({len(all_)})\t first {THRESHOLD.show} shown\n")
		for x in all_[:THRESHOLD.show]: summary.write(f"{str(x)}\n")
	
def fasta_id(filepath):
	""" get all ASV ids of fasta file """
	ids = []
	with open(filepath, "r") as file:
		for line in file:
			if line.startswith(">"): ids.append(line.strip().replace(">",""))
	return(ids)

def join_info(fasta,blast):
	""" combine blast output based on ASV label and apply given filter criteria """
	for asv in fasta:
		hit = " ".join(["NA" for x in range(7)]) #NOTE: defaults to NA, gets overwritten if ASV is valid (BLAST)
		if asv in blast.keys():	hit = f"{blast[asv]}"
		with open(OUTPUT.taxonomy, "a") as tax:
			tax.write(f"{str(asv)}\t{str(hit)}\n")

def main():
	""" output formatted table for phyloseq input """
	join_info(fasta_id(FASTA),filter_blast(BLAST))

if __name__ == "__main__":
	main()
