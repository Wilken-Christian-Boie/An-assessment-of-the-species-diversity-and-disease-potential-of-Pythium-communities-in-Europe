import os
import sys
import requests
from datetime import datetime
from collections import defaultdict

tax = {"Pythium" : "Sar;Stramenopiles;Oomycota;Pythiales;Pythiaceae;",
       "Albugo" : "Sar;Stramenopiles;Oomycota;Albuginales;Albuginaceae;",
       "Achlya" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;",
       "Aphanomyces" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;",
       "Apodachlya" : "Sar;Stramenopiles;Oomycota;Leptomitales;Leptomitales_incertea_sedis;",
       "Basidiophora" : "Sar;Stramenopiles;Oomycota;Peronosporales;Peronosporaceae;",
       "Brevilegnia" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;",
       "Eurychasma" : "Sar;Stramenopiles;Oomycota;Myzocytiopsidales;Eurychasmataceae;",
       "Halophytophthora" : "Sar;Stramenopiles;Oomycota;Pythiales;Pythiaceae;",
       "Hyaloperonospora" : "Sar;Stramenopiles;Oomycota;Peronosporales;Peronosporaceae;",
       "Peronospora" : "Sar;Stramenopiles;Oomycota;Peronosporales;Peronosporaceae;",
       "Lagenidium" : "Sar;Stramenopiles;Oomycota;Lagenidiales;Lagenidiaceae;",
       "Leptolegnia" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;",
       "Phytophthora" : "Sar;Stramenopiles;Oomycota;Lagenidiales;Lagenidiaceae;",
       "Phytophthora" : "Sar;Stramenopiles;Oomycota;Peronosporales;Peronosporaceae;",
       "Phytopythium" : "Sar;Stramenopiles;Oomycota;Pythiales;Pythiaceae;",
       "Plasmopara" : "Sar;Stramenopiles;Oomycota;Peronosporales;Peronosporaceae;",
       "Plasmoverna" : "Sar;Stramenopiles;Oomycota;Peronosporales;Peronosporaceae;",
       "Plectospira" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;",
       "Protoachlya" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;",
       "Pseudoperonospora" : "Sar;Stramenopiles;Oomycota;Peronosporales;Peronosporaceae;",
       "Pythiogeton" : "Sar;Stramenopiles;Oomycota;Pythiales;Pythiogetonaceae;",
       "Pythiopsis" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;",
       "Saprolegnia" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;",
       "Thraustotheca" : "Sar;Stramenopiles;Oomycota;Saprolegniales;Saprolegniaceae;"}

class FILE:
	ids = "../database/robideau_bold_ids"
	pre = "../database/Pythium_robideau_BOLD_pre.fa"
	name = "../database/Pythium_robideau_BOLD.fa"
	sample_id = "./figures/Data/new_taxonomy.txt"

class BOLD:
	url = "http://www.boldsystems.org/index.php/Public_DownloadData"
	obj = {"searchfield": "",
		   "downloadtype": "fasta:nuc",
		   "id_type": "processid",
		   "reporttype": "",
		   "cartToken": "general_445C1F6D-B260-4DAB-A53A-78D27B0EA255"}

def get_sample_ids(filepath):
	samples = {}
	with open(filepath, "r") as sfile:
		for x,line in enumerate(sfile, start=1):
			if x == 1: continue
			elements = line.split("\t")
			samples[elements[3]] = elements[2].replace(" ","_")
	return(samples)

def loadIds(filepath):
	return list(open(filepath, "r").read().splitlines())

def download(query, url, url_data): 
	url_data["searchfield"] = query
	x = requests.post(url, data = url_data)
	with open(FILE.pre, "w") as bold: bold.write(x.text)

def check(filepath):
	if os.path.exists(filepath):
		print(f"Remove {filepath}")
		os.remove(filepath)

def filter_and_save(d):
	sample_id = get_sample_ids(FILE.sample_id)
	with open(FILE.name, "a") as ffile:
		for header,sequence in d.items():
			if "ITS" in header:
				elements = header.split("|")
				bold_id = elements[0]
				taxonomy = elements[1].split(" ")
				if len(taxonomy) == 2: genus, species = taxonomy
				else: genus, species = taxonomy[0], "_".join(taxonomy[1:])
				
				#if "sp._nov" in species: taxonomy_name = f"{tax[genus]}{genus};{species}_{sample_id[bold_id]}"
				#else: taxonomy_name = f"{tax[genus]}{genus};{species}"
				taxonomy_name = f"{tax[genus]}{genus};{species}"
 
				ffile.write(f">{taxonomy_name}\n{sequence}\n") 

def read_fasta(filepath):
	entry = defaultdict()
	header = ""
	with open(filepath, "r") as fasta:
		for line in fasta:
			line = line.strip()
			if line == "": continue
			if line.startswith(">"):
				header = line.rstrip().replace(">","")
				entry[header] = ""
			else: entry[header] += line.replace("-","")
	return entry

def main():
	check(FILE.name)
	ids = loadIds(FILE.ids)
	query_string = " ".join(ids)
	print(f"Loaded {len(ids)} IDs for download")
	download(query_string, BOLD.url, BOLD.obj)
	filter_and_save(read_fasta(FILE.pre))

if __name__ == "__main__":
	main()
	os.system(f"echo 'download_bold.py\nDownload BOLD database successfully: {datetime.now()}\n' >> ../logs/pipeline_status.txt")
	os.system(f"rm {FILE.pre}")
