library(dplyr)
library(tidyr)

tax_label <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")

# read data
dat <- read.table("../output/ASVs_Pythium_counts.txt", sep = "\t", head = TRUE)
colnames(dat)[1] <- "OTU_ID" # set first colname to "OTU_ID
new_taxonomy <- read.table("./figures/Data/new_taxonomy_and_clades.txt", sep = "\t", header = TRUE) %>% 
					dplyr::select(Genus_neu, Species)

tax <- read.table("../output/ASVs_Pythium_taxonomy.txt", sep = "\t", head = FALSE)
tax <- tax %>%
		tidyr::separate(V2, into = tax_label, sep = " ") %>%
		dplyr::left_join(new_taxonomy, by = "Species") %>%
		dplyr::select(-Genus) %>%
		dplyr::rename(Genus = Genus_neu)

# create final form of count matrix
final <- dat %>%
            dplyr::left_join(tax, by = c("OTU_ID" = "V1")) %>%
            dplyr::mutate(tax = paste(Genus, Species, sep = " ")) %>%
            dplyr::select(-c("OTU_ID", all_of(tax_label)))

# get summed absolute read counts per GenusxSpecies pair
final <- aggregate(.~tax, final, sum)
names(final)[names(final) == "tax"] <- "OTU_ID"
final <- final %>%
			#dplyr::filter(!grepl("NA_NA|Saprolegnia", final$OTU_ID)) %>%
			dplyr::filter(grepl("Globisporangium|Elongisporangium|Phytopythium|Pilasporangium|Pythium", final$OTU_ID)) %>%
			dplyr::select(matches("OTU_ID|soil", perl = TRUE, ignore.case = TRUE))
			
# output summarised count table ready for fastspar
write.table(final, file = "../3_cooccurrence/ASVs_Pythium_counts_summarized.tsv", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(final$OTU_ID, file = "../3_cooccurrence/ASVs_Pythium_counts_summarized_name.tsv", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

