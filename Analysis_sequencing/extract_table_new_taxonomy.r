#!/usr/bin/env Rscript

library(tabulapdf)

col_names <- c("Genus", "Isolate", "Isolate_ID", "BOLD_ID", "Genbank_acc")

# read and properly format table of Nguyen et al. 2022 (Supplementary information)
taxonomy_table <- tabulapdf::extract_tables("./figures/Data/supportinginfo.pdf", pages = c(1:19))
taxonomy_table <- do.call(rbind, lapply(taxonomy_table, function(x) {setNames(rbind(names(x), x), col_names)}))
cat("Found ", nrow(taxonomy_table), "entries in table\n", sep = " ")

# export table and BOLD_ID
write.table(taxonomy_table, file = "./figures/Data/new_taxonomy.txt", col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(na.omit(taxonomy_table$BOLD_ID), file = "../database/robideau_bold_ids", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
