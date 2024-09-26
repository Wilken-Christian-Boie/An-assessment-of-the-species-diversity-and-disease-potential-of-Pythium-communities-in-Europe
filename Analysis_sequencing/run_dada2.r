# Pipeline for analyzing the geographic distribution of pythium sp. in europe
# based on several tutorials:
# https://benjjneb.github.io/dada2/tutorial.html#learn-the-error-rates
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

# load libraries
library(dada2)		# run main microbiome analysis
library(ggplot2)	# nice plots
library(cowplot)	# arrange multiple plots
library(crayon)		# allow colored console output
library(optparse)	# handle command line options

# get command line options
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Path of input folder"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Path to output folder"),
  make_option(c("--trimmf"), type = "numeric", default = 200, help = "Minimal length of forward read. Shorter reads will be discarded"),
  make_option(c("--trimmr"), type = "numeric", default = 200, help = "Minimal length of reverse read. Shorter reads will be discarded"),
  make_option(c("-p", "--project"), type = "character", default = "myproject" , help = "Uniqe name to label resulting files"),
  make_option(c("-d", "--database"), type = "character", default = NULL , help = "Path to database file (fasta) required as database")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied\n", call.=FALSE)
}

# Set path of sequencing reads
path <- opt$input
fnFs <- sort(list.files(path, pattern = "R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_filtered"), `[`, 1)

# Place filtered files in filtered/subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim sequences
cat(paste("\n[>]",green("Filter and trim sequences\n"), sep = " "))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(opt$trimmf,opt$trimmr), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
write.table(out, paste("../logs/filter_stats_dada2_",opt$project,".txt", sep = ""), quote = FALSE, col.names = TRUE, row.names = TRUE)
head(out)

# Learn error rates
cat(paste("\n[>]", green("Learn error rates\n"), sep = " "))
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Determine inference
cat(paste("\n[>]", green("Determine inference\n"), sep = " "))
## Empty files will crash the script, so we have to 
## remove empty samples of filtFs and filtRs
## https://github.com/benjjneb/dada2/issues/1279 
cat("Missing samples:\n")
cat(filtFs[!file.exists(filtFs)], sep = "\n")
cat(filtRs[!file.exists(filtRs)], sep = "\n")
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Merge pairs
cat(paste("\n[>]", green("Merge pairs\n"), sep = " "))
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
head(mergers[[1]])

# Construct count matrix
cat(paste("\n[>]", green("Construct sequence table\n"), sep = " "))
seqtab <- makeSequenceTable(mergers)
cat("[>] Samples:", dim(seqtab)[1], "\n", "ASVs:", dim(seqtab)[2], "\n")

# Remove chimera and display frequency
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
cat("[>] Samples:", dim(seqtab.nochim)[1], "\n", "ASVs:", dim(seqtab.nochim)[2], "\n")
cat(paste("[>] Ratio:", sum(seqtab.nochim)/sum(seqtab)))

# Summary of previous pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
tmp_track <- data.frame(track)
tmp_track$samples <- rownames(track)
write.table(tmp_track, paste(opt$output,"/ASVs_",opt$project,"_filter_statistics.txt",sep = ""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Assign taxonomy
cat(paste("\n[>]", green("Assign taxonomy\n"), sep = " "))
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa <- assignTaxonomy(seqtab.nochim, opt$database, multithread = 15, minBoot = 80, taxLevels = ranks)
taxa <- data.frame(taxa)
taxa$Taxonomy <- paste(taxa$Kingdom, taxa$Phylum, taxa$Class, taxa$Order, taxa$Family, taxa$Genus, taxa$Species)


# EXPORT OBJECTS #
cat(paste("\n[>]", green("Export results\n"), sep = " "))

# Edit ASV names
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# Save ASVs to fasta file
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(opt$output,"/ASVs_",opt$project,".fa",sep = ""))

# Create and save count matrix (ASV table)
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste(opt$output,"/ASVs_",opt$project,"_counts.txt",sep = ""), sep = "\t", quote = FALSE, col.names = NA)

# Save assigned taxonomy to file
#colnames(taxa) <- "Taxonomy"
rownames(taxa) <- gsub(pattern = ">", replacement = "", x = asv_headers)
write.table(taxa["Taxonomy"], paste(opt$output,"/ASVs_",opt$project,"_taxonomy.txt",sep = ""), sep = "\t", quote = FALSE, col.names = NA)
