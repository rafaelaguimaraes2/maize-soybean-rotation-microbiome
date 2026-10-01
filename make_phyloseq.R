# Load libraries:
library(phyloseq)
library(Biostrings)
library(RDPutils)
library(ape)

# Read in the ASV file:
otu <- read.table("otu_table.tsv", header = TRUE, row.names = 1, sep = "\t")


# Get the representative sequences:
rep.seqs <- Biostrings::readDNAStringSet("dna-sequences.fasta", format = "fasta")

# Import classification table
my.taxa <- make_tax_table("classification_table.tsv", confidence = 0.5)

# Import rooted tree file
my.tree <- ape::read.tree("rooted_tree.nwk")

# Read in sample data
sam <- read.table("sample_data.txt", header = TRUE, row.names = 1, sep = "\t")
sam <- sample_data(sam)

# Import into phyloseq:
otu.table <- otu_table(otu, taxa_are_rows = TRUE)
expt <- phyloseq(otu.table, rep.seqs, my.taxa, my.tree, sam)
expt
save(expt, file = "expt.rda")


