#!/usr/bin/env R

library(optparse)
library(tidyverse)

option_list <- list(make_option(c("-n", "--data_name"), type = "character", default = NULL, help = "Name of dataset"), 
                    make_option(c("-p", "--pseudobulk_psi"),
  type = "character", default = NULL, help = "Combined pseudobulk psi matrix as input, tab-deliminated"), 
                    make_option(c("-o", "--output"), type = "character",
  default = NULL, help = "Directory where multiple output RDS file will be written"), 
                    make_option(c("-r", "--num_reads_min"), type = "numeric", default = 10,
  help = "Minimum number of total reads covering node, which will be included in the output. This is for ensure meaningful psi quantification, default = 10"),
                    make_option(c("-d", "--psi_diff_cutoff"), type = "numeric", default = 0.2, help = "Minimum PSI difference from dataset average which will lead to the node being kept in the index, default = 0.2"))

# parse input
opt <- parse_args(OptionParser(option_list = option_list))

NAME <- opt$data_name
INPUT <- opt$pseudobulk_psi
OUTPUT <- opt$output
NUM_READS_MIN <- opt$num_reads_min
PSI_DIFF_CUTOFF <- opt$psi_diff_cutoff

######################## Functions build various inputs for scfindME


# temporary, hard-coded nodes info file will be included in final package

ni <- readRDS("/nfs/production/irene/ma/ysong/DATA/SCFIND/VASA-seq/data/nodes_info.rds")

# build original matrix function
buildMatrix.original <- function(file, num_reads_min) {

  data <- readr::read_tsv(file, col_names = TRUE, progress = show_progress())

  message("all values collected, generating matrix...")

  matrix.original <- data %>%
    tidyr::unite("Gene_node", Gene, Node, sep = "_") %>%
    dplyr::group_by(Sample) %>%
    dplyr::distinct(Gene_node, .keep_all = TRUE) %>%
    dplyr::select(Sample, Gene_node, Total_Reads, Psi) %>%
    dplyr::mutate(Filter = case_when(Total_Reads < num_reads_min ~ "DROP", Total_Reads >= num_reads_min ~ "KEEP", TRUE ~ NA_character_)) %>%
    dplyr::filter(Filter == "KEEP") %>%
    dplyr::select(Sample, Gene_node, Psi) %>%
    dplyr::ungroup(Sample) %>%
    dplyr::group_by(Gene_node) %>%
    tidyr::pivot_wider(names_from = Sample, values_from = Psi)

  message("original matrix constructed, ready to be scaled")

  return(matrix.original)
}

## save matrix.original

matrix_original <- buildMatrix.original(INPUT, NUM_READS_MIN)
head(matrix_original)

saveRDS(matrix_original, paste(OUTPUT, "/", NAME, "_matrix_original.rds", sep = ""))

print("Finish original matrix building")


# scale matrix to get contrast
scaleMatrix.diff <- function(matrix.original, psi.diff.cutoff) {
  df <- data.frame(matrix.original, row.names = matrix.original$Gene_node)
  dm <- as.matrix(df[, -1])
  head(dm)

  mean <- rowMeans(dm, na.rm = TRUE)

  matrix.scaled_diff <- dm - mean

  matrix.scaled_diff <- matrix.scaled_diff * 100

  # drop all-na rows
  matrix.scaled_diff <- matrix.scaled_diff[which(rowSums(is.na(matrix.scaled_diff)) < ncol(matrix.scaled_diff)), ]

  matrix.scaled_diff_selected <- data.frame(row.names = rownames(matrix.scaled_diff))

  for (cell in seq(1, ncol(matrix.scaled_diff))) {
    tv <- matrix.scaled_diff[, cell]
    temp <- which(abs(tv) < psi.diff.cutoff * 100)
    tv[temp] <- NA
    matrix.scaled_diff_selected[[colnames(matrix.scaled_diff)[cell]]] <- tv
  }

  # drop again all-NA rows
  matrix.scaled_diff_selected <- matrix.scaled_diff_selected[which(rowSums(is.na(matrix.scaled_diff_selected)) < ncol(matrix.scaled_diff_selected)), ]

  return(matrix.scaled_diff_selected)
}


# Select the desired type of splicing events nodes from the scaled matrix
buildMatrix.above <- function(matrix.scaled) {

  # above matrix for above index
  matrix.above <- data.frame(row.names = rownames(matrix.scaled))
  # set na and below ones to zero
  for (cell in seq(1, ncol(matrix.scaled))) {
    tv <- matrix.scaled[, cell]
    temp <- which(is.na(tv) | tv < 0)
    tv[temp] <- 0
    matrix.above[[colnames(matrix.scaled)[cell]]] <- tv
  }
  return(matrix.above)
}


buildMatrix.below <- function(matrix.scaled) {

  # below matrix for below index
  matrix.below <- data.frame(row.names = rownames(matrix.scaled))
  # set na and above ones to zero
  for (cell in seq(1, ncol(matrix.scaled))) {
    tv <- matrix.scaled[, cell]
    temp <- which(is.na(tv) | tv > 0)
    tv[temp] <- 0
    matrix.below[[colnames(matrix.scaled)[cell]]] <- tv
  }
  matrix.below <- matrix.below * (-1)

  return(matrix.below)
}

matrix_scaled <- scaleMatrix.diff(matrix_original, PSI_DIFF_CUTOFF)
matrix_above <- buildMatrix.above(matrix_scaled)
matrix_below <- buildMatrix.below(matrix_scaled)

saveRDS(matrix_scaled, paste(OUTPUT, "/", NAME, "_matrix_scaled.rds", sep = ""))
saveRDS(matrix_above, paste(OUTPUT, "/", NAME, "_matrix_above.rds", sep = ""))
saveRDS(matrix_below, paste(OUTPUT, "/", NAME, "_matrix_below.rds", sep = ""))

print("Finish scaled matrix building")

buildMatrix.stats <- function(matrix.original, matrix.scaled) {
  df <- data.frame(matrix.original, row.names = rownames(matrix.original))
  dm <- as.matrix(df)
  # create @metadata$stats calculate node means and SD value and store in index
  message("calculating mean PSI across dataset...")

  mean <- transform(dm, mean = apply(dm, 1, mean, na.rm = TRUE))

  message("calculating SD...")
  sd <- transform(dm, SD = apply(dm, 1, sd, na.rm = TRUE))

  mean$SD <- sd$SD
  mean <- mean[order(mean$SD), ]
  stats <- mean[, c("mean", "SD")]
  stats <- stats[which(rownames(stats) %in% rownames(matrix.scaled)), ]
  return(stats)
}


buildMatrix.node_list <- function(matrix.scaled, nodes_details_data) {
  node_list <- rownames(matrix.scaled)
  node_list_all <- nodes_details_data[which(nodes_details_data$Gene_node %in% node_list), ]

  # install.packages('XML', repos = 'http://www.omegahat.net/R') BiocManager::install('biomaRt')
  library("biomaRt")
  # listMarts()
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
  node_list_all$Gene_num <- gsub("\\.\\d+$", "", node_list_all$Gene)

  # takes a few minuites to match gene to name
  gene_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = node_list_all$Gene_num, mart = ensembl)

  # de-duplicate and match gene name using gene id
  gene_name <- gene_name[!duplicated(gene_name$ensembl_gene_id), ]

  gene_node_all <- merge(node_list_all, gene_name, by.x = "Gene_num", by.y = "ensembl_gene_id", all.x = TRUE)

  return(gene_node_all)

}


stats <- buildMatrix.stats(matrix_original, matrix_scaled)
nodes_detail <- buildMatrix.node_list(matrix_scaled, ni)

saveRDS(stats, paste(OUTPUT, "/", NAME, "_metadata_stats.rds", sep = ""))
saveRDS(nodes_detail, paste(OUTPUT, "/", NAME, "_metadata_nodes_detail.rds", sep = ""))

print("Finish stats and node list matrix building")
print("Succesfully completed")
