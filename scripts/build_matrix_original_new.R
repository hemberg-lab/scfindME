#!/usr/bin/env Rscript

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

# build original matrix function


data <- readr::read_tsv(INPUT, col_names = TRUE, progress = show_progress())
  message("all values collected, generating matrix...")

  matrix.original <- data %>%
    tidyr::unite("Gene_node", Gene, Node, sep = "_") %>%
    dplyr::group_by(Sample) %>%
    dplyr::distinct(Gene_node, .keep_all = TRUE) %>%
    dplyr::select(Sample, Gene_node, Total_Reads, Psi) %>%
    dplyr::mutate(Filter = case_when(Total_Reads < NUM_READS_MIN ~ "DROP", Total_Reads >= NUM_READS_MIN ~ "KEEP", TRUE ~ NA_character_)) %>%
    dplyr::filter(Filter == "KEEP") %>%
    dplyr::select(Sample, Gene_node, Psi) %>%
    dplyr::ungroup(Sample) %>%
    dplyr::group_by(Gene_node) %>%
    tidyr::pivot_wider(names_from = Sample, values_from = Psi)

  message("original matrix constructed, ready to be scaled")
  df <- data.frame(matrix.original, row.names = matrix.original$Gene_node)
  dm <- as.matrix(df[, -1])
  mean <- rowMeans(dm, na.rm = TRUE)
  matrix.scaled_diff <- dm - mean
  matrix.scaled_diff <- matrix.scaled_diff * 100
  # drop all-na rows
  matrix.scaled_diff <- matrix.scaled_diff[which(rowSums(is.na(matrix.scaled_diff)) < ncol(matrix.scaled_diff)), ]
  matrix.scaled_diff_selected <- data.frame(row.names = rownames(matrix.scaled_diff))

diff_cut <- matrix(0, nrow = nrow(matrix.scaled_diff), ncol = ncol(matrix.scaled_diff))

  for (cell in seq(1, ncol(matrix.scaled_diff))) {
    tv <- matrix.scaled_diff[, cell]
    tvd <- diff_cut[, cell]
    temp <- which(abs(tv) < 0.2 * 100)
    tv[temp] <- NA
    tvd[temp] <- 1
    matrix.scaled_diff_selected[[colnames(matrix.scaled_diff)[cell]]] <- tv
    diff_cut[, cell] <- tvd
  }
  
rownames(diff_cut) <- rownames(matrix.scaled_diff_selected)
colnames(diff_cut) <- colnames(matrix.scaled_diff_selected)

matrix.scaled_diff_selected <- matrix.scaled_diff_selected %>% filter(if_any(everything(), ~ !is.na(.)))

diff_cut <- diff_cut[rownames(matrix.scaled_diff_selected), ]

matrix.above <- data.frame(row.names = rownames(matrix.scaled_diff_selected))
  # set na and below ones to zero
  for (cell in seq(1, ncol(matrix.scaled_diff_selected))) {
    tv <- matrix.scaled_diff_selected[, cell]
    temp <- which(is.na(tv) | tv < 0)
    tv[temp] <- 0
    matrix.above[[colnames(matrix.scaled_diff_selected)[cell]]] <- tv
  }

matrix.below <- data.frame(row.names = rownames(matrix.scaled_diff_selected))
  # set na and above ones to zero
  for (cell in seq(1, ncol(matrix.scaled_diff_selected))) {
    tv <- matrix.scaled_diff_selected[, cell]
    temp <- which(is.na(tv) | tv > 0)
    tv[temp] <- 0
    matrix.below[[colnames(matrix.scaled_diff_selected)[cell]]] <- tv
  }
matrix.below <- matrix.below * (-1)



  df <- data.frame(matrix.original, row.names = matrix.original$Gene_node)
  dm <- as.matrix(df[, -1])


  mean <- transform(dm, mean = apply(dm, 1, mean, na.rm = TRUE))
  sd <- transform(dm, SD = apply(dm, 1, sd, na.rm = TRUE))
  mean$SD <- sd$SD
  mean <- mean[order(mean$SD), ]
  stats <- mean[, c("mean", "SD")]

stats <- stats[which(rownames(stats) %in% rownames(matrix.scaled_diff_selected)), ]


ni <- readRDS("/nfs/production/irene/ma/ysong/DATA/SCFIND/VASA-seq/data/node_info_vasa-seq.rds")

node_list <- rownames(matrix.scaled_diff_selected)

ni$Gene_node <- paste(ni$Gene, ni$Node, sep = "_")

node_list_all <- ni[which(ni$Gene_node %in% node_list), ]
node_list_all$Gene_num <- gsub("\\..*$", "", node_list_all$Gene)

  # install.packages('XML', repos = 'http://www.omegahat.net/R') BiocManager::install('biomaRt')
  library("biomaRt")
  # listMarts()
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

  # takes a few minuites to match gene to name
  gene_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = node_list_all$Gene_num, mart = ensembl)

  # de-duplicate and match gene name using gene id
  gene_name <- gene_name[!duplicated(gene_name$ensembl_gene_id), ]

  gene_node_all <- merge(node_list_all, gene_name, by.x = "Gene_num", by.y = "ensembl_gene_id", all.x = TRUE)


saveRDS(matrix.scaled_diff_selected, "matrix_scaled_diff_selected.rds")
saveRDS(matrix.above, "matrix_above.rds")
saveRDS(matrix.below, "matrix_below.rds")
saveRDS(diff_cut, "diff_cut.rds")
saveRDS(stats, "stats.rds")
saveRDS(gene_node_all, "gene_node_all.rds")


