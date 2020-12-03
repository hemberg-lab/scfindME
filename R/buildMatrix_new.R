#' Builds an original \code{matrix} from a given Whippet Collapsed Quant .psi.gz files containing directory and a metadata 
#'
#'
#' @param file the merged Whippet .psi.tsv file, can be an object or a url
#' @name buildMatrix.original
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv show_progress
#' @importFrom dplyr select mutate distinct filter group_by ungroup case_when
#' @importFrom tidyr unite pivot_wider
#' @return a tibble object for further scailing
#'
buildMatrix.original <- function(file){
  
  data <- readr::read_tsv(file, col_names=TRUE, progress=show_progress())
  
  message("all values collected, generating matrix...")
  
  matrix.original <- data %>% 
    tidyr::unite("Gene_node", Gene, Node, sep="_") %>%
    dplyr::group_by(Sample) %>%
    dplyr::distinct(Gene_node, .keep_all = TRUE) %>%
    dplyr::select(Sample, Gene_node, Total_Reads, Psi) %>%
    dplyr::mutate(Filter = case_when(Total_Reads < 10 ~ "DROP",
                                     Total_Reads >= 10 ~ "KEEP",
                                     TRUE ~ NA_character_
    )) %>%
    dplyr::filter(Filter == "KEEP") %>%
    dplyr::select(Sample, Gene_node, Psi) %>%
    dplyr::ungroup(Sample) %>%
    dplyr::group_by(Gene_node) %>%
    tidyr::pivot_wider(names_from = Sample, values_from= Psi)
  
  message("original matrix constructed, to be scaled")
  
  return(matrix.original)
}

#' Scale the original \code{matrix} with Z-score normalization
#'
#'
#' @param matrix.original the original matrix from buildMatrix.original
#' @name scaleMatrix.z
#' @return a scaled matrix object
#'

scaleMatrix.z <- function(matrix.original){
  df <- data.frame(matrix.original, row.names = rownames(matrix.original))
  dm <- as.matrix(df)
  
  matrix.scaled_z <- t(scale(t(dm)))
  
  matrix.scaled_z <- matrix.scaled_z*100
  
  matrix.scaled_z <- matrix.scaled_z[which(rowSums(is.na(matrix.scaled_z)) < ncol(matrix.scaled_z)), ]
  return(matrix.scaled_z)
}


#' Scale the original \code{matrix} with difference from dataset mean
#' This version is stringent that it requires the difference of PSI to mean is higher than 0.2
#' @param matrix.original the original matrix from buildMatrix.original
#' @name scaleMatrix.diff
#' @return a scaled matrix object
#'
scaleMatrix.diff <- function(matrix.original){
  df <- data.frame(matrix.original, row.names = rownames(matrix.original))
  dm <- as.matrix(df)
  
  mean <- rowMeans(dm, na.rm = TRUE)
  
  matrix.scaled_diff <- dm - mean
  
  matrix.scaled_diff <- matrix.scaled_diff*100
  
  matrix.scaled_diff <- matrix.scaled_diff[which(rowSums(is.na(matrix.scaled_diff)) < ncol(matrix.scaled_diff)), ]
  
  matrix.scaled_diff_selected  <- data.frame(row.names = rownames(matrix.scaled_diff))
  
  for (cell in seq(1, ncol(matrix.scaled_diff))){
    tv <- matrix.scaled_diff[, cell]
    temp <- which(abs(tv) < 20)
    tv[temp] = 0
    matrix.scaled_diff_selected[[colnames(matrix.scaled_diff)[cell]]] <- tv
  }
  
  return(matrix.scaled_diff_selected)
}


#' Select the desired type of splicing events nodes from the scaled matrix
#' Above dataset mean inclusion events
#' @param matrix.scaled the output of either scaleMatrix.diff or scaleMatrix.z
#' @name buildMatrix.above
#' @return a matrix object
#' 
buildMatrix.above <- function(matrix.scaled){
  
  # above matrix for above index
  matrix.above <- data.frame(row.names = rownames(matrix.scaled))
  # set na and below ones to zero
  for (cell in seq(1, ncol(matrix.scaled))){
    tv <- matrix.scaled[, cell]
    temp <- which(is.na(tv) | tv < 0)
    tv[temp] = 0
    matrix.above[[colnames(matrix.scaled)[cell]]] <- tv
  }
  return(matrix.above)
}


#' Select the desired type of splicing events nodes from the scaled matrix
#' Below dataset mean inclusion events
#' @param matrix.scaled the output of either scaleMatrix.diff or scaleMatrix.z
#' @name buildMatrix.below
#' @return a matrix object
#' 
buildMatrix.below <- function(matrix.scaled){
  
  # below matrix for below index
  matrix.below <- data.frame(row.names = rownames(matrix.scaled))
  
  # set na and above ones to zero
  for (cell in seq(1, ncol(matrix.scaled))){
    tv <- matrix.scaled[, cell]
    temp <- which(is.na(tv) | tv > 0)
    tv[temp] = 0
    matrix.below[[colnames(matrix.scaled)[cell]]] <- tv
  }
  matrix.below <- matrix.below*(-1)
  
  return(matrix.below)
}


#' Calculate dataset-wide mean and SD of psi values of all nodes
#' The stats matrix would be stored in the alternative splicing index metadata slot
#' @param matrix.original the output original matrix from the function "buildMatrix.original"
#' @param matrix.scaled the output scaled matrix from either the function "buildMatrix.scaled.above" or "buildMatrix.scaled.below"
#' @name buildMatrix.stats
#'
#' @return a matrix object(stats)
#' 
buildMatrix.stats <- function(matrix.original, matrix.scaled){
  df <- data.frame(matrix.original, row.names = rownames(matrix.original))
  matrix.original <- as.matrix(df)
  # create @metadata$stats
  # calculate node means and SD value and store in index
  message("calculating mean PSI across dataset...")
  
  mean <- transform(matrix.original, mean=apply(matrix.original, 1, mean, na.rm = TRUE))
  
  message("calculating SD...")
  sd <- transform(matrix.original, SD=apply(matrix.original, 1, sd, na.rm = TRUE))
  
  mean$SD <- sd$SD
  mean <- mean[order(mean$SD), ]
  stats <- mean[, c("mean", "SD")]
  stats <- stats[which(rownames(stats)%in%rownames(matrix.scaled)), ]
  return(stats)
}

#' Collect details of all nodes, including gene id, gene name, node type, strand, coordinates etc.
#' The node_list matrix would be stored in the alternative splicing index metadata slot
#' @param dir a directory containing .psi.gz files from Whippet quantification
#' @param matrix.scaled the output scaled matrix from either the function "buildMatrix.scaled.above" or "buildMatrix.scaled.below"
#' @name buildMatrix.node_list
#'
#' @return a matrix object(node_list)
#' 
#' @importFrom biomaRt useMart useDataset getBM
#' @importFrom utils read.table
#' 
buildMatrix.node_list <- function(matrix.scaled, nodes_details_data){
  node_list <- rownames(matrix.scaled)
  node_list_all <- nodes_details_data[which(nodes_details_data$Gene_node %in% node_list), ]
  
  # install.packages("XML", repos = "http://www.omegahat.net/R")
  # BiocManager::install("biomaRt")
  library("biomaRt")
  # listMarts()
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
  node_list_all$Gene_num <- gsub("\\.\\d+$", "", node_list_all$Gene)
  
  # takes a few minuites to match gene to name
  gene_name <- getBM(attributes= c("ensembl_gene_id","external_gene_name"),
                     filters= "ensembl_gene_id",
                     values = node_list_all$Gene_num,
                     mart= ensembl)
  
  # de-duplicate and match gene name using gene id
  gene_name <- gene_name[!duplicated(gene_name$ensembl_gene_id), ]
  
  gene_node_all <- merge(node_list_all, gene_name, by.x = "Gene_num", by.y = "ensembl_gene_id", all.x=TRUE)
  
  return(gene_node_all)
  
}

