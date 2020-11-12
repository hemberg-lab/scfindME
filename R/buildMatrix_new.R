#' Builds an original \code{matrix} from a given Whippet Collapsed Quant .psi.gz files containing directory and a metadata 
#'
#' @param dir a directory containing .psi.gz files from Whippet quantification, looks like "/lustre/scratch117/cellgen/team218/ys6/TM_FACS/MicroExonator/Whippet/Quant/Single_Cell/Unpooled/"
#' @param metadata metadata of the alternative splicing dataset, must include a column named "SRA" and has SRA ids for all Whippet quantified samples
#' @param column_label the column name which indicates cell id and .pai.gz file name in metadata, usually use "SRA" and use SRA ids
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


scaleMatrix.z <- function(matrix.original){
  df <- data.frame(matrix.original, row.names = "Gene_node")
  dm <- as.matrix(df)
  
  matrix.scaled_z <- t(scale(t(dm)))
  
  matrix.scaled_z <- matrix.scaled_z*100
  
  matrix.scaled_z <- matrix.scaled_z[which(rowSums(is.na(matrix.scaled_z)) < ncol(matrix.scaled_z)), ]
  return(matrix.scaled_z)
}

scaleMatrix.diff <- function(matrix.original){
  df <- data.frame(matrix.original, row.names = "Gene_node")
  dm <- as.matrix(df)
  
  mean <- rowMeans(dm, na.rm = TRUE)
  
  matrix.scaled_diff <- dm - mean
  
  matrix.scaled_diff <- matrix.scaled_diff*100
  
  matrix.scaled_diff <- matrix.scaled_diff[which(rowSums(is.na(matrix.scaled_diff)) < ncol(matrix.scaled_diff)), ]
  return(matrix.scaled_diff)
}


#' Scale the original matrix and build the input matrix of an above index
#' The output of this function correspond to the input param "psival" in the buildAltSpliceIndex method
#' @param matrix.original the output original matrix from the function "buildMatrix.original"
#'
#' @name buildMatrix.above
#'
#' @return a matrix object(scaled above)
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

#' Scale the original matrix and build the input matrix of a below index
#' The output of this function correspond to the input param "psival" in the buildAltSpliceIndex method
#' @param matrix.original the output original matrix from the function "buildMatrix.original"
#'
#' @name buildMatrix.below
#'
#' @return a matrix object(scaled below)
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
  df <- data.frame(matrix.original, row.names = "Gene_node")
  matrix.original <- as.matrix(df)
  # create @metadata$stats
  # calculate node means and SD value and store in index
  
  sd <- transform(matrix.original, SD=apply(matrix.original, 1, sd, na.rm = TRUE))
  mean <- transform(matrix.original, mean=apply(matrix.original, 1, mean, na.rm = TRUE))
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
buildMatrix.node_list <- function(dir, matrix.scaled){
  url <- paste(as.character(dir), as.character(colnames(matrix.scaled)[1]), ".psi.gz", sep = "")
  t1 = read.table(url, header=TRUE)
  t1$node_names <- as.factor(paste(as.character(t1[, "Gene"]), as.character(t1[, "Node"]), sep = "_"))
  # create @metadata$node_list
  test1 <- t1[!duplicated(t1$node_names), ]
  
  gene_id <- as.character(test1$Gene[which(test1$node_names%in%rownames(matrix.scaled))])
  node_num <- as.character(test1$Node[which(test1$node_names%in%rownames(matrix.scaled))])
  
  # get gene is to search for names (i.e. external gene name) in ensembl
  gene <- gsub("\\.\\d+$", "", gene_id)
  coord <- as.character(test1$Coord[which(test1$node_names%in%rownames(matrix.scaled))])
  strand <- as.character(test1$Strand[which(test1$node_names%in%rownames(matrix.scaled))])
  type <- as.character(test1$Type[which(test1$node_names%in%rownames(matrix.scaled))])
  
  
  gene.node <- data.frame(
    "gene_id" = gene_id,
    "gene" = gene,
    "node_id" = rownames(matrix.scaled),
    "node_num" = node_num,
    "coord" = coord,
    "strand" = strand,
    "type" = type
  )
  
  # install.packages("XML", repos = "http://www.omegahat.net/R")
  # BiocManager::install("biomaRt")
  # library("biomaRt")
  # listMarts()
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
  
  # takes a few minuites to match gene to name
  gene_name <- getBM(attributes= c("ensembl_gene_id","external_gene_name"),
                     filters= "ensembl_gene_id",
                     values = gene,
                     mart= ensembl)
  
  # de-duplicate and match gene name using gene id
  gene_name <- gene_name[!duplicated(gene_name$ensembl_gene_id), ]
  
  gene.node <- merge(gene.node, gene_name, by.x = "gene", by.y = "ensembl_gene_id", all.x=TRUE)
  
  gene.node <- setNames(gene.node,c("gene", "gene_id", "node_id", "node_num", "coord", "strand", "type", "gene_name"))
  
  return(gene.node)
  
}

