#' Builds an original \code{matrix} from a given Whippet Quant .psi.gz files containing directory and a metadata 
#'
#' @param dir a directory containing .psi.gz files from Whippet quantification, looks like "/lustre/scratch117/cellgen/team218/ys6/TM_FACS/MicroExonator/Whippet/Quant/Single_Cell/Unpooled/"
#' @param metadata metadata of the alternative splicing dataset, must include a column named "SRA" and has SRA ids for all Whippet quantified samples
#' @param column_label the column name which indicates cell id and .pai.gz file name in metadata, usually use "SRA" and use SRA ids
#' @name buildMatrix.original
#'
#' @return a matrix object(original)
#'
buildMatrix.original <- function(dir, metadata, column_label){
  url <- paste(as.character(dir), as.character(metadata[1, column_label]), ".psi.gz", sep = "")
  t1 = read.table(url, header=TRUE)
  t1$node_names <- as.factor(paste(as.character(t1[, "Gene"]), as.character(t1[, "Node"]), sep = "_"))
  node.list <- t1$node_names
  
  matrix <-  data.frame(node_names = as.factor(node.list))
  read.count <- data.frame(node_names = as.factor(node.list))
  
  # collect PSI values from input: SRA_id.psi.gz files from Whippet output
  for (run in metadata[[column_label]]){
    temp.run <- as.character(run)
    temp.url <- paste(as.character(dir), run, ".psi.gz", sep="")
    temp.matrix <- read.table(temp.url, header=TRUE)
    matrix[, temp.run] <- temp.matrix[, "Psi"]  
    read.count[, temp.run] <- temp.matrix[, "Total_Reads"]
    
  } 
  
  
  matrix<- matrix[!duplicated(matrix$node_names), ]
  
  read.count <- read.count[!duplicated(read.count$node_names), ]
  
  # assign node_names as rownames of the matrix
  rownames(matrix) <- matrix$node_names 
  matrix <- matrix[, -1]
  
  rownames(read.count) <- read.count$node_names 
  read.count <- read.count[, -1]
  
  
  # then remove nodes that only have NA values across all samples
  matrix.1 <- matrix[which(rowSums(is.na(matrix)) < ncol(matrix)), ]
  read.count.1 <- read.count[which(rownames(read.count)%in%rownames(matrix.1)), ]
  
  
  # if read count < 10, CI for PSI is too wide, so set to NA
  
  matrix.2 <- data.frame(row.names = rownames(matrix.1))
  
  for (cell in seq(1, ncol(matrix.1))){
    tv <- matrix.1[, cell]
    tc <- read.count.1[, cell]
    temp <- which(tc <= 10)
    tv[temp] = NA
    matrix.2[[as.character(colnames(matrix.1)[cell])]] <- tv
  }
  
  return(matrix.2)
}

#' Builds a read_count \code{matrix} from a given Whippet Quant .psi.gz files containing directory and a metadata 
#' The read_count matrix would be stored in the alternative splicing index metadata slot
#' @param dir a directory containing .psi.gz files from Whippet quantification
#' @param metadata metadata of the alternative splicing dataset, must include a column named SRA and has SRA ids for all Whippet quantified samples
#' @param column_label the column name which indicates cell id and .pai.gz file name in metadata, usually use "SRA" and use SRA ids
#' @name buildMatrix.read_count
#'
#' @return a matrix object(read count)
#'
#' @importFrom Matrix Matrix
#' 
buildMatrix.read_count <- function(dir, metadata, column_label){
  url <- paste(as.character(dir), as.character(metadata[1, column_label]), ".psi.gz", sep = "")
  t1 = read.table(url, header=TRUE)
  t1$node_names <- as.factor(paste(as.character(t1[, "Gene"]), as.character(t1[, "Node"]), sep = "_"))
  node.list <- t1$node_names
  
  matrix <-  data.frame(node_names = as.factor(node.list))
  read.count <- data.frame(node_names = as.factor(node.list))
  
  # collect PSI values from input: SRA_id.psi.gz files from Whippet output
  for (run in metadata[[column_label]]){
    temp.run <- as.character(run)
    temp.url <- paste(as.character(dir), run, ".psi.gz", sep="")
    temp.matrix <- read.table(temp.url, header=TRUE)
    matrix[, temp.run] <- temp.matrix[, "Psi"]  
    read.count[, temp.run] <- temp.matrix[, "Total_Reads"]
    
  } 
  
  
  matrix<- matrix[!duplicated(matrix$node_names), ]
  
  read.count <- read.count[!duplicated(read.count$node_names), ]
  
  # assign node_names as rownames of the matrix
  rownames(matrix) <- matrix$node_names 
  matrix <- matrix[, -1]
  
  rownames(read.count) <- read.count$node_names 
  read.count <- read.count[, -1]
  
  
  # then remove nodes that only have NA values across all samples
  matrix.1 <- matrix[which(rowSums(is.na(matrix)) < ncol(matrix)), ]
  read.count.1 <- read.count[which(rownames(read.count)%in%rownames(matrix.1)), ]
  
  
  # if read count < 10, CI for PSI is too wide, so set to NA
  
  matrix.2 <- data.frame(row.names = rownames(matrix.1))
  
  for (cell in seq(1, ncol(matrix.1))){
    tv <- matrix.1[, cell]
    tc <- read.count.1[, cell]
    temp <- which(tc <= 10)
    tv[temp] = NA
    matrix.2[[as.character(colnames(matrix.1)[cell])]] <- tv
  }
  
  
  # scailing by z-score normalization 
  # this could take several mins
  matrix.scaled <- t(scale(t(matrix.2)))
  
  matrix.scaled <- matrix.scaled*100
  
  # remove all NA rows 
  matrix.scaled <- matrix.scaled[which(rowSums(is.na(matrix.scaled)) < ncol(matrix.scaled)), ]
  read.count <- read.count.1[which(rownames(read.count.1)%in%rownames(matrix.scaled)), ]
  
  # create @metadata$read_count
  # make the read.count matrix sparse matrix 
  read.count$mean_count <- rowMeans(read.count, na.rm = TRUE)
  
  read.count[which(is.na(read.count))] <- 0
  
  
  # library(Matrix)
  read.count <- Matrix::Matrix(as.matrix(read.count), sparse = TRUE)
  
  return(read.count)
}
  



#' Scale the original matrix and build the input matrix of an above index
#' The output of this function correspond to the input param "psival" in the buildAltSpliceIndex method
#' @param matrix.original the output original matrix from the function "buildMatrix.original"
#'
#' @name buildMatrix.scaled.above
#'
#' @return a matrix object(scaled above)
#' 
buildMatrix.scaled.above <- function(matrix.original){
  
  # scailing by z-score normalization 
  # this could take several mins
  matrix.scaled <- t(scale(t(matrix.original)))
  
  matrix.scaled <- matrix.scaled*100
  
  # remove all NA rows 
  matrix.scaled <- matrix.scaled[which(rowSums(is.na(matrix.scaled)) < ncol(matrix.scaled)), ]

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
#' @name buildMatrix.scaled.below
#'
#' @return a matrix object(scaled below)
#' 
buildMatrix.scaled.below <- function(matrix.original){
  
  # scailing by z-score normalization 
  # this could take several mins
  matrix.scaled <- t(scale(t(matrix.original)))
  
  matrix.scaled <- matrix.scaled*100
  
  # remove all NA rows 
  matrix.scaled <- matrix.scaled[which(rowSums(is.na(matrix.scaled)) < ncol(matrix.scaled)), ]
  
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


