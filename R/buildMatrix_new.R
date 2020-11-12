#' Builds an original \code{matrix} from a given Whippet Collapsed Quant .psi.gz files containing directory and a metadata 
#'
#' @param dir a directory containing .psi.gz files from Whippet quantification, looks like "/lustre/scratch117/cellgen/team218/ys6/TM_FACS/MicroExonator/Whippet/Quant/Single_Cell/Unpooled/"
#' @param metadata metadata of the alternative splicing dataset, must include a column named "SRA" and has SRA ids for all Whippet quantified samples
#' @param column_label the column name which indicates cell id and .pai.gz file name in metadata, usually use "SRA" and use SRA ids
#' @name buildMatrix.original
#' @importFrom readr read_tsv
#' @importFrom dplyr select mutate distinct filter group ungroup 
#' @importFrom tidyr unite pivot_wider
#' @return a tibble object for further scailing
#'
buildMatrix.original <- function(file, Sample){
  
  data <- readr::read_tsv(file, col_names=TRUE, progress=show_progress())
  
  message("all values collected, generating matrix...")
  
  data_Psi <- data %>% 
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
  
  return(data_Psi)
}


scaleMatrix.z <- function(data_Psi, metadata){
  df <- data.frame(data_Psi, row.names = "Gene_node")
  dm <- as.matrix(df)
  data_scaled_z <- t(scale(t(dm)))
  return(data_scaled_z)
}

scaleMatrix.diff <- function(data_Psi, metadata){
  df <- data.frame(data_Psi, row.names = "Gene_node")
  dm <- as.matrix(df)
  mean <- rowMeans(dm, na.rm = TRUE)
  data_scaled_diff <- dm - mean
  return(data_scaled_diff)
}

