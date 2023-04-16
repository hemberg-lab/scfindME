#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)

option_list <- list(
  make_option(c("-n", "--data_name"),
    type = "character", default = NULL,
    help = "Name of dataset"
  ),
  make_option(c("-i", "--index"),
    type = "character", default = NULL,
    help = "Path to a complete scfindME index"
  ),
    make_option(c("-t", "--node_types", default='CE,AD,AA,NA,RI'),
    type = "character", default = NULL,
    help = "Node types to consider when detecting blocks, split by comma"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "Directory where discovered node blocks in all genes will be written as an .rds file"
  )
)

# parse input
opt <- parse_args(OptionParser(option_list = option_list))

name <- opt$data_name
index_path <- opt$index
output <- opt$output
types <- opt$node_types

# library(scfindME)
devtools::load_all("/nfs/research/irene/ysong/DATA/SCFIND/scfindME_package/scfindME")

index <- loadObject(index_path)

object <- index

stats <- object@metadata$stats

stats$node_id <- rownames(stats)

node.list <- object@metadata$node_list

all_genes <- as.character(levels(factor(node.list$Gene_name)))

types <- str_split(types, ",") %>% flatten_chr()

a <- merge(stats, node.list,
  by.x = "node_id",
  by.y = "Node_id",
  all.x = FALSE,
  all.y = FALSE
) %>% unique()


d <- data.frame()

node_types <- "node_types_all"

message("start processing all genes to find MXEs")

node_num <- 1

for (gene in all_genes) {
  nodes <- geneNodes(object, gene, "Gene_name")
  nodes$Type = as.character(nodes$Type)

  if (nrow(nodes) == 0) {
    break
  }

  if (nrow(nodes) > 1) {
    nodes_check <- nodes[which(nodes$Type %in% types), "Node_id"]

    if (length(nodes_check) >= 2) {
      pairs <- as.data.frame(combn(nodes_check, 2))

      for (i in seq(1, ncol(pairs))) {
        skip_to_next <- FALSE
        node_1 <- pairs[1, i]
        node_2 <- pairs[2, i]

        a_1 <- a[which(a$node_id %in% node_1), ]

        a_2 <- a[which(a$node_id %in% node_2), ]

        test_comb <- c(node_1, paste("-", node_2, sep = ""))

        test_comb_2 <- c(node_2, paste("-", node_1, sep = ""))

        if (0.9 <= (a_1[1, "mean"] + a_2[1, "mean"]) &
          (a_1[1, "mean"] + a_2[1, "mean"]) <= 1.1 &
          abs(a_1[1, "SD"] - a_2[1, "SD"]) < 0.1) {

          # at least in one cell typs it is specific

          condition <- tryCatch(
            {
              suppressMessages(sum(hyperQueryCellTypes(object, test_comb, datasets = "above")$pval < 0.05) >= 1 | sum(hyperQueryCellTypes(object, test_comb_2, datasets = "above")$pval < 0.05) >= 1)
            },
            error = function(e) {
              skip_to_next <<- TRUE
            }
          )

          if (skip_to_next) {
            next
          } else if (condition == TRUE) {
            # this is a promising mutually exclusive exon
            message("find a mutually exclusive exon pair that is cell type specific")
            
              if (sum(hyperQueryCellTypes(object, test_comb, datasets = "above")$pval < 0.05) >= 1 & sum(hyperQueryCellTypes(object, test_comb_2, datasets = "above")$pval < 0.05) >= 1) {
                  
                  sig_cell_types = suppressMessages(rbind(hyperQueryCellTypes(object, test_comb, datasets = "above")  %>% filter(pval < 0.05),  hyperQueryCellTypes(object, test_comb_2, datasets = "above")  %>% filter(pval < 0.05)))
                  
              } else if (sum(hyperQueryCellTypes(object, test_comb, datasets = "above")$pval < 0.05) >= 1){
                  
                  sig_cell_types = suppressMessages(hyperQueryCellTypes(object, test_comb, datasets = "above")  %>% filter(pval < 0.05))
                  
              } else if (sum(hyperQueryCellTypes(object, test_comb_2, datasets = "above")$pval < 0.05) >= 1){
                  
                  sig_cell_types = suppressMessages(hyperQueryCellTypes(object, test_comb_2, datasets = "above")  %>% filter(pval < 0.05))
                  
              }
           
            add <- rbind(a_1, a_2)
            add <- merge(add, sig_cell_types)
            add$node_num <- node_num
            node_num <- node_num + 1
            d <- rbind(d, add)
          }
        }
      }
    }
  }
}


saveRDS(d, paste(output, "/", name, "_mutually_exclusive_exons.rds", sep = ""))
