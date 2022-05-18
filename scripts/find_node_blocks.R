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
    make_option(c("-t", "--node_types"),
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

all_nodes <- scfindNodes(index)

all_genes <- levels(factor(nodeDetails(index, node.list = all_nodes)$Gene_name))

all_blocks <- data.frame()

node_types <- str_split(types, ",") %>% flatten_chr()

block_num <- 1

for (gene in all_genes) {
  message(gene)

  nodes <- geneNodes(index, gene, "Gene_name")

  tbl <- index@metadata$stats[which(rownames(index@metadata$stats) %in% nodes$Node_id), ] %>%
    tibble::rownames_to_column("Node_id") %>%
    mutate(node_num = as.numeric(gsub("^.*_", "", Node_id))) %>%
    arrange(node_num) %>%
    mutate(node_num_diff = ave(node_num, FUN = function(x) c(0, diff(x))))

  details <- nodeDetails(index, tbl$Node_id)

  tbl <- merge(tbl, details, by = "Node_id") %>%
    arrange(node_num) %>%
    filter(Type %in% node_types)

  new_block <- data.frame()
  block_num <- 1

  if (nrow(tbl) > 2) {
    
      # first node in block
    i <- 1

    mean <- tbl[i, "mean"]

    SD <- tbl[i, "SD"]

    new_block <- data.frame()

    while (i < nrow(tbl)) {
        
      if (abs(tbl[i + 1, "mean"] - mean) < 0.2 &
        abs(tbl[i + 1, "SD"] - SD < 0.1)) {
        
        message(i)
          
        add_block <- tbl[i + 1, ]

        add_block$block_num <- block_num
        
        mean <- mean(mean, add_block[, "mean"])
    
          
        SD <- mean(SD, add_block[, "SD"])


        if (nrow(new_block) == 0) {
            
          first_in_block <- tbl[i, ]
          first_in_block$block_num <- block_num
          new_block <- rbind(first_in_block, add_block)
            
        } else {
            
          new_block <- rbind(new_block, add_block)
            
        }
          
        i <- i + 1 ## keep adding nodes to this block
          
      } else { # nothing more to add for this group

        if (nrow(new_block) == 0) {
            
          i <- i + 1
          mean <- tbl[i, "mean"]
          SD <- tbl[i, "SD"]

            
            
        } else if (nrow(new_block) > 0) {

          # potential new block to be add
          # check cell types specificities
          # message('examining')
          # message(block_num)

          block_now <- block_num

          test_comb <- new_block %>%
            dplyr::filter(block_num == block_now) %>%
            select(Node_id)

          # print(test_comb)

          skip_to_next <- FALSE

          condition <- tryCatch(
            {
              sum(hyperQueryCellTypes(index, test_comb$Node_id)$pval < 0.05) >= 1
            },
            error = function(e) {
              skip_to_next <<- TRUE
            }
          ) # if the block is significant in some cell types

          if (condition == TRUE) {
              
            
            message(paste("add block No. ", block_num, sep = ""))
              
            sig_cell_types <- hyperQueryCellTypes(index, test_comb$Node_id) %>% filter(pval < 0.05)

            add_new_block <- merge(new_block, sig_cell_types, all = TRUE)

            all_blocks <- rbind(all_blocks, add_new_block)

            message(paste("finish block ", block_num, " in gene ", gene, sep = ""))

            block_num <- block_num + 1

            new_block <- data.frame()

            mean <- tbl[i + 1, "mean"]

            SD <- tbl[i + 1, "SD"]

            i <- i + 1
          } else {

            # message("skip block")
            # message(block_num)

            i <- i + 1

            mean <- tbl[i, "mean"]

            SD <- tbl[i, "SD"]

            # no need to update block num

            new_block <- data.frame()


            next
          }
        }
      }
    }
  }
                               
}

all_blocks <- unique(all_blocks)


saveRDS(all_blocks, paste(output, "/", name, "_all_node_blocks.rds", sep = ""))
