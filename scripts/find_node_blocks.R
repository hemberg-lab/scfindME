#!/usr/bin/env R

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


# library(scfindME)
devtools::load_all("~/scfindME")

index <- loadObject(index_path)


nodes <- scfindGenes(index)

all_genes <- levels(factor(nodeDetails(index, node.list = nodes)$external_gene_name))

all_blocks <- data.frame()

for (gene in all_genes) {
  message(gene)


  nodes <- geneNodes(index, gene, "external_gene_name")


  tbl <- index@metadata$stats[which(rownames(index@metadata$stats) %in% nodes$Gene_node), ] %>%
    tibble::rownames_to_column("Gene_node") %>%
    mutate(node_num = as.numeric(gsub("^.*_", "", Gene_node))) %>%
    arrange(node_num) %>%
    mutate(node_num_diff = ave(node_num, FUN = function(x) c(0, diff(x))))

  details <- nodeDetails(index, tbl$Gene_node)

  tbl <- merge(tbl, details, by = "Gene_node") %>%
    arrange(node_num) %>%
    filter(Type %in% c("CE", "RI", "AA", "AD"))

  new_block <- data.frame()
  block_num <- 1



  if (nrow(tbl) > 2) {
    i <- 1

    block_num <- 1

    mean <- tbl[i, "mean"]

    SD <- tbl[i, "SD"]

    new_block <- data.frame()

    while (i < nrow(tbl)) {
      if (abs(tbl[i + 1, "mean"] - mean) < 0.2 &
        abs(tbl[i + 1, "SD"] - SD < 0.1)) {
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

        # print(new_block)

        i <- i + 1
      } else { # nothing more to add for this group

        if (nrow(new_block) == 0) {
          i <- i + 1

          mean <- tbl[i + 1, "mean"]

          SD <- tbl[i + 1, "SD"]
        } else if (nrow(new_block) > 0) {

          # potential new block to be add
          # check cell types specificities
          # message('examining')
          # message(block_num)


          block_now <- block_num


          test_comb <- new_block %>%
            dplyr::filter(block_num == block_now) %>%
            select(Gene_node)

          # print(test_comb)

          skip_to_next <- FALSE

          condition <- tryCatch(
            {
              sum(hyperQueryCellTypes(object, test_comb$Gene_node)$pval < 0.05) >= 1
            },
            error = function(e) {
              skip_to_next <<- TRUE
            }
          ) # if the block is significant in some cell types

          if (condition == TRUE) {
            message("add block No.")
            message(block_num)

            sig_cell_types <- hyperQueryCellTypes(object, test_comb$Gene_node) %>% filter(pval < 0.05)

            add_new_block <- merge(new_block, sig_cell_types, all = TRUE)

            all_blocks <- rbind(all_blocks, add_new_block)

            message(paste("finish", block_num))

            block_num <- block_num + 1

            new_block <- data.frame()

            mean <- tbl[i + 1, "mean"]

            SD <- tbl[i + 1, "SD"]

            i <- i + 1
          } else {

            # message("skip block")
            # message(block_num)

            i <- i + 1

            mean <- tbl[i + 1, "mean"]

            SD <- tbl[i + 1, "SD"]

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


saveRDS(all_blocks, paste(output, name, "all_node_blocks.rds", sep = "_"))
