#!/usr/bin/R

library(optparse)
library(tidyverse)

option_list <- list( 
    make_option(c("-n", "--data_name"), type = "character", default=NULL,
        help="Name of dataset"),
    make_option(c("-i", "--index"), type = "character", default=NULL,
        help="Path to a complete scfindME index"),
    make_option(c("-o", "--output"), type = "character", default=NULL,
        help="Directory where discovered node blocks in all genes will be written as an .rds file")
    )

# parse input
opt <- parse_args(OptionParser(option_list=option_list))

name=opt$data_name
index_path=opt$index
output=opt$output


#library(scfindME)
devtools::load_all("~/scfindME")

index <- loadObject(index_path)
object <- index

stats <- object@metadata$stats

stats$node_id <- rownames(stats)
    
node.list <- object@metadata$node_list

all_genes <- as.character(levels(factor(node.list$external_gene_name)))
    
  a <- merge(stats, node.list, by.x = "node_id", 
             by.y = "Gene_node", 
             all.x = FALSE, 
             all.y = FALSE) %>% unique()
  
  b <- data.frame(row.names = a$node_id)
  d <- data.frame(row.names = a$node_id)

node_types <- 'node_types_all'

for(gene in all_genes){
    
    
    nodes <- geneNodes(index, gene,"external_gene_name")
    
    if(nrow(nodes) == 0) {break}
    
    if(nrow(nodes) > 1){
    
    nodes_check <- nodes[which(nodes$Type %in% c("CE", "AA", "AD","RI","AF","AL")), 'Gene_node']
    
    
    pairs <- as.data.frame(combn(nodes_check,2))
    
    
    for(i in seq(1, ncol(pairs))){
    
    node_1 <- pairs[1, i]
    node_2 <- pairs[2, i]
    
    a_1 <- a[which(a$node_id %in% node_1), ]
    
    a_2 <- a[which(a$node_id %in% node_2), ]

    test_comb <- c(node_1, paste("-", node_2, sep = ""))
    
    test_comb_2 <- c(node_2, paste("-", node_1, sep = ""))
    
    if(0.9 <= (a_1[1, "mean"] + a_2[1, "mean"]) & 
       (a_1[1, "mean"] + a_2[1, "mean"]) <= 1.1 & 
       abs(a_1[1, "SD"]- a_2[1, "SD"])<0.1){
        
            # at least in one cell typs it is specific
        
            if(sum(hyperQueryCellTypesAS(object, test_comb)$pval < 0.1) > 1 | sum(hyperQueryCellTypesAS(object, test_comb_2)$pval < 0.1) > 1){
                
            
                # this is a promising mutually exclusive exon
                d <-  rbind(d, a_1, a_2)
                
            }

  
        }
    }
    
 }
}

saveRDS(d,  paste(output, name, "mutually_exclusive_exons.rds", sep = "_"))