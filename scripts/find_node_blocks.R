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


library(scfindME)
# devtools::load_all("~/scfindME")

index <- loadObject(index_path)


nodes <- scfindGenes(index)

all_genes <- levels(factor(nodeDetails(index, node.list = nodes)$external_gene_name))

all_blocks <- data.frame()

for (gene in all_genes){
    
    nodes <- geneNodes(index, gene, 'external_gene_name')
    
    
    tbl <- index@metadata$stats[which(rownames(index@metadata$stats) %in% nodes$Gene_node), ] %>% 
        tibble::rownames_to_column("Gene_node") %>% 
        mutate(node_num = as.numeric(gsub("^.*_", "",  Gene_node))) %>% 
        arrange(node_num)%>% 
        mutate(node_num_diff = ave(node_num, FUN=function(x) c(0, diff(x))))
                           
    details <- nodeDetails(index, tbl$Gene_node)
                           
    tbl <- merge(tbl, details, by = 'Gene_node') %>% arrange(node_num) %>% filter(Type %in% c("CE", "RI", "AA", "AD"))
    
    new_block <- data.frame()
    block_num <- 1
                                   
    if(nrow(tbl) > 2){                               
    for(i in seq(2, nrow(tbl))){
    if(abs(tbl[i, 'mean'] - tbl[i-1, 'mean']) < 0.1 &
       abs(tbl[i, 'SD'] - tbl[i-1, 'SD']) < 0.08){
        
        add_block <- rbind(tbl[i-1, ], tbl[i, ])
        add_block$block_num <- block_num
        new_block <- rbind(new_block, add_block)
        
        }
    
    else{
        
        block_num <- ifelse(nrow(new_block) == 0, 1, max(new_block$block_num) + 1)
    }

}

    new_block <- unique(new_block)
    
    all_blocks <- rbind(all_blocks, new_block) 
    
}
                                   }
                                   
    
                                   
    
saveRDS(all_blocks, paste(output, name, "all_node_blocks.rds", sep = "_"))