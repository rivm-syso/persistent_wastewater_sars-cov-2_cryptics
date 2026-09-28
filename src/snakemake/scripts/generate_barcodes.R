# Auke Haver
# 2026-09-14
# NRS / Z&O / I&V / RIVM

# Load packages -----------------------------------------------------------

library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(yaml)
library(arrow)
library(ape)
library(jsonlite)
library(optparse)

# Load arguments ----------------------------------------------------------

option_list <- list(
    make_option(
        c("--ancestral"),
        help="ancestral",
        default = "data/snakemake_output/B.1.221_default/ancestral_data.json"
    ),
    make_option(
        c("--tree"),
        help="tree file",
        default = "data/snakemake_output/B.1.221_default/refined_tree.nwk"
    ),
    make_option(
        c("--node"),
        help="node data",
        default = "data/snakemake_output/B.1.221_default/refined_node_data.json"
    ),
    make_option(
        c("--barcodes"),
        help="barcodes",
        default = "data/input/usher_barcodes.feather"
    ),
    make_option(
        c("--problematic_sites"),
        help="problematic_sites",
        default = "data/reference/problematic_sites_sarsCov2.vcf"
    ),
    make_option(
        c("--new_barcodes_full"),
        help="new_barcodes_full",
        default = "data/snakemake_output/B.1.221_default/barcodes_full.feather"
    ),
    make_option(
        c("--new_barcodes_known"),
        help="new_barcodes_known",
        default = "data/snakemake_output/B.1.221_default/barcodes_known.feather"
    ),
    make_option(
        c("--selected_nodes"),
        help = "selected_nodes",
        default = "data/snakemake_output/B.1.221_default/selected_nodes.txt"
    ),
    make_option(
        c("--cryptic_mutations"),
        help = "cryptic_mutations",
        default = "data/snakemake_output/B.1.221_default/cryptic_mutations.tsv"
    )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

args <- commandArgs(trailingOnly = TRUE)




# Load data ---------------------------------------------------------------

ls.ancestral <- read_json(opt$ancestral)
ls.tree      <- read.tree(opt$tree)
ls.node      <- read_json(opt$node)
df.barcodes  <- read_feather(opt$barcodes) %>% pivot_longer(
    cols = !index, names_to = "muts", values_to = "value"
)

gc()

df.problematic_sites <- 
    read_tsv(
        file = opt$problematic_sites,
        skip = 89,
        col_names = c(
            "REGION",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO"
        ),
        col_types = "cicccccc"
    )
ls.wuhan_bases <- ls.ancestral$reference[[1]] %>% str_split(pattern = "", simplify = TRUE)

ls.known_sites <- unique(parse_number(df.barcodes$muts))


# functions ---------------------------------------------------------------

fun.get_tree_df <- function(
        tree, 
        data
){
    
    df = as.data.frame(tree$edge)
    names(df) <- c("parent", "node")
    df = full_join(
        df,
        tibble(
            label = c(tree$tip.label,tree$node.label),
            node =  seq(1, length(label)),
            muts = lapply(label, function(x) data$nodes[[x]]$muts)
        ),
        by = "node"
    ) %>%
        (
            function(x) mutate(
                x,
                parent = sapply(parent, function(y) ifelse(is.na(y), NA, x$label[x$node == y]))
            )
        )(.) %>%
        select(-node)
    return(df)
}

add_mutations <- function(old, new){
    if(length(old)==0){
        return(new)
    } else if(
        length(new)==0
    ){
        return(old)
    }
    newsites = parse_number(new)
    oldsites = parse_number(old)
    
    return(
        sort(unlist(c(old, new[!newsites %in% oldsites])))
    )
}
add_mutation_recursively <- function(dataframe, current_label, muts=list()){
    
    current_position = dataframe$label == current_label
    current_muts     = dataframe$muts[current_position] %>% unlist()
    current_parent   = dataframe$parent[current_position]
    
    muts = unlist(add_mutations(muts, current_muts))
    
    if(is.na(current_parent)){
        # End here
        return(muts)
    } else {
        # Continue one level down
        return(add_mutation_recursively(dataframe, current_parent, muts))
    }
}

fun.check_mutations <- function(mutation_list, sequence = ls.wuhan_bases){
    if(length(mutation_list) ==0){
        return(c())
    } else {
        output = c()
        for(mutation in mutation_list){
            pos = parse_number(mutation)
            if(!pos %in% df.problematic_sites$POS[df.problematic_sites$FILTER == "mask"]){
                
                alt = gsub(".*[0-9]+", "", mutation)
                if(alt != sequence[pos]){
                    output = c(output, paste0(sequence[pos], pos, alt))
                }
            }
        }
        return(unlist(as.character(reorder(output, parse_number(output)))))
    }
}

# Load data ---------------------------------------------------------------

df.tree <- fun.get_tree_df(
    ls.tree,
    ls.ancestral
) %>%
    (function(dataframe){
        mutate(
            dataframe,
            muts = lapply(
                label,
                function(nodetip){
                    add_mutation_recursively(dataframe, nodetip)
                }
            )) %>%
            mutate(
                muts = lapply(
                    muts,
                    FUN = fun.check_mutations
                ),
                sample_date = sapply(
                    X = label,
                    FUN = function(x) ls.node$nodes[[x]]$date
                ))
    }
    )(.) 

# Make sure from the latest nodes, all nodes are one lineage
df.tree_nodes <- filter(df.tree, str_detect(df.tree$label, "NODE")) %>%
    (function(x){
        
        included_nodes <- c(x$label[x$sample_date == max(x$sample_date)])
        parent = x$parent[x$label == included_nodes[[1]]]
        
        while(!is.na(parent)){
            included_nodes = c(parent, included_nodes)
            
            parent = x$parent[x$label == included_nodes[[1]]]
        }
        
        return(
            filter(x, label %in% included_nodes)
        )
    })(.) 

selected_nodes = c()

for(date in seq(
    from = as.Date(min(df.tree$sample_date)),
    to = as.Date(max(df.tree$sample_date)),
    length.out = 8
)){
    selected_nodes = c(
        selected_nodes,
        df.tree_nodes %>%
            mutate(diff = abs(as.integer(as.Date(sample_date) - as.Date(date)))) %>%
            filter(diff == min(diff)) %>%
            head(1) %>%
            pull("label")
    )
    df.tree_nodes <- filter(df.tree_nodes, !label %in% selected_nodes)
}


df.barcodes_new <- 
df.tree %>%
    filter(label %in% selected_nodes) %>%
    arrange(sample_date) %>%
    reframe(
        index = paste0("CRYPT.",  row_number()),
        muts  = muts,
        value = 1
    ) %>%
    unnest(muts)

write_tsv(
    x = df.barcodes_new,
    file = opt$cryptic_mutations
)


# Write to file -----------------------------------------------------------

write_lines(
    x = selected_nodes,
    file = opt$selected_nodes
)


write_feather(
    x = rbind(df.barcodes, df.barcodes_new) %>%
        complete(index, muts, fill = list(value = 0)) %>%
        pivot_wider(names_from = muts, values_from = value),
    sink = opt$new_barcodes_full
)
gc()
write_feather(
    x = rbind(df.barcodes, filter(df.barcodes_new, parse_number(muts) %in% ls.known_sites)) %>%
        complete(index, muts, fill = list(value = 0)) %>%
        pivot_wider(names_from = muts, values_from = value),
    sink = opt$new_barcodes_known
)

