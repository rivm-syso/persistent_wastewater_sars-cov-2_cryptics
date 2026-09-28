# Auke Haver
# 2026-09-14
# NRS / Z&O / I&V / RIVM

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# General IO --------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#' Read IVAR output TSV file
#'
#' Reads an IVAR-formatted TSV file containing variant information and returns
#' a tibble with predefined column names and types.
#'
#' @param filename Character. Path to the IVAR TSV file.
#'
#' @return A tibble with columns such as REGION, POS, REF, ALT, depth metrics,
#'   quality metrics, codon/AA information, and annotation fields.
#'
#' @seealso \code{\link[readr]{read_tsv}}
fun.read_ivar_file <- function(filename){
  read_tsv(
    file = filename,
    col_names = c(
      "REGION",
      "POS",
      "REF",
      "ALT",
      "REF_DP",
      "REF_RV",
      "REF_QUAL",
      "ALT_DP",
      "ALT_RV",
      "ALT_QUAL",
      "ALT_FREQ",
      "TOTAL_DP",
      "PVAL",
      "PASS",
      "GFF_FEATURE",
      "REF_CODON",
      "REF_AA",
      "ALT_CODON",
      "ALT_AA",
      "POS_AA"
    ),
    col_types =  "ciccinninnninlcccccc",
    skip = 1
  )
}

#' Read depth TSV file
#'
#' Reads a TSV file with per-position depth information and returns it as a
#' tibble.
#'
#' @param filename Character. Path to the depth TSV file.
#'
#' @return A tibble with columns REGION, POS, REF, and TOTAL_DP.
#'
#' @seealso \code{\link[readr]{read_tsv}}
fun.read_depth_file <- function(filename){
  read_tsv(
    file = filename,
    col_types = "cici",
    col_names = c("REGION", "POS", "REF", "TOTAL_DP")
  )
}

#' Write analysis inputs to disk
#'
#' Creates an analysis-specific directory and writes metadata and sequence
#' files derived from a data frame.
#'
#' @param df Data frame containing at least \code{location}, \code{sample_date},
#'   and \code{sequence} columns.
#' @param inputdir Character. Base directory where analysis subdirectories are
#'   created.
#' @param lineage Character. Lineage identifier used in the directory name.
#' @param analysis_name Character. Analysis name used in the directory name.
#'
#' @return Invisibly, the path to the created analysis directory.
#'
#' @seealso \code{\link[readr]{write_tsv}}, \code{\link[base]{writeLines}},
#'   \code{\link[base]{dir.create}}
write_analysis_inputs <- function(df, inputdir, lineage, analysis_name) {
  path <- paste0(inputdir, "/", lineage, "_", analysis_name)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

    # Metadata
  df %>%
    reframe(
      name = paste0(location, "_", sample_date),
      date = lubridate::decimal_date(sample_date)
    ) %>%
    write_tsv(file = file.path(path, "metadata.tsv"))
  
  # Sequences
  writeLines(
    text = paste0(
      ">", df$location, "_", df$sample_date, "\n", df$sequence,
      collapse = "\n"
    ),
    con = file.path(path, "sequences.fasta")
  )
}

#' Expand lineage alias to full lineage name
#'
#' Converts a possibly aliased lineage string to its full designation using a
#' provided alias list.
#'
#' @param string Character. Lineage string that may contain an alias prefix.
#' @param alias_list Named list or environment mapping alias prefixes to full
#'   lineage roots.
#'
#' @return Character. The full lineage string if an alias was found; otherwise
#'   the original string.
fun.return_full_lineage <- function(string, alias_list){
  prefix = gsub("\\..*", "", string)
  
  if(is.na(prefix) | startsWith(prefix, "X") | prefix == "B" | prefix == "A" | prefix == "not_reported"){
    return(string)
  } else{
    return(
      gsub(prefix, alias_list[[prefix]], string)
    )
  }
}

#' Match lineages including sublineages
#'
#' Tests whether a query lineage is identical to or a sublineage of a subject
#' lineage.
#'
#' @param query Character. Query lineage (e.g. \code{"BA.5.1"}).
#' @param subject Character. Subject lineage (e.g. \code{"BA.5"}).
#'
#' @return Logical. \code{TRUE} if \code{query} equals \code{subject} or is a
#'   sublineage (has \code{subject} as prefix followed by a dot), otherwise
#'   \code{FALSE}.
fun.match_lineages <- function(query, subject){
  # Matching needs to be either complete (e.g. BA.5=BA.5 or with sublineage BA.5 - BA.5.1, not BA.50)
  return(query == subject | startsWith(query, fixed(paste0(subject, "."))))
  
}


#' Alias lineage designation
#'
#' Converts full lineage names to their aliased short forms using a provided
#' alias mapping.
#'
#' @param string Character. Full lineage designation.
#' @param alias_list Named list or environment with lineage aliases, typically
#'   from \code{fun.load_lineage_alias()}.
#'
#' @return Character. Aliased lineage designation if a match is found; otherwise
#'   the original string.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' alias_info <- fun.load_lineage_alias()
#' alias_lineage("B.1.1.7", alias_info)
#' }
alias_lineage <- function(string, alias_list) {
  for (i in rev(names(alias_list))) {
    if (startsWith(string, fixed(rev(alias_list)[[i]])) && string != rev(alias_list)[[i]]) {
      return(
        sub("\\.$", "", sub("\\.\\.", ".", sub(rev(alias_list)[[i]], paste0(i, "."), string)))
      )
    }
  }
  return(string)
}

#' Bin lineage into predefined groups
#'
#' Assigns a lineage string to one of a set of predefined lineage bins, with
#' optional handling of recombinant lineages.
#'
#' @param string Character. Lineage string to be binned.
#' @param lineage_bins Character vector. Lineage bin prefixes to match against.
#' @param value Numeric. Value or weight assigned to the lineage within its bin
#'   (default \code{1}).
#' @param lineage_info List. Lineage alias information, typically from
#'   \code{fun.load_lineage_alias()}.
#' @param include_recombinants Logical. If \code{TRUE}, recombinant lineages
#'   (starting with \code{"X"}) are decomposed into constituent sublineages and
#'   distributed across bins.
#'
#' @return A tibble with columns \code{lineage_dealiased}, \code{bin}, and
#'   \code{frac}, possibly containing multiple rows for recombinant lineages.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bins <- c("BA.5", "BQ.1", "XBB")
#' info <- fun.load_lineage_alias()
#' bin_lineage("BA.5.1", bins, lineage_info = info)
#' }
bin_lineage <- function(string, lineage_bins, value = 1, lineage_info = fun.load_lineage_alias(), include_recombinants = TRUE){
  
  lineage_bins <- c((function(y){y[order(-nchar(y))]})(lineage_bins) )
  
  string = fun.return_full_lineage(string, lineage_info)
  
  for(i in lineage_bins){
    if(startsWith(string, fixed(i))){
      if(string == fixed(i) | startsWith(string, fixed(paste0(i, ".")))){
        return(
          tibble(
            lineage_dealiased = string,
            bin = i,
            frac = value
          ))}}
  }
  # If we are dealing with a recombinant
  if(include_recombinants & startsWith(string, "X")){
    recombinant_sublineages = gsub("\\*", "", unlist(lineage_info[[gsub("\\..*", "", string)]]))
    
    return(
      do.call(
        rbind,
        lapply(
          recombinant_sublineages,
          FUN = function(x) bin_lineage(x, lineage_bins, value = value/length(recombinant_sublineages), lineage_info = lineage_info)
        )
      ) %>%
        mutate(lineage_dealiased = string) %>%
        group_by(lineage_dealiased, bin)  %>%
        reframe(frac = sum(frac)))
  }
}


#' Read Freyja-Format TSV File and Parse Lineage Abundances
#'
#' Reads a Freyja output TSV file containing lineage abundance estimates and related metrics, and returns a tidy data frame with one row per lineage per file.
#'
#' This function parses the Freyja TSV format, extracts lineage names, abundance estimates, and other relevant columns, and returns a data frame suitable for downstream analysis.
#'
#' @param filename Character. Path to the Freyja-format TSV file to read.
#' @param add_filename Logical. If TRUE (default), include a \code{filename} column in the output.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{filename}{File name (if \code{add_filename} is TRUE)}
#'     \item{lineage}{Lineage or strain name}
#'     \item{abundance}{Estimated relative abundance (numeric)}
#'     \item{resid}{Residual error (numeric)}
#'     \item{coverage}{Fraction of variant sites passing depth cutoff (numeric)}
#'   }
#'
#' @details
#' The function removes the "summarized" attribute if present, cleans whitespace and special characters from the parsed values, splits multi-lineage fields into separate rows, and coerces abundance, residual, and coverage columns to numeric types.
#'
#' @seealso \code{\link[readr]{read_tsv}}
#' @examples
#' \dontrun{
#' freyja_data <- read_freyja_tsv("sample_freyja_output.tsv")
#' head(freyja_data)
#' }
#' @export
read_freyja_tsv <- function(filename, add_filename=TRUE){
    suppressMessages(
        data <- read_tsv(
            file = gsub("", "", filename),
            col_names = c("attribute", "value"),
            skip = 1,
            col_types = "cc"
        ) %>%
            filter(attribute != "summarized") %>%
            mutate(
                filename = filename,
                value = gsub("\n|\\[|\\]|\'| $", "", value),
                value = gsub(" +", " ", value),
                value = trimws(value, which = "both")
            ) %>%
            pivot_wider(
                names_from = attribute,
                values_from = value
            ) %>%
            mutate(
                lineage_alias = str_split(lineages, " "),
                abundances = str_split(abundances, " ")
            ) %>%
            unnest(c(lineage_alias, abundances)) %>%
            reframe(
                filename,
                lineage = lineage_alias,
                abundance = as.numeric(abundances),
                resid = as.numeric(resid),
                coverage = as.numeric(coverage)
            )) %>%
        distinct()
    if(!add_filename){
        data <- data %>%
            select(-filename)
    }
    return(data)
}

#' Load lineage alias mapping
#'
#' Downloads or reads from disk the Pango lineage alias key JSON and returns
#' it as an R list.
#'
#' If the local file \code{"data/processed/alias_key.json"} does not exist, it
#' is downloaded from the cov-lineages GitHub repository and cached.
#'
#' @return List containing lineage alias mappings.
#'
#' @seealso \code{\link[jsonlite]{read_json}}, \code{\link[jsonlite]{write_json}}
fun.load_lineage_alias <- function(){
  if(!file.exists("data/processed/alias_key.json")){
    jsondata <- read_json("https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json")
    write_json(x = jsondata, path = "data/processed/alias_key.json")
    return(jsondata)
  } else{
    return(read_json("data/processed/alias_key.json"))
  }
}

#' Read barcodes feather file
#'
#' Reads a barcodes file in Feather format and reshapes it to long format,
#' keeping only non-zero entries.
#'
#' @param filepath Character. Path to the Feather file.
#'
#' @return A tibble with columns \code{index}, \code{muts}, and \code{vals}
#'   for non-zero entries.
#'
#' @seealso \code{\link[arrow]{read_feather}}, \code{\link[tidyr]{pivot_longer}}
fun.read_barcodes_file <- function(filepath){
  read_feather(filepath) %>%
    pivot_longer(
      cols = !index,
      names_to = "muts",
      values_to = "vals"
    ) %>%
    filter(vals != 0)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE 1 ------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#' Create abundance legend grob
#'
#' Builds a standalone legend for abundances using a viridis color scale,
#' suitable for inclusion in composite plots.
#'
#' @param legend_direction Character. Legend direction
#'   (\code{"horizontal"} or \code{"vertical"}).
#' @param tick_direction Integer. Direction argument passed to
#'   \code{scale_fill_viridis_c()}.
#' @param key_width Numeric. Legend key width in points.
#' @param key_height Numeric. Legend key height in points.
#'
#' @return A \code{ggdraw} object containing only the legend.
#'
#' @seealso \code{\link[cowplot]{get_legend}}, \code{\link[cowplot]{ggdraw}}
fun.abundance_legend <- function(legend_direction = "horizontal", tick_direction = 1, key_width = 48, key_height = 12){
  ggdraw(
    get_legend(
      ggplot()+
        geom_point(aes(x= 1,y=1,fill=.2), pch = 21) + 
        scale_fill_viridis_c(limits= c(0.001,1), name = NULL,direction = tick_direction, breaks = c(0.001,seq(.2,1,.2)))+
        theme(legend.direction = legend_direction, legend.key.width = unit(key_width, "pt"), legend.key.height = unit(key_height, "pt"))
    )
  )
}

#' Plot viral load ratio versus lineage abundance
#'
#' Produces a scatter plot of log10 viral load ratios over time, overlaid with
#' a ribbon indicating confidence intervals and colored by lineage abundance.
#'
#' @param df.log10_ratio_confint Data frame with columns including
#'   \code{Date_measurement} and \code{maxratio} for the confidence ribbon.
#' @param df.wwtp_ratioes Data frame with viral load ratios per measurement
#'   (must contain \code{RWZI_AWZI_name}, \code{Date_measurement}, and
#'   \code{ratio}).
#' @param df.lineage_estimates Data frame with lineage abundance estimates
#'   (must contain \code{location}, \code{sample_date}, \code{bin},
#'   \code{abundance}, \code{sample_hash}).
#' @param target_location Character. Target location name to subset data.
#' @param target_lineage Character. Target lineage bin to highlight.
#' @param layout_commands A \code{ggplot2} theme or layout component appended
#'   to the plot.
#'
#' @return A \code{ggplot} object showing viral load ratio and lineage abundance
#'   over time.
fun.plot_vl_ratio <- function(
    df.log10_ratio_confint,
    df.wwtp_ratioes,
    df.lineage_estimates,
    target_location,
    target_lineage,
    layout_commands){
  ggplot() + 
    geom_ribbon(
      data = df.log10_ratio_confint,
      mapping = aes(
        x = Date_measurement,
        ymin = -1.1,
        ymax = maxratio
      ),
      fill = "grey"
    )+
    geom_point(
      data = left_join(
        filter(df.wwtp_ratioes, RWZI_AWZI_name == c(target_location)),
        df.lineage_estimates %>%
          filter(location == target_location) %>%
          complete(nesting(sample_hash, location, sample_date), bin, fill = list(abundance = 0)) %>%
          filter(bin==target_lineage) %>%
          group_by(sample_hash,location, sample_date, bin) %>% 
          reframe(abundance = sum(abundance)),
        by = c("RWZI_AWZI_name" = "location", "Date_measurement" = "sample_date")
      ),
      mapping = aes(x = Date_measurement, y = ifelse(log10(ratio) < -1, -1, log10(ratio)), fill = abundance), pch = 21) +
    scale_fill_viridis_c() +
    scale_y_continuous(
      limits = c(-1.1, 2.6),
      breaks = seq(-1,2.5, .5),
      expand = c(0,0)
    ) +
    layout_commands
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE 2 ------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#' Create cryptic mutation alluvial plot
#'
#' Builds an alluvial plot of cryptic mutations for a given lineage, separating
#' known and novel sites and visualizing presence/absence patterns across index
#' positions.
#'
#' @param target_lineage Character. Target lineage identifier used to locate
#'   input files.
#'
#' @return A \code{ggplot} object representing the alluvial structure of cryptic
#'   mutations.
#'
#' @seealso \code{\link[tidyverse]{pivot_longer}}, \code{\link[ggalluvial]{stat_flow}}
fun.make_cryptic_alluvial <- function(target_lineage){
    df.alluvial <- 
        read_tsv(file = file.path("data/snakemake_output",paste0(target_lineage, "_default"), "cryptic_mutations.tsv"), col_types = "cci") %>%
        mutate(sites = ifelse(parse_number(muts) %in% ls.known_sites, "known", "novel")) %>%
        pivot_wider(names_from = "index", values_from = value, values_fill = 0) %>%
        group_by(across(-c(muts))) %>%
        reframe(count = n()) %>%
        mutate(category = row_number()) %>%
        pivot_longer(
            cols = !c(count, category, sites),
            names_to = "node",
            values_to = "status"
        )%>%
        mutate(status = factor(ifelse(status == 1, "Present", "Absent"), levels=  c("Present", "Absent"))) 
    
    ggplot(
        data = df.alluvial %>%
            group_by(node, status,  sites) %>%
            mutate(sum = sum(count)) %>% ungroup() ,
        mapping = aes(
            x = node,
            y = count,
            stratum = status,
            alluvium = category,
            fill = status,
            label = sum
        )
    )+
        stat_flow( color = "black", alpha = .5) +
        stat_stratum(color = "black", alpha = .5) +
        geom_text(stat = 'stratum', angle = 0, size = 3)+
        scale_fill_viridis_d(begin = .2, end = .8, direction = -1)+
        facet_grid(rows = vars(sites), scales = "free", switch="both")   +
        theme_void()+
        theme(
            legend.position = "none",
            axis.text.x = element_blank(),
            strip.text = element_text(angle = 90),
            strip.background = element_blank(),
            plot.margin =unit(rep(10,4), unit = "mm")
        )
    
}

#' Plot cryptic tree with highlighted path
#'
#' Generates a time-scaled phylogenetic tree for a target lineage, highlighting
#' selected nodes and the path between them.
#'
#' @param target_lineage Character. Target lineage identifier used to locate
#'   tree and node data files.
#' @param daterange Date vector. Range of dates used to set and label the
#'   x-axis.
#'
#' @return A \code{ggplot} object of the tree with highlighted path and nodes.
#'
#' @seealso \code{\link[ape]{read.tree}}, \code{\link[ggtree]{ggtree}}
fun.make_cryptic_tree <- function(target_lineage, daterange){
    
    ls.selected_nodes <- readLines(
        con  = file.path("data/snakemake_output",paste0(target_lineage, "_default"), "selected_nodes.txt")
    )
    
    ls.tree <- read.tree(
        file = file.path("data/snakemake_output",paste0(target_lineage, "_default"), "refined_tree.nwk")
    )
    
    ls.node_data <- read_json(
        path = file.path("data/snakemake_output",paste0(target_lineage, "_default"), "refined_node_data.json")
    )
    
    # Converting edge lengths into time (days)
    ls.tree$edge.length <- ls.tree$edge.length / ls.node_data$clock$rate*365.25
    
    ls.dates <- sort(unique(sapply(names(ls.node_data$nodes), function(x){as.Date(ls.node_data$nodes[[x]]$date)})))
    
    baseplot <- ggtree(ls.tree, color = NA)
    
    tree_tibble <- as_tibble(ls.tree) 
    
    node <- ls.selected_nodes[[8]] # there are always only eight, so this selects the last
    path = c(node)
    
    while(!node == ls.selected_nodes[[1]]){
        node = tree_tibble$label[tree_tibble$node == tree_tibble$parent[tree_tibble$label == node]]
        path = c(path, node)
    }
    
    return(
        baseplot+
            geom_tree(
                data = baseplot$data,
                color = "black",
                linetype = "solid"
            ) +
            geom_tree(
                data = baseplot$data %>% filter(label %in% path),
                color = "#21918c",
                linetype = "solid"
            ) +
            geom_point(
                data = baseplot$data %>% filter(label %in% ls.selected_nodes),
                mapping = aes(
                    x = x,
                    y = y
                ),
                pch = 21,
                fill = "#21918c",
                size = 2
            )  +
            scale_x_continuous(
                limits = c(
                    min(as.numeric(daterange - min(ls.dates))),
                    max(as.numeric(daterange - min(ls.dates)))
                ),
                breaks = as.numeric(daterange - min(ls.dates)),
                labels = gsub("-", "\n", strftime(daterange, "%b-%Y")),
                expand = c(0,0),
                name = NULL
            )+
            
            theme_bw() +
            theme(
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                plot.margin = margin(t=10,b=10,l=10,r=10, unit = "mm")
            )
    )
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE 3 ------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#' Create cryptic lineage heatmap
#'
#' Creates a heatmap of lineage abundances over time for a given lineage and
#' location, distinguishing cryptic and other bins.
#%
#' @param target_lineage Character. Target lineage identifier used to locate
#%   input files.
#' @param target_location Character. Target location name used to filter
#'   samples.
#'
#' @return A \code{ggplot} object showing abundance heatmap across dates and
#'   lineage bins.
#'
#' @seealso \code{\link[readr]{read_tsv}}, \code{\link[ggplot2]{geom_tile}}
fun.make_cryptic_heatmap <- function(target_lineage, target_location){
    freyja_file <- file.path(ls.parameters$output_dir, paste0(target_lineage, "_default"), "barcodes_known.tsv")
    
    freyja_data <- read_tsv(
        file = freyja_file,
        col_types = "ccnnnc"
    )
    
    return(
        freyja_data %>%
            mutate(
                bin = sapply(lineage, function(x){
                    if(x == target_lineage | str_detect(x, paste0(target_lineage, ".")) | str_detect(x, "CRYPT")){
                        return(x)
                    } else {
                        return("other")
                    }
                })
            ) %>%
            group_by(sample_hash, bin) %>%
            reframe(abundance = sum(abundance)) %>%
            left_join(df.metadata %>% distinct(sample_hash, location, sample_date), by = "sample_hash") %>%
            filter(location == target_location) %>%
            complete(
                sample_date,
                bin,
                fill = list(abundance = 0.001)
            )%>%
            mutate(
                bingroup = ifelse(str_detect(bin, "CRYPT"), "a", "b"),
                bin = gsub("CRYPT.", "", bin),
                bin = factor(bin, levels = c("other", unique(bin)[unique(bin) != "other"]))
            ) %>%
            ggplot(aes(
                x = as.factor(sample_date),
                y = bin,
                fill = abundance
            )) +
            geom_tile(
                color = "black",
                alpha = .8
            ) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 3),
                legend.position = "none",
                strip.background = element_blank(),
                strip.text = element_blank(),
                axis.title = element_blank()
            ) +
            facet_grid(rows = vars(bingroup), scales = "free_y", space = "free_y") +
            scale_fill_viridis_c(
                breaks = c(0.001, seq(.2, 1, .2)), 
                labels = format(c(0.001, seq(.2, 1, .2)), nsmall = 3), limits = c(0.001,1), name = NULL) +
            scale_x_discrete(expand=  c(0,0)) +
            scale_y_discrete(expand = c(0,0)) 
    )
}

#' Plot bootstrap values on phylogenetic tree
#'
#' Reads a bootstrap-annotated tree file and plots bootstrap support values on
#' internal nodes using a color scale.
#'
#' @param treefile Character. Path to the tree file in Newick format.
#' @param width Numeric. Horizontal plot limit (x-axis maximum) controlling
#'   tree width.
#'
#' @return A \code{ggplot} object visualizing bootstrap values on the tree.
#'
#' @seealso \code{\link[ape]{read.tree}}, \code{\link[ggtree]{ggtree}}
fun.plot_bootstrap_tree <- function(treefile, width = 40){
  ls.tree <- read.tree(treefile) 
  df.tree <- ls.tree%>% fortify() %>%
    mutate(bootstrap = ifelse(!str_detect(label, "[a-z]"), label, NA)) %>%
    filter(!is.na(bootstrap))
  
  ls.tree$tip.label <- gsub(".*_", "",ls.tree$tip.label)  
  
  ggtree(ls.tree, branch.length = "none")%<+% df.tree +
    geom_label2(aes(subset = !is.na(bootstrap), label = as.numeric(bootstrap), fill =as.numeric(bootstrap)), size = 2.9)+#, size = .5)+
    scale_fill_viridis_c( limits=  c(0,100), name = "Bootstrap value") +
    geom_tiplab( angle  =0)+
    theme(
      legend.position = "none",
      legend.position.inside = c(0.1, 0.8)
    )  +
    xlim(0, width) 
}