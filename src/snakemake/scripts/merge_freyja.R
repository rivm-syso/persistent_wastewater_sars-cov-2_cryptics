# Auke Haver
# 2026-09-14
# NRS / Z&O / I&V / RIVM

# Load packages -----------------------------------------------------------

library(stringr)
library(dplyr)
library(tidyr)
library(readr)

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

# Load data ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

tibble(
    data  = lapply(
        args[-1],
        read_freyja_tsv
    )
) %>%
    unnest(data) %>% 
    mutate(
        sample_hash = gsub(".*/|_lineages.tsv", "", filename)
    ) %>%
    distinct() %>%
    write_tsv(
        file = args[1]
    )


