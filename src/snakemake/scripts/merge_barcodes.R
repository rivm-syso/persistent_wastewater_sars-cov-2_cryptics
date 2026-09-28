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
    c("--B.11_full"),
    help="B.11_full",
    default = "data/snakemake_output/B.11_default/barcodes_full.feather"
  ),
  make_option(
    c("--B.11_known"),
    help="B.11_known",
    default = "data/snakemake_output/B.11_default/barcodes_known.feather"
  ),
  make_option(
    c("--B.1.221_full"),
    help="B.1.221_full",
    default = "data/snakemake_output/B.1.221_default/barcodes_full.feather"
  ),
  make_option(
    c("--B.1.221_known"),
    help="B.1.221_known",
    default = "data/snakemake_output/B.1.221_default/barcodes_known.feather"
  ),
  make_option(
    c("--merged_full"),
    help="merged_full",
    default = "data/snakemake_output/merged_default/barcodes_full.feather"
  ),
  make_option(
    c("--merged_known"),
    help="merged_known",
    default = "data/snakemake_output/merged_default/barcodes_known.feather"
  )
)



opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

args <- commandArgs(trailingOnly = TRUE)


read_barcodes_file <- function(filepath){
  read_feather(filepath) %>%
    pivot_longer(
      cols = !index,
      names_to = "muts",
      values_to = "vals"
    ) %>%
    filter(vals != 0)
}

write_barcodes_file <- function(long_df, filepath){
  pivot_wider(
    long_df,
    names_from = "muts",
    values_from = "vals",
    values_fill = 0
  ) %>%
    write_feather(filepath)
}

write_barcodes_file(
  long_df = rbind(
    mutate(read_barcodes_file(opt$B.11_full),    index = gsub("CRYPT", "B.11.x", index)),
    mutate(read_barcodes_file(opt$B.1.221_full), index = gsub("CRYPT", "B.1.221.x", index))
  ) %>%
    distinct(),
  filepath = opt$merged_full
)

gc()

write_barcodes_file(
  long_df = rbind(
    mutate(read_barcodes_file(opt$B.11_known),    index = gsub("CRYPT", "B.11.x", index)),
    mutate(read_barcodes_file(opt$B.1.221_known), index = gsub("CRYPT", "B.1.221.x", index))
  ) %>%
    distinct(),
  filepath = opt$merged_known
)


