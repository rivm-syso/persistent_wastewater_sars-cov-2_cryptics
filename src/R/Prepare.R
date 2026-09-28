# Auke Haver
# 2026-09-14
# NRS / Z&O / I&V / RIVM

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Setup ------------------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

source("src/R/functions.R")
source("src/R/packages.R")

# Load parameters
ls.parameters <- read_yaml("data/parameters.yaml")
df.metadata   <- read_tsv("data/input/metadata.tsv")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Download problematic sites ---------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

fp.problematic_sites <- "data/reference/problematic_sites_sarsCov2.vcf"
if(!file.exists(fp.problematic_sites)){
    download.file(
        url = "https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/refs/heads/master/problematic_sites_sarsCov2.vcf",
        destfile = fp.problematic_sites
    )
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Viral load data --------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if(!file.exists("data/processed/viral_load_sites.tsv.gz")){
    df.vl_site <- 
        read_delim(
            delim = ";",
            file = "https://data.rivm.nl/covid-19/COVID-19_rioolwaterdata.csv"
        )%>%
        filter(
            between(
                as.Date(Date_measurement), 
                as.Date(ls.parameters$study_period[[1]]), 
                as.Date(ls.parameters$study_period[[2]])
            )
        )
    
    write_tsv(
        df.vl_site,
        file = "data/processed/viral_load_sites.tsv.gz"
    )
} else {
    df.vl_site <- 
        read_tsv(
            file = "data/processed/viral_load_sites.tsv.gz",
            col_types = "nDDicn"
        )
}

if(!file.exists("data/processed/viral_load_country.tsv.gz")){
    df.vl_country <- 
        read_delim(
            delim = ";",
            file = "https://data.rivm.nl/covid-19/COVID-19_rioolwaterdata_landelijk.csv"
        ) %>%
        filter(
            between(
                as.Date(Date_measurement), 
                as.Date(ls.parameters$study_period[[1]]), 
                as.Date(ls.parameters$study_period[[2]])
            )
        )
    
    
    write_tsv(
        df.vl_country,
        file = "data/processed/viral_load_country.tsv.gz"
    )
} else {
    df.vl_country <- 
        read_tsv(
            file = "data/processed/viral_load_country.tsv.gz",
            col_types = "nDDn"
        )  %>%
        filter(
            between(
                as.Date(Date_measurement), 
                as.Date(ls.parameters$study_period[[1]]), 
                as.Date(ls.parameters$study_period[[2]])
            )
        )
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Select high-abundance, high-coverage samples ---------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


df.samples_selection <- 
    df.metadata%>%
    filter(location %in% c("Woerden", "Katwoude")) %>% left_join(
      df.vl_site,
      by = c("sample_date" = "Date_measurement", "location" = "RWZI_AWZI_name")
    ) %>%
    filter(
        sapply(
            X = sample_hash,
            FUN = function(x){
                sum(read_tsv(
                    file = file.path(ls.parameters$pipeline_dir, "mpileup_depths/", paste0( x, ".tsv")),
                    col_names = c("index", "pos", "ref", "depth"),
                    col_types = "cici"
                )$depth >= 10)/29903 >= 0.80
            }
        )
    ) %>%
    mutate(
        freyja_data = lapply(
            X = sample_hash,
            FUN = function(x){
                read_freyja_tsv(
                    filename = file.path(ls.parameters$pipeline_dir, "freyja_demix/sars-cov-2/", paste0(x, "_lineages.tsv"))
                )
            }
        )
    ) %>%
    unnest(freyja_data) %>%group_by(lineage) %>%
    group_nest() %>%
    mutate(bin = sapply(lineage, function(x){
        for(lineage in c("B.11", "B.1.221")){
            if(x == lineage | str_detect(x, paste0(lineage,".[0-9]"))){
                return(lineage)
            }
        }
        return("other")
        
        
    })) %>% unnest(data) %>%
    group_by(location, sample_date, sample_hash, bin) %>%
    reframe(abundance = sum(abundance)) %>%
    filter(bin %in% c("B.11", "B.1.221")) %>%
    filter(abundance >= 0.85) %>%
    mutate(
        sequence = sapply(
            X = sample_hash,
            FUN = function(x){
                paste0(read_lines(paste0(file.path(ls.parameters$pipeline_dir,"min_depth_10",paste0(x,".fasta"))), skip = 1))
            }
        )
    )

df.samples_selection$sample_hash %>%
  writeLines(
    "data/processed/selected_samples.txt"
  )

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Load Lineage table -----------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

df.barcodes  <- read_feather(file = file.path(ls.parameters$pipeline_dir, "usher_barcodes.feather")) %>% 
  pivot_longer(cols = !index, names_to = "muts", values_to = "value")

df.lineage_table <- distinct(df.barcodes, index) %>%
  reframe(
    lineage = index,
    bin = lapply(index,  function(x) bin_lineage(x, names(ls.parameters$lineage_bins)[-1], include_recombinants = TRUE))
  ) %>%
  unnest(bin)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Load Freyja data -------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Load original Freyja data
df.freyja_original <- 
  df.metadata  %>%
  select(sample_hash, location, sample_date) %>%
  mutate(freyja_data = lapply(
    sample_hash, 
    function(x) read_freyja_tsv(file.path(ls.parameters$pipeline_dir,"freyja_demix/sars-cov-2", paste0(x, "_lineages.tsv"))))
    ) %>%
  unnest(freyja_data) %>%
  left_join(
    df.lineage_table,
    by = "lineage",
    relationship = "many-to-many"
    ) %>%
    group_by(sample_hash, location, sample_date, lineage, bin, resid, coverage) %>%
    reframe(abundance = sum(abundance * frac)) 

write_tsv(
  x = df.freyja_original,
  file = "data/processed/freyja_data.tsv.gz"
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Prepare Snakemake analsyis ---------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

dir.create("data/snakemake_input", showWarnings = FALSE)
dir.create("data/snakemake_output", showWarnings = FALSE)

# initialize configs
config <- list(
    inputdir = normalizePath("data/snakemake_input/"),
    outputdir = normalizePath("data/snakemake_output/"),
    outputfiles = c(),
    containers = list(
        augur  = normalizePath("environment/augur.sif"),
        freyja = normalizePath("environment/freyja.sif"),
        renv   = normalizePath("environment/r.sif")
    ),
    original_barcodes    = file.path(normalizePath(ls.parameters$pipeline_dir),"usher_barcodes.feather"),
    sample_depths_dir    = file.path(normalizePath(ls.parameters$pipeline_dir),"mpileup_depths/"),
    sample_variants_dir  = file.path(normalizePath(ls.parameters$pipeline_dir),"snv/"),
    reference_genome     = normalizePath("data/reference/reference.fasta"),
    problematic_sites    = normalizePath("data/reference/problematic_sites_sarsCov2.vcf"),
    barcode_gen_script   = normalizePath("src/snakemake/scripts/generate_barcodes.R"),
    barcode_merge_script = normalizePath("src/snakemake/scripts/merge_barcodes.R"),
    freyja_merge_script  = normalizePath("src/snakemake/scripts/merge_freyja.R"),
    analyses = list(
        freyja_samples = list(
            default    = list(
                B.11    = df.metadata$sample_hash[df.metadata$location == "Woerden"],
                B.1.221 = df.metadata$sample_hash[df.metadata$location == "Katwoude"],
                merged  = df.metadata$sample_hash
            )
        )
    )
)


# Loop -------------------------------------------------------------------------

for(file in c("barcodes_full.tsv", "barcodes_known.tsv")){
  config$outputfiles = c(
    config$outputfiles, 
    file.path(config$outputdir, "merged_default", file)
  )
}

for(lineage in c("B.11", "B.1.221")){
    
  for(analysis in c("default", "gisaid")){
    dir.create(
      recursive = TRUE,
      showWarnings = FALSE,
      path = paste0(config$inputdir,"/", lineage, "_", analysis)
    )
  }
    samples = df.samples_selection$sample_hash[df.samples_selection$bin== lineage]
    
    df.default <- df.samples_selection %>% filter(sample_hash %in% samples)
    
    write_analysis_inputs(df.default, config$inputdir, lineage, "default")
    
    
    for(file in c("barcodes_full.tsv", "barcodes_known.tsv")){
        config$outputfiles = c(
            config$outputfiles, 
            file.path(config$outputdir, paste0(lineage, "_default"), file)
        )
    }
    
    for(file in c(
      "tree.nwk",
      "refined_tree.nwk",
      "refined_node_data.json",
      "ancestral_data.json"
    )){
      config$outputfiles = c(
        config$outputfiles, 
        file.path(config$outputdir, paste0(lineage, "_gisaid"), file)
      )
    }
    
    for(iteration in seq(1,50)){
        
        ## VALIDATION ANALYSIS
        analysis = paste0("validate", iteration)
        set.seed(iteration)
        
        test_samples = sample(samples, size = ceiling(length(samples)*0.2))
        train_samples= samples[!samples %in% test_samples]
        
        config[["analyses"]][["freyja_samples"]][[analysis]][[lineage]] = test_samples
        
        df.train <- df.samples_selection %>% filter(sample_hash %in% train_samples)
        write_analysis_inputs(df.train, config$inputdir, lineage, analysis)
        
        for(file in c("barcodes_full.tsv", "barcodes_known.tsv")){
            config$outputfiles = c(
                config$outputfiles, 
                file.path(config$outputdir, paste0(lineage, "_", analysis), file)
            )
        }
        
        ## RANDOMIZE DATE ANALYSIS
        set.seed(iteration)
        df.random <- 
            filter(df.samples_selection, bin == lineage) %>%
            mutate(sample_date = sample(sample_date, size = nrow(.), replace = FALSE))
        
        analysis = paste0("randomdates", iteration)
        
        write_analysis_inputs(
            df           = df.random,
            inputdir     = config$inputdir,
            lineage      = lineage,
            analysis_name = analysis
        )
        
        for(filetype in c(
            "refined_tree.nwk",
            "refined_node_data.json",
            "ancestral_data.json",
            "barcodes_full.feather",
            "barcodes_known.feather"
        )){
            config$outputfiles = c(
                config$outputfiles, 
                file.path(config$outputdir, paste0(lineage, "_", analysis), filetype)
            )
        }
        
    }
}
write_yaml(
    config,
    paste0(config$inputdir, "/config.yaml")
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Load genome file -------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if(!file.exists("data/processed/genome.gff")){
  read_tsv(
    file = "https://github.com/nextstrain/nextclade_data/raw/refs/heads/master/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs/genome_annotation.gff3",
    skip = 2,
    col_names = c("1", "2", "type", "start", "stop", "3", "strand", "4", "name"),
    col_types = "ccciicccc"
  ) %>%
    reframe(
      gene = gsub("gene_name=", "", name),
      start,
      stop
    ) %>%
    arrange(start) %>%
    write_tsv("data/processed/genome.gff")
  
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# GISAID -----------------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sapply(c("B.11", "B.1.221"), function(x){
  write_analysis_inputs(
    df =  rbind(
        filter(df.samples_selection, bin == x) %>% select(location, sample_date, sequence),
        reframe(read_tsv(file = file.path("data/gisaid/", paste0(x, ".tsv.gz")), col_types = "ccDcc"), location = accession_id, sample_date, sequence =seq)
      ),
    inputdir = config$inputdir,
    lineage = x,
    analysis_name = "gisaid"
  )
})

