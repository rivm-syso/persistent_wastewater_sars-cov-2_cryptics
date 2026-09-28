# Auke Haver
# 2026-09-14
# NRS / Z&O / I&V / RIVM

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 0. Setup environment ---------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Load functions and packages
source("src/R/functions.R")
source("src/R/packages.R")

# Set options
options(ignore.negative.edge=TRUE)

# Load parameters
ls.parameters <- read_yaml("data/parameters.yaml")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 1. PLOTTING FIGURE 1: Viral load ---------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Load data ---------------------------------------------------------------

# Metadata
df.metadata   <- read_tsv("data/input/metadata.tsv")
df.vl_country <- read_tsv("data/processed/viral_load_country.tsv.gz", col_types = "cDDn")
df.vl_sites   <- read_tsv("data/processed/viral_load_sites.tsv.gz", col_types = "nDDccn")

# Join instances of viral loads observed on specific dates with the National
# average for that day; calculate the ratio site / average
df.ratios <-
  inner_join(
    x = reframe(
      df.vl_sites, 
      RWZI_AWZI_name, 
      Date_measurement, 
      site_RNA_flow_per_100000 = RNA_flow_per_100000
    ),
    y = df.vl_country,
    by = join_by(Date_measurement)
  ) %>%
  mutate(ratio = site_RNA_flow_per_100000 / RNA_flow_per_100000)

# Generate a continuous 7-day average confidence interval of viral load ratioes
df.confint_log10_vl_ratio <-
  tibble(
    Date_measurement = seq(
      from = as.Date(ls.parameters$study_period[[1]]), 
      to   = as.Date(ls.parameters$study_period[[2]]), 
      by   = "1 day"
    ),
    # Only the upper 0.975 quantile is relevant, because the ratio can be -inf
    maxratio = sapply(Date_measurement, function(x){
      filter(
        df.ratios,
        between(df.ratios$Date_measurement, x-3, x+3)
      ) %>%
        pull("ratio") %>%
        log10() %>% 
        quantile(.,.95) %>%
        return()
    })
  )

# Generate a dataframe where the viral load is divided into fractions corresponding
# to the 7-day average lineage bin abundance
df.lineage_fraction_vl <- (function(freyja_data){
  df.freyja_binned <- freyja_data %>%
    filter(location %in% ls.parameters$reference_sites) %>%
    group_by(sample_hash, location, sample_date, bin) %>%
    reframe(abundance = sum(abundance)) %>%
    complete(nesting(sample_hash, location, sample_date), bin, fill = list(abundance = 0))
  
  df.confint_log10_vl_ratio %>%
    mutate(
      lineage_data = lapply(
        X = Date_measurement,
        FUN = function(x){
          df.freyja_binned %>%
            filter(between(sample_date, x-3, x+3)) %>%
            group_by(bin) %>% reframe(perc = sum(abundance) / n_distinct(sample_hash))
        })
    ) %>%
    unnest(lineage_data, keep_empty = TRUE)%>%
    inner_join(df.vl_country,by = join_by(Date_measurement)) %>%
    mutate(perc = replace_na(perc, 1), bin = replace_na(bin, "Unknown")) %>%
    complete(nesting(Date_measurement,RNA_flow_per_100000), bin, fill = list(perc =0))%>%
    mutate(vl = RNA_flow_per_100000*perc) %>%
    return()
})(
  read_tsv("data/processed/freyja_data.tsv.gz", col_types = "ccDccnnn")
)


## 1.2 BASE FIGURE LAYOUT ------------------------------------------------------
# This is the base layout for figure 1, which is shared in subfigure A-C
ls.fig1_layout <- 
  list(
    theme_bw(),
    scale_x_date(
      name = "Date",
      breaks = seq(
        as.Date("2020-09-01"), 
        ls.parameters$study_period[[2]], 
        by = "3 months"
      ),
      date_labels = "%b-%y",
      limits = c(
        as.Date(ls.parameters$study_period[[1]]),
        as.Date(ls.parameters$study_period[[2]])
      ),
      expand = c(0,0),
      minor_breaks = seq(as.Date("2020-10-01"), as.Date("2025-12-01"), by = "1 months")
    ),
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank()
    )
  )


## PLOTTING THE FIGURE ---------------------------------------------------------
# Here we initialize the plot
ggsave(
  filename  = "output/Figure_1.pdf",
  plot      = plot_grid(
    plotlist = list(
      ### SUBFIGURE A ----------------------------------------------------------
      ggplot() + 
        geom_area(
          data = df.lineage_fraction_vl,
          mapping = aes(
            x = Date_measurement,
            y = vl/1e14,
            fill = bin
          ), 
          color = "black", alpha =.6
        ) +
        scale_fill_manual(values = ls.parameters$lineage_bins)+
        scale_y_continuous(
          limits = c(0, 4.5),
          breaks = seq(0,6, 1),
          expand = c(0,0),
          name = "National Average\n1e14 GC per 1e5 inhabitants"
        ) +
        theme(legend.position = "none") + 
        ls.fig1_layout,
      ### SUBFIGURE B ----------------------------------------------------------
      fun.plot_vl_ratio(
        df.log10_ratio_confint = df.confint_log10_vl_ratio,
        df.wwtp_ratioes = df.ratios,
        df.lineage_estimates = read_tsv("data/processed/freyja_data.tsv.gz", col_types = "ccDccnnn"),
        target_location = "Woerden",
        target_lineage = "B.11",
        layout_commands = ls.fig1_layout
      ),
      ### SUBFIGURE C ----------------------------------------------------------
      fun.plot_vl_ratio(
        df.log10_ratio_confint = df.confint_log10_vl_ratio,
        df.wwtp_ratioes = df.ratios,
        df.lineage_estimates = read_tsv("data/processed/freyja_data.tsv.gz", col_types = "ccDccnnn"),
        target_location = "Katwoude",
        target_lineage = "B.1.221",
        layout_commands = ls.fig1_layout
      ),
      fun.abundance_legend()
      
    ),
    nrow = 4, ncol = 1, rel_heights = c(5,5,5,1)
  ),
  width = 10.5,
  height = 9.45
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE 2 ------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ls.known_sites <- colnames(read_feather(file = file.path(ls.parameters$pipeline_dir, "usher_barcodes.feather"))) %>%
  .[. != "index"] %>%
  parse_number() %>%
  unique() %>% sort()

plt.timescaled_trees <- 
  plot_grid(
    fun.make_cryptic_tree("B.11",    daterange = seq(as.Date("2021-03-01"), as.Date("2023-03-01"), by = "2 months")),
    fun.make_cryptic_tree("B.1.221", daterange = seq(as.Date("2022-03-01"), as.Date("2026-03-01"), by = "4 months")),
    fun.make_cryptic_alluvial("B.11"),
    fun.make_cryptic_alluvial("B.1.221"),
    ncol = 2,
    align = "hv",
    axis = "trbl"
  )

ggsave(
  filename = "output/Figure_2.pdf",
  plot = plt.timescaled_trees,
  width = 13,
  height = 10
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE 3 ------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ggsave(
  filename = "output/Figure_3.pdf",
  plot     = plot_grid(
    fun.make_cryptic_heatmap(target_lineage = "B.11",    target_location = "Woerden"),
    fun.make_cryptic_heatmap(target_lineage = "B.1.221", target_location = "Katwoude"),
    fun.abundance_legend(),
    align = "v",
    ncol = 1,
    rel_heights = c(3,3,1),
    axis = "tblr"
  ),
  height = 4.1,
  width = 12
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE S1: Log10 ratio - abundance ----------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ggsave(
  filename = "output/Figure_S1.pdf", 
  width = 7, height = 4, 
  plot = inner_join(
    read_tsv("data/processed/freyja_data.tsv.gz", col_types = "ccDccnnn"),
    filter(df.ratios, RWZI_AWZI_name %in% c("Woerden", "Katwoude")),
    by = c("location"  = "RWZI_AWZI_name", "sample_date" = "Date_measurement")
  ) %>%
    complete(nesting(location, sample_date, ratio), bin, fill = list(abundance = 0)) %>%
    filter((location == "Katwoude" & bin == "B.1.221") | (location == "Woerden" & bin == "B.11")) %>%
    group_by(location, sample_date, ratio, bin) %>%
    reframe(abundance = sum(abundance)) %>% 
    (function(x){
      ggplot()  +
        geom_point(
          data = x,
          mapping = aes(
            x = log10(ratio),
            y = abundance,
            fill = location
          ),
          pch = 21) +
        geom_line(
          data = mutate(x, est = predict(glm(abundance~log10(ratio)*location, "binomial",x), type = "response")),
          mapping = aes(
            x = log10(ratio),
            y = est,
            color = location
          )
        )+
        theme_bw() +
        labs(x = "Log10 GC Ratio\nWWTP/National Average", y = "Cryptic Lineage Abundance") +
        scale_y_continuous(
          limits = c(0.00,1), breaks = c(0.001, seq(.2,1,.2)),
          labels = c("<=0.001*", "0.200", "0.400", "0.600", "0.800","1.000")
        ) +
        lims(x = c(-1, 2.5)) +
        scale_fill_manual(
          values = c("Katwoude" = "#5ec962", "Woerden"  = "#3b528b"),
          labels = c("Katwoude" = "WWTP-KW", "Woerden"  = "WWTP-WR"),
          name   = NULL
        ) +
        scale_color_manual(
          values = c("Katwoude" = "#5ec962", "Woerden"  = "#3b528b"),
          labels = c("Katwoude" = "WWTP-KW", "Woerden"  = "WWTP-WR"),
          name = NULL
        )
    })(.)
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE S2: WWTP_UT B.11 ---------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


ggsave(
  filename = "output/Figure_S2.pdf",
  width = 7.5,
  height = 6,
  plot = filter(
    df.ratios,
    between(Date_measurement, as.Date("2023-01-01"), as.Date("2023-02-28")),
    RWZI_AWZI_name %in% c("Woerden", "Utrecht")
  ) %>%
    left_join(
      distinct(df.metadata, sample_hash, location, sample_date),
      by = c( "RWZI_AWZI_name" = "location", "Date_measurement" = "sample_date")
    ) %>%
    left_join(
      read_tsv(
        file = file.path(ls.parameters$output_dir,"merged_default/barcodes_full.tsv"),
        col_types = "ccnnnc"
      ) %>%
        group_by(lineage) %>%
        group_nest() %>%
        mutate(bin = sapply(lineage, function(x){
          if(x == "B.11" | str_detect(x, "x")){return(x)} else{return("other")}
        })) %>%
        unnest(data) %>%
        group_by(sample_hash, bin) %>%
        reframe(abundance = sum(abundance)),
      by = "sample_hash"
    ) %>%
    mutate(
      bin = factor(gsub(".*x.", "", bin), levels = c( "other", "B.11", NA,seq(1,8))),
      abundance = replace_na(abundance,1), 
      logratio = log10(ratio),
      # Only the absolute height of ratio is relevant and the composition of lineages within that bar. 
      binratio = ifelse(logratio >0, abs(abundance*logratio), -abs(abundance*logratio))
    )%>%
    (function(x){
      ggplot() +
        geom_col(
          data = x,
          mapping = aes(x = Date_measurement,y = binratio, fill = bin
          ),
          color = "black",
          width = 1
        )+
        facet_grid(rows = vars(RWZI_AWZI_name))+
        geom_hline(
          data = tibble(location = c("Woerden", "Utrecht")),
          mapping = aes(yintercept = 0),
          linetype = "dashed"
        ) +
        scale_x_date(
          date_minor_breaks = "1 day",
          breaks = seq(as.Date("2023-01-01"), as.Date("2023-02-28"), by = "7 days"),
          limits = c(as.Date("2023-01-01"), as.Date("2023-02-28")),
          expand = c(0,0)
        ) +
        scale_y_continuous(
          expand = c(0,0),
          limits = c(-.6, 2.1),
          breaks = seq(-.5, 2, .5)
        ) +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          panel.spacing = unit(1, "cm"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_blank()
        ) +
        scale_fill_viridis_d(na.value = "grey",  name = NULL)
    })(.)
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE S3: Node support for original trees  -------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ggsave(
  filename = "output/Figure_S3.pdf",
  plot = plot_grid(
    fun.plot_bootstrap_tree(
      treefile = file.path(ls.parameters$output_dir, "B.11_default",    "tree.nwk"), 
      width    = 21
    ),
    fun.plot_bootstrap_tree(
      treefile = file.path(ls.parameters$output_dir, "B.1.221_default", "tree.nwk"), 
      width    = 42
    ),
    rel_widths = c(2,3)
  ),
  width = 10, 
  height = 15
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE S4: Substitution rates validation ----------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

df.rates <- 
  tibble(
    file = list.files(ls.parameters$output_dir, recursive = TRUE, pattern = "refined_node_data.json", full.names = TRUE),
    data = lapply(file, function(x)jsonlite::read_json(x))
  ) %>% 
  filter(!str_detect(file, "gisaid")) %>%
  reframe(
    lineage = gsub(".*snakemake_output/+|_.*", "", file),
    analysis = gsub(".*B.11_|.*B.1.221_|/.*", "", file),    
    intercept = sapply(data, function(x) x$clock$intercept),
    rate      = sapply(data, function(x) x$clock$rate),
    rate_std  = sapply(data, function(x) x$clock$rate_std),
    rtt_Tmrca = sapply(data, function(x) x$clock$rtt_Tmrca),
    type = sapply(file, function(x){
      for(i in c("default", "validate", "random")){
        if(str_detect(x, i)){
          return(i)
        }
      }
    })
  )

df.rates %>%
  group_by(lineage, type) %>%
  reframe(count = n())


ggsave(
  plot = ggplot(
    data = mutate(
      filter(df.rates, type != "default"),
      iteration = sapply(analysis, function(x){as.numeric(gsub("[a-z]+", "",x))}))
  ) +
    geom_rect(
      data = filter(df.rates, type == "default"),
      mapping = aes(
        xmin = -0.5, xmax = 50.5,
        ymin = rate- rate_std,
        ymax = rate+rate_std,
        fill = "Analysis"
      ),
      alpha = .2
    ) +
    geom_hline(mapping = aes(yintercept = 1.8*10^-4, color = "NextStrain rate"), linetype = "dashed")+
    geom_point(
      mapping = aes(
        x = iteration, 
        y = rate, 
        color = ifelse(str_detect(analysis, "random"), "Randomized dates", "Resampling")
      )
    )  +
    geom_errorbar(
      mapping = aes(
        x = iteration, 
        color = ifelse(str_detect(analysis, "random"), "Randomized dates", "Resampling"), 
        ymin = rate- rate_std, 
        ymax = rate+ rate_std
      )
    )+ 
    scale_color_manual(
      name = NULL,
      values = list(
        "Randomized dates" = "#21918c", 
        "NextStrain rate" = "black",
        "Resampling" = "#3b528b"
      )
    )+
    scale_fill_manual(
      name = NULL,
      labels = list("Analysis" = "Full dataset estimate" ),
      values = list("Analysis" = "#440154")
    )+
    facet_grid(cols = vars(lineage))+ 
    theme_bw() +
    theme(
      strip.background = element_blank()
    )+
    scale_y_continuous(name = "Substitutions/site/year") +
    scale_x_continuous(limits = c(-0.5, 50.5), expand = c(0,0), name = "Bootstrap iteration") , 
  filename = "output/Figure_S4.pdf", 
  width = 10, 
  height = 5
)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE S5: Reference datasets, sub/ins/del frequencies --------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Load Depth and SNV data
df.selection_depth_snv <- 
  filter(df.metadata, location %in% c(ls.parameters$reference_sites, "Katwoude", "Woerden")) %>%
  select(sample_hash, location) %>%
  mutate(
    snv_data = lapply(
      X =  file.path(ls.parameters$pipeline_dir, "snv",paste0(sample_hash,"_ivar.tsv")),
      FUN = fun.read_ivar_file
    ),
    depth_data = lapply(
      X = file.path(ls.parameters$pipeline_dir, "mpileup_depths",paste0( sample_hash,".tsv")),
      FUN = fun.read_depth_file
    )
  )

# Load genome ORF file
df.genome <- read_tsv(
  file = "data/processed/genome.gff",
  col_types = "cii"
)

df.indel_selection <- 
  df.selection_depth_snv %>% 
  filter(sample_hash %in% readLines("data/processed/selected_samples.txt")) %>% # The 93 selected samples
  select(-depth_data) %>%
  unnest(snv_data) %>%  select( location, POS, ALT, ALT_FREQ, ALT_DP, sample_hash) %>%
  group_by(POS, ALT, location) %>%
  filter(str_detect(ALT, "-|\\+")) %>%
  filter(ALT_DP >= 10, ALT_FREQ >= .5) %>%
  filter(n()>2) %>%
  reframe(
    residue = mapply(
      x = POS,
      y = ALT,
      SIMPLIFY = TRUE,
      FUN = function(x,y){
        genes = df.genome$gene[df.genome$start <= x & df.genome$stop >= x]
        
        if(length(genes)==0){
          # Intergenic
          if(str_detect(y, "-")){
            # Intergenic deletion
            n_bases = sum(strsplit(y, "")[[1]] == "N")
            
            return(paste0(x, "-", x+n_bases -1, "del"))
          } else{
            # Intergenic insertion
            
            return(paste0(x), "ins", gsub("\\+", "", y))
          }
          
        } else {
          # Not intergenic
          output = c()
          for(gene in genes){
            residue =  floor((x -df.genome$start[df.genome$gene == gene])/3)+1
            output = c(output, paste0(gene, "_", residue))
            
          }
          return(paste0(output, collapse = ";"))
        }
      }
    )
  )

df.mutation_freqs_regular_wwtp <- 
  df.selection_depth_snv %>%
  mutate(
    data = lapply(sample_hash, function(x){
      # Add S:D215G and S:L828F substitutions and possible forms of the 3' UTR stemp-loop deletion
      tibble(
        POS = c(22206, 24044, 29735, 29737, 29733),
        ALT = c("G",     "T",     "-NNNNNNNNNNNNNNNNNNNNNNNN", "-NNNNNNNNNNNNNNNNNNNNNNN", "-NNNNNNNNNNNNNNNNNNNNNNNNNN")
      ) %>%
        rbind(
          distinct(df.indel_selection, POS, ALT)
        ) %>%
        distinct()
    })
  )%>%
  mutate(data = mapply(x = data, y = snv_data, z = depth_data, SIMPLIFY =FALSE, FUN=  function(x,y,z){
    left_join(x,z, by = c("POS")) %>%
      left_join(select(y,-TOTAL_DP),   by = c("POS", "ALT", "REF")) %>%
      select(REF,POS, ALT, ALT_FREQ, ALT_DP, TOTAL_DP)
  })) %>%
  select(-snv_data, -depth_data) %>%
  unnest(data)%>%
  mutate(ALT_FREQ = ifelse(TOTAL_DP < 10, NA, replace_na(ALT_FREQ,0)),sub = paste0(REF, POS, ALT)) %>%
  distinct()%>%
  left_join(distinct(df.metadata, location, sample_date, sample_hash)) %>%
  filter(location %in% ls.parameters$reference_sites)

plt.mutation_frequencies_req <-
  df.mutation_freqs_regular_wwtp %>%
  filter(!is.na(ALT_FREQ)) %>%
  mutate(
    sub = sapply(sub, function(y){
      if(str_detect(y, "-N+$")){
        gsub("\\-N+$", paste0("-", sum(strsplit(y, "")[[1]] == "N")), y)
      } else {
        return(y)
      }
    }
    )
  ) %>%
  arrange(POS) %>%
  mutate(sub = factor(sub, levels = unique(sub))) %>%
  ggplot(aes(x = sample_date, y = ALT_FREQ)) +
  geom_bin2d(color = "black")+
  theme_bw()+
  scale_x_date(
    breaks = seq(as.Date("2021-01-01"), as.Date("2026-01-01"), by = "1 year"),
    date_labels = "%Y",
    limits=  c(as.Date("2021-01-01"), as.Date("2026-01-01")),
    expand = c(0,0)
  )  +
  theme(
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(.9,0.05),
    legend.key.height=unit(.375,"cm"),
    legend.key.width =unit(.20,"cm")
  ) +
  scale_fill_viridis_c(limits = c(1,120), trans = "log2")+
  facet_wrap(~sub, ncol = 6)

ggsave(
  "output/Figure_S5.pdf",
  plot = plt.mutation_frequencies_req,
  width = 6,
  height= 10
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE S6: Model residuals ------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

df.residuals <- full_join(
  # Original residuals
  read_tsv("data/processed/freyja_data.tsv.gz", col_types = "ccDccnnn")  %>%
    mutate(bin = sapply(lineage, function(x){
      if(x == "B.11" | x == "B.1.221" | str_detect(x, fixed("B.1.221."))){
        return("crypt")
      } else {
        return("other")
      }
    })) %>%
    complete(nesting(sample_hash, resid), bin , fill = list(abundance = 0)) %>%
    filter(bin == "crypt") %>% group_by(sample_hash, resid) %>% 
    reframe(crypt_frac = sum(abundance)),
  # Updated residuals
  tibble(
    type = c("known", "full"),
    data = lapply(type, function(y){
      read_tsv(file = file.path(ls.parameters$output_dir, "merged_default", paste0("barcodes_", y,".tsv")), col_types = "ccnnnc") %>%
        mutate(bin = sapply(lineage, function(x){
          if(x == "B.11" | x == "B.1.221" | str_detect(x, fixed("B.1.221."))| str_detect(x, "x") | str_detect(x, "CRYPT") ){
            return("crypt")
          } else {
            return("other")
          }
        })) %>%
        complete(nesting(sample_hash, resid), bin , fill = list(abundance = 0)) %>%
        filter(bin == "crypt") %>% group_by(sample_hash) %>% 
        reframe(new_resid = resid, new_crypt_frac = sum(abundance))
    })
  ) %>% unnest(data),
  by = "sample_hash"
) 

# Load residuals obtained using new barcodes
df.residuals_validate  <- 
  tibble(
    files = list.files(
      path       = ls.parameters$output_dir, 
      pattern    = "barcodes_[a-z]+.tsv", 
      recursive  = TRUE, 
      full.names = TRUE) %>% .[str_detect(., "validate")],
    data  = lapply(
      X = files,
      FUN = function(x)
        read_tsv(
          file = x,
          col_types = "ccnnnc"
        )
    )
  ) %>%
  unnest(data) %>%
  reframe(
    lineage = gsub(".*/", "", gsub("_validate.*", "", files)),
    validate =  gsub(".*_validate|/barcodes.*", "", files),
    barcodes = gsub(".*barcodes_|.tsv", "", files),
    sample_hash,
    resid
  ) %>% arrange(sample_hash) %>%
  distinct()

ls.fig6_layout <- list(
  theme_bw() +
    theme(
      legend.position = "right",
      legend.title.position = "top",
      legend.key.height=unit(2,"cm")
    ),
  scale_fill_viridis_c(limits = c(1,2049), name = "# Samples", trans = "log2", breaks = c(1, 4, 16, 64, 256, 1024))
)

ggsave(
  plot = plot_grid(
    ggpubr::ggarrange(
      ggplot() +
        geom_abline(slope = 1, linetype = "dashed")  +
        geom_bin2d(
          data = filter(df.residuals,type == "full"),
          mapping = aes( x = resid, y = new_resid),
          color = "black"
        ) + 
        scale_x_continuous(limits = c(0,70), expand = c(0,0)) +
        scale_y_continuous(limits = c(0,70), expand = c(0,0)) +
        labs(x = "Orignal residuals", y = "Updated residuals (all-sites)") + 
        ls.fig6_layout,
      ggplot() +
        geom_abline(slope = 1, linetype = "dashed") +
        geom_bin2d(
          data = filter(df.residuals,type == "known"),
          mapping = aes( x = resid, y = new_resid),
          color = "black"
        ) + 
        scale_x_continuous(limits = c(0,70), expand = c(0,0)) +
        scale_y_continuous(limits = c(0,70), expand = c(0,0)) +
        labs(x = "Orignal residuals", y = "Updated residuals (known-sites only)") + 
        ls.fig6_layout,
      ggplot() +
        geom_hline(mapping = aes(yintercept = 100), linetype = "dashed") +
        geom_bin2d(
          data = filter(df.residuals,type == "known"),
          aes(
            x = new_crypt_frac,
            y = new_resid / resid*100
          ),
          color = "black"
        ) + 
        labs(x = "Cryptic lineage abundance ", y = "Residuals") + 
        scale_x_continuous(
          breaks = c(0.001, seq(0.2,1,0.2)),
          labels = c("<=0.001*","0.200","0.400","0.600","0.800","1.000")
        )+
        ls.fig6_layout,
      common.legend = TRUE,
      legend = "right",
      ncol = 3
    ),
    ggplot() +
      geom_hline(mapping = aes(yintercept = 100), linetype = "dashed") +
      geom_boxplot(
        data = left_join(
          df.residuals_validate, 
          reframe(df.residuals, barcodes = type, sample_hash, old_resid = resid), 
          by = c("sample_hash" ,"barcodes"),
          relationship = "many-to-many") %>%
          left_join(select(df.metadata, sample_hash, sample_date, location),by = join_by(sample_hash)),
        mapping = aes(
          x= as.factor(sample_date),
          y = resid / old_resid *100,
          fill = barcodes
        )
      )+
      facet_grid(cols = vars(location), space = "free_x", scales = "free_x") +
      scale_y_continuous(name = "% change in residuals") + 
      theme_bw()+
      theme(
        axis.text.x = element_text(angle = 90, hjust = .5, vjust = .5),
        panel.spacing = unit(1, "cm"),
        strip.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_blank(),
        strip.text = element_blank(),
        legend.position = "right"
      )+
      scale_fill_viridis_d(begin = .2, name = "Barcodes type", direction = -1, labels = list("full" = "All-sites", "known" = "Known-sites only")),
    rel_heights = c(3,2),
    ncol = 1
  ), 
  filename = "output/Figure_S6.pdf", 
  width = 16, 
  height = 9
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE S7 ----------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

df.hq_selection_freq_indels <- 
  df.indel_selection %>%
  group_by(location) %>%
  group_nest() %>%
  left_join(
    select(df.selection_depth_snv, sample_hash, location,  snv_data, depth_data),
    by = c("location")
  ) %>%
  mutate(data = mapply(x = data, y = snv_data, z = depth_data, SIMPLIFY =FALSE, FUN=  function(x,y,z){
    left_join(x,z, by = c("POS")) %>%
      left_join(select(y,-TOTAL_DP),   by = c("POS", "ALT")) %>%
      select(residue,POS, ALT, ALT_FREQ, ALT_DP, TOTAL_DP)
  })) %>%
  select(-snv_data, -depth_data) %>%
  unnest(data)%>%
  mutate(ALT_FREQ = ifelse(TOTAL_DP < 10, NA, replace_na(ALT_FREQ,0))) %>%
  distinct()

plt.indels <- 
  ggpubr::ggarrange(
    plotlist = lapply(
      c("Woerden", "Katwoude"),
      FUN = function(x){
        df.hq_selection_freq_indels %>%
          left_join(select(df.metadata, sample_hash, sample_date)) %>%
          filter(location == x) %>%
          ggplot(
            aes(
              y = reorder(paste0(POS, ALT), POS), 
              x =as.factor(sample_date), 
              fill = ALT_FREQ)
          ) +
          geom_tile(color = "black")+
          theme_bw() +
          scale_x_discrete(expand = c(0,0))+
          scale_y_discrete(expand = c(0,0))+
          scale_fill_viridis_c(
            limits=  c(0,1)
          ) +
          theme(
            axis.text.x = element_text(angle = 90, size = 4),
            axis.text.y = element_text(size = 4),
            legend.direction = "horizontal",
            legend.title.position = "top",
            legend.position = "bottom",legend.key.width=unit(1,"cm")
          )
      }
    ),
    ncol = 1,
    align = "hv",
    heights = c(15, 30),
    common.legend = TRUE
  )
ggsave(
  "output/figure_S7.pdf",
  plot = plt.indels,
  width = 10,
  height=  6
)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# PLOTTING FIGURE S8 -----------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

plt.circular_trees <-
  plot_grid(
    (function(filepath){
      tree <- read.tree(file = filepath) %>% drop.tip("EPI_ISL_425051") 
      p<- ggtree(tree, layout = "circular", aes(color = branch.length), branch.length = "none")
      p +
        geom_tiplab(color = "black")+
        theme(legend.key.height=unit(1,"cm"), legend.position = "inside", legend.position.inside = c(0.9,0.9))+
        scale_color_viridis_c(direction = 1, end = .9, limits = c(0, 66)/29903, breaks = seq(0, 66, 6)/29903, labels = seq(0, 66, 6))
    })("data/snakemake_output/B.11_gisaid/tree.nwk"),
    (function(filepath){
      tree <- read.tree(file = filepath) 
      p<- ggtree(tree, layout = "circular", aes(color = branch.length), branch.length = "none")
      p +
        geom_tiplab(color = "black")+
        theme(legend.key.height=unit(1,"cm"), legend.position = "inside", legend.position.inside = c(0.9,0.9))+
        scale_color_viridis_c(direction = 1, end = .9, limits = c(0, 44)/29903, breaks = seq(0, 44, 4)/29903, labels = seq(0, 44, 4))
    })("data/snakemake_output/B.1.221_gisaid/tree.nwk"),
    ncol = 2
  )

ggsave(plot = plt.circular_trees,"output/Figure_S8.pdf", width = 40, height =20)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Supplementary Table 3 --------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

df.metadata %>%
  reframe(
    `Date of sampling` = sample_date,
    `WWTP name` = `location`,
    `Enrichment protocol` = NA, # Information from other source
    `PCR Settings` = NA,
    `Primerset` = primerset,
  ) %>%
  writexl::write_xlsx("output/table_S3.xlsx")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Zenodo export ----------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

dir.create("output/Zenodo", showWarnings = FALSE)

df.reconstructed_sequences <- 
  tibble(
    lineage = c("B.11", "B.1.221"),
    sequences = lapply(
      lineage,
      function(x){
        ls.ancestral <- read_json(file.path(ls.parameters$output_dir, paste0(x, "_default/ancestral_data.json")))
        
        return(
          tibble(
            node = readLines(file.path(ls.parameters$output_dir, paste0(x, "_default/selected_nodes.txt"))),
            seq   = sapply(node, function(x){ls.ancestral$nodes[[x]]$seq})
          ) %>%
            mutate(cryptic_sublineage = paste0(x, ".x.", row_number()))
        )
      }
    )
  ) %>%
  unnest(sequences)

writeLines(
  text = paste0(">", df.reconstructed_sequences$cryptic_sublineage, "\n", df.reconstructed_sequences$seq),
  con = "output/Zenodo/reconstructed_sequences.fasta"
)

file.copy(
  overwrite = TRUE,
  from = file.path(ls.parameters$output_dir, "merged_default/barcodes_known.feather"),
  to   = "output/Zenodo/barcodes_known_sites.feather"
)
file.copy(
  overwrite = TRUE,
  from = file.path(ls.parameters$output_dir, "merged_default/barcodes_full.feather"),
  to   = "output/Zenodo/barcodes_all_sites.feather"
)
tibble(
  files = list.files(
    path = file.path(ls.parameters$output_dir), pattern = "refined_node_data.json", recursive = TRUE, full.names = TRUE
  ) %>% .[str_detect(., "default")],
  data = lapply(
    files, function(x){
      
    }
  )
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Statistics -------------------------------------------------------------------
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

df.muts <- 
  tibble(lineage = c("B.11", "B.1.221")) %>%
  mutate(data = lapply(lineage, function(x) read_tsv(file.path("data/snakemake_output/", paste0(x, "_default"), "cryptic_mutations.tsv")))) %>%
  unnest(data)

# # mutations per gene
df.muts %>% distinct(lineage, muts) %>% 
  mutate(gene = sapply(muts, function(x){
    return(paste0(df.genome$gene[df.genome$start <= parse_number(x) & df.genome$stop >= parse_number(x)], collapse = "_"))
  })) %>% group_by(gene) %>%
  reframe(count = n()) %>%
  mutate(perc = count / sum(count) * 100) %>% arrange(-count)%>% View()

# Observed in all nodes
df.muts %>% group_by(muts) %>%
  filter(n()==16) %>%
  distinct(muts)%>% View()

# Muts within cryptic sublineages
df.muts %>% group_by(muts, lineage) %>%
  filter(n()==8) %>%
  distinct(lineage, muts) %>% View()



# # Sublineages which contain a related mutation
df.barcodes <- fun.read_barcodes_file(file.path(ls.parameters$pipeline_dir,"usher_barcodes.feather"))

tibble(
  group = c("B.11", "B.1.221"),
  sublineages = list("B.11", c("B.1.221", "B.1.221.1", "B.1.221.2", "B.1.221.3", "B.1.221.4"))
) %>%
  unnest(sublineages) %>%
  mutate(muts = lapply(sublineages, function(x) df.barcodes$muts[df.barcodes$index == x & df.barcodes$vals == 1])) %>%
  unnest(muts) %>%
  mutate(
    cryptic_sublineages = mapply(x = group, y = muts, function(x,y){df.muts$index[df.muts$lineage== x & df.muts$muts == y]}, SIMPLIFY = TRUE),
    n_cryptic_sublineages = sapply(cryptic_sublineages, length),
    cryptic_sublineages = sapply(cryptic_sublineages, function(x) paste0(gsub("CRYPT.", "", x), collapse = ";"))
  ) %>% View()


