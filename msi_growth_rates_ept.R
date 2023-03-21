# script to calculate growth rates in the first and second half of the time-series

# get occupancy model predictions - using ones processed by rob - see 'fresh_eng.R'

library(remotes) # remotes: package installation
library(sf) # sf: spatial manipulation
library(tidyr) # tidyr: data manipulation
library(dplyr) # dplyr: data manipulation
library(rjags) # rjags: occupancy models
library(R2jags) # R2jags: occupancy models
library(HDInterval) # HDInterval: credible intervals
library(effsize) # effsize: effect sizes
library(ggplot2) # ggplot2: plotting
library(ggrepel) # ggrepel: plotting
library(cowplot) # cowplot: plotting

# # wrappeR package from github
# remotes::install_github("https://github.com/03rcooke/wrappeR", ref = "main")
library(wrappeR) # wrappeR: multi-species indicators

# # commit packages to renv
# renv::snapshot()

## small helper functions
ggplot2::theme_set(cowplot::theme_cowplot() +
                     ggplot2::theme(plot.background = element_rect(fill = "white", colour = "white")))

# log of 0 is undefined
nudgeOcc <- function(x, nudgeFac = 0.0001) {
  x[x == 0] <- nudgeFac
  return(x)
}

# function to calculate geometric mean
geomean <- function(x) exp(mean(log(nudgeOcc(x))))

# load_rdata function
# loads an RData file, and assigns it to an object name
load_rdata <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## global parameters
startyear <- 1971
endyear <- 2020
ci <- 0.95
ggplot2::theme_set(cowplot::theme_cowplot())

## compile models

# taxonomic group(s)
taxa <- c("Trichoptera", "Ephemeroptera", "Plecoptera")
ntax <- length(taxa)

regions <- c("Anglian", "Severn", "Humber", "Solway.Tweed", "Northumbria", "South.East", 
             "South.West", "Thames", "North.West")
nregions <- length(regions)

# # create roster - run for multiple groups & multiple regions
# roster <- wrappeR::createRoster(index = 1:(nregions * ntax),
#                        modPath = "/data-s3/occmods/",
#                        metaPath = "/data-s3/metadata/",
#                        ver = "2021_Ellcur",
#                        indicator = "all",
#                        region = rep(regions, times = ntax),
#                        nSamps = 999,
#                        minObs = 1,
#                        scaleObs = "region",
#                        write = TRUE,
#                        outPath = "~/fresh_rob/data/filtered/",
#                        group = rep(taxa, each = nregions),
#                        t0 = startyear,
#                        tn = endyear)
# 
# # filtered data
# filt_tax <- lapply(roster, wrappeR::applySamp, parallel = FALSE)

## prep occupancy data per region

# focal years 1971 to 2020
foc_yrs <- paste0("year_", startyear:endyear)

# read in filtered data
occ_read <- function(x, reg) {
  
  load_rdata(paste0("data/filtered/", x, "_all_", reg, "_samp.rdata")) %>% 
    .[[1]] %>% 
    dplyr::mutate(grp = x)
  
}

occ_prep <- function(taxa, reg, foc_yrs) {
  
  # observation metadata
  meta_out <- lapply(taxa, function(x) {
    
    load_rdata(paste0("data/filtered/", x, "_all_", reg, "_samp.rdata")) %>% 
      .[[2]] %>% 
      dplyr::mutate(tax_grp = x)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  # species models to keep
  pass_keep <- meta_out %>% 
    dplyr::filter(!!sym(paste0("rot_EqualWt_r_", reg)) == TRUE) %>% 
    .[,paste0("Species_r_", reg)]
  
  # prepare occupancy data for pass_keep species
  occ_prep <- lapply(taxa, function(x) occ_read(x = x, reg = reg)) %>%
    dplyr::bind_rows(.) %>%
    # trim to focal years 1971 to 2020
    dplyr::select(dplyr::all_of(foc_yrs), iteration, species, grp) %>%
    # filter to pass_keep species
    dplyr::filter(species %in% pass_keep) %>% 
    # add region aggregate name
    dplyr::mutate(agg = reg)
  
}

occ_reg <- lapply(regions, function(x) occ_prep(taxa = taxa, reg = x, foc_yrs = foc_yrs)) %>% 
  dplyr::bind_rows()

spp_agg <- occ_reg %>% 
  dplyr::distinct(species, grp, agg) %>% 
  dplyr::count(agg, grp)

spp_tot <- occ_reg %>% 
  dplyr::distinct(species, grp) %>% 
  dplyr::count(grp)

## annual multispecies index - geometric mean occupancy

## per catchment

msi_occ <- occ_reg %>%
  tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
  dplyr::group_by(year, iteration, agg) %>%
  # geometric mean
  dplyr::summarise(gm = geomean(occ)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(year = as.numeric(gsub("year_", "", year)))

msi_sum_occ <- msi_occ %>% 
  dplyr::group_by(agg, year) %>% 
  dplyr::summarise(dplyr::across(c(gm), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::mutate(dplyr::across(c(gm_median), ~ .x / (ifelse(!is.na(dplyr::first(.x)), dplyr::first(.x), dplyr::nth(.x, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::mutate(dplyr::across(c(gm_low, gm_upp), ~ .x / (ifelse(!is.na(dplyr::first(gm_median)), dplyr::first(gm_median), dplyr::nth(gm_median, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = dplyr::recode(agg, North.West = "North West", Solway.Tweed = "Solway Tweed", Northumbria = "Northumbria", Humber = "Humber", Anglian = "Anglian", Thames = "Thames", South.East = "South East", South.West = "South West", Severn = "Severn")) %>% 
  dplyr::mutate(agg = factor(agg, levels = c("Solway Tweed", "Northumbria", "North West", "Humber", "Anglian", "Severn", "Thames", "South West", "South East")))

## overall

msi_overall <- occ_reg %>% 
  tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
  dplyr::group_by(year, iteration) %>%
  # geometric mean
  dplyr::summarise(gm = geomean(occ)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(year = as.numeric(gsub("year_", "", year)))

msi_sum_overall <- msi_overall %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(dplyr::across(c(gm), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::mutate(dplyr::across(c(gm_median), ~ .x / (ifelse(!is.na(dplyr::first(.x)), dplyr::first(.x), dplyr::nth(.x, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::mutate(dplyr::across(c(gm_low, gm_upp), ~ .x / (ifelse(!is.na(dplyr::first(gm_median)), dplyr::first(gm_median), dplyr::nth(gm_median, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::ungroup() 

## identify trough for each group using the overall MSI

msi_sum_overall

msi_sum_overall %>%
  filter(year > 1980) %>%
  summarise(minYear = year[which.min(ind_gm_median)])

# 1993

# msi growth rate for each period

# first period
msi_sum_occ_first <- msi_occ %>%
  filter(year %in% c(1971, 1992)) %>%
  pivot_wider(everything(), names_from = year, values_from = "gm") %>%
  janitor::clean_names() %>%
  dplyr::mutate(gr_rate =  ((x1992 / x1971) ^ (1 / length(1971:1992)) - 1) * 100) %>%
  dplyr::group_by(agg) %>% 
  dplyr::summarise(dplyr::across(c(gr_rate), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = dplyr::recode(agg, North.West = "North West", Solway.Tweed = "Solway Tweed", Northumbria = "Northumbria", Humber = "Humber", Anglian = "Anglian", Thames = "Thames", South.East = "South East", South.West = "South West", Severn = "Severn")) %>% 
  dplyr::mutate(agg = factor(agg, levels = c("Solway Tweed", "Northumbria", "North West", "Humber", "Anglian", "Severn", "Thames", "South West", "South East"))) %>% 
  tibble::add_column(period = "First (1971-1992)")

# second period
msi_sum_occ_second <- msi_occ %>%
  filter(year %in% c(1993, 2020)) %>%
  pivot_wider(everything(), names_from = year, values_from = "gm") %>%
  janitor::clean_names() %>%
  dplyr::mutate(gr_rate =  ((x2020 /x1993) ^ (1 / length(1993:2020)) - 1) * 100) %>%
  dplyr::group_by(agg) %>% 
  dplyr::summarise(dplyr::across(c(gr_rate), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = dplyr::recode(agg, North.West = "North West", Solway.Tweed = "Solway Tweed", Northumbria = "Northumbria", Humber = "Humber", Anglian = "Anglian", Thames = "Thames", South.East = "South East", South.West = "South West", Severn = "Severn")) %>% 
  dplyr::mutate(agg = factor(agg, levels = c("Solway Tweed", "Northumbria", "North West", "Humber", "Anglian", "Severn", "Thames", "South West", "South East"))) %>% 
  tibble::add_column(period = "Second (1993-2020)")

# plot
msi_growth_rates <- dplyr::bind_rows(msi_sum_occ_first, msi_sum_occ_second) %>%
  dplyr::mutate(period = factor(period, levels = rev(c("First (1971-1992)", "Second (1993-2020)")))) %>% 
  ggplot2::ggplot() +
  ggplot2::geom_pointrange(ggplot2::aes(x = agg, 
                                        y = gr_rate_median, ymin = gr_rate_low, ymax = gr_rate_upp,
                                        colour = period),
                           size = rel(0.25), position = ggplot2::position_dodge(0.3)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::coord_flip() +
  ggplot2::scale_colour_brewer(name = "Period", palette = "Set1", breaks = c("First (1971-1992)", "Second (1993-2020)")) +
  ggplot2::scale_x_discrete(limits = rev) +
  ggplot2::scale_y_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  ggplot2::labs(y = "Growth rate") +
  ggplot2::theme(legend.position = "bottom", 
                 axis.text = ggplot2::element_text(size = rel(0.65)),
                 axis.title.y = ggplot2::element_blank(),
                 legend.justification = "centre",
                 plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

cowplot::save_plot("outputs/fig_2_msi_growth_rates.png", msi_growth_rates, base_width = 6, base_height = 5)

# relate trends to data availability

# see data_availability_ept.R
dataAvailability <- readRDS("data/period_df.rds") %>%
  dplyr::mutate(agg = catchment)

# add onto trend data
all_data <- dplyr::bind_rows(msi_sum_occ_first, msi_sum_occ_second) %>%
  dplyr::inner_join(., dataAvailability, by = c("period", "agg"))

ggplot2::ggplot(all_data, ggplot2::aes(x = nuSY, y = gr_rate_median)) +
  ggplot2::geom_pointrange(ggplot2::aes(x = nuSY, y = gr_rate_median, ymin = gr_rate_low, ymax = gr_rate_upp, colour = period)) +
  ggplot2::geom_text(ggplot2::aes(label = agg), size = 2.75, angle = 90) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  ggplot2::scale_colour_brewer(name = "Period", palette = "Set1", direction = -1, breaks = c("First (1971-1992)", "Second (1993-2020)")) +
  ggplot2::labs(x = "Number of site-years", y = "Growth rate") +
  ggplot2::theme(legend.position = "bottom",
                 plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

ggplot2::ggsave("outputs/fig_s3_dataAvailability_impact.png", width = 5, height = 4)

# explore effect on precision

sd_occ <- occ_reg %>%
  tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
  dplyr::group_by(year, grp, species, agg) %>%
  dplyr::summarise(occ_sd = sd(occ)) %>% 
  dplyr::group_by(year, agg) %>%
  dplyr::summarise(median_occ_sd = median(occ_sd)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(year = as.numeric(gsub("year_", "", year))) %>%
  dplyr::filter(year > 1970) %>%
  dplyr::filter(year %in% c(1971, 1992, 1993, 2020)) %>%
  dplyr::mutate(period = ifelse(year %in% 1971:1992, "First (1971-1992)", "Second (1993-2020)")) %>%
  dplyr::group_by(period, agg) %>%
  dplyr::summarise(median_occ_sd = median(median_occ_sd)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(agg = gsub("\\."," ", agg))

# merge with data available
all_data_sd <- dplyr::inner_join(sd_occ, all_data, by = c("period","agg"))

g1 <- ggplot2::ggplot(all_data_sd, ggplot2::aes(x = nuSY, y = median_occ_sd)) +
  ggplot2::geom_point(ggplot2::aes(colour = period)) +
  ggplot2::geom_text(ggplot2::aes(label = catchment), size = 2.75, angle = 90) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_colour_brewer(name = "Period", palette = "Set1", direction = -1, breaks = c("First (1971-1992)", "Second (1993-2020)")) +
  ggplot2::labs(x = "Number of site-years", y = "Median SD of species\nannual occupancy") +
  ggplot2::theme(legend.position = "none",
                 plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

# relationship between mean growth rate and sd of growth rate
g2 <- ggplot2::ggplot(all_data_sd, ggplot2::aes(y = abs(gr_rate_median), x = median_occ_sd)) +
  ggplot2::geom_point(ggplot2::aes(colour = period)) +
  ggplot2::geom_text(ggplot2::aes(label = catchment), size = 2.75, angle = 90) +
  ggplot2::scale_colour_brewer(name = "Period", palette = "Set1", direction = -1, breaks = c("First (1971-1992)", "Second (1993-2020)")) +
  ggplot2::labs(x = "Median SD of species annual occupancy", y = "\n|Median Growth rate|") +
  ggplot2::theme(legend.position = "bottom",
                 plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

cowplot::plot_grid(g1, g2, nrow = 2)

ggplot2::ggsave("outputs/fig_s4_dataAvailability_rel.png", width = 6, height = 7)
