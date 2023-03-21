## Rob freshwater

# install JAGS in terminal
# sudo apt install jags

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
library(ggh4x) # ggh4x: plotting

# # wrappeR package from github
# remotes::install_github("https://github.com/03rcooke/wrappeR", ref = "main")
library(wrappeR) # wrappeR: multi-species indicators

# # commit packages to renv
# renv::snapshot()

## small helper functions

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
ggplot2::theme_set(cowplot::theme_cowplot() +
                     ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = "white")))

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

## Map of catchments

catch <- sf::read_sf(dsn = "data", layer = "WFD_River_Basin_Districts_Cycle_2")
UK <- sf::read_sf(dsn = "data", layer = "Countries__December_2019__Boundaries_UK_BFC")

england <- dplyr::filter(UK, objectid == "1")

catch_eng <- sf::st_intersection(england, catch) %>% 
  # drop Dee catchment
  dplyr::filter(rbd_name != "Dee") %>% 
  dplyr::mutate(rbd_name = factor(rbd_name, levels = c("North West", "Solway Tweed", "Northumbria", "Humber", "Anglian", "Thames", "South East", "South West", "Severn"))) %>% 
  dplyr::mutate(rbd_name2 = as.factor(as.character(as.numeric(rbd_name))))

set.seed(123)
tol_catch <- sample(c('#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00', '#DDDDDD'))

catch_map <- ggplot2::ggplot(data = dplyr::filter(UK, objectid %in% c(1,3,4))) +
  ggplot2::geom_sf(size = 0.2, fill = "white", colour = "black") +
  ggplot2::geom_sf(data = catch_eng, ggplot2::aes(fill = rbd_name), col = "black", alpha = 0.5, lwd = 0.5) +
  ggplot2::scale_fill_manual(values = tol_catch) +
  ggplot2::theme(axis.line = ggplot2::element_blank(),
                 axis.title = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 axis.text = ggplot2::element_blank(),
                 legend.position = "none")

## multispecies index - geometric mean occupancy

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
  dplyr::mutate(agg = factor(agg, levels = c("Solway Tweed", "Northumbria", "North West", "Humber", "Anglian", "Severn", "Thames", "South West", "South East"))) %>% 
  dplyr::mutate(dplyr::across(c(ind_gm_median, ind_gm_low, ind_gm_upp), ~ ((.x - 100) / (2020 - 1971)), .names = "py_{.col}"))

## Overall

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
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = "Overall")

tol <- c(tol_catch[[2]], tol_catch[[3]], tol_catch[[1]], tol_catch[[4]], tol_catch[[5]], tol_catch[[9]], tol_catch[[6]], tol_catch[[8]], tol_catch[[7]])

strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = tol))

## ind_plot

ind_plot <- ggplot2::ggplot(msi_sum_occ, ggplot2::aes(x = year, y = ind_gm_median)) +
  ggplot2::geom_hline(yintercept = 100, colour = "black", lty = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ind_gm_low, ymax = ind_gm_upp, fill = agg), alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = ind_gm_median, colour = agg)) +
  ggplot2::geom_point(ggplot2::aes(colour = agg), size = 0.4) +
  ggplot2::scale_fill_manual(values = tol) +
  ggplot2::scale_colour_manual(values = tol) +
  ggplot2::scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggh4x::facet_wrap2(~ agg, strip = strip) +
  ggplot2::labs(x = "Year", y = "Index of occupancy") +
  ggplot2::theme(legend.position = "none",
                 panel.spacing = unit(1.5, "lines"))

ind_overall <- ggplot2::ggplot(msi_sum_overall, ggplot2::aes(x = year, y = ind_gm_median)) +
  ggplot2::geom_hline(yintercept = 100, colour = "black", lty = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ind_gm_low, ymax = ind_gm_upp), fill = "darkgrey", alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = ind_gm_median), colour = "darkgrey") +
  ggplot2::geom_point(colour = "darkgrey", size = 0.4) +
  ggplot2::scale_x_continuous(name = "Year", breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggplot2::facet_wrap(~ agg) +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "darkgrey"),
                 legend.position = "none",
                 axis.title.y = ggplot2::element_blank())

comb_plot <- cowplot::ggdraw(cowplot::plot_grid(ind_plot, cowplot::plot_grid(NULL, ind_overall, NULL, nrow = 3, rel_heights = c(0.9, 0.8, 0.1)), rel_widths = c(1, 0.5))) +
  cowplot::draw_plot(catch_map, 0.55, 0.49, 0.6, 0.6) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

cowplot::save_plot("outputs/fig_1_ind.png", comb_plot, base_width = 16, base_height = 9)

## occ_plot

occ_plot <- ggplot2::ggplot(msi_sum_occ, ggplot2::aes(x = year, y = gm_median)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = gm_low, ymax = gm_upp, fill = agg), alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = gm_median, colour = agg)) +
  ggplot2::geom_point(ggplot2::aes(colour = agg), size = 0.4) +
  ggplot2::scale_fill_manual(values = tol) +
  ggplot2::scale_colour_manual(values = tol) +
  ggplot2::scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggh4x::facet_wrap2(~ agg, strip = strip) +
  ggplot2::labs(x = "Year", y = "Geometric mean occupancy") +
  ggplot2::theme(legend.position = "none",
                 panel.spacing = unit(1.5, "lines"))

occ_overall <- ggplot2::ggplot(msi_sum_overall, ggplot2::aes(x = year, y = gm_median)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = gm_low, ymax = gm_upp), fill = "darkgrey", alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = gm_median), colour = "darkgrey") +
  ggplot2::geom_point(size = 0.4, colour = "darkgrey") +
  ggplot2::scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggplot2::facet_wrap(~ agg) +
  ggplot2::labs(x = "Year") +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "darkgrey"),
                 legend.position = "none",
                 axis.title.y = ggplot2::element_blank())

occ_comb_plot <- cowplot::ggdraw(cowplot::plot_grid(occ_plot, cowplot::plot_grid(NULL, occ_overall, NULL, nrow = 3, rel_heights = c(0.9, 0.8, 0.1)), rel_widths = c(1, 0.5))) +
  cowplot::draw_plot(catch_map, 0.55, 0.49, 0.6, 0.6) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

cowplot::save_plot("outputs/fig_s1_raw_occ.png", occ_comb_plot, base_width = 16, base_height = 9)

## Freshwater riverflies

msi_fresh <- occ_reg %>% 
  tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
  dplyr::group_by(year, iteration) %>%
  # geometric mean
  dplyr::summarise(gm = geomean(occ)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(year = as.numeric(gsub("year_", "", year)))

msi_sum_fresh <- msi_fresh %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(dplyr::across(c(gm), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::mutate(dplyr::across(c(gm_median), ~ .x / (ifelse(!is.na(dplyr::first(.x)), dplyr::first(.x), dplyr::nth(.x, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::mutate(dplyr::across(c(gm_low, gm_upp), ~ .x / (ifelse(!is.na(dplyr::first(gm_median)), dplyr::first(gm_median), dplyr::nth(gm_median, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = "Freshwater")

# when is trough
msi_sum_fresh$year[which.min(msi_sum_fresh$ind_gm_median)]

## number of records

meta_prep <- function(taxa, reg) {
  
  # observation metadata
  meta_out <- lapply(taxa, function(x) {
    
    load_rdata(paste0("data/filtered/", x, "_all_", reg, "_samp.rdata")) %>% 
      .[[2]] %>% 
      dplyr::mutate(tax_grp = x)
    
  }) %>% 
    dplyr::bind_rows(.) %>% 
    # add region aggregate name
    dplyr::mutate(agg = reg)
  
  colnames(meta_out) <- gsub(x = colnames(meta_out), pattern = reg, replacement = " ")
  
  return(meta_out)
  
}

meta_reg <- lapply(regions, function(x) meta_prep(taxa = taxa, reg = x)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::filter(`Species_r_ ` %in% unique(occ_reg$species))

rec_tot <- meta_reg %>% 
  dplyr::group_by(tax_grp) %>% 
  dplyr::summarise(recs = sum(`n_obs_regional_r_ `),
                   recs_med = median(`n_obs_regional_r_ `)) %>% 
  dplyr::rename(grp = tax_grp)

spp_tot <- occ_reg %>% 
  dplyr::distinct(species, grp) %>% 
  dplyr::count(grp)

rep_tot <- meta_reg %>% 
  dplyr::group_by(tax_grp) %>% 
  dplyr::mutate(reps = `rot_prop_repeats_grp_r_ ` * `n_obs_regional_r_ `) %>% 
  dplyr::summarise(reps_med = median(reps)) %>% 
  dplyr::rename(grp = tax_grp)

info_tot <- dplyr::left_join(spp_tot, rec_tot, by = c("grp")) %>% 
  dplyr::left_join(rep_tot, by = c("grp"))

sum(rec_tot$recs)

rec_reg <- meta_reg %>% 
  dplyr::group_by(tax_grp, agg) %>% 
  dplyr::summarise(recs = sum(`n_obs_regional_r_ `),
                   recs_med = median(`n_obs_regional_r_ `)) %>% 
  dplyr::rename(grp = tax_grp)

spp_reg <- occ_reg %>% 
  dplyr::distinct(species, grp, agg) %>% 
  dplyr::count(grp, agg)

rep_reg <- meta_reg %>% 
  dplyr::group_by(tax_grp, agg) %>% 
  dplyr::mutate(reps = `rot_prop_repeats_grp_r_ ` * `n_obs_regional_r_ `) %>% 
  dplyr::summarise(reps_med = median(reps)) %>% 
  dplyr::rename(grp = tax_grp)

info_reg <- dplyr::left_join(spp_reg, rec_reg, by = c("grp", "agg")) %>% 
  dplyr::left_join(rep_reg, by = c("grp", "agg")) %>% 
  dplyr::mutate(grp = dplyr::recode(grp, Ephemeroptera = "Mayflies", Plecoptera = "Stoneflies", Trichoptera = "Caddisflies"))

## multispecies index - geometric mean occupancy per TAXONOMIC GROUP

msi_occ_tax <- occ_reg %>%
  tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
  dplyr::group_by(year, grp, iteration, agg) %>%
  # geometric mean
  dplyr::summarise(gm = geomean(occ)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(year = as.numeric(gsub("year_", "", year)))

msi_sum_occ_tax <- msi_occ_tax %>% 
  dplyr::group_by(agg, grp, year) %>% 
  dplyr::summarise(dplyr::across(c(gm), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::mutate(dplyr::across(c(gm_median), ~ .x / (ifelse(!is.na(dplyr::first(.x)), dplyr::first(.x), dplyr::nth(.x, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::mutate(dplyr::across(c(gm_low, gm_upp), ~ .x / (ifelse(!is.na(dplyr::first(gm_median)), dplyr::first(gm_median), dplyr::nth(gm_median, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = dplyr::recode(agg, North.West = "North West", Solway.Tweed = "Solway Tweed", Northumbria = "Northumbria", Humber = "Humber", Anglian = "Anglian", Thames = "Thames", South.East = "South East", South.West = "South West", Severn = "Severn")) %>% 
  dplyr::mutate(agg = factor(agg, levels = c("Solway Tweed", "Northumbria", "North West", "Humber", "Anglian", "Severn", "Thames", "South West", "South East"))) %>% 
  dplyr::mutate(grp = dplyr::recode(grp, Ephemeroptera = "Mayflies (Ephemeroptera)", Plecoptera = "Stoneflies (Plecoptera)", Trichoptera = "Caddisflies (Trichoptera)")) %>% 
  dplyr::mutate(dplyr::across(c(ind_gm_median, ind_gm_low, ind_gm_upp), ~ ((.x - 100) / (2020 - 1971)), .names = "py_{.col}"))

## Overall

msi_overall_tax <- occ_reg %>% 
  tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
  dplyr::group_by(year, grp, iteration) %>%
  # geometric mean
  dplyr::summarise(gm = geomean(occ)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(year = as.numeric(gsub("year_", "", year)))

msi_sum_overall_tax <- msi_overall_tax %>% 
  dplyr::group_by(grp, year) %>% 
  dplyr::summarise(dplyr::across(c(gm), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::mutate(dplyr::across(c(gm_median), ~ .x / (ifelse(!is.na(dplyr::first(.x)), dplyr::first(.x), dplyr::nth(.x, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::mutate(dplyr::across(c(gm_low, gm_upp), ~ .x / (ifelse(!is.na(dplyr::first(gm_median)), dplyr::first(gm_median), dplyr::nth(gm_median, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = "Overall") %>% 
  dplyr::mutate(grp = dplyr::recode(grp, Ephemeroptera = "Mayflies (Ephemeroptera)", Plecoptera = "Stoneflies (Plecoptera)", Trichoptera = "Caddisflies (Trichoptera)")) %>% 
  dplyr::mutate(dplyr::across(c(ind_gm_median, ind_gm_low, ind_gm_upp), ~ ((.x - 100) / (2020 - 1971)), .names = "py_{.col}"))

tax_col <- c('#DDAA33', '#BB5566', '#004488')

## ind_plot_tax

ind_plot_tax <- ggplot2::ggplot(msi_sum_occ_tax, ggplot2::aes(x = year, y = ind_gm_median)) +
  ggplot2::geom_hline(yintercept = 100, colour = "black", lty = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ind_gm_low, ymax = ind_gm_upp, fill = grp), alpha = 0.15) +
  ggplot2::geom_line(ggplot2::aes(y = ind_gm_median, colour = grp)) +
  ggplot2::geom_point(ggplot2::aes(colour = grp), size = 0.4) +
  ggplot2::scale_fill_manual(values = tax_col) +
  ggplot2::scale_colour_manual(values = tax_col) +
  ggplot2::scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggh4x::facet_wrap2(~ agg, strip = strip) +
  ggplot2::labs(x = "Year", y = "Index of occupancy") +
  ggplot2::theme(legend.position = "none",
                 panel.spacing = unit(1.5, "lines"))

ind_overall_tax <- ggplot2::ggplot(msi_sum_overall_tax, ggplot2::aes(x = year, y = ind_gm_median)) +
  ggplot2::geom_hline(yintercept = 100, colour = "black", lty = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ind_gm_low, ymax = ind_gm_upp, fill = grp), alpha = 0.15) +
  ggplot2::geom_line(ggplot2::aes(y = ind_gm_median, colour = grp)) +
  ggplot2::geom_point(ggplot2::aes(colour = grp), size = 0.4) +
  ggplot2::scale_fill_manual(values = tax_col) +
  ggplot2::scale_colour_manual(values = tax_col) +
  ggplot2::scale_x_continuous(name = "Year", breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggplot2::facet_wrap(~ agg) +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "darkgrey"),
                 legend.position = "bottom", 
                 legend.direction = "vertical",
                 legend.justification = "centre",
                 legend.title = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank())

comb_plot_tax <- cowplot::ggdraw(cowplot::plot_grid(ind_plot_tax, cowplot::plot_grid(NULL, ind_overall_tax, NULL, nrow = 3, rel_heights = c(1.2, 1, 0.1)), rel_widths = c(1, 0.5))) +
  cowplot::draw_plot(catch_map, 0.55, 0.49, 0.6, 0.6) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

cowplot::save_plot("outputs/fig_s5_ind_tax.png", comb_plot_tax, base_width = 16, base_height = 9)

## occ_plot

occ_plot_tax <- ggplot2::ggplot(msi_sum_occ_tax, ggplot2::aes(x = year, y = gm_median)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = gm_low, ymax = gm_upp, fill = grp), alpha = 0.15) +
  ggplot2::geom_line(ggplot2::aes(y = gm_median, colour = grp)) +
  ggplot2::geom_point(ggplot2::aes(colour = grp), size = 0.4) +
  ggplot2::scale_fill_manual(values = tax_col) +
  ggplot2::scale_colour_manual(values = tax_col) +
  ggplot2::scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggh4x::facet_wrap2(~ agg, strip = strip) +
  ggplot2::labs(x = "Year", y = "Geometric mean occupancy") +
  ggplot2::theme(legend.position = "none",
                 panel.spacing = unit(1.5, "lines"))

occ_overall_tax <- ggplot2::ggplot(msi_sum_overall_tax, ggplot2::aes(x = year, y = gm_median)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = gm_low, ymax = gm_upp, fill = grp), alpha = 0.15) +
  ggplot2::geom_line(ggplot2::aes(y = gm_median, colour = grp)) +
  ggplot2::geom_point(ggplot2::aes(colour = grp), size = 0.4) +
  ggplot2::scale_fill_manual(values = tax_col) +
  ggplot2::scale_colour_manual(values = tax_col) +
  ggplot2::scale_x_continuous(name = "Year", breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggplot2::facet_wrap(~ agg) +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "darkgrey"),
                 legend.position = "bottom", 
                 legend.direction = "vertical",
                 legend.justification = "centre",
                 legend.title = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank())

occ_comb_plot_tax <- cowplot::ggdraw(cowplot::plot_grid(occ_plot_tax, cowplot::plot_grid(NULL, occ_overall_tax, NULL, nrow = 3, rel_heights = c(1.2, 1, 0.1)), rel_widths = c(1, 0.5))) +
  cowplot::draw_plot(catch_map, 0.55, 0.49, 0.6, 0.6) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

cowplot::save_plot("outputs/fig_s6_raw_occ_tax.png", occ_comb_plot_tax, base_width = 16, base_height = 9)

## Only present in all catchments

## Disaggregated by TAXONOMIC GROUP

spp_unif <- occ_reg %>% 
  dplyr::distinct(species, grp, agg) %>% 
  dplyr::count(species, grp) %>% 
  dplyr::filter(n >= 9)

table(spp_unif$grp)

msi_unif <- occ_reg %>%
  dplyr::filter(species %in% spp_unif$species) %>% 
  tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
  dplyr::group_by(year, iteration, agg) %>%
  # geometric mean
  dplyr::summarise(gm = geomean(occ)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(year = as.numeric(gsub("year_", "", year)))

msi_sum_unif <- msi_unif %>% 
  dplyr::group_by(agg, year) %>% 
  dplyr::summarise(dplyr::across(c(gm), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::mutate(dplyr::across(c(gm_median), ~ .x / (ifelse(!is.na(dplyr::first(.x)), dplyr::first(.x), dplyr::nth(.x, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::mutate(dplyr::across(c(gm_low, gm_upp), ~ .x / (ifelse(!is.na(dplyr::first(gm_median)), dplyr::first(gm_median), dplyr::nth(gm_median, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = dplyr::recode(agg, North.West = "North West", Solway.Tweed = "Solway Tweed", Northumbria = "Northumbria", Humber = "Humber", Anglian = "Anglian", Thames = "Thames", South.East = "South East", South.West = "South West", Severn = "Severn")) %>% 
  dplyr::mutate(agg = factor(agg, levels = c("Solway Tweed", "Northumbria", "North West", "Humber", "Anglian", "Severn", "Thames", "South West", "South East"))) %>% 
  dplyr::mutate(dplyr::across(c(ind_gm_median, ind_gm_low, ind_gm_upp), ~ ((.x - 100) / (2020 - 1971)), .names = "py_{.col}"))

msi_plot_unif <- ggplot2::ggplot(msi_sum_unif, ggplot2::aes(x = year, y = gm_median)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = gm_low, ymax = gm_upp, fill = agg), alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = gm_median, colour = agg)) +
  ggplot2::geom_point(ggplot2::aes(colour = agg), size = 0.4) +
  ggplot2::scale_fill_manual(values = tol) +
  ggplot2::scale_colour_manual(values = tol) +
  ggplot2::scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggh4x::facet_wrap2(~ agg, strip = strip) +
  ggplot2::labs(x = "Year", y = "Geometric mean occupancy") +
  ggplot2::theme(legend.position = "none",
                 panel.spacing = unit(1.5, "lines"))

ind_plot_unif <- ggplot2::ggplot(msi_sum_unif, ggplot2::aes(x = year, y = ind_gm_median)) +
  ggplot2::geom_hline(yintercept = 100, colour = "black", lty = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ind_gm_low, ymax = ind_gm_upp, fill = agg), alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = ind_gm_median, colour = agg)) +
  ggplot2::geom_point(ggplot2::aes(colour = agg), size = 0.4) +
  ggplot2::scale_fill_manual(values = tol) +
  ggplot2::scale_colour_manual(values = tol) +
  ggplot2::scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggh4x::facet_wrap2(~ agg, strip = strip) +
  ggplot2::labs(x = "Year", y = "Index of occupancy") +
  ggplot2::theme(legend.position = "none",
                 panel.spacing = unit(1.5, "lines"))

## Overall

msi_overall_unif <- occ_reg %>% 
  dplyr::filter(species %in% spp_unif$species) %>% 
  tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
  dplyr::group_by(year, iteration) %>%
  # geometric mean
  dplyr::summarise(gm = geomean(occ)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(year = as.numeric(gsub("year_", "", year)))

msi_sum_overall_unif <- msi_overall_unif %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(dplyr::across(c(gm), list(
    median = ~median(.x, na.rm = TRUE), 
    low = ~HDInterval::hdi(.x, credMass = ci)[[1]],
    upp = ~HDInterval::hdi(.x, credMass = ci)[[2]]))) %>% 
  dplyr::mutate(dplyr::across(c(gm_median), ~ .x / (ifelse(!is.na(dplyr::first(.x)), dplyr::first(.x), dplyr::nth(.x, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::mutate(dplyr::across(c(gm_low, gm_upp), ~ .x / (ifelse(!is.na(dplyr::first(gm_median)), dplyr::first(gm_median), dplyr::nth(gm_median, 2)) / 100), .names = "ind_{.col}")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(agg = "Overall")

ind_overall_unif <- ggplot2::ggplot(msi_sum_overall_unif, ggplot2::aes(x = year, y = ind_gm_median)) +
  ggplot2::geom_hline(yintercept = 100, colour = "black", lty = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ind_gm_low, ymax = ind_gm_upp), fill = "darkgrey", alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = ind_gm_median), colour = "darkgrey") +
  ggplot2::geom_point(colour = "darkgrey", size = 0.4) +
  ggplot2::scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015), expand = c(0, 0)) +
  ggplot2::facet_wrap(~ agg) +
  ggplot2::labs(x = "Year") +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "darkgrey"),
                 legend.position = "none",
                 axis.title.y = ggplot2::element_blank())

comb_plot_unif <- cowplot::ggdraw(cowplot::plot_grid(ind_plot_unif, cowplot::plot_grid(NULL, ind_overall_unif, NULL, nrow = 3, rel_heights = c(0.9, 0.8, 0.1)), rel_widths = c(1, 0.5))) +
  cowplot::draw_plot(catch_map, 0.55, 0.49, 0.6, 0.6) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = "white"))

cowplot::save_plot("outputs/fig_s7_ind_unif.png", comb_plot_unif, base_width = 16, base_height = 9)

# overall
cor(msi_sum_overall$ind_gm_median, msi_sum_overall_unif$ind_gm_median, use = "complete.obs")

catch_names <- as.character(unique(msi_sum_occ$agg))

cor_catch <- lapply(1:length(catch_names), function(x) {
  
  catch_name <- catch_names[[x]]
  
  cor_out <- cor(dplyr::filter(msi_sum_occ, agg == !!catch_name)$ind_gm_median, 
                 dplyr::filter(msi_sum_unif, agg == !!catch_name)$ind_gm_median, use = "complete.obs")
  
  data.frame(agg = catch_name, r = cor_out)
  
}) %>% 
  dplyr::bind_rows()
