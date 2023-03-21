# get each visit dataset underlying each model

# data owned by recording schemes

# df_Trichoptera <- readRDS("/data-s3/occmods/Trichoptera/input_data/2021_Ellcur/visitData_caddisflies.rds")$occDetdata %>%
#   tibble::add_column(Taxa = "Trichoptera")
# 
# df_Aquatic <- readRDS("/data-s3/occmods/AquaticBugs/input_data/2021_Ellcur/visitData_aquatic_bugs_20210219.rds")$occDetdata %>%
#   tibble::add_column(Taxa = "Aquatic bugs")
# 
# df_Plecoptera <- readRDS("/data-s3/occmods/Plecoptera/input_data/2021_Ellcur/visitData_stoneflies.rds")$occDetdata %>%
#   tibble::add_column(Taxa = "Plecoptera")
# 
# df_Ephemeroptera <- readRDS("/data-s3/occmods/Ephemeroptera/input_data/2021_Ellcur/visitData_mayflies.rds")$occDetdata %>%
#   tibble::add_column(Taxa = "Ephemeroptera")

# combine all
all_df <- dplyr::bind_rows(df_Trichoptera,
                           df_Aquatic,
                           df_Plecoptera,
                           df_Ephemeroptera)

# map sites to basins
catchmentData <- read.csv("data/districts_ALL_FINAL.csv") %>%
  dplyr::rename(site = "SQ1_SQUARE",
                catchment = "rbd_name")
unique(catchmentData$catchment)

# for plots of data availability
sum_df <- all_df %>%
  dplyr::inner_join(., catchmentData, by="site") %>%
  dplyr::filter(TP > 1969 & TP <2021) %>%
  dplyr::filter(catchment != "Dee") %>%
  dplyr::group_by(TP, catchment) %>%
  dplyr::mutate(SY = paste(site,TP)) %>%
  dplyr::summarise(nuSY = length(unique(SY)),
                   nuVisits = length(unique(visit))) %>% 
  dplyr::mutate(catchment = factor(catchment, levels = c("Solway Tweed", "Northumbria", "North West", "Humber", "Anglian", "Severn", "Thames", "South West", "South East")))

set.seed(123)
tol_catch <- sample(c('#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00', '#DDDDDD'))
tol <- c(tol_catch[[2]], tol_catch[[3]], tol_catch[[1]], tol_catch[[4]], tol_catch[[5]], tol_catch[[9]], tol_catch[[6]], tol_catch[[8]], tol_catch[[7]])

g1 <- ggplot2::ggplot(sum_df, ggplot2::aes(x = TP, y = nuVisits, group = catchment, colour = catchment)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(values = tol) +
  ggplot2::labs(x = "Year", y = "Number of visits", colour = "") +
  ggplot2::theme_bw()

g2 <- ggplot2::ggplot(sum_df, ggplot2::aes(x = TP, y = nuSY, group = catchment, colour = catchment)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(values = tol) +
  ggplot2::labs(x = "Year", y = "Number of site-years", colour = "") +
  ggplot2::theme_bw()

cowplot::plot_grid(g1, g2, nrow = 2)

ggplot2::ggsave("outputs/fig_s2_dataAvail_ts.png", width = 6, height = 6)

# for models to check relationship model predictions
period_df <- all_df %>%
  dplyr::inner_join(., catchmentData, by = "site") %>%
  dplyr::filter(TP > 1969 & TP <2021) %>%
  dplyr::mutate(SY = paste(site, TP),
                period = ifelse(TP %in% 1970:1992, "First (1971-1992)", "Second (1993-2020)")) %>%
  dplyr::group_by(period, catchment) %>%
  dplyr::summarise(nuSY = length(unique(SY)),
                   nuVisits = length(unique(visit))) 

saveRDS(period_df, file = "data/period_df.rds")
