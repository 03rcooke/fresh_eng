#get each visit dataset underlying each model

df_Trichoptera <- readRDS("/data-s3/occmods/Trichoptera/input_data/2021_Ellcur/visitData_caddisflies.rds")$occDetdata %>%
  tibble::add_column(Taxa = "Trichoptera")

df_Aquatic <- readRDS("/data-s3/occmods/AquaticBugs/input_data/2021_Ellcur/visitData_aquatic_bugs_20210219.rds")$occDetdata %>%
  tibble::add_column(Taxa = "Aquatic bugs")

df_Plecoptera <- readRDS("/data-s3/occmods/Plecoptera/input_data/2021_Ellcur/visitData_stoneflies.rds")$occDetdata %>%
  tibble::add_column(Taxa = "Plecoptera")

df_Ephemeroptera <- readRDS("/data-s3/occmods/Ephemeroptera/input_data/2021_Ellcur/visitData_mayflies.rds")$occDetdata %>%
  tibble::add_column(Taxa = "Ephemeroptera")


#combine all
all_df <- bind_rows(df_Trichoptera,
                    df_Aquatic,
                    df_Plecoptera,
                    df_Ephemeroptera)

#map sites to catchments
catchmentData <- read.csv("districts_ALL_FINAL.csv") %>%
                  dplyr::rename(site = "SQ1_SQUARE",
                                catchment = "rbd_name")
unique(catchmentData$catchment)

#for plots of data availability
sum_df <- all_df %>%
  inner_join(., catchmentData, by="site") %>%
  filter(TP > 1969 & TP <2021) %>%
  filter(catchment!="Dee") %>%
  group_by(TP, catchment) %>%
  mutate(SY = paste(site,TP)) %>%
  summarise(nuSY = length(unique(SY)),
            nuVisits = length(unique(visit))) 

g1 <- ggplot(sum_df, aes(x=TP, y = nuVisits, group=catchment, colour=catchment)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_bw() + 
  xlab("Year") + ylab("Number of visits")


g2 <- ggplot(sum_df, aes(x=TP, y = nuSY, group=catchment, colour=catchment)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_bw() + 
  xlab("Year") + ylab("Number of site-year data points")

cowplot::plot_grid(g1,g2,nrow=2)

ggsave("dataAvail_catchment_ts.png", width=5, height=6)

#for models to check relationship model predictions
period_df <- all_df %>%
            inner_join(., catchmentData, by="site") %>%
            filter(TP > 1969 & TP <2021) %>%
            mutate(SY = paste(site,TP),
                   Period = ifelse(TP %in% 1970:1992, "first (1971-1992)", "second(1993-2020)")) %>%
            group_by(Period, catchment) %>%
            summarise(nuSY = length(unique(SY)),
                      nuVisits = length(unique(visit))) 

g1 <- ggplot(period_df) +
  geom_point(aes(x = catchment, y = nuSY, colour = Period))+
  ylab("number of site-years") +
  coord_flip() + theme_bw() + 
  theme(legend.position = "top")

g2 <- ggplot(period_df) +
  geom_point(aes(x = catchment, y = nuVisits, colour = Period))+
  ylab("number of visits") +
  coord_flip() + theme_bw() + 
  theme(legend.position = "top")

cowplot::plot_grid(g1, g2, nrow=1)
    
saveRDS(period_df, file = "/data/notebooks/rstudio-fresh/fresh_diana/dataAvailabilty_ept.rds")
            

