################################################################################
####            MAPPING BENTHIC BIODIVERSITY TO FACILITATE FUTURE           ####
####              SUSTAINABLE DEVELOPMENT (ECOSPHERE JOURNAL)               ####
################################################################################

## This script relates to the work in Cooper, K.M., Thompson, M.S.A., Bolam, S.G.,
# Peach, C.M., Webb, T.J., Downie, A-L. Mapping benthic biodiversity to facilitate
# future sustainable development. Ecosphere.

# Data used in the script is sourced from the OneBenthic 
# (https://rconnect.cefas.co.uk/onebenthic_portal/) database using sql 
# queries. For users without direct access to this database, data
# can be sourced using either OneBenthic APIs:
# Faunal data: https://rconnect.cefas.co.uk/onebenthic_api_1/__docs__/)
# Sediment data: https://rconnect.cefas.co.uk/onebenthic_api_3/__docs__/
# or using the OneBenthic Data Extraction tool: Grab/Core
# (https://rconnect.cefas.co.uk/onebenthic_dataextractiongrabcore/)

# This script includes code for running Random Forest models (a quick looksee),
# but the accompanying files: 
#  ClusterModel_2025.R and
#  ContinuousVariablesModel_2025.R 
# should be used for final modelling.

# Input data files for modelling are produced in this script and include:
# biodiv_metrics_4_modelling.csv (BIODIV METRICS)
# biodiv_metric_rare_4_modelling.csv (RARE TAXA)
# biodiv_cluster_4_modelling.csv (BIODIV CLUSTERS)

# Noye, the biodiversity metric calculations were implemented using R code adapted from Thompson’s open-source script:
#https://github.com/MurraySAThompson/biodiversity-estimation-across-spatial-scales-and-Hill-numbers
#_______________________________________________________________________________
#### GET DATA ####

## Set working directory
setwd("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/R")
getwd()

## Load packages
library(pool)
library(DBI)
library(RPostgres)
library(dplyr)
library(ggplot2)

## Create connection to OneBenthic
Sys.setenv(R_CONFIG_ACTIVE = "one_benthic")

dw <- config::get()

pool <- dbPool(drv = dbDriver(dw$driver),
               dbname = dw$database,
               host = dw$server,
               port =  dw$port,
               user = dw$uid,
               password = dw$pwd)

## SQL select query.
data = dbGetQuery(pool,"SELECT
su.surveyname,
s.samplecode,
s.samplelat,
s.samplelong,
w.validname,
w.validaphiaid,
w.family,
w.genus,
ts.abund,
s.date,
s.year,
s.month,
s.gear_gearcode,
su.datapubliclyavailable,
w.rank

FROM 
associations.survey as su
INNER JOIN associations.surveysample as ss ON ss.survey_surveyname = su.surveyname 
INNER JOIN samples.sample as s ON ss.sample_samplecode = s.samplecode
INNER JOIN faunal_data.taxasample as ts ON s.samplecode= ts.sample_samplecode 
--LEFT JOIN faunal_data.taxaqual as tq ON ts.taxaqual_qualifier = tq.qualifier 
INNER JOIN faunal_data.worrms as w ON w.aphiaid = ts.worrms_aphiaid 
INNER JOIN associations.sampleowner as so ON so.sample_samplecode = s.samplecode
INNER JOIN associations.owner as o ON so.owner_ownername = o.ownername

WHERE (s.gear_gearcode = 'MHN' OR
s.gear_gearcode = 'DG' OR
s.gear_gearcode = 'VV' OR
s.gear_gearcode = 'SM' OR
s.gear_gearcode = 'NIOZ' OR
s.gear_gearcode = 'BC_0.1' OR
s.gear_gearcode = 'C/VV'OR
s.gear_gearcode = 'BC' OR
s.gear_gearcode = 'DVV' OR
s.gear_gearcode = 'DG/VV')
AND (ts.taxaqual_qualifier NOT IN ('J','E','EP', 'L', 'MEGA', 'PR', 'PU', 'TAIL', 'Z', 'SP', 'FRAG','FUR', 'DEAD') OR ts.taxaqual_qualifier IS NULL)
AND (s.treatment = 'R' or s.treatment IS NULL)
AND s.macrosieve= 1
AND w.include = TRUE
AND s.samplelat > 47.92938
AND ts.abund IS NOT NULL
AND s.id <= 52951
AND s.year NOT IN (1976,1989)
AND su.surveyname NOT IN ('Long-Term Monitoring Program Muddy-sandy intertidal flats in Kandalaksha Bay (White Sea)')
--AND su.datapubliclyavailable = TRUE
ORDER by su.surveyname, s.samplecode,ts.abund desc;")

## Save file
#write.csv(data, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\data.csv", row.names=FALSE)
#_______________________________________________________________________________
#### PREPARE DATA: GENERATE COUNTS (ALL TAXA) ####

## Inspect data
head(data)

## Total count by sample
counts = data %>% group_by(samplecode,samplelat, samplelong,year, month, date) %>% summarise(count = sum(abund))

## Change column names
colnames(counts)[2] <- 'latitude'
colnames(counts)[3] <- 'longitude'
counts
dim(counts) #37909     7

## Plot samples by year
p= ggplot()+
  geom_point(data=counts,aes(longitude,latitude,col="blue"), size=0.15,show.legend = FALSE)+
  coord_map(xlim = c(-10.7, 4),ylim = c(48, 62)) +#set x,y limits of plot 
  theme_bw(base_size=24)+ 
   guides(colour = guide_legend(override.aes = list(size=5)))+ # Change size of legend dots #(were too small) 
  labs(x="Longitude",y="Latitude")+ 
  facet_wrap(~year) 
p

## Check the number of samples
length(unique(data$samplecode))# 37909

## Number of records
dim(data)# 1194234      15

## Taxonomic identification level summary
data %>%
  group_by(rank) %>%
  summarize(count=n())%>%
  mutate(per= prop.table(count) * 100)%>%
  arrange(desc(per))

## Subset for only Species, Genus and Family
data_ss <- data[which(data$rank=='Species'|data$rank=='Genus'|data$rank=='Family'),]
dim(data_ss)# 1133535      15

## How much of the data are we dropping? Answer ~5%
dim(data_ss)/dim(data)*100# 94.9%

## Re-establish 'data' name for df 'data_ss'
data <- data_ss
dim(data)# 1133535      15
#_______________________________________________________________________________
#### PREPARE DATA: ADDRESS MISSING DATE INFO (REQUIRED FOR DIVERSITY METRIC CALCS) ####

## DF for rows without NA (dates etc)
data_noNA <- data[complete.cases(data$date), ]
dim(data_noNA)# 1080003      15

## DF for rows with NA (dates etc)
data_withNA <-data[!complete.cases(data$date), ]
dim(data_withNA)# 53532    15

## Make sure you have all rows
tot <- nrow(data_noNA)+nrow(data_withNA)
tot# 1133535

## Make month '6' where value is NA
data_withNA$month[is.na(data_withNA$month)] <- 6

## Check month present in all rows
head(data_withNA$month[is.na(data_withNA$month)])# it is

## Check all rows have a year
head(data_withNA$year[is.na(data_withNA$year)]) # they do

## Create missing dates with nominal day of 30
data_withNA$date=paste(data_withNA$month,30, data_withNA$year, sep="/")

## Get date into correct format
data_withNA$date  <- as.Date(data_withNA$date , format = "%m/%d/%Y")

## Bind together the two dfs
data2 <- rbind(data_noNA,data_withNA)

## Check dimesions
dim(data2)# 1133535      15
head(data2)

## Re-establish df name
data <- data2

## Save file
write.csv(data, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\MURRAY\\counts_fam_genus_speciesv2.csv", row.names=FALSE)
#_______________________________________________________________________________
#### SAMPLE LOCATIONS (FIGURE 1)  ####

## Load packages
library(sf)
library(rasterVis)# to use raster in ggplot
library(raster)
library(ggnewscale)
library(scales)
library(ggplot2)
library(dplyr)
library(shadowtext)

## Load countries polygon and Northern Ireland border
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))
ni_border <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\ni_border.shp"))

## Load DEM
dem <- raster("C:\\Users\\kmc00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\GEBCO_15_Nov_2023_be52dc9c9b2d\\gebco_2023_n61.0_s48.0_w-11.0_e11.0.tif")

## Reduce size of DEM (optional)
#dem2 <- aggregate(dem, fact=5)
#dem <- dem2

## Remove elevations above sea level.
dem[dem>0] <- 0

## Create slope and hillshade
slope = terrain(dem, opt='slope')
aspect = terrain(dem, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)

dem_spdf <- as(dem, "SpatialPixelsDataFrame")
dem_spdf <- as.data.frame(dem_spdf)
colnames(dem_spdf) <- c("value", "x", "y")

hill_spdf <- as(hill, "SpatialPixelsDataFrame")
hill_spdf <- as.data.frame(hill_spdf)
colnames(hill_spdf) <- c("value", "x", "y")

## Get unique position coordinates
unique_pos <- data %>% dplyr::select(samplelong, samplelat)%>%unique

## Convert latitude and longitude into geometries using st_as_sf(). 
points <- unique_pos %>%
  st_as_sf(coords = c("samplelong", "samplelat"), crs = 4326)

## st_coordinates() extracts the lon/lat values as a data frame with X and Y columns so you can use geom_point (ability to resize points etc)
points2 <- st_coordinates(points)
points2 <- as.data.frame(st_coordinates(points))

## Create plot (map of sample locations, place names and background bathymetry)
# For bathy colours and breakpoints: https://stackoverflow.com/questions/70739780/is-there-a-scale-function-with-which-i-can-use-4-breaks-points

PSam2=ggplot()+ 
  geom_raster(data = hill_spdf, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_raster(data = dem_spdf, aes(x = x, y = y, fill = value), alpha=0.7)+
  scale_fill_gradientn(
    colours = c("#051650","#02367b","#006ca5","#0496c7", "#04bade","#55e2e9"),
    limits  = c(-3800,0),
    values  = scales::rescale(c(-3800,-1000,-100,-50,-25,0), from = c(-3800,0)))+
  geom_point(data = points2,aes(x = X, y = Y), fill="yellow",alpha = 0.7,colour="yellow",size=0.35)+
  geom_sf(data=countries, fill ="black",col ="black")+ 
  geom_sf(data=ni_border, ,col ="grey",linewidth = 0.2)+ 
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw()+#text = element_text(size=40)
  xlab("Longitude") +
  ylab("Latitude")+
  annotate("text",x=c(-1.3),y=c(52.5),label=c("UNITED \nKINGDOM"),color="white", size=5)+#3
  annotate("text",x=c(-6.7),y=c(51.2),label=c("Celtic \nSea"),color="white", size=7)+
  annotate("text",x=c(-1.1),y=c(50.15),label=c("English Channel"),color="white", size=7)+
  annotate("text",x=c(-4.7),y=c(53.8),label=c("Irish Sea"),color="white", size=6)+
  annotate("text",x=c(3),y=c(56.5),label=c("North Sea"),color="white", size=7)+
  annotate("text",x=c(6.3),y=c(54.5),label=c("German Bight"),color="white", size=5)+
  annotate("text",x=c(4.3),y=c(54.3),label=c("Oyster Ground"),color="white", size=5)+
  annotate("text",x=c(2.9),y=c(51.9),label=c("Southern \nBight"),color="white", size=5)+
  annotate("text",x=c(0),y=c(58.5),label=c("Fladen Ground"),color="white", size=5)+
  annotate("text",x=c(2.33),y=c(54.9),label=c("Dogger Bank"),color="white", size=5)+
  annotate("text",x=c(5.6),y=c(58.3),label=c("Norwegian Trench"),color="white", size=5, angle = -40)+
  annotate("text",x=c(8.8),y=c(58),label=c("Skagerrak"),color="white", size=5, angle = 40)+
  annotate("text",x=c(-1),y=c(56.5),label=c("Scalp \nBank"),color="white", size=5)+
  annotate("text",x=c(3),y=c(58),label=c("Ling Bank"),color="white", size=5)+
  annotate("text",x=c(2.5),y=c(59),label=c("Utsira \nHigh"),color="white", size=5)+
  annotate("text",x=c(5.1),y=c(57),label=c("Fisher \nBanks"),color="white", size=5)+
  annotate("text",x=c(7.5),y=c(57),label=c("Jutland \nBank"),color="white", size=5)+
  annotate("text",x=c(-5.8),y=c(52),label=c("St George's \nChannel"),color="white", size=4.5)+
  annotate("text",x=c(-7.6),y=c(53.2),label=c("IRELAND"),color="white", size=5)+
  annotate("text",x=c(1),y=c(49),label=c("FRANCE"),color="white", size=5)+
  annotate("text",x=c(6.5),y=c(52.5),label=c("NETHERLANDS"),color="white", size=5)+
  annotate("text",x=c(9.05),y=c(56),label=c("DENMARK"),color="white", size=5)+
  annotate("text",x=c(3.9),y=c(51),label=c("BELGIUM"),color="white", size=5)+
  annotate("text",x=c(7.6),y=c(59),label=c("NORWAY"),color="white", size=5)+
  annotate("text",x=c(9),y=c(53),label=c("GERMANY"),color="white", size=5)+
  annotate("text",x=c(-4.5),y=c(54.3),label=c("Isle \nof \nMan"),color="white", size=4)+
  annotate("text",x=c(-5.5),y=c(54.97),label=c("North\nChannel"),color="white", size=3)+
  annotate("text",x=c(-6.7),y=c(56.9),label=c("Hebrides"),color="white", size=5)+
  annotate("text",x=c(0.32),y=c(56.5),label=c("Devil's \nHole"),color="white", size=5)+
  annotate("text",x=c(0.1666),y=c(55.883),label=c("Swallow \nHole"),color="white", size=5)+
  annotate("text",x=c(-2.4),y=c(56.12),label=c("Firth of Forth"),color="white", size=4)+
  annotate("text",x=c(-9),y=c(50),label=c("SW \nApproaches"),color="white", size=6)+
  shadowtext::geom_shadowtext(aes(label = 'Inner \nSilver \nPit'),x=1,y=53.5, size=5,bg.color="#0496c7", bg.r=0.1, color = "white")+#hjust=1,
  shadowtext::geom_shadowtext(aes(label = 'Bristol Channel'),x=-4.8,y=51.4, size=5,bg.color="#0496c7", bg.r=0.1, color = "white")+#hjust=1,
  shadowtext::geom_shadowtext(aes(label = 'Outher \nThames'),x=1.35,y=51.57,angle = 0, size=4,bg.color="#0496c7", bg.r=0.1, color = "white")+#hjust=1,
  shadowtext::geom_shadowtext(aes(label = 'Strait of Dover'),x=1.5,y=51,angle = 45, size=5,bg.color="#0496c7", bg.r=0.1, color = "white")+#hjust=1,
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 20),legend.title = element_text(color = "white", size = 20))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.16))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.key.size = unit(1.5, "cm"))+
  theme(text = element_text(size = 22))+
  labs(fill = "Bathymetry (m)\n")

## Save Figure 1
ggsave(plot = PSam2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_1.png"),
       width = 41,height = 40,units = "cm", pointsize = 48,
       device = "png",limitsize = FALSE,bg="white")
#_______________________________________________________________________________
#### CALCULATE BIODIV METRICS: SCRIPT SETUP #####

## Specify required packages
pkgs = c("geosphere", "lubridate", "sf", "tidyverse", 
         "mapplots", "purrr", "iNEXT", "mapplots", 
         "RColorBrewer", "ggpubr", "parallel", 
         "janitor")

## Loading specified packages
invisible(lapply(pkgs, library, character.only = TRUE))

## Path to store results
outpath = 'C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/OUTPUTS/diversity_estimates/'

## Converting dataframe to tibble and rename
sppdf <- as_tibble(data)

## Set parameters for diversity assessment
n_samp = 6# number of samples to use for gamma estimation
sub_sample = 'Y'# whether to restrict to n_samp for gamma
n_days = 182# the number of days between samples (i.e. use samples up to 6 months either side for gamma estimation) 
radius = 75000# spatial radius in meters used to estimate gamma 
years=sort(unique(sppdf$year))# years to be assessed
no_cores <- detectCores() - 4# Calculate the number of cores for parallel processing
file_name='one_benthic_10_24'# file name for saving

## Load functions
source("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/R/diversity_functions_MT.R")

# Running the data tidying function for desired years 
# Adjust firstyear and lastyear variables if required
tidy_sppdf <- tidy_benthic_dat(dat = sppdf, firstyear = min(years), lastyear = max(years)) #firstyear = 1981, lastyear = 2020
#_______________________________________________________________________________
#### CALCULATE BIODIV METRICS: DO CALCS AND SAVE ####

# Running the parallel_process function using the tidy benthic data
# Aggregates all the .RData files for processing before providing the .RData output 
all_df <- parallel_process()

#_______________________________________________________________________________
#### CALCULATE BIODIV METRICS: LOAD RESULTS ####

## Load results 
load(paste0(outpath, 'biodiversity_estimates_1985_2023_one_benthic_10_24.Rdata'))
#_______________________________________________________________________________
#### BIODIV METRICS: MAP RESULTS (QUICK LOOK) ####

## Load packages
library(mapplots)
library(RColorBrewer)
library(ggpubr)
library(dplyr)

## World map
world_shp = sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))

## Average gamma at ices rectangle level just for look see 
av_spatial_div = all_df %>% 
  mutate(log10_av_count = log10(av_count),
         log10_tot_count = log10(tot_count),) %>%
  dplyr::select(sample, latitude, longitude, year, month, date,
                log10_av_count, cv_count, log10_tot_count, 
                sample_a_q0, sample_a_q1, sample_a_q2,
                #a_q0, a_q1, a_q2,
                b_q0, b_q1, b_q2,
                g_q0, g_q1, g_q2) %>%
  pivot_longer(cols= -c(sample:date),
               names_to = "metric",
               values_to = "value") %>%
  mutate(ices = ices.rect2(longitude, latitude),
         metric = factor(metric, levels = c('log10_av_count', 'cv_count', 'log10_tot_count', 
                                            'sample_a_q0', 'sample_a_q1', 'sample_a_q2',
                                            'b_q0', 'b_q1', 'b_q2',
                                            'g_q0', 'g_q1', 'g_q2'))) %>%
  group_by(ices, metric) %>%
  summarise(av_metric = mean(value, na.rm = T))
coords = ices.rect(av_spatial_div$ices)
av_spatial_div = cbind(av_spatial_div, coords)

plots_sp_div = list()
for(met in levels(av_spatial_div$metric)) { # met=levels(av_spatial_div$metric)[1] 
  
  p = av_spatial_div %>%
    filter(metric == met,
           !is.na(av_metric),
           !is.na(lat)) %>%
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = av_metric)) +
    scale_fill_gradientn(colours = c('white', brewer.pal(9,"YlOrRd")), 
                         name = '') +
    labs(x='Longitude', y='Latitude', title=met) +
    guides(fill = guide_colourbar(barwidth = 1)) +
    geom_sf(data = world_shp, 
            fill = 'black', 
            color = 'black',
            size = 0.1) +
    coord_sf(xlim = c(-16, 25), ylim = c(36, 65))+
    theme(panel.background = element_rect(fill = 'grey80'),
          panel.border = element_rect(colour='black', fill=NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  plots_sp_div[[met]] = p
}

sp_div_plts = do.call(ggarrange, c(plots_sp_div, align='hv', ncol=3, nrow=4)) %>%
  annotate_figure(fig.lab = 'Spatial diversity', fig.lab.size = 14) 

# Final plot for all metrics
sp_div_plts
#_______________________________________________________________________________
#### BIODIV METRICS: ADD COUNTS DATA AND SAVE ####

## Add total count information. Sub-setting (outlier removal) code not used as issue addressed later in code 
# justification for subsetting:
# alpha cutoff close to alpha_n (i.e. where all individuals are singletons), yielding unreliable estimates (these could be plot in supporting material?)
# beta and gamma estimate only used where n_samp == 6
# beta diversity >14 would represent >100% turnover, hence removed
# hill nos should be 0>1>2 

final_df = counts %>%
  rename(sample=samplecode) %>%
  left_join(all_df, by =c('sample', 'latitude', 'longitude', 'year', 'month', 'date') ) %>%
  mutate(b_q0 = case_when(b_q0 <= 14 & n_samp == 6 ~  b_q0, TRUE ~ NA),
         b_q1 = case_when(b_q1 <= 14 & n_samp == 6 ~  b_q1, TRUE ~ NA),
         b_q2 = case_when(b_q2 <= 14 & n_samp == 6 ~  b_q2, TRUE ~ NA),
         sample_a_q0 = case_when(sample_a_q0 < alpha_n-1 ~ sample_a_q0, TRUE ~ NA),
         sample_a_q1 = case_when(sample_a_q1 < alpha_n-1 & sample_a_q1 < sample_a_q0 ~ sample_a_q1,
                                 TRUE ~ NA),
         sample_a_q2 = case_when(sample_a_q2 < alpha_n-1 &
                                   sample_a_q2 < sample_a_q1 ~ sample_a_q2,
                                 TRUE ~ NA),
         g_q0 = case_when(is.na(sample_a_q0) | is.na(b_q0) | n_samp <6 ~ NA, TRUE ~ g_q0),
         g_q1 = case_when(is.na(sample_a_q1) | is.na(b_q1) | n_samp <6 ~ NA, TRUE ~ g_q1),
         g_q2 = case_when(is.na(sample_a_q2) | is.na(b_q2) | n_samp <6 ~ NA, TRUE ~ g_q2),
         assemblage = file_name)
 
## Check max values
max(final_df$sample_a_q0, na.rm=T)
colnames(final_df)

## Pivot data
library(tidyr)
final_df_long = final_df %>% 
  ungroup()%>%
  dplyr::select(longitude, latitude, sample,
                count, cv_count, tot_count, 
                #av_count, cv_count, tot_count, 
                sample_a_q0, sample_a_q1, sample_a_q2,
                b_q0, b_q1, b_q2,
                g_q0, g_q1, g_q2) %>%
  tidyr::pivot_longer(cols= -c(longitude:sample),
               names_to = "metric",
               values_to = "measurement") %>%
  rename(X=longitude, Y=latitude, Sample=sample)

## Check data
head(final_df_long)
dim(final_df_long)# 454908      5

## Find max values by metric
final_df_long %>%
  group_by(metric) %>%
  summarise(max = max(measurement, na.rm=TRUE))

## Examine data distribution for each metric
library(ggplot2)
ggplot(data=final_df_long,aes(x=measurement))+
  geom_histogram()+
  facet_wrap(~metric, scales='free')

## Save data
write.csv(final_df_long,"C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/DATA/biodiv_data.csv", row.names = FALSE)
#_______________________________________________________________________________
#### BIODIV METRICS: LOAD BIODIVERSITY/COUNTS DATA ####

## Load biodiversity metric data in long format (cols: 'X','Y','Sample','metric','measurement')
biodiv8 <- read.csv("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/DATA/biodiv_data.csv",header=TRUE, stringsAsFactors=FALSE)#

## Make 'metric' a factor
biodiv8$metric <- as.factor(biodiv8$metric)

## Inspect data
head(biodiv8)
#View(biodiv8)

## Create a tally of all values in the 'metric' column
biodiv8 %>%
  filter(!is.na(measurement)) %>%
  count(metric)

#        metric     n
#1         b_q0 35348
#2         b_q1 32195
#3         b_q2 28094
#4        count 37909
#5     cv_count 35519
#6         g_q0 35348
#7         g_q1 29277
#8         g_q2 25770
#9  sample_a_q0 35519
#10 sample_a_q1 32193
#11 sample_a_q2 32178
#12   tot_count 35519

#_______________________________________________________________________________
#### PREPARE BIODIV METRICS DATA: REMOVE 'REPLICATES' TO ADDRESS SPATIAL AUTOCORRELATION ####

## Load packages
library(tidyr)
library(sp)

## Long to wide format
biodiv8_wide <- biodiv8 %>%
pivot_wider(names_from = metric, values_from = measurement)
head(biodiv8_wide)

## Set coordinates
coordinates(biodiv8_wide) <- c("X", "Y")

## Work out 50m distance in decimal degrees. 1 degree of latitude =111,000m
50/111000# degrees for 50m #0.0004504505

## Set distance within which to remove replicates (zero)
zd <- zerodist(biodiv8_wide,zero = 0.0004504505)

## Drop replicates
biodiv8_wide_norep <- biodiv8_wide[-zd[,2], ]
dim(biodiv8_wide_norep)# 22793   13
class(biodiv8_wide_norep)

## Change class to df
biodiv8_wide_norep_2=data.frame(biodiv8_wide_norep)
class(biodiv8_wide_norep_2)# df
names(biodiv8_wide_norep_2)

## Drop col 'optional'
biodiv8_wide_norep_3=biodiv8_wide_norep_2[,1:(ncol(biodiv8_wide_norep_2)-1)]
head(biodiv8_wide_norep_3)
max(biodiv8_wide_norep_3$sample_a_q2,na.rm=T)

## Turn back to long format
biodiv8 <- gather(biodiv8_wide_norep_3, metric, measurement, count:g_q2, factor_key=TRUE)
head(biodiv8)
dim(biodiv8)# 273516      5
#_______________________________________________________________________________
#### PREPARE BIODIV METRICS DATA: REMOVE OUTLIERS ####

## Load package
library(dplyr)

## Calculate medians for each biodiversity metric
biodiv9 = biodiv8
str(biodiv9)
biodiv9$metric <- as.character(biodiv9$metric)
biodiv9 %>%
  na.omit() %>%
  group_by(metric)%>% 
  summarise(median=median(measurement))

## Add medians (see results from previous step)
biodiv9 <- biodiv9 %>%
    mutate(median = case_when(
      metric == 'sample_a_q0' ~ 33.8,
      metric == 'sample_a_q1' ~ 17.7,
      metric == 'sample_a_q2' ~ 9.67,
      metric == 'b_q0' ~ 4.88,
      metric == 'b_q1' ~ 6.62,
      metric == 'b_q2' ~ 6.89,
      metric == 'g_q0' ~ 168,
      metric == 'g_q1' ~ 136,
      metric == 'g_q2' ~ 98.6,
      metric == 'tot_count' ~ 739,
      metric == 'count' ~ 70,
      metric == 'cv_count' ~0.858
      ))

## Find absolute difference between each value and the median
biodiv9$abs_diff_med <-abs( biodiv9$median -biodiv9$measurement)

## Find the Median Absolute Deviation (MAD)
biodiv9 %>%
  na.omit() %>%
  group_by(metric)%>% 
  summarise(mad=median(abs_diff_med))

## Add MAD column
biodiv9 <- biodiv9 %>%
    mutate(mad = case_when(
      metric == 'sample_a_q0' ~ 14.1,
      metric == 'sample_a_q1' ~ 8.94,
      metric == 'sample_a_q2' ~ 5.16,
      metric == 'b_q0' ~ 1.25,
      metric == 'b_q1' ~ 2.01,
      metric == 'b_q2' ~ 2.41,
      metric == 'g_q0' ~ 53.6,
      metric == 'g_q1' ~ 42.9,
      metric == 'g_q2' ~ 32.6,
      metric == 'tot_count' ~ 371,
      metric == 'count' ~ 51,
      metric == 'cv_count' ~ 0.30
      ))
head(biodiv9)

## Find Modified Z-Score for each data value
biodiv9$mod_z_score = 0.6745*(biodiv9$abs_diff_med) / biodiv9$mad
head(biodiv9)

## Remove values with a modified Z-score less than -3.5 or greater than 3.5
biodiv9 = biodiv9 %>%
  na.omit() %>%
  group_by(metric) %>%
  filter(mod_z_score > -3.5 & mod_z_score < 3.5)# %>%
head(biodiv9)
#View(biodiv9)

## Remove unwanted cols
biodiv9 <- biodiv9[,1:5]

## Create df of biodiv metrics (long format) for modelling (outliers excluded)
biodiv9_mod <-as.data.frame(biodiv9)
head(biodiv9_mod)

## Examine data distribution for each metric
library(ggplot2)
ggplot(data=biodiv9_mod,aes(x=measurement))+
  geom_histogram()+
  facet_wrap(~metric, scales='free')
#_______________________________________________________________________________
## DO THIS IF NOT REMOVING OUTLIERS IN STEP ABOVE (ie FOLLOWING APPROACH IN SECTION 'BIODIV METRICS: ADD COUNTS DATA AND SAVE')
#biodiv9 <- biodiv8
#biodiv9_mod <- biodiv9
#_______________________________________________________________________________
#### SAVE BIODIV METRICS DATA FOR MODELLING ####

## Reorder columns
biodiv9_mod2 <- biodiv9_mod[,c(3,1,2,4,5)]

## Change type
biodiv9_mod2$type <- 'numeric'

## Add paper(to help manage modelling tasks)
biodiv9_mod2$paper <- 'biodiversity'

## Update column names
colnames(biodiv9_mod2) <- c('sample','x','y','metric','value','type','paper')

## Reorder columns
biodiv9_mod2 <- biodiv9_mod2[,c(7,1:4,6,5)]
head(biodiv9_mod2)

## How many records do we have for each metric?
biodiv9_mod2 %>%
  filter(!is.na(value)) %>%
  count(metric)

#        metric     n
#1         b_q0 21562
#2         b_q1 19526
#3         b_q2 17171
#4        count 20508
#5     cv_count 21679
#6         g_q0 21754
#7         g_q1 17648
#8         g_q2 15616
#9  sample_a_q0 21829
#10 sample_a_q1 19469
#11 sample_a_q2 19288
#12   tot_count 20322

dim(biodiv9_mod2)# 236372      7

## Save metric data for use with RF modelling script (see xx.R)
write.csv(biodiv9_mod2, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\biodiv_metrics_4_modelling.csv", row.names=FALSE)
#_______________________________________________________________________________
#### PREPARE BIODIV METRICS DATA FOR CLUSTERING: REMOVE ROWS WITH MISSING DATA ####

## Load package
library(tidyr)

## Change data into wide format
biodiv9_wide <- spread(biodiv9, metric, measurement)
head(biodiv9_wide)
dim(biodiv9_wide)#22661    15

## Remove rows with NA
biodiv9_wide <- na.omit(biodiv9_wide)
head(biodiv9_wide)
dim(biodiv9_wide)# 13654    15
#_______________________________________________________________________________
#### PREPARE BIODIV METRICS DATA FOR CLUSTERING: ASSESS SKEWNESS AND TRANSFORM AS NECESSARY ####

## Load packages
library(moments)
library(dplyr)
library(ggplot2)

## Return skewness values
skew <- biodiv9_wide[,4:15] %>% 
  dplyr::mutate(across(1:12, moments::skewness))%>% 
  distinct()
skew

## Histograms
hist <- ggplot(biodiv9, aes(measurement)) +
  geom_histogram() +
  facet_wrap(~metric, scales='free')
hist

# Note all metrics positively skewed (i.e. values >0)
# Symmetric: Values between -0.5 to 0.5
# Moderated Skewed data: Values between -1 and -0.5 or between 0.5 and 1
# Highly Skewed data: Values less than -1 or greater than 1

## Apply appropriate transformations
biodiv9_wide$sample_a_q0_none=biodiv9_wide$sample_a_q0#0.227
biodiv9_wide$sample_a_q1_sqrt=sqrt(biodiv9_wide$sample_a_q1)#0.627
biodiv9_wide$sample_a_q2_sqrt=sqrt(biodiv9_wide$sample_a_q2)#0.778
biodiv9_wide$b_q0_sqrt=sqrt(biodiv9_wide$b_q0)#0.570
biodiv9_wide$b_q1_none=biodiv9_wide$b_q1#0.465
biodiv9_wide$b_q2_none=biodiv9_wide$b_q2#0.160
biodiv9_wide$g_q0_sqrt=sqrt(biodiv9_wide$g_q0)#0.521
biodiv9_wide$g_q1_sqrt=sqrt(biodiv9_wide$g_q1)#0.501
biodiv9_wide$g_q2_sqrt=sqrt(biodiv9_wide$g_q2)#0.530
biodiv9_wide$tot_count_log=log(biodiv9_wide$tot_count)#1.42
biodiv9_wide$count_log=log(biodiv9_wide$count)#1.22
biodiv9_wide$cv_count_log=log(biodiv9_wide$cv_count)#1.07

## Inspect data
head(biodiv9_wide)
dim(biodiv9_wide)# 13654    27

## Take only transformed data cols
names(biodiv9_wide)
biodiv9_trans <- biodiv9_wide[,c(16:24,26,27,25)]
head(biodiv9_trans)

## Convert to long format
biodiv9_trans_long <- gather(biodiv9_trans, metric, measurement, sample_a_q0_none:tot_count_log, factor_key=TRUE)
head(biodiv9_trans_long)

## Check histograms following transformation
hist2 <- ggplot(biodiv9_trans_long, aes(measurement)) +
  geom_histogram() +
  facet_wrap(~metric, scales='free')
hist2

## Check skewness values following transformation
biodiv9_trans2 <- biodiv9_trans %>% 
  mutate(across(1:12, moments::skewness))%>% 
  distinct()
head(biodiv9_trans2) # Skewness following transformation (now all =<0.5 )
#_______________________________________________________________________________
#### PREPARE BIODIV METRICS DATA FOR CLUSTERING: NORMALISATION ####

## Load  packages
library(scales)

## Create a df for scaling
biodiv9_trans_scale <- biodiv9_trans

## Scale Values Between 0 and 1 Using scales Package
biodiv9_trans_scale$sample_a_q0_none <- scales::rescale(biodiv9_trans_scale$sample_a_q0_none)
biodiv9_trans_scale$sample_a_q1_sqrt <- scales::rescale(biodiv9_trans_scale$sample_a_q1_sqrt)
biodiv9_trans_scale$sample_a_q2_sqrt <- scales::rescale(biodiv9_trans_scale$sample_a_q2_sqrt)
biodiv9_trans_scale$b_q0_sqrt <- scales::rescale(biodiv9_trans_scale$b_q0_sqrt)
biodiv9_trans_scale$b_q1_none <- scales::rescale(biodiv9_trans_scale$b_q1_none)
biodiv9_trans_scale$b_q2_none <- scales::rescale(biodiv9_trans_scale$b_q2_none)
biodiv9_trans_scale$g_q0_sqrt <- scales::rescale(biodiv9_trans_scale$g_q0_sqrt)
biodiv9_trans_scale$g_q1_sqrt <- scales::rescale(biodiv9_trans_scale$g_q1_sqrt)
biodiv9_trans_scale$g_q2_sqrt <- scales::rescale(biodiv9_trans_scale$g_q2_sqrt)
biodiv9_trans_scale$tot_count_log <- scales::rescale(biodiv9_trans_scale$tot_count_log)
biodiv9_trans_scale$count_log <- scales::rescale(biodiv9_trans_scale$count_log)
biodiv9_trans_scale$cv_count_log <- scales::rescale(biodiv9_trans_scale$cv_count_log)

## Check scaling results
head(biodiv9_trans_scale)
#_______________________________________________________________________________
#### PREPARE BIODIV METRICS DATA FOR CLUSTERING: IDENTIFY COVARIATES ####

## Load packages
library(ggplot2)       
library(GGally)    

## Apply ggpairs function 
p <- ggpairs(as.data.frame(biodiv9_trans_scale), columnLabels = c(
"''^0*D[alpha]",
"''^1*D[alpha]",
"''^2*D[alpha]", 
"''^0*D[beta]",
"''^1*D[beta]", 
"''^2*D[beta]",
"''^0*D[gamma]", 
"''^1*D[gamma]",
"''^2*D[gamma]",
"N",
"N[cv]",
"N[tot]"
), labeller = 
          label_parsed)

figure_s1 <- p + theme(text=element_text(size=14))+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Save pairs plot
ggsave(plot = figure_s1,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_S1.png"),
       height = 300, width =300, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285

## Following variables highly correlated:
#1Dα	vs	0Dα	0.892
#2Dα	vs	0Dα	0.731
#2Dα	vs	1Dα	0.945
#1Dβ	vs	0Dβ	0.906
#2Dβ	vs	1Dβ	0.723
#1Dγ	vs	0Dγ	0.987
#2Dγ	vs	0Dγ	0.964
#2Dγ	vs	1Dγ	0.992

## Therefore drop 1Dα, 2Dα, 1Dβ, 1Dγ, 2Dγ
## ie keep 0Dα, 0Dβ, 2Dβ, 0Dγ, N, Nav, Ntot

## Convert to long format
clus_data <- as.data.frame(biodiv9_trans_scale[,c(1,4,6,7,10,11,12)])
head(clus_data)

## Bring in coordinates and samplecode
clus_data1 <- cbind(biodiv9_wide[,1:3],clus_data)
head(clus_data1)
dim(clus_data1)#  13654    10

## Save data (ready for clustering)
#write.csv(clus_data1, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\data_4_clustering.csv", row.names=FALSE)

## Load data
#data_4_clustering <- read.csv("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\data_4_clustering.csv",header=TRUE, stringsAsFactors=FALSE)
#head(data_4_clustering)
#_______________________________________________________________________________
#### BIODIVERSITY CLUSTERS: ELBOW PLOT (FIGURE 3A) ####

## Subset clustering metrics only
data_4_clustering2 <- clus_data1[,4:10]

## Change class of df data_4_clustering2 to a matrix 
clus4=data.matrix(data_4_clustering2)

## Load package
library(factoextra)

## Generate plot (clus4 comes from below step)
elbow <- fviz_nbclust(clus4,kmeans, method = "wss",linecolor = "black",k.max = 20)+
  geom_vline(xintercept = 8, linetype = 2)+theme_classic(base_size = 16)
plot(elbow)
#_______________________________________________________________________________
#### BIODIVERSITY CLUSTERS: DO CLUSTERING ####

## Perform K-means clustering of data. Results (cluster group) to the object 'results' 
set.seed(1234) 
results=kmeans(clus4,8,algorithm="MacQueen",iter.max=100,nstart=25) 

## Number of samples belonging to each cluster group 
results$size # 1322 1948 1245 2096 1346 2569 1461 1667
#_______________________________________________________________________________
#### BIODIVERSITY CLUSTERS: DENDROGRAM (FIGURE 3B) ####

## Load packages
library(ggplot2)
library(ggdendro)
library(ggplot2)
library(dplyr)
library(dendextend)

## Function to calculate absolute differences between cluster centres over all variables.
nclusters = 8
absdiff = matrix(0, nrow=nclusters, ncol=nclusters)
centers = results$centers
for (j in 1:nclusters) {
  for (k in 1:nclusters) {
    absdiff[j,k] = sum(abs(centers[j,] - centers[k,]))
  }
}
d=round(absdiff, 1)

## Find distance matrix
d1 <- dist(as.matrix(d))

## Produce dendrogram
test2 <- d1%>% hclust %>% as.dendrogram %>%set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 9) %>%  # node point size
  set("leaves_col", c( 
"#C0C1BC",#3
"#6FD326",#8
"#37C331",#2
"#A8E21B",#5
"#F3F223",#4
"#E2E256",#6
"#E0F210",#1
"#D1D189"#7
)) %>%
  set("labels_cex",0.9)%>% 
  #set("labels", c('     3', '     8', '     2', '     5', '     4', '     6', '     1', '     7'))%>% 
  set("labels", c('     Bio-H', '     Bio-B', '     Bio-A', '     Bio-C', '     Bio-E', '     Bio-F', '     Bio-D', '     Bio-G'))%>%
  set("branches_lwd", 0.7)

## Change dendrogram into a ggplot
ggd1 <- as.ggdend(test2)
dendrogram <-  ggplot(ggd1, horiz = T)+theme_classic(base_size = 16)+theme(axis.title.y=element_blank(),
                                                               axis.text.y=element_blank(),
                                                               axis.ticks.y=element_blank(),
                                                              axis.line.y=element_blank())+labs(y='Height')
dendrogram
#_______________________________________________________________________________
#### BIODIVERSITY CLUSTERS: ELBOW & DENDROGRAM (FIGURE 3) ####

## Combined elbow plot and dendrogram plots
png("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_3.png", width = 33, height = 20, units = "cm", res = 500,pointsize = 12)
ggpubr::ggarrange(elbow,NULL,dendrogram,labels = c("a)","", "b)"),nrow=1,widths = c(1, 0.05, 1))
dev.off()
#_______________________________________________________________________________
#### BIODIVERSITY CLUSTERS: PREPARE CLUSTER RESULTS FOR PLOTTING ####

## Add cluster group from k-means results file to df 'clus_data1' which includes 'Sample', 'Latitude_WGS84' and 'Longitude_WGS84' 
faunal.cluster=cbind(clus_data1[,1:3],results$cluster)
head(faunal.cluster)

## Change name of col 'results$cluster' to 'FaunalCluster' 
names(faunal.cluster)[4]<-paste("FaunalCluster") 

## Make FaunalCluster a factor
faunal.cluster$FaunalCluster <- as.factor(faunal.cluster$FaunalCluster)
#str(faunal.cluster)

## Load packages
library(sf)
library(dplyr)

## Select coordinates
faunal.cluster2<- faunal.cluster %>% dplyr::select(X, Y)
head(faunal.cluster2)

## Use st_as_sf() to convert latitude and longitude into the magic geometry column. 
faunal.cluster3 <-faunal.cluster2 %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326)

## st_coordinates() extracts the lon/lat values as a data frame with X and Y columns so you can use geom_point (ability to resize points etc)
faunal.cluster4 <- st_coordinates(faunal.cluster3)

## Change to a dataframe
faunal.cluster5 <- as.data.frame(faunal.cluster4)
head(faunal.cluster5)

## Add Samplecode and FaunalCluster to coordinates
faunal.cluster6 <- cbind(faunal.cluster5,faunal.cluster[,c(3,4)])
head(faunal.cluster6)

## Change order of cols
faunal.cluster7 <- faunal.cluster6[,c(3,1,2,4)]
head(faunal.cluster7)

## Save results
#write.csv(faunal.cluster7, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\faunal.cluster6.csv", row.names=FALSE)
#_______________________________________________________________________________
#### BIODIVERSITY CLUSTERS: GROUP CHARACTERISTICS ####

## Load package
library(tidyr)

## Cluster centers
center <-results$centers
head(center)

## Create dataset with the cluster number
cluster <- c(1: 8)
center_df <- data.frame(cluster, center)
head(center_df)

## Reshape the data to long format
center_reshape <- gather(center_df, features, values, sample_a_q0_none: tot_count_log)#
head(center_reshape)

## Change data from long to wide format
zscores4gt <- center_reshape %>% spread(cluster, values, fill = NA, convert = FALSE)
head(zscores4gt)

## Bring together metrics and cluster result
metric.cluster <- base::merge(biodiv9_wide,faunal.cluster,by="Sample")
head(metric.cluster)
names(metric.cluster)

## Remove unwanted columns
metric.cluster2 <- metric.cluster[,c(30,4:15)]
head(metric.cluster2)

## Create col 'ClusterNum'
metric.cluster2$ClusterNum <- as.factor(metric.cluster2$FaunalCluster)
head(metric.cluster2)

## Output means by cluster group
metric.cluster3 <- metric.cluster2 %>%
  group_by(FaunalCluster) %>%
  dplyr::summarize(
    a_q0 = mean(sample_a_q0, na.rm=TRUE),
    a_q1 = mean(sample_a_q1, na.rm=TRUE),
    a_q2 = mean(sample_a_q2, na.rm=TRUE),
    b_q0 = mean(b_q0, na.rm=TRUE),
    b_q1 = mean(b_q1, na.rm=TRUE),
    b_q2 = mean(b_q2, na.rm=TRUE),
    g_q0 = mean(g_q0, na.rm=TRUE),
    g_q1 = mean(g_q1, na.rm=TRUE),
    g_q2 = mean(g_q2, na.rm=TRUE),
    count = mean(count, na.rm=TRUE),
    cv_count = mean(cv_count, na.rm=TRUE),
    tot_count = mean(tot_count, na.rm=TRUE))
head(metric.cluster3)

## Update column name and change to character
colnames(metric.cluster3)[1] <- "ClusterNum"
metric.cluster3$ClusterNum <- as.character(metric.cluster3$ClusterNum)
head(metric.cluster3)

## Output S.D. by cluster group
metric.sd.cluster3 <- metric.cluster2 %>%
  group_by(ClusterNum) %>%
  dplyr::summarize(
    a_q0 = sd(sample_a_q0, na.rm=TRUE),
    a_q1 = sd(sample_a_q1, na.rm=TRUE),
    a_q2 = sd(sample_a_q2, na.rm=TRUE),
    b_q0 = sd(b_q0, na.rm=TRUE),
    b_q1 = sd(b_q1, na.rm=TRUE),
    b_q2 = sd(b_q2, na.rm=TRUE),
    g_q0 = sd(g_q0, na.rm=TRUE),
    g_q1 = sd(g_q1, na.rm=TRUE),
    g_q2 = sd(g_q2, na.rm=TRUE),
    count = sd(count, na.rm=TRUE),
    cv_count = sd(cv_count, na.rm=TRUE),
    tot_count = sd(tot_count, na.rm=TRUE))
head(metric.sd.cluster3)

## Change from tibble to df
metric.sd.cluster3 <- as.data.frame(metric.sd.cluster3)

## Transpose so cols are clusters
zscores4gt2 <- as.data.frame(t(zscores4gt))
head(zscores4gt2)

## Make 1st row the column names
colnames(zscores4gt2) <- zscores4gt2[1,]
zscores4gt2 <- zscores4gt2[-1, ] 
head(zscores4gt2)

## Make row names the first column
zscores4gt3 <- cbind(names = rownames(zscores4gt2), zscores4gt2)
row.names(zscores4gt3) <- NULL
head(zscores4gt3)

## Update column names
colnames(zscores4gt3)[1] <- "ClusterNum"
colnames(zscores4gt3)[2] <- "b_q0"
colnames(zscores4gt3)[3] <- "b_q2"
colnames(zscores4gt3)[4] <- "count"
colnames(zscores4gt3)[5] <- "cv_count"
colnames(zscores4gt3)[6] <- "g_q0"
colnames(zscores4gt3)[7] <- "a_q0"
colnames(zscores4gt3)[8] <- "tot_count"

## Add in missing columns (metrics left out of the clustering for which centres to not exist)
zscores4gt3$a_q1 <- NA
zscores4gt3$a_q2 <- NA
zscores4gt3$b_q1 <- NA
zscores4gt3$g_q1 <- NA
zscores4gt3$g_q2 <- NA

## Update column order
names(zscores4gt3)
zscores4gt3 <- zscores4gt3[,c(1,7,9,10,2,11,3,6,12,13,4,5,8)]#a_q0,a_q1,a_q2,b_q0,b_q1,b_q2,count,cv_count,tot_count
head(zscores4gt3)

## Merge with zscore
metric.cluster6 <- rbind(zscores4gt3,metric.cluster3,metric.sd.cluster3)
head(metric.cluster6)

## Change values to numeric
str(metric.cluster6)
metric.cluster6$a_q0 <- as.numeric(metric.cluster6$a_q0)
metric.cluster6$a_q1 <- as.numeric(metric.cluster6$a_q1)
metric.cluster6$a_q2 <- as.numeric(metric.cluster6$a_q2)
metric.cluster6$b_q0 <- as.numeric(metric.cluster6$b_q0)
metric.cluster6$b_q1 <- as.numeric(metric.cluster6$b_q1)
metric.cluster6$b_q2 <- as.numeric(metric.cluster6$b_q2)
metric.cluster6$g_q0 <- as.numeric(metric.cluster6$g_q0)
metric.cluster6$g_q1 <- as.numeric(metric.cluster6$g_q1)
metric.cluster6$g_q2 <- as.numeric(metric.cluster6$g_q2)
metric.cluster6$count <- as.numeric(metric.cluster6$count)
metric.cluster6$cv_count <- as.numeric(metric.cluster6$cv_count)
metric.cluster6$tot_count <- as.numeric(metric.cluster6$tot_count)
head(metric.cluster6)

## Add row group col
metric.cluster6$group <- c('Centres','Centres','Centres','Centres','Centres','Centres','Centres','Centres','Mean','Mean','Mean','Mean','Mean','Mean','Mean','Mean','SD','SD','SD','SD','SD','SD','SD','SD')
head(metric.cluster6)

## Arrange data in long format
data_long <- tidyr::gather(metric.cluster6, metric, value, a_q0:tot_count, factor_key=TRUE)
head(data_long)

## Arrange data in wide format 
data_wide <- spread(data_long, ClusterNum, value)
head(data_wide)

## Change order of rows
data_wide_order <- data_wide
head(data_wide_order)

## Make metric a character
str(data_wide_order)
data_wide_order$metric <- as.character(data_wide_order$metric)

## Change names to symbols
data_wide_order$metric[data_wide_order$metric == 'a_q0'] <- 'α q=0'
data_wide_order$metric[data_wide_order$metric == 'a_q1'] <- 'α q=1'
data_wide_order$metric[data_wide_order$metric == 'a_q2'] <- 'α q=2'
data_wide_order$metric[data_wide_order$metric == 'b_q0'] <- 'β q=0'
data_wide_order$metric[data_wide_order$metric == 'b_q1'] <- 'β q=1'
data_wide_order$metric[data_wide_order$metric == 'b_q2'] <- 'β q=2'
data_wide_order$metric[data_wide_order$metric == 'g_q0'] <- 'γ q=0'
data_wide_order$metric[data_wide_order$metric == 'g_q1'] <- 'γ q=1'
data_wide_order$metric[data_wide_order$metric == 'g_q2'] <- 'γ q=2'
head(data_wide_order)
#View(data_wide_order)

## Remove unwanted row labels
data_wide_order[c(2:12,14:24,26:36),1] <- ''

## Load packages
library(reporter)
library(magrittr)

## Create metric symbols
data_wide_order[1,2] <-supsc("0")  %p% "D α"
data_wide_order[2,2] <-supsc("1")  %p% "D α"
data_wide_order[3,2] <-supsc("2")  %p% "D α"
data_wide_order[4,2] <-supsc("0")  %p% "D β"
data_wide_order[5,2] <-supsc("1")  %p% "D β"
data_wide_order[6,2] <-supsc("2")  %p% "D β"
data_wide_order[7,2] <-supsc("0")  %p% "D γ"
data_wide_order[8,2] <-supsc("1")  %p% "D γ"
data_wide_order[9,2] <-supsc("2")  %p% "D γ"
data_wide_order[10,2] <-"N"
data_wide_order[11,2] <-"N cv"
data_wide_order[12,2] <-"N tot"

data_wide_order[13,2] <-supsc("0")  %p% "D α"
data_wide_order[14,2] <-supsc("1")  %p% "D α"
data_wide_order[15,2] <-supsc("2")  %p% "D α"
data_wide_order[16,2] <-supsc("0")  %p% "D β"
data_wide_order[17,2] <-supsc("1")  %p% "D β"
data_wide_order[18,2] <-supsc("2")  %p% "D β"
data_wide_order[19,2] <-supsc("0")  %p% "D γ"
data_wide_order[20,2] <-supsc("1")  %p% "D γ"
data_wide_order[21,2] <-supsc("2")  %p% "D γ"
data_wide_order[22,2] <-"N"
data_wide_order[23,2] <-"N cv"
data_wide_order[24,2] <-"N tot"

data_wide_order[25,2] <-supsc("0")  %p% "D α"
data_wide_order[26,2] <-supsc("1")  %p% "D α"
data_wide_order[27,2] <-supsc("2")  %p% "D α"
data_wide_order[28,2] <-supsc("0")  %p% "D β"
data_wide_order[29,2] <-supsc("1")  %p% "D β"
data_wide_order[30,2] <-supsc("2")  %p% "D β"
data_wide_order[31,2] <-supsc("0")  %p% "D γ"
data_wide_order[32,2] <-supsc("1")  %p% "D γ"
data_wide_order[33,2] <-supsc("2")  %p% "D γ"
data_wide_order[34,2] <-"N"
data_wide_order[35,2] <-"N cv"
data_wide_order[36,2] <-"N tot"
head(data_wide_order)
#View(data_wide_order)

## Remove the SD values
data_wide_order <- data_wide_order[1:24,]
View(data_wide_order)
#_______________________________________________________________________________
#### BIODIVERSITY CLUSTERS: ASSIGN COLOURS ####

## Use Biodiversity Stripes (see https://findingnature.org.uk/2022/08/10/biodiversity-stripes/)

## Load package
library(janitor)

## Example colours - set up 8 equally spaced groups
values <- c(10, 20, 30, 40, 50, 60, 70, 80)

## Scale your values to range between 0 and 1
rr <- range(values)
svals <- (values-rr[1])/diff(rr)

## Colour range
f <- colorRamp(c( "#C0C1BC","#FCFA0A","#37C331" )) #"#00B050"
colors <- rgb(f(svals)/255)

## Check that it works
image(seq_along(svals), 1, as.matrix(seq_along(svals)), col=colors,
      axes=FALSE, xlab="", ylab="")

## Return colour hex codes
colors
"#C0C1BC" "#D1D189" "#E2E256" "#F3F223" "#E0F210" "#A8E21B" "#6FD326" "#37C331"
#_______________________________________________________________________________
#### BIODIVERSITY CLUSTERS: GROUP CHARACTERISTICS CENTRES ONLY (TABLE 4) ####

## Take just the centres for variables used in clustering
col_scale_data <- data_wide_order[c(1,4,6,7,10,11,12),2:10]
str(col_scale_data)
 
# Add total row (centres summed). Note these values used to set cluster colours.
col_scale_data_new <- col_scale_data %>%
  bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~'Total')))

# View new dataframe
print(col_scale_data_new)

# View new dataframe
#colnames(col_scale_data_new)[1] <- 'Metric'

## Update order of columns
col_scale_data_new2 <- col_scale_data_new[,c(1,3,9,6,2,5,7,8,4)]
print(col_scale_data_new2)

## Load package
library(gt)

## Create table
metric.cluster8 <-  col_scale_data_new2%>%gt()%>%
   sub_missing()%>%
  cols_align(
    align = c("center"),
    columns = c(2:9))%>%
  fmt_number(columns = c('2','8','5','1','4','6','7','3'), decimals = 2)%>%
  fmt_number(rows = c(8), decimals = 2)%>%# change number of dp on rows ?
  tab_spanner(
    label = md("**Cluster**"),
    columns = c( '2','8','5','1','4','6','7','3'))%>%
   data_color(
    direction = "row",
    columns = c(2:9),
    rows = c(1:8),# rows to have white to gray shading (SD excluded)
    method = "numeric",
    palette = c(
      "white","#737373"),
    na_color = "white")%>%
 tab_style(
    style = list(
      cell_fill(color =  "#37C331")),
    location = list(
      cells_column_labels(columns = c('2'))))%>%
  tab_style(
    style = list(
      cell_fill(color = "#6FD326")),
    location = list(
      cells_column_labels(columns = c('8'))))%>% 
  tab_style(
    style = list(
      cell_fill(color ="#A8E21B")),
    location = list(
      cells_column_labels(columns = c('5'))))%>%
  tab_style(
    style = list(
      cell_fill(color ="#E0F210")),
    location = list(
      cells_column_labels(columns = c('1'))))%>%
  tab_style(
    style = list(
      cell_fill(color ="#F3F223")),
    location = list(
      cells_column_labels(columns = c('4'))))%>%
  tab_style(
    style = list(
      cell_fill(color = "#E2E256")),
    location = list(
      cells_column_labels(columns = c('6'))))%>%
      tab_style(
    style = list(
      cell_fill(color =  "#D1D189")),
    location = list(
      cells_column_labels(columns = c('7'))))%>%
   tab_style(
    style = list(
      cell_fill(color =  "#C0C1BC")),
    location = list(
      cells_column_labels(columns = c('3'))))%>%
     cols_label(
       metric = md("**Metric**"),
    '1' = md("**Bio-D**"),
    '2' = md("**Bio-A**"),
    '3' = md("**Bio-H**"),
    '4' = md("**Bio-E**"),
    '5' = md("**Bio-C**"),
    '6' = md("**Bio-F**"),
    '7' = md("**Bio-G**"),
    '8' = md("**Bio-B**"))%>%
    tab_style(
    style = cell_borders(
      sides = c( "bottom"),
      color = "#d3d3d3",
      weight = px(2),
      style = "solid"),
    locations = cells_body(
      columns = 2:8,
      rows = 7))%>% 
  tab_style(
    style = cell_borders(
      sides = c( "bottom"),
      color = "#d3d3d3",
      weight = px(2),
      style = "solid"),
    locations = cells_body(
      columns = 2:8,
      rows = 8))%>%
  cols_width(
    everything() ~ px(75)  # Set all columns to 150 pixels wide
  )%>%
  tab_style(
    style = list(
      cell_text(style = "italic")),
    locations = cells_body(
      columns = metric))%>%
  fmt_markdown(columns = metric)#%>% 
#tab_caption(caption = md("**"))
##Table 6.** Biodiversity cluster group centres for metrics used in clustering.Shading  provides a  visual means of identifying highest (dark grey) and lowest (no shading) values by row. Column totals in bottom row are used for colouring.
## View table
metric.cluster8

## Save Response traits table in html
metric.cluster8 %>%
 # gtsave("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Table_6.png", expand = 10)# save table as .png
 gtsave("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Table_6.html")# save in html (tabular) format
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: PREPARE RASTER DATA ####

# Load packages
library(raster)

## Load biodiv data used for clustering (i.e. outliers removed, transformed, standardized, covariates removed - see step 2.3.1.5)
bio <- clus_data1

## Load factor data (faunal.cluster7 from step 2.3.2.5)
fac <- faunal.cluster7

## Change names of cols
colnames(fac)=c("Sample","lon","lat","cluster")

## Make cluster a factor
fac$cluster=as.factor(fac$cluster)

## Get PHY data from rasters using variables used in RF modelling. Note sediment data comes from actual samples
phyc_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/phyc_mean.tif")# PHYTOPLANKTON
thetao_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/thetao_mean.tif")# BOTTOM TEMP
no3_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/no3_mean.tif")# NITRATE
Current_Sp<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Current_Sp.tif")# CURRENT SPEED
vd<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/VD.tif")# VALLEY DEPTH
so_range<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/so_range.tif")# SALINITY RANGE
dfe_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/dfe_mean.tif")# DISSOLVED IRON
Wave_veloc<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/wave_veloc.tif")# WAVE VELOCITY
CND<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/CND.tif")#CHANNEL NETWORK DISTANCE
LSF<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/LSF.tif")#LS-FACTOR
CDP0<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/CDP0.tif")#CLOSED DEPRESSIONS
#si_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/si_mean.tif")# SILICATE
#so_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/so_mean.tif")# SALINITY MEAN
#Bathymetry<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Bathymetry.tif")
#chl_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/chl_mean.tif")
#CNBL<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/CNBL.tif")
#RSP<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/RSP.tif")
#gravel<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Predicted_Gravel_Fraction.tif")
#mud<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Predicted_Mud_Fraction.tif")
#SPM_MEAN<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/SPM_MEAN.tif")
#SPM_SUMMER<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/SPM_SUMMER.tif")
#SPM_WINTER<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/SPM_WINTER.tif")
#o2_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/o2_mean.tif")
#ph_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/ph_mean.tif")
#ph_range<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/ph_range.tif")
#po4_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/po4_mean.tif")
#chl_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/chl_mean.tif")
#KDPAR_mean_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/KDPAR_mean_mean.tif")

## Create raster stack
predictors <- stack(
  phyc_mean,
  thetao_mean,
  no3_mean,
  Current_Sp,
  vd,
  so_range,
  dfe_mean,
  Wave_veloc,
  CND,
  LSF,
  CDP0
)

## Plot raster stack
plot(predictors)

## Update names for predictor variables
names(predictors)=c(
'Phytoplankton',
'Bottom temp.',
'Nitrate',
'Current speed',
'Valley depth',
'Salinity range',
'Diss. Iron',
'Wave velocity',
'Ch. network distance',
'LS-factor',
'Closed depressions'

#'Bathymetry',
#'Chlorophyll',
#'CNBL',
#'RSP',
#'gravel',
#'mud',
#'Mean SPM',
#'Summer SPM',
#'Winter SPM',
#'Diss. Oxygen',
#'Nitrate',
#'pH',
#'Bottom temp. range',
#'ph_range',
#'Phosphate',
#'chl_mean',
#'KDPAR_mean_mean',
)#"Sand",

## Unload extract function from tidyr package (otherwise it won't work)
#.rs.unloadPackage("tidyr")

## Extract predictor variables from raster stack and store in df
sdata <- raster::extract(predictors, fac[,2:3])

##Change from matrix to df
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(fac$Sample,fac$cluster,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "Cluster"
#head(sdata2)
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: PREPARE SEDIMENT DATA ####

## Load packages
library(pool)
library(DBI)
library (RPostgres)
library(dplyr)

## Connect to OneBenthic DB
Sys.setenv(R_CONFIG_ACTIVE = "one_benthic")

dw <- config::get()

pool <- dbPool(drv = dbDriver(dw$driver),
               dbname = dw$database,
               host = dw$server,
               port =  dw$port,
               user = dw$uid,
               password = dw$pwd)

## Load sediment sieve data
sed_data = dbGetQuery(pool,
"
  select 
svs.sample_samplecode,
sv.wentworth_id,
w.wentworthclass,
svs.percentage
from sediment_data.sedvarsample as svs
inner join sediment_data.sedvar as sv on svs.sedvar_sievesize = sv.sievesize
inner join sediment_data.wentworth as w on w.wwid = sv.wentworth_id;")
#head(sed_data)
#dim(sed_data)

## Subset for samples present in bio data (df 'fac')
sed_data2 <- subset(sed_data, sample_samplecode %in% fac$Sample)

## Change percentage into numeric
sed_data2$percentage <- as.numeric(sed_data2$percentage)

## Add sed categories
sed_data2$wentworth_id[sed_data2$wentworth_id == 1] <- "Mud"
sed_data2$wentworth_id[sed_data2$wentworth_id == 2] <- "Sand"
sed_data2$wentworth_id[sed_data2$wentworth_id == 3] <- "Sand"
sed_data2$wentworth_id[sed_data2$wentworth_id == 4] <- "Sand"
sed_data2$wentworth_id[sed_data2$wentworth_id == 5] <- "Gravel"
sed_data2$wentworth_id[sed_data2$wentworth_id == 6] <- "Gravel"
sed_data2$wentworth_id[sed_data2$wentworth_id == 7] <- "Gravel"
sed_data2$wentworth_id[sed_data2$wentworth_id == 8] <- "Cobbles"
head(sed_data2)

# Group by sum using dplyr
library(dplyr)
sed_data3 <- sed_data2 %>% group_by(sample_samplecode,wentworth_id)%>%summarise(sum=sum(percentage))
head(sed_data3)                                                                                

## Change from long to wide format
library(tidyr)
sed_data4 <- as.data.frame(spread(sed_data3, wentworth_id, sum))

## Update column name
colnames(sed_data4)[1] <- 'Sample'

## Calculate total percentage
sed_data4$total <- rowSums(sed_data4[,3:5])

## Merge response and predictor variables (df 'sdata2') with sediment data
phy <- merge(sdata2,sed_data4,by="Sample")
head(phy)

## Drop the cobbles column
phy <- within(phy, rm(Cobbles))
phy <- within(phy, rm(total))
phy <- within(phy, rm(Sand))

## Drop rows with missing data
phy <-  phy[complete.cases(phy), ]

## Merge the bio and phy data (already includes factor)
phy_bio <- left_join(phy, bio, by = "Sample")
head(phy_bio)

## Check the number of available samples by faunal cluster group
table(phy_bio$Cluster)# Answer is 1421  (8)
#   2    8    5    1    4    6    7    3 
#1788 1619 1244 1088 1628 2226 1185  909 

## Set up objects for number of samples by cluster
one <-1088
two<- 1788
three<- 909
four<- 1628
five<- 1244
six<- 2226
seven <- 1185
eight<- 1619

## Create a subset of the data for cluster group 1
SS1= subset(phy_bio,Cluster=="1")

## Check number of samples in df 'SS1' tallies with the above table
#dim(SS1)

## Generate random numbers (without replacement) from 1 to n
set.seed(1)# to stop random numbers changing
SS1num=sample(1:one,one)

## Add random numbers to new col 'SS1NUM' in df SS1
SS1$SSNumber <-SS1num
#head(SS1)

## Now do the same for cluster group 2
SS2= subset(phy_bio,Cluster=="2")
#dim(SS2)
set.seed(2)# to stop random numbers changing
SS2num=sample(1:two,two)
SS2$SSNumber <-SS2num

## Now do the same for cluster group 3
SS3= subset(phy_bio,Cluster=="3")
#dim(SS3)
set.seed(2)# to stop random numbers changing
SS3num=sample(1:three,three)
SS3$SSNumber <-SS3num

## Now do the same for cluster group 4
SS4= subset(phy_bio,Cluster=="4")
#dim(SS4)
set.seed(2)# to stop random numbers changing
SS4num=sample(1:four,four)
SS4$SSNumber <-SS4num

## Now do the same for cluster group 5
SS5= subset(phy_bio,Cluster=="5")
#dim(SS5)
set.seed(2)# to stop random numbers changing
SS5num=sample(1:five,five)
SS5$SSNumber <-SS5num

## Now do the same for cluster group 6
SS6= subset(phy_bio,Cluster=="6")
#dim(SS6)
set.seed(2)# to stop random numbers changing
SS6num=sample(1:six,six)
SS6$SSNumber <-SS6num

## Now do the same for cluster group 7
SS7= subset(phy_bio,Cluster=="7")
#dim(SS6)
set.seed(2)# to stop random numbers changing
SS7num=sample(1:seven,seven)
SS7$SSNumber <-SS7num

## Now do the same for cluster group 8
SS8= subset(phy_bio,Cluster=="8")
#dim(SS6)
set.seed(2)# to stop random numbers changing
SS8num=sample(1:eight,eight)
SS8$SSNumber <-SS8num

## Stitch together all the data subsets
best.data.ss=rbind(SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8)

## Now select only samples with random no. < 201
best.data.ss.final=subset(best.data.ss, SSNumber<201)
#dim(best.data.ss.final)# should be 1600 (8 x 200) -it is
#names(best.data.ss.final)

## Create a df 'bestBIO' for biodiversity data (metrics used for clustering) and save
bestBIO=best.data.ss.final[,c(18:24)]
head(bestBIO)
dim(bestBIO)#1600 7
#write.csv(bestBIO,file = "OUTPUTS/bestBIO.csv",row.names=TRUE)

## Create a df 'bestPHY' for physical variables and save.
bestPHY=best.data.ss.final[,c(3:15)]
#bestPHY$Bathymetry <- abs(bestPHY$Bathymetry)# Make bathymetric values positive
head(bestPHY)
dim(bestPHY)#1600   13

sum(is.na(bestPHY))

#View(bestPHY)
#write.csv(bestPHY,file = "OUTPUTS/bestPHY.csv",row.names=TRUE)

## Create a df 'bestFAC' for factor cluster and save
bestFAC=as.data.frame(best.data.ss.final[,2])
head(bestFAC)
dim(bestFAC)#1600 1
sum(is.na(bestFAC))
colnames(bestFAC)[1] <- 'Cluster'
#write.csv(bestFAC,file = "OUTPUTS/bestFAC.csv",row.names=TRUE)
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: DETECT COLLINEARITY IN ENVIRONMENTAL VARIABLES ####

#Using Variance Inflation Factors (VIFs) from usdm package

## Load package
library(usdm)

## Get correlation coefficients for environmental variables
cor(bestPHY, use="all.obs", method="spearman")

## Use vifstep function to identify set of variables which are not highly correlated (i.e. Variance Inflation Factor VIF <2.5)
vifstep(bestPHY, th=2.5)

## Create new df for variable with VIF SCORES <2.5
bestPHY2=subset(bestPHY, select = c(
Phytoplankton,
Bottom.temp.,
Current.speed,
Valley.depth,
Salinity.range,
Diss..Iron,
Ch..network.distance,
LS.factor,
Closed.depressions,
Gravel,
Mud

))
dim(bestPHY2)#1600 11
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: CHECK FOR SKEWNESS AND TRANSFORM AS NECESSARY ####

# Check for any skewness in the data. Negative skewness indicates that the mean of the data
# values is less than the median,and the data distribution is left-skewed. Positive skewness
# would indicate that the mean of the data values is larger than the median, and the data
# distribution is right-skewed.
#bestPHY2 <- bestPHY
#View(bestPHY2)
summary(bestPHY2)# Positive skewness (mean>median) for SPM_SUMMER (1), gravel (2), Mud (4), VD (7), CND (10), LSF (12)

## Return skewness values
bestPHY2[,1:11] %>% 
  dplyr::mutate(across(1:11, moments::skewness))%>% 
  distinct()



## Transform relevant columns 
#bestPHY2[,c(1,2,4,5,7,8,9,10,11)]=log(bestPHY2[,c(1,2,4,5,7,8,9,10,11)]+0.1)
bestPHY2[,c(7,8,9,10,11)]=log(bestPHY2[,c(7,8,9,10,11)]+0.1)
View(bestPHY2)

## Find row with NA
sum(is.na(bestPHY2))
which(is.na(bestPHY2), arr.ind = TRUE)

## Need to remove row with NaN
bestPHY2 <- bestPHY2[-1278, ]

## Check dim of bestPHY2
dim(bestPHY2)#1599

## Need to update bestBIO and bestFAC too
bestBIO <- bestBIO[-1278, ]
bestFAC <- bestFAC[-1278, ]

## Check dim of bestBIO
dim(bestBIO)#1599

## Check length of bestFAC
length(bestFAC)#1599
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: BIOENV (TABLE 7) ####

## Run bioenv with the transformed faunal (df 'bestBIO') and the selected env variables
#(df 'bestPHY2').

# First make sure no NAs in data. Start by combining dataframes
df_combined <- cbind(bestBIO, bestPHY2)
view(df_combined)
dim(df_combined)#1599

## Remove rows with na
df_combined <- na.omit(df_combined)
names(df_combined)

## Reinstate objects
bestBIO <- df_combined[,1:7]
bestPHY2 <- df_combined[,8:18]

## Check dfs same length
dim(bestBIO)# 1599 7
dim(bestPHY2)# 1599 11
length(bestFAC)# 1599

# Normalize the data
bestPHY2 <- scale(bestPHY2)


## Perform BIOENV analysis
library(vegan) # Load library
res<-bioenv(bestBIO, bestPHY2) 
res # correlation = 0.1664933 

## See all best results
summary(res)

## Output results as a dataframe
res2 <-data.frame(unclass(summary(res)), check.names = FALSE, stringsAsFactors = FALSE)
res3 <- res2[,c(1,3,2)]
class(res3)

## Reduce number of decimal places in correlation column
op = function(x, d=2) sprintf(paste0("%1.",d,"f"), x) 
res3$correlation <- op(res3$correlation, 4)
colnames(res3) <- c('Size','Variables','Correlation (ρ)')

## Update variable names
res3 <-res3%>% mutate(Variables = gsub("\\Gravel", "Gravel %", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Mud", "Mud %", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\LS.factor", "LS factor", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Silicate", "Silicate", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Salinity.range", "Salinity range (ppt)", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Ch..network.distance", "Ch. network distance (m)", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Current.speed", "Current speed", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Summer.SPM", "Summer SPM", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Bottom.temp.", "ottom temp.", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Valley.depth", "Valley depth", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Closed.depressions", "Closed depressions", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Diss..Iron", "Diss. Iron (mmol m-3)", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Phytoplankton", "Phytoplankton (mmol m-3)", Variables))
res3 <-res3%>% mutate(Variables = gsub("\\Current speed", "Current speed (m s-1)", Variables))
res3

## Load packages
library(knitr)
library(dplyr)
library(kableExtra)

## Create table as kable output
kable_table <- kable(res3[1:8,], escape=FALSE,format = "html",caption = "")%>%
  kable_styling()%>%
  row_spec(6,bold=T,hline_after = T)
  
# Save the kable table as an HTML file
save_kable(kable_table, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Table_7.html")

## Create results table as .png

## Load packages
#library(knitr)
#library(kableExtra)
#library(magrittr)
#library(webshot2)
#library(magick)
#library(dplyr)
#library(kable)
#webshot::install_phantomjs(force=TRUE)#force=TRUE

## Set up file patch for html output
#html_file <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/OUTPUTS/html_file.html"

## Save kable. Do this first to html and then to png (won't save to png directly)
#kable(res3[1:10,], escape=FALSE,format = "html",caption = "Table 8. Results of a 'best' analysis identifying the subset of environmental variables which are most correlated with the biodiversity data.")%>%
#  kable_styling()%>%
#  row_spec(3,bold=T,hline_after = T)%>%
#  save_kable(file = html_file)
 
## Capture the screenshot of the HTML file
#png_file <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/OUTPUTS/Table_8.png"
#webshot(html_file, file = png_file)

## Read the PNG file
#img <- image_read("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/OUTPUTS/Table_8.png")

## Trim the white space
#img_trimmed <- image_trim(img)

## Save the trimmed image
#image_write(img_trimmed, "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/OUTPUTS/Table_8.png")

#_______________________________________________________________________________
#### EXPLAINING PATTERNS: ADONIS ####

## Use 'adonis' to quantify the variation in biodiversity data explained by the environmental predictors identified in 'best'. Requires normalised phy data in df format

## Load package
library(vegan)

## Normalise env data prior - already done

## bestPHY2 is a matrix so need to convert back to a df for use with adonis
bestPHY2 <- as.data.frame(bestPHY2)

## Use Adonis to see how much of the variation is explained by the different variables. Enter
# phy variable which are important (see BEST results) http://www.talkstats.com/showthread.php/15912-Permanova-Adonis
adonis.res=adonis2(formula = bestBIO ~ Phytoplankton + Current.speed + Ch..network.distance + LS.factor + Gravel + Mud, data = bestPHY2, permutations = 999,method = "bray",by = "terms")
adonis.res

# Results
# 79.2% of the variation remains unexplained by the model
# variables explain 20.8 % pf variation (gravel = 8.0%, phytoplankton = 6.7%, , mud = 3.1%, current speed =2.2%, LS-factor = 0.6%
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: dbRDA ORDINATION (FIGURE 6) ####

## Load libraries
library(ggplot2)
library(vegan)
library(gridExtra)
library(ggrepel)
# Use distance based redundancy analysis (dbRDA) ordination to visualise the  
# relationship between macrofaunal data and predictor variables https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/capscale

## Perform dbRDA
vare.cap <- capscale(bestBIO ~  Phytoplankton + Current.speed + Ch..network.distance + LS.factor + Gravel + Mud, bestPHY2,dist="bray")
vare.cap

## Extract site scores and explanatory variables
site_scores <- scores(vare.cap, display = "sites")
species_scores <- scores(vare.cap, display = "species")
explanatory_scores <- scores(vare.cap, display = "bp")

## Scale the explanatory scores to make arrows longer
explanatory_scores <- explanatory_scores *  1.8# Adjust the scaling factor as needed

## Convert to data frames for ggplot
site_scores_df <- as.data.frame(site_scores)
species_scores_df <- as.data.frame(species_scores)
explanatory_scores_df <- as.data.frame(explanatory_scores)

# Add cluster information to site scores
site_scores_df$Cluster <-bestFAC#$Cluster

## Check types
class(site_scores_df)
class(bestFAC)

## Check dims
length( bestFAC)#1599
dim(site_scores_df)#1593 3

## Define custom labels
custom_labels <- c('Phytoplankton','Current speed', 'Ch. network distance', 'LS-factor', 'Gravel' , 'Mud')

## Define colors
colors <- c("#37C331", "#6FD326", "#A8E21B", "#E0F210", "#F3F223", "#E2E256", "#D1D189", "#C0C1BC")

## Define the order of legend items
legend_order <- c( '2','8','5','1','4','6','7','3')


# Plot
p <- ggplot() +
  geom_point(data = site_scores_df, aes(x = CAP1, y = -CAP2, color = as.factor(Cluster)), size = 2) +
  geom_segment(data = explanatory_scores_df, aes(x = 0, y = 0, xend = CAP1, yend = -CAP2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = explanatory_scores_df, aes(x = CAP1, y = -CAP2, label = custom_labels), color = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
   #scale_color_manual(values = colors, breaks = legend_order) +
  ####
  

scale_colour_manual(
  breaks = c('2','8','5','1','4','6','7','3'),
  values = c(
    '2' = "#37C331",
    '8' = "#6FD326",
    '5' = "#A8E21B",
    '1' = "#E0F210",
    '4' = "#F3F223",
    '6' = "#E2E256",
    '7' = "#D1D189",
    '3' = "#C0C1BC"
  ),
  labels = c(
    '2' = "Bio-A",
    '8' = "Bio-B",
    '5' = "Bio-C",
    '1' = "Bio-D",
    '4' = "Bio-E",
    '6' = "Bio-F",
    '7' = "Bio-G",
    '3' = "Bio-H"
  ),
  name = "Cluster",
  na.value = "transparent",
  na.translate = FALSE)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
  legend.text = element_text(color = "black", size = 14),  # Increase text size
  legend.title = element_text(color = "black", size = 16), # Optional: increase title siz
) +
  labs(color = "Cluster") +
  xlim(-2.5, 2.5) +
  ylim(-3, 2.5)+
  guides(color = guide_legend(override.aes = list(size = 4)))  # Increase legend point size

p
# Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicHotSpots/OUTPUTS/Figure_6.png", plot = p, width = 20, height = 20, units = "cm", dpi = 800)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF INDIVIDUAL METRICS: RASTER PREDICTOR VARIABLES (QUICK LOOK) ####

# This code is intended for a quick look-see. For final models and associated RF outputs (inc associated model confidence maps), use R script in file xx

## Load packages
library(raster)

## Load environmental predictor rasters
bathy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/bathy3.tif")
cur <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/cur3.tif")
gravel <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Predicted_Gravel_Fraction.tif")
light <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/light3.tif")
mud <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Predicted_Mud_Fraction.tif")
oxy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/oxy3.tif")
phyto <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/phyto3.tif")
sal <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/sal3.tif")
sil <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/sil3.tif")
spm <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/SPM_MEAN.tif")
temp <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/temp3.tif")
wov <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Wave_veloc.tif")

## Update crs
crs(bathy) <- "+proj=longlat +datum=WGS84 +no_defs" 

## Create raster stack
predictors <- stack(bathy,cur,gravel,light,mud,oxy,phyto,sal,sil,spm,temp,wov)
#predictors

# Simple names for predictor variables
names(predictors)=c("Bathymetry","Current","Gravel","Light","Mud","Oxygen","Phyto","Salinity","Silicate","SPM","Temperature","WOV")#"Sand",
#names(predictors)

## Plot raster stack
#plot(predictors)
#plot(bathy)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF INDIVIDUAL METRICS: METRIC SELECTION ####

## Select data by metric (do metrics one at a time)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='sample_a_q0', na.rm = TRUE)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='sample_a_q1', na.rm = TRUE)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='sample_a_q2', na.rm = TRUE)

metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='b_q0', na.rm = TRUE)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='b_q1', na.rm = TRUE)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='b_q2', na.rm = TRUE)

metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='g_q0', na.rm = TRUE)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='g_q1', na.rm = TRUE)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='g_q2', na.rm = TRUE)

metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='count', na.rm = TRUE)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='cv_count', na.rm = TRUE)
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='tot_count', na.rm = TRUE)
head(metric)
#_______________________________________________________________________________
## Calculate rare taxa metric (gamma 2 / gamma 0 * 100)

## Select 
metric <- biodiv9_mod %>% group_by(metric) %>% filter(metric=='g_q0' |metric=='g_q2' , na.rm = TRUE)
dim(metric)#37370     5
#View(metric )

## Wide data format 
library(tidyr)
metric_wide <- spread(metric, metric, measurement)
dim(metric_wide)# 21766     5
#View(metric_wide)

## Remove rows with NA
metric_wide_complete <- metric_wide[complete.cases(metric_wide), ]
dim(metric_wide_complete)#22741     5

## Create new col for 'g_rare'
metric_wide_complete$g_rare <- (1-metric_wide_complete$g_q2/ metric_wide_complete$g_q0)*100
metric_wide_complete
dim(metric_wide_complete)#15604     6

## Take just columns of interest with col names like other 'metric' objects
metric_wide_complete2 <- metric_wide_complete[,c(1:3,6)]
colnames(metric_wide_complete2)[4] <- 'measurement'
metric_wide_complete2$metric <- 'g_rare'
metric_wide_complete3 <- metric_wide_complete2[,c(1:3,5,4)]#re-order cols
metric <- as_tibble(metric_wide_complete3)
#View(metric)
str(metric)
class(metric)
head(metric)

## Save metric 'rare' for modelling
metric_rare <- metric[,c(3,1,2,4,5)]
metric_rare$type <- 'numeric'
metric_rare$paper <- 'biodiversity'
head(metric_rare)
colnames(metric_rare) <- c('sample','x','y','metric','value','type','paper')
metric_rare <- metric_rare[,c(7,1:4,6,5)]
head(metric_rare)
dim(metric_rare)

## Save metric rare data for use with RF modelling script (see xx.R)
write.csv(metric_rare, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\biodiv_metric_rare_4_modelling.csv", row.names=FALSE)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF INDIVIDUAL METRICS: EXTRACT PREDICTOR VARIABLE VALUES FROM RASTER STACK ####

## Extract
sdata <- raster::extract(predictors, metric[,1:2])

## Change from matrix to df
class(sdata)
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(metric$Sample,metric$measurement,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "value"
head(sdata2)

## Change cols to appropriate type
str(sdata2)
#sdata2$Cluster=as.factor(sdata2$Cluster)
#sdata2$Sample=as.character(sdata2$Sample)
#sdata2$AvCur=as.numeric(sdata2$AvCur)
#sdata2$Chla=as.numeric(sdata2$Chla)
#sdata2$Depth=as.numeric(sdata2$Depth)
#$Gravel=as.numeric(sdata2$Gravel)
#sdata2$Mud=as.numeric(sdata2$Mud)
#sdata2$Sal=as.numeric(sdata2$Sal)
#sdata2$Sand=as.numeric(sdata2$Sand)
#sdata2$SPM=as.numeric(sdata2$SPM)
#sdata2$Stress=as.numeric(sdata2$Stress)
#sdata2$Temp=as.numeric(sdata2$Temp)
#sdata2$WOV=as.numeric(sdata2$WOV)

## First check cols of correct type
str(sdata2)
dim(sdata2)# 21829    14
#____________________________________________________________________________________________________________________
#### SPATIAL MODELLING OF INDIVIDUAL METRICS: MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Call library
#install.packages("caTools")
library(caTools)

## Vector for response variable
Y <- sdata2[,2]

## Take a random 90% of samples in proportion across the 11 groups
set.seed(2)
msk= sample.split( Y, SplitRatio = 9/10, group = NULL )

## Check it was 90% 10% split

## The training set
train = sdata2[ msk,] 
dim(train)# 19646    14
head(train)

## Remove station labels for train
train2 =train[,2:14]#was 11
head(train2)

## The test set
test  = sdata2[ !msk,]
dim(test) # 2183   14
head(test)

## Remove station labels for test
test2 =test[,2:14]#was 11
head(test2)
str(test2)

## Check number of observations for train (TRUE) and test (FALSE) sets
print(table(Y, msk)) 

## Check number of samples in train and test sets is equal to total
dim(sdata) # 21829    12
dim(train)+dim(test)# 21829    28
#____________________________________________________________________________________________________________________
#### SPATIAL MODELLING OF INDIVIDUAL METRICS: DO MODELLING ####

## Call library
library(randomForest)

## Model
model <- value ~Bathymetry+Current+Gravel+Light+Mud+Oxygen+Phyto+Salinity+Silicate+SPM+Temperature+WOV

## Prepare training dta
train3 <- train2[complete.cases(train2), ]
train3$value <- as.numeric(train3$value)
str(train3)
head(train3)

## Run model
rf2 <- randomForest(model, data=train3)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF INDIVIDUAL METRICS: EXAMINE HOW VARIABLES ARE INFLUENCING THE MODEL ####

## Produce plots 
#varImpPlot(rf2)

#png('OUTPUTS/STRUCTURE/variables_affecting_model_structure.png') # height=nrow(pr), width=ncol(pr) EFFECTS TRAITS
#varImpPlot(rf2)
#dev.off()

#____________________________________________________________________________________________________________________
#### SPATIAL MODELLING OF INDIVIDUAL METRICS: PRODUCE FULL COVERAGE RASTER FOR EACH METRIC ####

## Use model to predict cluster group for each raster cell
pr <- predict(predictors, rf2)
#cols <-  terrain.colors(255)
#plot(pr,col= cols)
plot(pr)

## Save as .tiff
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\sample_a_q0.tif',overwrite=TRUE,format = "GTiff")
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\sample_a_q1.tif',overwrite=TRUE,format = "GTiff")
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\sample_a_q2.tif',overwrite=TRUE,format = "GTiff")

writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\b_q0.tif',overwrite=TRUE,format = "GTiff")
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\b_q1.tif',overwrite=TRUE,format = "GTiff")
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\b_q2.tif',overwrite=TRUE,format = "GTiff")

writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\g_q0.tif',overwrite=TRUE,format = "GTiff")
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\g_q1.tif',overwrite=TRUE,format = "GTiff")
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\g_q2.tif',overwrite=TRUE,format = "GTiff")

writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\count.tif',overwrite=TRUE,format = "GTiff")
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\cv_count.tif',overwrite=TRUE,format = "GTiff")
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\tot_count.tif',overwrite=TRUE,format = "GTiff")

writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\g_rare.tif',overwrite=TRUE,format = "GTiff")

#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: SETUP ####

## Load packages
library(sf)
library(dplyr)
library(RColorBrewer)
library(terra)
library(raster)
library(stars)
library(ggplot2)
library(colorRamps)
library(terra)
library(tidyterra)
library(raster)
library(dplyr)
library(stars)
library(ggplot2)
library(colorRamps)
library(climateStability)
library(ggpubr)

## Bring in countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))

## Data required for modelling biodiversity metrics (long format, outliers removed, all metrics)
head(biodiv9_mod)
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D0 Alpha ####

## Load model raster
a_q0_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/0Dα (Hill 0, Alpha)/Model/sample_a_q0_Mean_SQ_2025.tif')

## Reduce size of raster
a_q0_model_agg <- aggregate(a_q0_model, fact=7,fun = modal)

## Remove raster to save space
rm(a_q0_model)

## Annotate label
label1 =  "{}^{0}~italic(D)[italic(α)]" #

## Model plot
D0_alpha <-  ggplot() +
  geom_spatraster(data = a_q0_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    panel.background = element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D0 Beta ####

## Load raster layer
b_q0_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/0Dβ (Hill 0, Beta)/Model/b_q0_Mean_SQ_2025.tif')

## Reduce size of raster
b_q0_model_agg <- aggregate(b_q0_model, fact=7,fun = modal)

## Remove raster to save space
rm(b_q0_model)

## Annotate label
label1 = "{}^{0}~italic(D)[italic(β)]"# "~beta~ q == 0"

## Model plot
D0_beta <- ggplot() +
  geom_spatraster(data = b_q0_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
   x = -8.3,
           y = 59.7)+
theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D0 Gamma ####

## Load raster layer
g_q0_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/0Dγ (Hill 0, Gamma)/Model/g_q0_Mean_SQ_2025.tif')

## Reduce size of raster
g_q0_model_agg <- aggregate(g_q0_model, fact=7,fun = modal)

## Remove raster to save space
rm(g_q0_model)

## Annotate label
label1 ="{}^{0}~italic(D)[italic(γ)]" # "~gamma~ q == 0"

## Model plot
D0_gamma <- ggplot() +
  geom_spatraster(data = g_q0_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D1 alpha ####

## Load raster layer
a_q1_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/1Dα (Hill 1, Alpha)/Model/sample_a_q1_Mean_SQ_2025.tif')

## Reduce size of raster
a_q1_model_agg <- aggregate(a_q1_model, fact=7,fun = modal)

## Remove raster to save space
rm(a_q1_model)

## Annotate label
label1 = "{}^{1}~italic(D)[italic(α)]"#"~alpha~ q == 1"

## Model plot
D1_alpha <- ggplot() +
  geom_spatraster(data = a_q1_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D1 beta ####

## Load model raster D1 beta
b_q1_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/1Dβ (Hill 1, Beta)/Model/b_q1_Mean_SQ_2025.tif')

## Reduce size of raster
b_q1_model_agg <- aggregate(b_q1_model, fact=7,fun = modal)

## Remove raster to save space
rm(b_q1_model)

## Annotate label
label1 = "{}^{1}~italic(D)[italic(β)]" # "~beta~ q == 1"

## Model plot
D1_beta <- ggplot() +
  geom_spatraster(data = b_q1_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D1 gamma ####

## Load raster
g_q1_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/1Dγ (Hill 1, Gamma)/Model/g_q1_Mean_SQ_2025.tif')

## Reduce size of raster
g_q1_model_agg <- aggregate(g_q1_model, fact=7,fun = modal)

## Remove raster to save space
rm(g_q1_model)

## Annotate label
label1 = "{}^{1}~italic(D)[italic(γ)]"# "~gamma~ q == 0"

## Model plot
D1_gamma <- ggplot() +
  geom_spatraster(data = g_q1_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D2 alpha ####

## Load raster
a_q2_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/2Dα (Hill 2, Alpha)/Model/sample_a_q2_Mean_SQ_2025.tif')

## Reduce size of raster
a_q2_model_agg <- aggregate(a_q2_model, fact=7,fun = modal)

## Remove raster to save space
rm(a_q2_model)

## Annotate label
label1 = "{}^{2}~italic(D)[italic(α)]" #"~alpha~ q == 2"

## Model plot
D2_alpha <- ggplot() +
  geom_spatraster(data = a_q2_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D2 beta ####

## Load model raster D2 beta
b_q2_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/2Dβ (Hill 2, Beta)/Model/b_q2_Mean_SQ_2025.tif')

## Reduce size of raster
b_q2_model_agg <- aggregate(b_q2_model, fact=7,fun = modal)

## Remove raster to save space
rm(b_q2_model)

## Annotate label
label1 ="{}^{2}~italic(D)[italic(β)]" # "~beta~ q == 2"

## Model plot
D2_beta <- ggplot() +
  geom_spatraster(data = b_q2_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: D2 gamma ####

## Load raster
g_q2_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/2Dγ (Hill 2, Gamma)/Model/g_q2_Mean_SQ_2025.tif')

## Reduce size of raster
g_q2_model_agg <- aggregate(g_q2_model, fact=7,fun = modal)

## Remove raster to save space
rm(g_q2_model)

## Annotate label
label1 = "{}^{2}~italic(D)[italic(γ)]" #"~gamma~ q == 2"

## Model plot
D2_gamma <- ggplot() +
  geom_spatraster(data = g_q2_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: N ####

## Load raster
count_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/N (Abund)/Model/av_count_Mean_SQ_2025.tif')

## Reduce size of raster
count_model_agg <- aggregate(count_model, fact=7,fun = modal)

## Remove raster to save space
rm(count_model)

## Annotate label
label1 = "italic(N)"##"N~alpha" "~av~count"

## Model plot
count_map <- ggplot() +
  geom_spatraster(data = count_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: N CV ####

## Load raster
cv_count_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Ncv (Abund cv)/Model/cv_count_Mean_SQ_2025.tif')

## Reduce size of raster
cv_count_model_agg <- aggregate(cv_count_model, fact=7,fun = modal)

## Remove raster to save space
rm(cv_count_model)

## Annotate label
label1 = "italic(N[cv])"#"~cv~count"

## Model plot
cv_count_map <- ggplot() +
  geom_spatraster(data = cv_count_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: N TOT ####

## Load raster
tot_count_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Ntot (Abund total)/Model/tot_count_Mean_SQ.tif')

## Reduce size of raster
tot_count_model_agg <- aggregate(tot_count_model, fact=7,fun = modal)

## Remove raster to save space

## Annotate label
label1 = "italic(N[tot])"#"~tot~count"

## Model plot
tot_count_map <- ggplot() +
  geom_spatraster(data = tot_count_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS: COMBINED PLOT (FIGURE 2) ####

## Create each row
D0_stitch <- egg::ggarrange(D0_alpha, D0_beta, D0_gamma, labels = c("", "",""),nrow=1)#ggpubr
D1_stitch <- egg::ggarrange(D1_alpha, D1_beta, D1_gamma, labels = c("", "",""),nrow=1)#ggpubr
D2_stitch <- egg::ggarrange(D2_alpha, D2_beta, D2_gamma, labels = c("", "",""),nrow=1)#ggpubr
N_stitch <- egg::ggarrange(count_map,cv_count_map,tot_count_map, labels = c("","",""),nrow=1)#ggpubr

## Stitch rows together
figure2 <- ggpubr::ggarrange(D0_stitch,D1_stitch,D2_stitch,N_stitch,nrow=4,font.label=list(color="black",size=16,face='plain'),align="v",widths = c(0.5,0.5,0.5, 1))#font.label=list(color="black",size=6,face='plain'),align="v") ,labels = c("D0", "D1", "D2","N")

## Add x and y labels
fig2 <- annotate_figure(
  figure2, 
  bottom = text_grob("          Longitude", 
  color = "black", face = "plain", size = 16),#,
  left = text_grob("Latitude", 
  color = "black", face = "plain", size = 16,rot = 90)
)

## Save combined plot
ggsave(plot = fig2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_2.png"),
       height = 350, width =280, units = "mm", dpi = 500,#height = 400, width =320, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285
#_______________________________________________________________________________
#### SPATIAL MODELLING OF BIODIVERSITY CLUSTERS: RASTER PREDICTOR VARIABLES (QUICK LOOK) ####

# This code is intended for a quick look-see. For final models and associated RF outputs (inc associated model confidence maps), use R script in file xx

## Load packages
library(raster)

## Load environmental predictor rasters
bathy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/bathy3.tif")
cur <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/cur3.tif")
gravel <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Predicted_Gravel_Fraction.tif")
light <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/light3.tif")
mud <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Predicted_Mud_Fraction.tif")
oxy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/oxy3.tif")
phyto <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/phyto3.tif")
sal <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/sal3.tif")
sil <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/sil3.tif")
spm <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/SPM_MEAN.tif")
temp <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/temp3.tif")
wov <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Wave_veloc.tif")

## Update crs
crs(bathy) <- "+proj=longlat +datum=WGS84 +no_defs" 

## Create raster stack
predictors <- stack(bathy,cur,gravel,light,mud,oxy,phyto,sal,sil,spm,temp,wov)
#predictors

# Simple names for predictor variables
names(predictors)=c("Bathymetry","Current","Gravel","Light","Mud","Oxygen","Phyto","Salinity","Silicate","SPM","Temperature","WOV")#"Sand",
#names(predictors)

## Plot raster stack
#plot(predictors)
#plot(bathy)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF BIODIVERSITY CLUSTERS: PREPARE DATA ####

## Data used for modelling
BiodivCluster <- faunal.cluster7
head(BiodivCluster)

## Save biodiv cluster data for passing on to AD
BiodivCluster2 <- BiodivCluster
BiodivCluster2$type <- 'categorical'
BiodivCluster2$paper <- 'biodiversity'
BiodivCluster2$metric <- 'biodiversity_cluster'
head(BiodivCluster2)
colnames(BiodivCluster2) <- c('sample','x','y','value','type','paper','metric')
BiodivCluster2 <- BiodivCluster2[,c(6,1:3,7,5,4)]
dim(BiodivCluster2)

## Save biodiv clusters for use with RF modelling script (see xx.R)
write.csv(BiodivCluster2, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\biodiv_cluster_4_modelling.csv", row.names=FALSE)

## Change names of cols
colnames(BiodivCluster)=c("Sample","lon","lat","cluster")
head(BiodivCluster)
dim(BiodivCluster)# 13654     4

## Extract variables from raster stack
sdata <- raster::extract(predictors, BiodivCluster[,2:3])
head(sdata)
#class(sdata)

## Change from matrix to df
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(BiodivCluster$Sample,BiodivCluster$cluster,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "Cluster"
head(sdata2)
#____________________________________________________________________________________________________________________
#### SPATIAL MODELLING OF BIODIVERSITY CLUSTERS: MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Call library
#install.packages("caTools")
library(caTools)

## Vector for response variable
Y <- sdata2[,2]

## Take a random 90% of samples in proportion across the 11 groups
set.seed(2)
msk= sample.split( Y, SplitRatio = 9/10, group = NULL )

## Check it was 90% 10% split

## The training set
train = sdata2[ msk,] 
#dim(train)#12287    14
#View(train)

## Remove station labels for train
train2 =train[,2:14]#was 11
#View(train2)

## The test set
test  = sdata2[ !msk,]
#dim(test)#962  11
#View(test)

## Remove station labels for test
test2 =test[,2:14]#was 11
#class(test2)
#View(test2)
#str(test2)

## Check number of observations for train (TRUE) and test (FALSE) sets
print(table(Y, msk)) 

## Check number of samples in train and test sets is equal to total
dim(sdata) # 13654    12
dim(train)+dim(test)# 13654    28
#____________________________________________________________________________________________________________________
#### SPATIAL MODELLING OF BIODIVERSITY CLUSTERS: DO MODELLING ####

## Call library
#install.packages("randomForest")
library(randomForest)

## Model
model <- factor(Cluster)~Bathymetry+Current+Gravel+Light+Mud+Oxygen+Phyto+Salinity+Silicate+SPM+Temperature+WOV# cluster (multiple biodiv metrics)

## Run model
train3 <- train2[complete.cases(train2), ]
#str(train3)
train3$Cluster <- as.numeric(train3$Cluster)
train3$Cluster <- as.factor(train3$Cluster)
#str(train3)
#head(train3)
rf2 <- randomForest(model, data=train3)
#____________________________________________________________________________________________________________________
#### SPATIAL MODELLING OF BIODIVERSITY CLUSTERS: EXAMINE HOW VARIABLES ARE INFLUENCING THE MODEL ####

## Produce plots 
varImpPlot(rf2)

#png('OUTPUTS/BIODIVERSITY/variables_affecting_model_structure.png') # height=nrow(pr), width=ncol(pr) 
#varImpPlot(rf2)
#dev.off()
#_______________________________________________________________________________
#### SPATIAL MODELLING OF BIODIVERSITY CLUSTERS: EVALUATE THE MODEL PERFORMANCE ####

## Predict cluster group for test set
pred <- predict(rf2, newdata = test2)
#table(pred, test2$Cluster)

## We can test the accuracy as follows:
(46+108+61+103+38+159+54+93)/ nrow(test2)# 48%

## Confusion matrix plot
#https://stackoverflow.com/questions/37897252/plot-confusion-matrix-in-r-using-ggplot

confusion_matrix <- as.data.frame(table(pred, test2$Cluster))

cm <- ggplot(confusion_matrix, aes(pred,sort(Var2,decreasing = T), fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#737373") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c("1","2","3","4","5","6")) +
  scale_y_discrete(labels=c("6","5","4","3","2","1"))

## Save confusion matrix
#png('C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_S4.png')
#cm
#dev.off()
#_______________________________________________________________________________
#### SPATIAL MODELLING OF BIODIVERSITY CLUSTERS: PRODUCE FULL COVERAGE RASTER ####

##Use model to predict cluster group for each raster cell
pr <- predict(predictors, rf2)

## Basic plot
plot(pr)

## Save raster
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Biodiversity_new.tif',overwrite=TRUE,format = "GTiff")
#_______________________________________________________________________________
#### MAP BIODIVERSITY CLUSTERS: POINT MAP PLOT ####

## Load packages
library(sf)
library(RColorBrewer)
library(terra)
library(tidyterra)
library(raster)
library(stars)
library(ggplot2)
library(ggpubr)

## Bring in countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))

## Set CRS
st_crs(countries) = 4326

## Results from clustering
head(faunal.cluster7)
str(faunal.cluster7)

## Change order of levels from highest to lowest biodiversity
faunal.cluster7$FaunalCluster <- factor(faunal.cluster7$FaunalCluster, levels= c('2','8','5','1','4','6','7','3'))
levels(faunal.cluster7$FaunalCluster)

## Point sample map
biodiv_point= ggplot()+
  geom_sf(data=countries, fill ="black",col ="black")+ 
  geom_point(data = faunal.cluster7,aes(x = X, y = Y,col=FaunalCluster), size = 0.1)+
  scale_colour_manual(breaks = c('2','8','5','1','4','6','7','3'), values = c(
    "#37C331", "#6FD326" ,"#A8E21B", "#E0F210", "#F3F223","#E2E256",  "#D1D189",  "#C0C1BC" ),
    labels = c('Bio-A','Bio-B','Bio-C','Bio-D','Bio-E','Bio-F','Bio-G','Bio-H'),
    na.value="transparent")+
  guides(colour = guide_legend(override.aes = list(size=3)))+#NEW
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw(base_size = 14)+#base_size = 24
  theme(legend.key = element_rect(colour = NA,fill = NA),legend.background=element_blank(),legend.text = element_text(color="white",size= 8))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(
    axis.title.x=element_blank(),
    axis.title.y = element_text(colour = "white"))+
  theme(legend.title= element_blank())+
  labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.4, "cm"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

biodiv_point
#_______________________________________________________________________________
#### MAP BIODIVERSITY CLUSTERS: MAP PLOT (FIGURE 4) ####

## Load raster layer
biodiv <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Bio (Combined Biodiversity)/Model/BiodiversityClusterMaxClass_Aug25.tif')

## Reduce raster resolution
biodiv.agg  <- aggregate(biodiv, fact=4,fun = modal)

## Make cluster a factor
values(biodiv.agg) <- as.factor(values(biodiv.agg))

## Produce map
pbio <-ggplot() +
  geom_spatraster(data = biodiv.agg) +
   scale_fill_manual(breaks = c('2','8','5','1','4','6','7','3'), values = c(
     "#37C331", "#6FD326" ,"#A8E21B", "#E0F210", "#F3F223","#E2E256",  "#D1D189",  "#C0C1BC" ),
     labels = c('Bio-A','Bio-B','Bio-C','Bio-D','Bio-E','Bio-F','Bio-G','Bio-H'),
     na.value="transparent")+
   geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 10)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.21))+
  theme(axis.title.x = element_text(colour = "black"))+
  labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.8, "cm"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pbio

## Save plot
ggsave(plot = pbio,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_4.png"),
       height = 195, width =195,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
#### MAP BIODIVERSITY CLUSTERS CONFIDENCE: MAP PLOT (SUPPLEMENTARY FIGURE S5) ####

## Load biodiversity cluster confidence raster
biodiv_cluster_conf <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Bio (Combined Biodiversity)/Confidence/BiodiversityClusterConfidence_Aug25.tif')

## Reduce resolution of confidence layer
biodiv_cluster_conf.agg  <- aggregate(biodiv_cluster_conf, fact=5,fun = modal)

## Remove raster to save space
rm(biodiv_cluster_conf)

## Produce map
pbio_conf <- ggplot() +
  geom_spatraster(data = biodiv_cluster_conf.agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "Confidence")+#Oranges
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw(base_size = 10)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.2))+
  theme(legend.title = element_text(color = "white", size = 12))+
  labs(fill = "Confidence")+
  theme(axis.title.x = element_text(colour = "white"),
        axis.title.y = element_text(colour = "white"))+
  labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.8, "cm"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pbio_conf

## Save Figure S5
ggsave(plot = pbio_conf,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_S5.png"),
       height = 195, width =195,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
#### MAP RARE TAXA WITHIN 3 MOST DIVERSE CLUSTER GROUPS (SUPPLEMENTARY FIGURE S6) ####

## Load packages
library(terra)
library(tidyterra)
library(colorRamps)
library(ggplot2)
library(sf)
library(ggpubr)

## Bring in countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))

## Load biodiversity cluster raster
biodiv <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Bio (Combined Biodiversity)/Model/BiodiversityClusterMaxClass_Aug25.tif')

## Plot raster
plot(biodiv)
class(biodiv)# Spatraster

## Create a copy of biodiv raster
biodiv2 <- biodiv

## Remove unwanted groups. Top 3 groups are 2, 8 and 5
biodiv2[biodiv2 == 1] <- NA
biodiv2[biodiv2 == 4] <- NA
biodiv2[biodiv2 == 6] <- NA
biodiv2[biodiv2 == 7] <- NA
biodiv2[biodiv2 == 3] <- NA
plot(biodiv2)

## Convert raster to polygon
p2 = as.polygons(biodiv)
class(p2)

## Set colours
cols2 <- c( "#E0F210",
            "#37C331",
            "#C0C1BC",
            "#F3F223",
            "#A8E21B",
            "#E2E256",
            "#D1D189",
            "#6FD326"
           )#ordered 1:8#
plot(p2, col=cols2) 

## Subset for top 3 biodiv clusters
rare_1 <- subset(p2, p2$BiodiversityClusterMaxClass_Aug25 == 2)
rare_2 <- subset(p2, p2$BiodiversityClusterMaxClass_Aug25 == 8)
rare_6 <- subset(p2, p2$BiodiversityClusterMaxClass_Aug25 == 5)
rare_top3 <- subset(p2, p2$BiodiversityClusterMaxClass_Aug25 == 2 | p2$BiodiversityClusterMaxClass_Aug25 == 8 | p2$BiodiversityClusterMaxClass_Aug25 == 5)
plot(rare_top3)

## Load rare species layer
g_rare <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/γ rare/Model/g_rare_Mean_SQ_2025.tif')
plot(g_rare)

## Plot the rare taxa layer
ggplot() +
  geom_spatraster(data = g_rare) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 20)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(axis.title.x = element_text(colour = "black"))+
  labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.8, "cm"))+
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(colour = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor = element_blank())

## Cut out a geographic subset
rare_top3_crop <- terra::crop(g_rare, rare_top3, mask=TRUE)

## Plot rare taxa within top 3 biodiv clusters
rare_top3_crop_map<-ggplot() +
  geom_spatraster(data = rare_top3_crop) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 20)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 18))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(axis.title.x = element_text(colour = "black"))+
  labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0.2,0,1), "cm"),legend.key.size = unit(1, "cm"))+#t, r, b, l
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.x=element_text(colour = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor = element_blank())

rare_top3_crop_map

## Make cluster a factor
values(biodiv2) <- as.factor(values(biodiv2))

## Raster map of top 3 clusters (2,8, 5 or Bio-A, Bio-B, Bio-C)
rare_top3_map <-ggplot() +
  geom_spatraster(data = biodiv2) +
  scale_fill_manual(breaks = c('2','8','5','1','4','6','7','3'), values = c(
     "#37C331", "#6FD326" ,"#A8E21B", "#E0F210", "#F3F223","#E2E256",  "#D1D189",  "#C0C1BC" ),
     labels = c('Bio-A','Bio-B','Bio-C','Bio-D','Bio-E','Bio-F','Bio-G','Bio-H'),
     na.value="transparent")+ 
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 20)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 18))+
  theme(legend.position=c(0.9,0.21))+
  theme(axis.title.x = element_text(colour = "black"))+
  labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0.4,0,1), "cm"),legend.key.size = unit(1, "cm"))+#t, r, b, l
  theme(panel.grid.major = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        #axis.text.x=element_text(colour = "white"),
        panel.grid.minor = element_blank())

rare_top3_map

## Stitch above 2 plots into a single row
sep_cluster_stitch <- egg::ggarrange(rare_top3_map,rare_top3_crop_map, labels = c("a)","b)"),nrow=1,label.args = list(gp = grid::gpar(font = 4, cex =
2)))#ggpubr

## Add annotations
figureS5 <- annotate_figure(
  sep_cluster_stitch, 
  #top = text_grob("     Cluster                                                          % Rare", 
  #color = "black", face = "bold", size = 24),
  bottom = text_grob("          Longitude", 
  color = "black", face = "plain", size = 24),#,
  left = text_grob("Latitude", 
  color = "black", face = "plain", size = 24,rot = 90)
)

## Save final plot
ggsave(plot = figureS5,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_S6.png"),
       height = 260, width =520, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285

#_______________________________________________________________________________
#### IDENTIFY SPATIAL RESOLUTION OF MODELS ####

#Load packages
library(raster)
library(terra)
library(tidyterra)
library(viridis)
library(ggplot2)
library(colorRamps)
library(ggpubr)
require(grid)

## Load rasters
gravel<- rast("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Predicted_Gravel_Fraction.TIF")

## Get approximate longitude of raster center
lon_center <- (xmin(gravel) + xmax(gravel)) / 2

## Compute UTM zone automatically
utm_zone <- floor((lon_center + 180) / 6) + 1

## Build UTM CRS string
utm_crs <- paste0("+proj=utm +zone=", utm_zone, " +datum=WGS84 +units=m")

## Create a template raster in projected CRS
# We'll keep roughly the same resolution as the original (converted to meters)
# Approximate conversion: 1 degree ≈ 111,000 meters
res_m <- res(gravel) * 111000  # rough conversion for lon/lat to meters

## convert SpatRaster → RasterLayer
gravel_r <- raster(gravel)  # convert SpatRaster → RasterLayer
r_template <- raster(extent(gravel_r), crs=utm_crs)
res(r_template) <- res_m

## Reproject raster safely
# Use 'bilinear' for continuous data, 'ngb' for categorical
r_utm <- projectRaster(gravel_r, r_template, method="bilinear")

## Check the resolution in meters
res(r_utm)
#_______________________________________________________________________________
#### RASTER PREDICTORS 'BEST' EXPLINING BIODIVERSITY CLUSTER PATTERNS (SUPPLEMENTARY FIGURE S7) ####

#Load packages
library(raster)
library(terra)
library(tidyterra)
library(viridis)
library(ggplot2)
library(colorRamps)
library(ggpubr)
require(grid)

## Load rasters
Current_Sp<- rast("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Current_Sp.tif")
gravel<- rast("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Predicted_Gravel_Fraction.TIF")

## Update names
names(gravel)<-'Gravel'

## Make gravel data a percentage
gravel=gravel*100

## Define breaks
breaks <- c(0, 5, 10, 20, 30, 40,50,60,70,80,90, 100)

## Plot label
label1 = "Gravel"#"~alpha~ q == 1"

## Gravel raster plot
gravel_rast <- ggplot() +
  geom_spatraster(data = gravel) +
  scale_fill_viridis_b(breaks = breaks, option = "C", na.value = "transparent") +
  coord_sf(xlim = c(-10, 9), ylim = c(49, 60)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw(base_size = 30) +
  theme(legend.background=element_blank(),legend.text = element_text(color="black",size= 18))+#make legend background transparent and text white
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(1,"cm"))+#t, r, b, l
  theme(legend.position=c(0.87,0.2))+#  # Adjust the position as needed
  theme(panel.grid.major = element_blank(),
        axis.title.y=element_blank(),
                     panel.grid.minor = element_blank())+
   
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
      annotate(geom = "text",
           label = label1,
           #parse = TRUE,
           size=12,
          colour = "#707070",
          # fontface="bold",
           x = -9,
           y = 60.15
          )+
  guides(fill=guide_legend(title="%"))
gravel_rast

## Define breaks for Current speed
breaks3 <- c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6)

## Plot label
label3 = "Current speed"#"~alpha~ q == 1"

## Current speed raster plot
Current_Sp_rast <- ggplot() +
  geom_spatraster(data = Current_Sp) +
  scale_fill_viridis_b(breaks = breaks3, option = "C", na.value = "transparent") +
  coord_sf(xlim = c(-10, 9), ylim = c(49, 60)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw(base_size = 30) +
   theme(legend.background=element_blank(),legend.text = element_text(color="black",size= 18))+#make legend background transparent and text white
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(1,"cm"))+#t, r, b, l
  theme(legend.position=c(0.87,0.22))+#  # Adjust the position as needed
  theme(
    axis.title.y=element_blank(),
      axis.text.y=element_text(colour = "white"))+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  annotate(geom = "text",
           label = label3,
           #parse = TRUE,
           size=12,
          colour = "#707070",
          # fontface="bold",
           x = -8,
           y = 60.15
          )+
  guides(fill=guide_legend(title=expression(m~s^{-1})))
Current_Sp_rast

## Stitch plots together
var_rast <-ggpubr::ggarrange(gravel_rast+ rremove("xlab"),Current_Sp_rast+ rremove("xlab"), ncol = 2, nrow = 1,
          labels = c("", "",""),
          font.label = list(size = 30, color = "black"))

## Add annotations
var_rast2 <- annotate_figure(var_rast, 
                    bottom = textGrob("Longitude", gp = gpar(cex =2)),
                    left = textGrob("Latitude", rot = 90,gp = gpar(cex =2)))

## Save plot
ggsave(plot = var_rast2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_S7.png"),
       height = 350, width =700, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285
#_______________________________________________________________________________
####  CONFIDENCE MAPS FOR BIODIV METRIC PLOTS FROM FIGURE 2 (SUPPLEMENTARY FIGURE S6) 

## Load packages
library(sf)
library(dplyr)
library(RColorBrewer)
library(terra)
library(raster)
library(stars)
library(ggplot2)
library(colorRamps)
library(terra)
library(tidyterra)
library(raster)
library(dplyr)
library(stars)
library(ggplot2)
library(colorRamps)
library(climateStability)
library(ggpubr)

## Bring in countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D0 Alpha confidence ####

## Load model raster
a_q0_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/0Dα (Hill 0, Alpha)/Confidence/sample_a_q0_CV_SQ_2025.tif')

## Reduce size of raster
a_q0_model_agg <- aggregate(a_q0_model, fact=7,fun = modal)

## Remove raster to save space
rm(a_q0_model)

## Annotate label
label1 =  "{}^{0}~italic(D)[italic(α)]" #

## Model plot
D0_alpha <-  ggplot() +
  geom_spatraster(data = a_q0_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    panel.background = element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D0 Beta confidence ####

## Load raster layer
b_q0_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/0Dβ (Hill 0, Beta)/Confidence/b_q0_CV_SQ_2025.tif')

## Reduce size of raster
b_q0_model_agg <- aggregate(b_q0_model, fact=7,fun = modal)

## Remove raster to save space
rm(b_q0_model)

## Annotate label
label1 = "{}^{0}~italic(D)[italic(β)]"# "~beta~ q == 0"

## Model plot
D0_beta <- ggplot() +
  geom_spatraster(data = b_q0_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
   x = -8.3,
           y = 59.7)+
theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D0 Gamma confidence ####

## Load raster layer
g_q0_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/0Dγ (Hill 0, Gamma)/Confidence/g_q0_CV_SQ_2025.tif')

## Reduce size of raster
g_q0_model_agg <- aggregate(g_q0_model, fact=7,fun = modal)

## Remove raster to save space
rm(g_q0_model)

## Annotate label
label1 ="{}^{0}~italic(D)[italic(γ)]" # "~gamma~ q == 0"

## Model plot
D0_gamma <- ggplot() +
  geom_spatraster(data = g_q0_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D1 alpha confidence ####

## Load raster layer
a_q1_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/1Dα (Hill 1, Alpha)/Confidence/sample_a_q1_CV_SQ_2025.tif')

## Reduce size of raster
a_q1_model_agg <- aggregate(a_q1_model, fact=7,fun = modal)

## Remove raster to save space
rm(a_q1_model)

## Annotate label
label1 = "{}^{1}~italic(D)[italic(α)]"#"~alpha~ q == 1"

## Model plot
D1_alpha <- ggplot() +
  geom_spatraster(data = a_q1_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D1 beta confidence ####

## Load model raster D1 beta
b_q1_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/1Dβ (Hill 1, Beta)/Confidence/b_q1_CV_SQ_2025.tif')

## Reduce size of raster
b_q1_model_agg <- aggregate(b_q1_model, fact=7,fun = modal)

## Remove raster to save space
rm(b_q1_model)

## Annotate label
label1 = "{}^{1}~italic(D)[italic(β)]" # "~beta~ q == 1"

## Model plot
D1_beta <- ggplot() +
  geom_spatraster(data = b_q1_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D1 gamma confidence ####

## Load raster
g_q1_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/1Dγ (Hill 1, Gamma)/Confidence/g_q1_CV_SQ_2025.tif')

## Reduce size of raster
g_q1_model_agg <- aggregate(g_q1_model, fact=7,fun = modal)

## Remove raster to save space
rm(g_q1_model)

## Annotate label
label1 = "{}^{1}~italic(D)[italic(γ)]"# "~gamma~ q == 0"

## Model plot
D1_gamma <- ggplot() +
  geom_spatraster(data = g_q1_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D2 alpha confidence ####

## Load raster
a_q2_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/2Dα (Hill 2, Alpha)/Confidence/sample_a_q2_CV_SQ_2025.tif')

## Reduce size of raster
a_q2_model_agg <- aggregate(a_q2_model, fact=7,fun = modal)

## Remove raster to save space
rm(a_q2_model)

## Annotate label
label1 = "{}^{2}~italic(D)[italic(α)]" #"~alpha~ q == 2"

## Model plot
D2_alpha <- ggplot() +
  geom_spatraster(data = a_q2_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D2 beta confidence ####

## Load model raster D2 beta
b_q2_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/2Dβ (Hill 2, Beta)/Confidence/b_q2_CV_SQ_2025.tif')

## Reduce size of raster
b_q2_model_agg <- aggregate(b_q2_model, fact=7,fun = modal)

## Remove raster to save space
rm(b_q2_model)

## Annotate label
label1 ="{}^{2}~italic(D)[italic(β)]" # "~beta~ q == 2"

## Model plot
D2_beta <- ggplot() +
  geom_spatraster(data = b_q2_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: D2 gamma confidence ####

## Load raster
g_q2_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/2Dγ (Hill 2, Gamma)/Confidence/g_q2_CV_SQ_2025.tif')

## Reduce size of raster
g_q2_model_agg <- aggregate(g_q2_model, fact=7,fun = modal)

## Remove raster to save space
rm(g_q2_model)

## Annotate label
label1 = "{}^{2}~italic(D)[italic(γ)]" #"~gamma~ q == 2"

## Model plot
D2_gamma <- ggplot() +
  geom_spatraster(data = g_q2_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: N AV confidence ####

## Load raster
av_count_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/N (Abund)/Confidence/av_count_CV_SQ_2025.tif')

## Reduce size of raster
av_count_model_agg <- aggregate(av_count_model, fact=7,fun = modal)

## Remove raster to save space
rm(av_count_model)

## Annotate label
#label1 = "italic(N[av])"##"N~alpha" "~av~count"
label1 = "italic(N)"

## Model plot
av_count_map <- ggplot() +
  geom_spatraster(data = av_count_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: N CV confidence ####

## Load raster
cv_count_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Ncv (Abund cv)/Confidence/cv_count_CV_SQ_2025.tif')

## Reduce size of raster
cv_count_model_agg <- aggregate(cv_count_model, fact=7,fun = modal)

## Remove raster to save space
rm(cv_count_model)

## Annotate label
label1 = "italic(N[cv])"#"~cv~count"

## Model plot
cv_count_map <- ggplot() +
  geom_spatraster(data = cv_count_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: N TOT confidence ####

## Load raster
tot_count_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Ntot (Abund total)/Confidence/tot_count_CV_SQ.tif')

## Reduce size of raster
tot_count_model_agg <- aggregate(tot_count_model, fact=7,fun = modal)

## Remove raster to save space
rm(gamma_model)

## Annotate label
label1 = "italic(N[tot])"#"~tot~count"

## Model plot
tot_count_map <- ggplot() +
  geom_spatraster(data = tot_count_model_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "CV")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 14)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  #theme(legend.title= element_blank())+
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=6.3,
           colour = "#707070",
           fontface="bold",
           x = -8.3,
           y = 59.7)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP BIODIVERSITY METRICS CONFIDENCE: STITCH CONFIDENCE PLOTS TOGETHER (FIGURE S3) ####

## Create each row
D0_stitch <- egg::ggarrange(D0_alpha, D0_beta, D0_gamma, labels = c("", "",""),nrow=1)#ggpubr
D1_stitch <- egg::ggarrange(D1_alpha, D1_beta, D1_gamma, labels = c("", "",""),nrow=1)#ggpubr
D2_stitch <- egg::ggarrange(D2_alpha, D2_beta, D2_gamma, labels = c("", "",""),nrow=1)#ggpubr
N_stitch <- egg::ggarrange(av_count_map,cv_count_map,tot_count_map, labels = c("","",""),nrow=1)#ggpubr

## Combine rows
figures3 <- ggpubr::ggarrange(D0_stitch,D1_stitch,D2_stitch,N_stitch,nrow=4,font.label=list(color="black",size=16,face='plain'),align="v",widths = c(0.5,0.5,0.5, 1))

## Annotate
figs3 <- annotate_figure(
  figures3, 
  bottom = text_grob("          Longitude", 
  color = "black", face = "plain", size = 16),#,
  left = text_grob("Latitude", 
  color = "black", face = "plain", size = 16,rot = 90)
)

## Save plot
ggsave(plot = figs3,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\OUTPUTS\\Figure_S3.png"),
       height = 400, width =320, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285
