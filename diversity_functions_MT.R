# Loading required packages
pkgs = c("geosphere", "lubridate", "sf", "tidyverse", 
         "mapplots", "purrr", "iNEXT", 
         "RColorBrewer", "ggpubr", "parallel", 
         "janitor")

invisible(lapply(pkgs, library, character.only = TRUE))

## Loading the data tidying function
tidy_benthic_dat <- function(dat = sppdf, firstyear = min_yr, lastyear = max_yr){#firstyear = 1981, lastyear = 2020
  
  
  
  ## Tidying the data (renaming variables, grouping, filtering, adding dates)
  spp <- dat %>% 
    rename(sample = samplecode,
           latitude = samplelat,
           longitude = samplelong,
           species = validname,
           #species = family,
           count = abund) %>%
    clean_names() %>% 
    group_by(sample, latitude, longitude, year,month,date, species) %>%
    summarise(count = sum(count)) %>%
    filter(year %in% c(min(years-1): max(years+1))) %>%
    ungroup() %>% 
    mutate(date = as.Date(date, "%Y-%m-%d"))
}

get_samples <- function(tidy_sppdf, yr){
  
  # All samples to estimate from inc. years either side
  samples_yr  <- tidy_sppdf %>%
    filter(year >= yr - 1 & year <= yr + 1)
  
  # If data can be used
  if(nrow(samples_yr) > 0) {
    
    # Pulling samples 
    samp_samps <- samples_yr %>% 
      distinct(sample) %>%
      pull()
    
    # Defining focal locations for year of interest
    samp_locs <- tidy_sppdf %>% 
      filter(year == yr) %>%
      dplyr::select(longitude, latitude, sample, date) %>% 
      distinct(sample) %>%
      pull()
    
    # Define all points to generate info from
    comp_1_dist <- tidy_sppdf %>% 
      filter(sample %in% unique(c(samp_locs, samp_samps))) %>%
      dplyr::select(longitude, latitude, sample, date) %>% 
      distinct() 
    
    # Making a list 
    comp_1_lon_lat <- comp_1_dist %>%
      rownames_to_column() %>%  
      dplyr::select(-c(date)) %>%
      pivot_longer(c(-rowname, -sample)) %>% 
      dplyr::select(-rowname) %>% 
      pivot_wider(names_from = sample, values_from = value) %>%
      dplyr::select(-name) %>% 
      as.list()
    
    # Creating sample dates dataframe
    sample_dates <- comp_1_dist %>%
      distinct(sample, date) 
    
    # List of dataframes to return from the function
    ret <- list(
      focal_year = yr,
      sample_dates = sample_dates,
      samp_locs = samp_locs,
      comp_1_dist = comp_1_dist,
      comp_1_lon_lat = comp_1_lon_lat,
      samp_samps = samp_samps
    )
    
  } else {
    # Creating a NULL object to return if samples weren't found
    ret <- NULL
    print(paste(yr, 'missing temporal data', sep=' '))
  }
  
  ret
}

# Pairwise differences with sub-sampling option----
get_pairwise_diffs <- function(tidy_sppdf, yr, n_days, radius, n_samp, sub_sample){
  
  
  # Run the samples_yr function to get the necessary set-up for the focal year
  samples_yr <- get_samples(tidy_sppdf, yr = yr)
  
  # Create a distance matrix
  dist_days <- difftime(samples_yr$sample_dates$date,
                        as.Date("2000/01/01"), units = 'days') %>% 
    dist(diag = TRUE, upper = TRUE) %>% 
    as.matrix() %>% 
    as_tibble() %>% 
    setNames(samples_yr$sample_dates$sample) %>% 
    mutate(sample = samples_yr$sample_dates$sample) %>% 
    dplyr::select(sample, everything())
  
  # Removing duplicates 
  dist_days[upper.tri(dist_days)]  <- NA
  
  # Pivoting data
  samps_days  <- dist_days %>%
    pivot_longer(cols = -c(sample),
                 names_to = "sample2",
                 values_to = "time_days") %>%
    filter(!is.na(time_days))
  
  # Checking locs are in sample and samps in sample2
  samps_in_date1 <- samps_days %>%
    filter(sample %in% samples_yr$samp_locs &
             sample2 %in% samples_yr$samp_samps) %>%
    mutate(sample_locs = sample,
           sample_samps = sample2)
  samps_in_date2 <- samps_days %>%
    filter(sample %in% samples_yr$samp_samps &
             sample2 %in% samples_yr$samp_locs) %>%
    mutate(sample_locs = sample2,
           sample_samps = sample)
  
  # Identify samps in required timeframe
  samps_in_date <-  samps_in_date1 %>%
    bind_rows(samps_in_date2) %>%
    filter(time_days <= n_days) %>%
    dplyr::select(-c(sample, sample2))
  
  # List of all samps within timeframe to generate spatial info
  samps <- unique(c(samps_in_date$sample_locs, samps_in_date$sample_samps))
  
  sub_comp_1_lon_lat  <- samples_yr$comp_1_lon_lat[samps]
  
  # Pairwise distances of all samps within timeframe 
  pairwise_distances <- map_df(sub_comp_1_lon_lat, function(x) {
    purrr::map(sub_comp_1_lon_lat, distCosine, x)}) %>%
    mutate(sample = names(sub_comp_1_lon_lat)) %>%
    dplyr::select(sample, everything())
  
  # Removing duplicates 
  pairwise_distances[upper.tri(pairwise_distances)] <- NA
  
  # Pivoting data
  samps_dists <-  try(pairwise_distances %>% 
                        pivot_longer(cols= -c(sample),
                                     names_to = "sample2",
                                     values_to = "distance_m") %>%
                        filter(!is.na(distance_m)), silent=T)
  
  if(is(samps_dists, "try-error") == FALSE) {
    
    # Checking locs are in sample and samps in sample2
    samps_in_region1 <- samps_dists %>%
      filter(sample %in% samples_yr$samp_locs &
               sample2 %in% samples_yr$samp_samps) %>%
      mutate(sample_locs = sample,
             sample_samps = sample2)
    samps_in_region2 <- samps_dists %>%
      filter(sample %in% samples_yr$samp_samps &
               sample2 %in% samples_yr$samp_locs) %>%
      mutate(sample_locs = sample2,
             sample_samps = sample)
    
    
    # Subset samples within radius
    samps_in_dist <- samps_in_region1 %>%
      bind_rows(samps_in_region2) %>%
      filter(distance_m < radius,
             distance_m > 0) %>%
      dplyr::select(-c(sample, sample2)) 
    
    # Identifying all samps in time-space window
    samps_in_dist_days <- samps_in_dist %>%
      left_join(samps_in_date, c('sample_locs', 'sample_samps')) %>%
      arrange(sample_locs) %>%
      filter(!is.na(time_days))  %>%
      distinct()  
    
    ret <- samps_in_dist_days
    
    # To randomly sub-sample based on n_samp
    if(sub_sample == 'Y') {
      
      # Create data.frame for sub-sampling using n_samp
      sub_samps_in_dist_days <- data.frame(sample_locs=NA,
                                           sample_samps=NA,
                                           distance_m=NA,
                                           time_days=NA)
      
      for(loc in unique(samps_in_dist_days$sample_locs)){ #loc = unique(samps_in_dist_days$sample_locs)[1]
        
        df <-  subset(samps_in_dist_days, sample_locs == loc)
        
        samps <-  df$sample_samps
        
        set.seed(21)
        
        samps <- samps %>% sample(min(length(unique(df$sample_samps)), n_samp)) 
        
        sub_samps_in_dist_days <- sub_samps_in_dist_days %>%
          bind_rows(subset(df, sample_locs == loc & sample_samps %in% samps)) %>%
          filter(!is.na(sample_locs))
      }
      
      samps_in_dist_days <- sub_samps_in_dist_days
      
    } else{ 
      
      samps_in_dist_days <- samps_in_dist_days
      
    }
    
    # You may want to subset samples further based on n_samp in area using sub_samps_in_dist_days
    samps_rad_info <- samps_in_dist_days %>%
      group_by(sample_locs) %>%
      summarise(n_samps = length(unique(sample_samps))) %>%
      as.data.frame()
    
    # All samples to estimate from inc. years either side
    samples_yr <- tidy_sppdf %>%
      subset(year==yr | 
               (year==yr+1) |
               (year==yr-1)) %>%
      mutate(month = month(date)) %>%
      ungroup() 
    
    # Format data, and include loc to estimate gamma diversity
    divdf <- samples_yr %>%
      filter(sample %in% unique(c(samps_in_dist_days$sample_locs, samps_in_dist_days$sample_samps))) %>%
      group_by(sample, species) %>%
      summarise(count = sum(count)) %>%
      filter(count >0) %>%
      pivot_wider(names_from = sample,
                  values_from = count,
                  values_fn = sum,
                  values_fill = 0) %>%
      as.data.frame()
    rownames(divdf) = divdf[,1]
    divdf = divdf[,-1]
    
    # List of dataframes to return from function
    returnlist <- list("divdf" = divdf, 
                       "samples_yr" = samples_yr, 
                       "samps_rad_info" = samps_rad_info,
                       "samps_in_dist_days" = samps_in_dist_days) 
    
    return(returnlist)
    
    
  } else {
    ret <- NULL
  }
  
  ret
}


div_est <-function(tidy_sppdf, yr, n_days, radius, n_samp, sub_sample){
  
  ### Alpha diversity estimation for all samples 
  
  # Format data, and include loc to estimate gamma diversity
  a_divdf <- tidy_sppdf %>%
    filter(year %in% yr) %>%
    group_by(sample, species) %>%
    summarise(count = round(sum(count))) %>%
    #summarise(count = sum(count)) %>%
    filter(count >0) %>%
    pivot_wider(names_from = sample,
                values_from = count,
                values_fn = sum,
                values_fill = 0) %>%
    as.data.frame()
  rownames(a_divdf) = a_divdf[,1]
  a_divdf = a_divdf[,-1]
  
  # Sampling depth for rarefaction = 2*median sum of sample count across all data
  bt <- tidy_sppdf %>%
    group_by(sample) %>% 
    summarise(total = sum(count),
              nspp = length(unique(species))) %>% 
    #ungroup() %>% 
    summarise(b = round(median(total)*2),
              t = round(max(total)))
  
  # alpha diversity estimation  
  myest_yr <- iNEXT(a_divdf, q=c(0,1,2), size=c(bt$b,bt$t), 
                    datatype="abundance", nboot = 50) 
  
  # sample-level alpha estimates
  # formatted to average
  sample_div_yr <- myest_yr$iNextEst$size_based %>%
    filter(m == bt$b)
  
  # formatted for join as sample-based estimates
  a_samp_dat <- sample_div_yr %>%
    rename(sample=Assemblage) %>%
    mutate(hill_n = case_when(Order.q==0 ~ 'sample_a_q0',
                              Order.q==1 ~ 'sample_a_q1',
                              Order.q==2 ~ 'sample_a_q2')) %>%
    dplyr::select(sample, qD, hill_n) %>%
    pivot_wider(names_from = hill_n,
                values_from = qD)
  
  samples <- get_pairwise_diffs(tidy_sppdf, 
                                n_days=n_days, 
                                radius=radius, 
                                n_samp=n_samp, 
                                yr=yr, 
                                sub_sample=sub_sample)
  
  # All samples to estimate from inc. years either side
  samples_yr  <- tidy_sppdf %>%
    filter(year >= yr - 1 & year <= yr + 1)
  
  # Creating a dataframe to store the values
  samps_rad_info_yr <- tibble(year = yr)
  
  # n samples to extrapolate to when estimating gamma diversity
  gamma_n_samp <- (n_samp+1)*2
  
  # Gather samples within radius, calculate sum of distances and days
  for(i in 1:nrow(samples$samps_rad_info)){#i=159
    loc <- unique(samples$samps_rad_info$sample_locs)[i]
    samps_in_region <- unique(subset(samples$samps_in_dist_days, sample_locs == loc)$sample_samps)
    # Sample distance and time clustering
    dddat <- samples$samps_in_dist_days %>%
      filter(sample_locs %in% loc,
             sample_samps %in% samps_in_region)
    
    # Sample-level abundance
    dat <- samples_yr %>%
      filter(sample %in% c(loc, samps_in_region)) %>%
      group_by(sample) %>%
      summarise(count = sum(count)) %>%
      ungroup() %>%
      summarise(av_count = mean(count),
                sd_count = sd(count),
                # Coefficient of variation (a type of beta-diversity for abundance)
                cv_count = sd(count)/mean(count),
                tot_count = sum(count)) %>%
      mutate(sample = loc)
    
    # Alpha
    a_dat <- sample_div_yr %>%
      filter(Assemblage %in% c(loc, samps_in_region)) %>%
      group_by(Order.q) %>%
      summarise(av_alpha = mean(qD)) %>%
      mutate(hill_n = case_when(Order.q==0 ~ 'q0',
                                Order.q==1 ~ 'q1',
                                Order.q==2 ~ 'q2')) %>%
      dplyr::select(-Order.q) %>%
      pivot_wider(names_from = hill_n,
                  values_from = av_alpha) 
    
    # Add to created dataframe
    samples$samps_rad_info[i, 'sum_km_between_samps'] = round(sum(dddat$distance_m)/1000, 1)
    samples$samps_rad_info[i, 'sum_time_between_samps'] = sum(dddat$time_days)
    samples$samps_rad_info[i, 'av_count'] = dat$av_count
    samples$samps_rad_info[i, 'sd_count'] = dat$sd_count
    samples$samps_rad_info[i, 'cv_count'] = dat$cv_count 
    samples$samps_rad_info[i, 'tot_count'] = dat$tot_count
    samples$samps_rad_info[i, 'a_q0'] = a_dat$q0
    samples$samps_rad_info[i, 'a_q1'] = a_dat$q1
    samples$samps_rad_info[i, 'a_q2'] = a_dat$q2
    samples$samps_rad_info[i, 'alpha_n'] = bt$b ## new
    samples$samps_rad_info[i, 'gamma_n_samp'] = gamma_n_samp ## new
    if(length(samps_in_region)>1) {
      
      # Incidence matrix for sample-based rarefaction
      g_divdf <- samples$divdf[,c(loc, samps_in_region)]
      g_divdf <- g_divdf[rowSums(g_divdf) > 0, ]
      g_divdf[g_divdf>0] = 1
      
      samples$samps_rad_info[i, 'g_taxa_count'] = nrow(g_divdf) ##
      samples$samps_rad_info[i, 'year'] = yr ##
      
      # G diversity estimation using sample-based rarefaction
      gam_samp_yr <- try(iNEXT(list(g_divdf), q=c(0,1,2), datatype="incidence_raw", 
                               size=c(5, n_samp+1, gamma_n_samp), nboot = 50), silent = T)
      
      if(is(gam_samp_yr,"try-error")==F) {
        
        g_dat <- gam_samp_yr$iNextEst$size_based %>%
          filter(t == gamma_n_samp) %>%
          mutate(hill_n = case_when(Order.q == 0 ~ 'q0',
                                    Order.q == 1 ~ 'q1',
                                    Order.q == 2 ~ 'q2')) %>%
          dplyr::select(hill_n, qD) %>%
          pivot_wider(names_from = hill_n,
                      values_from = qD) 
        
        
        if(nrow(g_dat)>0) {
          
          samples$samps_rad_info[i, 'g_q0'] = g_dat$q0
          samples$samps_rad_info[i, 'g_q1'] = g_dat$q1
          samples$samps_rad_info[i, 'g_q2'] = g_dat$q2 
          
        } else {
          
          samples$samps_rad_info[i, 'g_q0'] = NA
          samples$samps_rad_info[i, 'g_q1'] = NA
          samples$samps_rad_info[i, 'g_q2'] = NA 
          
          print('missing spatial samples for gamma estimates')
          
        }
        
      } else {
        
        samples$samps_rad_info[i, 'g_q0'] = NA
        samples$samps_rad_info[i, 'g_q1'] = NA
        samples$samps_rad_info[i, 'g_q2'] = NA 
        
        print('insufficient species occurrences for gamma estimate')
        
      }
      
      # Beta diversity estimation
      samples$samps_rad_info <- samples$samps_rad_info %>%
        mutate(b_q0 = if_else(is.na(g_q0)| is.na(a_q0), NA, g_q0/a_q0),
               b_q1 = if_else(is.na(g_q1)| is.na(a_q1), NA, g_q1/a_q1),
               b_q2 = if_else(is.na(g_q2)| is.na(a_q2), NA, g_q2/a_q2)) 
      
      samps_rad_info_yr <- samples$samps_rad_info
      
    }
    
  }
  
  
  # Store additional data with info on sample coverage, confidence intervals etc that may be useful down the line
  save(gam_samp_yr, myest_yr, file=paste0(outpath, 'biodiversity_estimates_additional_info_',yr,'_',file_name,'.Rdata'))
  
  # Store the data we want to use
  save(samps_rad_info_yr, a_samp_dat, file=paste0(outpath, 'biodiversity_estimates_',yr,'_',file_name,'.Rdata')) 
  
  return(glimpse(samps_rad_info_yr))
  
}


all_func <- function(yr){
  
  # Biodiversity estimations function
  div_est(tidy_sppdf, yr=yr, n_samp=n_samp, n_days=n_days, radius=radius, sub_sample=sub_sample)
  
}


# Parallel processing function
parallel_process <- function(){
  
  # Calculate the number of cores
  no_cores <- no_cores
  
  # Initiate cluster:
  cl <- makeCluster(no_cores)
  # Export packages:
  clusterEvalQ(cl, {   library(iNEXT); library(tidyverse); library(purrr); library(ggpubr);
    library(dplyr); library(lubridate); library(geosphere); library(janitor)})  
  # Export variables:
  clusterExport(cl, varlist=c('tidy_sppdf', 'years', 'all_func', 'get_samples',  
                              'get_pairwise_diffs', 'div_est',
                              'n_samp', 'sub_sample', 'radius', 'n_days',
                              'no_cores', 'outpath', 'file_name', 'years'))  
  #ll <- parLapply(cl, min(years):max(years), function(yr) all_func(yr)) 
  ll <- parLapply(cl, years, function(yr) all_func(yr))
  stopCluster(cl)
  
  #### Now need to read all the annual files and rbind them together:
  setwd(outpath)
  
  # List all .RData files in the folder
  rdata_files <- list.files(pattern = "\\.Rdata$")
  
  # Exclude additional info files to load
  exclude_patterns <- c("additional_info", paste0(min(years), '_', max(years)))
  rdata_files_to_load <- rdata_files[!grepl(paste(exclude_patterns, collapse = "|"), basename(rdata_files))]
  rdata_files_to_load <- rdata_files_to_load[grepl(paste(file_name, collapse = "|"), basename(rdata_files_to_load))]
  
  #### Now need to read all the annual files and rbind them together:
  load(rdata_files_to_load[1])
  
  all_df <- a_samp_dat %>%
    left_join(samps_rad_info_yr, c('sample'='sample_locs'))
  
  for(df in rdata_files_to_load[-1]) {
    
    load(df)
    
    yr_df = a_samp_dat %>%
      left_join(samps_rad_info_yr, c('sample'='sample_locs'))
    
    all_df <- all_df %>%
      bind_rows(yr_df)
  }
  
  # append to sample info for plotting 
  samp_info <- tidy_sppdf %>%
    distinct(sample, latitude, longitude, year, month, date)
  
  all_df <- all_df %>%
    dplyr::select(-year) %>%
    left_join(samp_info, c('sample')) 
  
  save(all_df, file=paste0('biodiversity_estimates_', 
                           min(years), '_', max(years), '_', file_name, '.Rdata'))
  load(paste0('biodiversity_estimates_', 
              min(years), '_', max(years), '_', file_name, '.Rdata'))
  
  return(all_df)
  
}