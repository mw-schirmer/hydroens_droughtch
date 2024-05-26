# -----------------------------------------------------------------------------
# Script Name:   run_models.R
# Created By:    Michael Schirmer
# Description:   extracts historical input and monthly forecast for a set of catchments
#                runs the openQuarrel models
#                visualizes meteorological input and and hydroogical results as
#                  leaflet

suppressPackageStartupMessages(library(dplyr))
library(tidyr) %>% suppressPackageStartupMessages()
library(purrr) %>% suppressPackageStartupMessages()
library(stringr) %>% suppressPackageStartupMessages()
library(readr) %>% suppressPackageStartupMessages()
library(sf) %>% suppressPackageStartupMessages()
library(ggplot2) %>% suppressPackageStartupMessages()

library(airGR) %>% suppressPackageStartupMessages()


# user settings ---------------

# determine if you want to use parallel computing code
comp_parallel <- FALSE

# show progress bars
show_progress_bars <- FALSE

# clean data during run
clean_data <- TRUE

# output image at end
output_image <- FALSE

# define local or vilan run
local <- FALSE

if (local) {
  git_folder <- "D:/gitlabext/hydroens_droughtch/operational"
  # base folder for data stored not in git
  data_input_folder <- "D:/hydroens_droughtch/data/input/openQuarrel/drought_ch_oper"
  forecast_basefolder <- "D:/hydroens_droughtch/data/example_meteogrids/mch/"
  # if parellel computing, specify the number of cores
  num_cores <- 7

} else {
  git_folder <- "/home/schirmer/gitlabext/hydroens_droughtch/operational"
  # base folder for data stored not in git
  data_input_folder <- "/home/schirmer/hydroens_droughtch_oper/data/input"
  forecast_basefolder <- "/home/hydro/monthlyIFS/input_pp/ifs_grid_prevah"
  # in order to make png working for leaflet map on vilan
  options(bitmapType = "cairo")
  # if parellel computing, specify the number of cores
  num_cores <- 30
}

# base folder for data input stored in git
data_input_git_folder <- file.path(git_folder, "data", "input")

# shapefile folder with catchment information
shape_folder <- file.path(data_input_git_folder, "shapefile")

# basin information folder
basin_info_folder <- file.path(data_input_git_folder, "basin_info")
raster_name <- "dhm_25_l2"

# climatology folder
clim_folder <- file.path(data_input_git_folder, "climatology")

# model_parameter folder, in git
model_par_folder <- file.path(data_input_git_folder, "model_parameters")

# meteo data source, so far only used to load historic data
precip_data_type <- "RhiresD"
temp_data_type <- "TabsD"

# historic driving data settings, not in git
# todo: think about a daily update
historic_data_folder <- file.path(data_input_folder, "historic_data")
historic_data_filename <- "historic_data_3.rds"

# recent driving data settings, not in git
# todo: think about second hierarchy with 2 km lat long data, in case 1 km is not available
P_verified_folder <- file.path(data_input_folder, "griddata/verified/daily/RhiresD_daily_precipitation/swissgrid/netcdf")
T_verified_folder <- file.path(data_input_folder, "griddata/verified/daily/TabsD_daily_mean_temperature/swissgrid/netcdf")
P_operational_folder <- file.path(data_input_folder, "griddata/operational/daily/RprelimD_daily_precipitation/swissgrid/netcdf")
T_operational_folder <- file.path(data_input_folder, "griddata/operational/daily/TabsD_daily_mean_temperature/swissgrid/netcdf")

# specify which last forecast should be run, -1 is newest, -2 is second newest etc
last_forecast <- -1

# determine if script needs to run ---------------------------
# determine newest initial time, i.e. last in alphanumerical order, e.g. 20230724
forecast_init_str <- list.dirs(forecast_basefolder, recursive = FALSE) %>% nth(last_forecast) %>% basename()

# where the forecast data is
forecast_folder <- file.path(forecast_basefolder, forecast_init_str)

# store the forecast type, i.e. at third position of forecast_basefolder, e.g. meteo_ifs_pp_new
# forecast_type <- str_split_i(forecast_basefolder, pattern = "/", i = 3)
# this does not work for original folder, so it is hardcoded
forecast_type <- "meteo_ifs_pp_new"

# output folder, not in git
output_folder <- file.path(data_input_folder, "../results", forecast_type, forecast_init_str)

# # check if output folder exists, if yes, don't run
# if (dir.exists(output_folder)) {
#   quit(save = "no")
# }



# functions ----------------------
source(file.path(git_folder, "R/extract_rast.R"))
# todo: as this is the official version (except for not including stars and sf, this function has to be linked at vilan)
source(file.path(git_folder, "R/read_prevah.R"))
source(file.path(git_folder, "R/leaflet_weekly_drougthch_map.R"))
source(file.path(git_folder, "R/runoff_plotting_fns.R"))
source(file.path(git_folder, "R/class_dist_plot.R"))
source(file.path(git_folder, "R/util_fns.R"))
source(file.path(git_folder, "R/meteoswiss_app_plot.R"))
source(file.path(git_folder, "../../ensemblehydromodel/R/functions_and_settings.R"))

# x is a vector, returns the vector without trailing NAs, other NAs are kept
remove_trailing_NAs <- function(x) {
  rev(cumsum(!is.na(rev(x)))) != 0
}

# fill NAs linearly in column var of a data frame called df_name
# assume that all groups (e.g. catchments) have the same NA for warning message
fill_na <- function(df, var, df_name) {
  sum_na <- sum(is.na(df[[var]]))/n_groups(df)
  if (sum_na > 0) {
    warning(sprintf("%s NAs per group %s in variable %s of %s, apply linear interpolation ...",
                    sum_na, group_vars(df), var, df_name), call. = FALSE)
    df[[var]] <- zoo::na.approx(df[[var]])
  }
  return(df)
}

# remove duplicates in a data frame containing a date column and take the second one
rem_duplicates <- function(df) {
  df <- mutate(df, temp_order = seq_along(df$date)) %>% arrange(date, desc(temp_order)) %>%
    distinct(date, .keep_all = TRUE) %>% select(-temp_order)
}

# extract forecast helper function
extract_forecast <- function(forecast_var, forecast_folder, shape, y_ID_col, output_var = "value"){

  # list all files recursively in the forecast folder starting with forecast_var and ending with .2km
  forecast_files <- list.files(file.path(forecast_folder), pattern = paste0("^", forecast_var, ".*2km$"),
                               full.names = TRUE, recursive = TRUE)

  # member is considered to be the second last folder hierarchy
  member <- forecast_files %>% str_split_i(pattern = "/", i = -2)
  # date is considered to be first number in the file name
  forecast_dates <- forecast_files %>% str_split_i(pattern = "/", i = -1) %>%
    parse_number()
  # join forecast day and member to a new name
  names(forecast_files) <- paste(member, forecast_dates, sep = "_")

  # read forecast files and convert it to terra SpatRaster
  forecast_data_list <- map(forecast_files, \(x) read.prevah(x) %>% terra::rast())

  # extract data, needs
  extracted_data <- extract_rast(forecast_data_list, shape_lv03, y_ID_col = "mach_ID", output_var = output_var) %>%
    # separate member from date
    separate(col = name, into = c("member", "date"), sep = "_") %>%
    # convert date vector
    mutate(date = lubridate::ymd(date), source = "forecast", .keep = "unused")

  return(extracted_data)
}


# run for a forecast member a model (stemming from a calibration method and calibration criterion and Q transformation) stored in model_parameter_df,
# for a catchment ID with Basin_info for all catchments in a list
# driving data is stored in vectors and matrices as this is much faster as a data.frame (mach_ID_col, DatesR_col, member_col, driving_mat)
# the length of the forecast (days) is used to subset the simulated runoff to just the forecast period
run_single_model <- function(catchment_ID, model_comb_thisrun, driving_data_catchment,
                             cal_q_transfo,
                             model_parameter_catchment, Basin_info_catchment,
                             length_forecast, cal_par) {

  # output NA if run fails, numeric, otherwise cut_Q fails
  Qsim <- NA %>% as.numeric()

  # catch conditions if run fails
  # todo, record conditions somehow that is still fast
  conds_single_model_run <- catch_cnds({

    # select model parameters of very model and transfo combination
    is_single_model <- model_parameter_catchment$model_comb == model_comb_thisrun &
      model_parameter_catchment$transfo_new == cal_q_transfo
    model_parameter_df <- model_parameter_catchment[is_single_model,]

    # parameters stored in column data, which is a list which contains a df with the column parameters
    param <- model_parameter_df %>% pull(data) %>% first() %>% pull(parameters)

    if (nrow(model_parameter_df) != 1) stop(nrow(model_parameter_df))

    # create input for runoff_model
    runoff_input <- create_input(model_parameter_df$runoff_model, driving_data_catchment, Basin_info_catchment) %>%
      suppressMessages() %>% suppressWarnings()

    # retrieve snow module param and update precipitation input
    # create snow module and input
    snow_module <- model_parameter_df$snow_module
    if (!is.na(snow_module)) {
      # create snow input
      snow_input <- create_input(snow_module, driving_data_catchment, Basin_info_catchment) %>%
        suppressMessages() %>% suppressWarnings()

      nof_param_snow <- cal_par[[snow_module]][["nof_param"]]
      snow_param <- param[1:nof_param_snow]
      # remaining runoff model parameters
      param <- param[(nof_param_snow + 1):length(param)]

      # model solid precipitation
      # ensure that P is present in input
      if (!"P" %in% names(runoff_input)) stop("P is not an entry of a list like input")
      snow_module_results <- simulate_snow(snow_module, snow_param, snow_input) %>%
        suppressMessages() %>% suppressWarnings()

      # update precipitation with snow module surface water runoff
      runoff_input$P <- snow_module_results$surface_water_runoff

    }

    # simulate the last x years, but output only the forecast days
    length_sim <- nrow(driving_data_catchment)
    Qsim <- simulate_model(model_parameter_df$runoff_model, param, runoff_input
    )$Qsim[seq(length_sim-length_forecast+1, length_sim)] %>%
      suppressMessages() %>% suppressWarnings()

  })

  return(Qsim)
}


subset_driving_data <- function(catchment_ID, forecast_members, driving_data) {

  is_catch <- driving_data$mach_ID == catchment_ID
  is_member <- driving_data$member %in% c(forecast_members, "none")
  is_relevant <- is_catch & is_member
  driving_data <- driving_data[is_relevant,]

  return(driving_data)

}

# run for all models stored in model parameters with the same data for a forecast member
# for a catchment ID with Basin_info for all catchments in a list
# driving data is stored in vectors and matrices as this is much faster as a data.frame (mach_ID_col, DatesR_col, member_col, driving_mat)
# the length of the forecast (days) is used to subset the simulated runoff to just the forecast period
run_all_models <- function(catchment_ID, forecast_member,
                           model_parameter_df, Basin_info_list,
                           driving_data,
                           length_forecast, cal_par) {

  # subset driving data
  driving_data_catchment <- subset_driving_data(catchment_ID, forecast_member, driving_data)

  # subset spatial data
  Basin_info_catchment <- Basin_info_list[[catchment_ID]]

  # select model parameters available for this catchment
  model_parameter_catchment <- model_parameter_df %>% filter(mach_ID == catchment_ID) #%>%
  # todo_urgent: delete
  # filter(runoff_model %in% c("topmodel", "TUW", "sacramento", "CemaNeigeGR4J"))

  # loop over all model parameters in this catchment
  loop_var <- model_parameter_catchment %>% select(model_comb, transfo_new)

  # it is  required that the looping variables are not existent
  rm(model_comb, transfo_new) %>% suppressWarnings()

  # run the models over all looping variables
  Qsim_list <- pmap(loop_var,
                    \(model_comb, transfo_new) run_single_model(catchment_ID, model_comb, driving_data_catchment,
                                                                transfo_new,
                                                                model_parameter_catchment, Basin_info_catchment,
                                                                length_forecast, cal_par))

  # names list entries according to model_parameter_catchment
  # (without model parameters in column data (not needed) and mach_ID (in outer loop))
  names(Qsim_list) <- model_parameter_catchment %>%
    select(-all_of(c("data", "mach_ID"))) %>%
    # replace 0.2 in 0/2 as we need the . for unlist later
    mutate(transfo_new = str_replace(transfo_new, "\\.", "/")) %>%
    mutate(error_crit_transfo = str_replace(error_crit_transfo, "\\.", "/")) %>%
    mutate(exponent = str_replace(exponent, "\\.", "/")) %>%
    unite(combination, everything(), sep = ".") %>% pull(combination)

  return(Qsim_list)
}


# load data  ========================
# get get hash and report if there are no local changes
git_hash <- system(paste0("git -C ", git_folder, " rev-parse HEAD"), intern = TRUE)
git_hash_short <- system(paste0("git -C ", git_folder, " rev-parse --short HEAD"), intern = TRUE)

# change git hash when local changes are present
# todo: think about reporting changes in the file line by line
git_changes <- system(paste0("git -C ", git_folder, " status --porcelain"), intern = TRUE)
if (!identical(git_changes, character(0))) {
  git_hash = paste(git_hash, git_changes)
  git_hash_short <- paste(git_hash_short, git_changes)
}

message(Sys.time(), " >> start forecast for ", forecast_folder, " with git hash ", git_hash_short, " ...")

# load shape file to define catchments and to extract meteo data
# lv03 version for prevah
shape_lv03 <- readRDS(file.path(shape_folder, "shape_drought_ch.rds"))
# lv95 version for meteoswiss data
shape_lv95 <- readRDS(file.path(shape_folder, "shape_drought_ch_lv95.rds"))
shape_wgs84 <- readRDS(file.path(shape_folder, "shape_drought_ch_wgs84_simpl100.rds"))
catchment_IDs <- shape_lv03$mach_ID
names(catchment_IDs) <- catchment_IDs
catchment_owner_IDs <- shape_lv03$owner_ID

# load all model parameters
model_parameter_df <- readRDS(file.path(model_par_folder,
                                        "model_parameters_split_ensmos.rds"))

# load basin information
Basin_info_list <- map(catchment_IDs, \(x) read_rds(file.path(basin_info_folder, paste0("BasinInfo_", x, "_", raster_name, ".rds"))))

# load climatology
# todo_urgent: get new quantiles
climate_quantiles <- readRDS(file.path(clim_folder, "climate_quantiles_KGE_5transfo_11models_11quantiles_daily.rds"))


## load historic data (already extracted) ---------------
source_data_type <- paste(precip_data_type, temp_data_type, sep = "_") %>% tolower()
historic_data <- read_rds(file.path(historic_data_folder, source_data_type, historic_data_filename)) %>%
  mutate(source = "historic") %>%
  # just for safety, otherwise a grouping variable is added later on
  ungroup()

## extract recent historic data ---------------
message(Sys.time(), " > start extracting historic data ...")
recent_folder_list <- list(P_operational = list.files(P_operational_folder, full.names = TRUE),
                           T_operational = list.files(T_operational_folder, full.names = TRUE),
                           P_verified = list.files(P_verified_folder, full.names = TRUE),
                           T_verified = list.files(T_verified_folder, full.names = TRUE))

# loop over all folder and extract data
recent_extract <- map_dfr(recent_folder_list,
                          \(x) extract_rast(x, shape_lv95, y_ID_col = "mach_ID", name_time = TRUE),
                          .id = "source", .progress = show_progress_bars)  %>%
  # convert date vector
  mutate(date = as.Date(name), .keep = "unused")  %>%
  # arrange source verifies that operational is before verified, which is needed for the second remove duplicates
  arrange(source, date)  %>%
  # a first remove duplicates within sources, this is needed otherwise pivot_wider will not work
  group_by(mach_ID, source) %>% group_modify(~ rem_duplicates(.x)) %>%
  # separate variable from source and create P and T columns
  separate(source, into = c("var", "source"), sep = "_") %>% pivot_wider(names_from = var, values_from = value) %>%
  # remove duplicates and take second, this round takes into acount duplicates between sources
  group_by(mach_ID) %>% group_modify(~ rem_duplicates(.x)) %>%
  # remove trailing NAs as they cannot be filled with na.approx
  # (and I do not want them to be extrapolated, rather ignored, otherwise use na.spline in fill.na)
  filter(remove_trailing_NAs(P), remove_trailing_NAs(T)) %>%
  # linear interpolation of P and T
  fill_na("P", "recent_extract") %>%
  fill_na("T", "recent_extract") %>%
  # calculate evapotranspiration
  # todo: test if this grouping is needed
  group_by(mach_ID) %>%
  group_modify( ~ mutate(.x, E = airGR::PE_Oudin(JD = as.POSIXlt(.x$date)$yday + 1,
                                                 Temp = .x$T,
                                                 Lat = 46.8, LatUnit = "deg"))) %>%
  ungroup() %>%
  # round for two decimal digits
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

# merge with historic and keep only distinct with the  order of priority historic, recent_verified, recent_operational
historic_data_merged <- bind_rows(historic_data, recent_extract) %>%
  distinct(mach_ID, date, .keep_all = TRUE)

# provide info on time range
first_historic <- historic_data_merged %>% pull(date) %>% min
latest_verified <- historic_data_merged %>% filter(source == "verified" | source == "historic") %>% pull(date) %>% max
latest_operational <- historic_data_merged %>% filter(source == "operational") %>% pull(date) %>% max

message(Sys.time(), " > ", sprintf("Historic data extracted, available from %s to %s (with latest verified %s) ...",
                                   first_historic, latest_operational, latest_verified))

# overwrite historic_data updated with new verified data as one rds file
# only if there are any new verified time steps
if (latest_verified > historic_data$date %>% max) {
  output_historic_data <- historic_data_merged %>% filter(source %in% c("historic", "verified")) %>%
    # drop source folder as data is verified
    select(-source) %>%
    # # todo: exclude dates at the start in order to have only x years rather than a start year to not accumulate the data over the years
    # filter(date >= "2020-01-01") %>%
    write_rds(file.path(historic_data_folder, source_data_type, historic_data_filename))
}

if (clean_data) {
  rm(historic_data, recent_extract)
  invisible(gc())
}


## extract forecast data-------------
message(Sys.time(), " > start extracting forecast data ...")

input_var_vec <- c(P = "prec", T = "temp")
forecast_extract <- imap(input_var_vec,
                         \(x, y) extract_forecast(x, forecast_folder,
                                                  shape = shape_lv03, y_ID_col = "mach_ID", output_var = y),
                         .progress = show_progress_bars) %>%
  # merge temp and precip
  reduce(full_join) %>% suppressMessages() %>%
  # calculate evapotranspiration
  group_by(mach_ID, member) %>%
  group_modify(~ mutate(.x,E = airGR::PE_Oudin(JD = as.POSIXlt(.x$date)$yday + 1,
                                               Temp = .x$T,
                                               Lat = 46.8, LatUnit = "deg"))) %>%
  ungroup() %>%
  # round for two decimal digits
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

# store forecast members and forecast dates
forecast_member_vec <- forecast_extract$member %>% unique()
forecast_dates <- forecast_extract$date %>% unique()

# store forecast outliers
outlier_forecasts <- forecast_extract %>% filter(P > 5e2) %>% ungroup() %>%
  distinct(member, .keep_all = "TRUE")
outlier_member <- outlier_forecasts %>% pull(member)

# combine historic, recent and forecast data
driving_data <- historic_data_merged %>%
  # take always full time period of forecast, i.e. delete all historic time stamps within the forecast range
  filter(date < forecast_dates[1]) %>%
  # add a forecast member type none for historic data
  mutate(member = "none") %>%
  # bind with forecast
  bind_rows(forecast_extract) %>%
  # # this is not needed as all data should already be distinct, but kept as security:
  # # take only the first of double entries for each date, catchment and forecast member
  distinct(date, mach_ID, member, .keep_all = TRUE) %>%
  # convert to ensemblehydromodel needs:
  # change date into DatesR and as POSIXct
  rename(DatesR = date) %>%
  mutate(DatesR = as.POSIXct(DatesR)) %>%
  # just for safety, otherwise a grouping variable is added later on
  ungroup()
# set time zone, needs to be UTC in order to keep midnight in later steps (e.g. outputting)
attr(driving_data$DatesR, "tzone") <- "UTC"


# input leaflet ---------------------
message(Sys.time(), " > save input leaflet ...")
conds_input_leaflet <- catch_cnds({

  # create popup_plot_list for leaflet map
  # todo think about doing this with group_map and without subsetting in patch_both_runoff_plots
  meteoswiss_plot_list <- map2(catchment_IDs, catchment_owner_IDs,
                               \(x, y) suppressMessages(meteoswiss_app_plot(x, y, forecast_extract, transfer_add = 20)),
                               .progress = show_progress_bars)

  # create dummy input for no classification
  dummy_df <- forecast_extract %>% filter(member == "ctr", mach_ID == "CHZH-0042") %>% slice(1:7) %>%
    mutate(class = factor("1"))

  # create leaflet map
  leaflet_map <- leaflet_weekly_droughtch_map(shape_wgs84, dummy_df, mach_ID, meteoswiss_plot_list, class_colours = col2hex("grey"),
                                              fillOpacity = 0.2, pop_width = 500, pop_height = 300)

  # create directory
  meteo_leaflet_folder <- file.path(output_folder, "meteo_leaflet")
  dir.create(meteo_leaflet_folder, recursive = TRUE, showWarnings = FALSE)

  # save input leaflet
  htmlwidgets::saveWidget(widgetframe::frameableWidget(leaflet_map),
                          file.path(meteo_leaflet_folder, "meteo_leaflet.html"),
                          selfcontained = FALSE)

  # remove all subsequent steps
  if (clean_data) {
    rm(meteoswiss_plot_list, dummy_df, leaflet_map)
    invisible(gc())
  }

})

# remove all subsequent steps
if (clean_data) {
  rm(historic_data_merged, forecast_extract)
  invisible(gc())
}

# run models -----------------
message(Sys.time(), " > run models ...")

# get default calibration parameters
cal_par <- default_cal_par

# lambda routing for hydromad
hydromad_routing_vec <- c("bucket", "awbm", "cwi", "cmd", "snow")
names(hydromad_routing_vec) <- hydromad_routing_vec
hydromad_par <- purrr::map(hydromad_routing_vec, ~ set_hydromad_par(.x, routing = "lambda"))

# overwrite hydromad options
cal_par[names(hydromad_par)] <- hydromad_par

# run all all models, outer loop over all models which have not the same data, inner loop: same data
loop_var <- expand.grid(mach_ID = catchment_IDs, member = forecast_member_vec,
                        stringsAsFactors = FALSE)

# it is  required that the looping variables are not existent
rm(mach_ID, member) %>% suppressWarnings()


# helper function to apply both on parallel an non-parallel code
processing_function <- \(mach_ID, member) run_all_models(mach_ID, member,
                                                         model_parameter_df, Basin_info_list, driving_data,
                                                         length(forecast_dates), cal_par)

if (comp_parallel) {
  library(furrr)
  library(future)

  plan(multisession, workers = num_cores)
  Qsim_list <- future_pmap(loop_var, processing_function,
                           .options = furrr_options(packages = c("dplyr", "tidyr", "purrr", "stringr", "readr", "airGR")),
                           .progress = show_progress_bars)

} else {
  Qsim_list <- pmap(loop_var, processing_function, .progress = show_progress_bars)
}

message(Sys.time(), " > modify Qsim ...")
# names list entries according to loop_var
names(Qsim_list) <- loop_var %>% unite(combination, everything(), sep = ".") %>% pull(combination)

# unlist Qsim_list, works only with "."
Qsim_list <- unlist(Qsim_list, recursive = FALSE)

# define which function should be used
mapping_function <- if (comp_parallel) future_map_dfr else map_dfr

# create a long data frame
# long table with a combined type
Qsim_long <- mapping_function(Qsim_list, \(x) tibble(date = forecast_dates, Qsim = x), .id = "type", .progress = show_progress_bars) %>%
  # clean long table separated in multiple columns
  separate_wider_delim(type, names = c("mach_ID", "member", colnames(model_parameter_df %>% select(-all_of(c("data", "mach_ID"))))),
                       delim = ".", too_few = "align_start", cols_remove = FALSE) %>%
  # replace "/" in "." in all columns with exponent again to have 0.2 instead of 0/2
  mutate(transfo_new = str_replace(transfo_new, "/", "\\.")) %>%
  mutate(error_crit_transfo = str_replace(error_crit_transfo, "/", "\\.")) %>%
  mutate(exponent = str_replace(exponent, "/", "\\.")) %>% ungroup() %>%
  # add package of the runoff model
  # not clear why rowwise is needed
  rowwise %>% mutate(family = get_family(runoff_model)) %>% ungroup()

if (clean_data) {
  rm(Qsim_list)
  invisible(gc())
}


# classify data ----------
message(Sys.time(), " > classify data ...")
# define relevant columns
relevant_cols <- c("yday", "mach_ID", "model_comb", "transfo_new")
# merge week based climatology with forecast for all groups
Qsim_long <- Qsim_long %>%
  # add yday column need to join with climatology
  mutate(yday = lubridate::yday(date)) %>%
  # now join using relevant columns
  left_join(climate_quantiles %>% select(c(all_of(relevant_cols), starts_with("quantile"))), by = relevant_cols) %>%
  # delete all entries with NA in limits
  filter(if_all(starts_with("quantile"), ~ !is.na(.))) %>%
  # loop over all group combinations and classify data
  group_by(mach_ID, model_comb, transfo_new) %>% group_modify(~ cut_Q(
    .x, Qsim,
    breaks_expr = c(-Inf, quantile2.5, quantile10, quantile25, quantile33, quantile66, quantile75, quantile90, quantile97.5, Inf),
    labels_expr = c("< q2.5", "< q10", "< q25", "< q33", "< q66", "< q75", "< q90", "< q97.5", ">= q97.5")
  )) %>%
  ungroup()


# save data --------------
message(Sys.time(), " > save data ...")
conds_save <- catch_cnds({

  data_folder <- file.path(output_folder, "data")
  dir.create(data_folder, recursive = TRUE, showWarnings = FALSE)

  # save git hash in a text file
  write(git_hash, file.path(data_folder, "git_hash.txt"))

  # save forecast data frame
  write_rds(Qsim_long, file.path(data_folder, "Qsim_long.rds"), compress = "gz")

  # save outlier data
  if (!identical(outlier_member, character(0))) {
    write_rds(outlier_forecasts, file.path(data_folder, "outlier_forecasts.rds"))
    write_csv(outlier_forecasts, file.path(data_folder, "outlier_forecasts.csv"))
    message(Sys.time(), " > outlier present in forecast and information stored ...")
  }

  # store runs with NAs in class as ratio of whole run
  na_affected <- Qsim_long %>% group_by(mach_ID, member, model_comb, transfo_new) %>%
    summarise(NA_ratio_class = sum(is.na(class))/n(),
              NA_ratio_Qsim = sum(is.na(Qsim))/n()) %>%
    filter(NA_ratio_class > 0) %>%
    arrange(desc(NA_ratio_class))
  if (nrow(na_affected) != 0) {
    write_rds(na_affected, file.path(data_folder, "na_forecasts.rds"))
    write_csv(na_affected, file.path(data_folder, "na_forecasts.csv"))
    message(Sys.time(), " > NAs present in forecast and information stored ...")
  }

}, "cond_save")


# leaflet map --------------
message(Sys.time(), " > save leaflet map ...")
# filter out NAs in Qsim and class for later plots
Qsim_long <- filter(Qsim_long, !is.na(Qsim), !is.na(class))

# choose class colours
batlow <- khroma::colour("batlow")
class_colours <- batlow(Qsim_long$class %>% levels %>% length, range = c(0, .83)) %>% rev()

# create popup_plot_list for leaflet map
# todo think about doing this with group_map and without subsetting in patch_both_runoff_plots
runoff_plot_list <- map2(catchment_IDs, catchment_owner_IDs,
                         \(x, y) suppressMessages(get_runoff_plot(x, y, Qsim_long)),
                         .progress = show_progress_bars)

class_dist_plot_list <- map2(catchment_IDs, catchment_owner_IDs,
                             \(x, y) suppressMessages(get_class_dist_plot(x, y, Qsim_long, class_colours)),
                             .progress = show_progress_bars)

# determine forecast weeks for weekly highlighting
forecast_start_date <- forecast_dates %>% min
forecast_end_date <- forecast_dates %>% max
week_lead_max_ind <- difftime(forecast_end_date + 1, forecast_start_date, unit = "weeks") %>% as.double() %>% ceiling()
week_lead_ind <- 1:week_lead_max_ind

# create looping variables
loop_var <- expand.grid(class_dist_plots = class_dist_plot_list, week_highlight = week_lead_ind,
                        stringsAsFactors = FALSE)
# add also the runoff plot list
loop_var$runoff_plots <- runoff_plot_list

# remove all subsequent steps
if (clean_data) {
  rm(runoff_plot_list, class_dist_plot_list)
  invisible(gc())
}

# it is  required that the looping variables are not existent
rm(class_dist_plots, runoff_plots, week_highlight) %>% suppressWarnings()

# loop (simultaneous over the two) plots and the to be highlighted weeks
popup_plot_list <- pmap(loop_var, \(class_dist_plots, runoff_plots, week_highlight)
                        highlight_and_patch_both_runoff_plots(class_dist_plots, runoff_plots, week_highlight),
                        .progress = show_progress_bars)

# remove all subsequent steps
if (clean_data) {
  rm(loop_var)
  invisible(gc())
}

# create leaflet map
leaflet_map <- leaflet_weekly_droughtch_map(shape_wgs84, Qsim_long, mach_ID, popup_plot_list, class_colours,
                                            pop_width = 500, pop_height = 400)

# save html
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
htmlwidgets::saveWidget(widgetframe::frameableWidget(leaflet_map),
                        file.path(output_folder, "leaflet.html"),
                        selfcontained = FALSE)

# save data image for release runs
if (output_image) {
  rm(popup_plot_list, leaflet_map, Qsim_long)
  save.image(file.path(output_folder, "debugging.Rdata"))
}

# save all caught warnings and errors ----------------

# combine all subconditions in a list
subconds <- mget(ls(pattern = "conds_"))

# if there are conditions caught, save them
if (length(subconds) != 0) {

  # re-order conditions to errors and others
  append_error_elements <- function(conds, subconds) {

    conds$errors <- append(conds$errors, subconds$errors)
    conds$others <- append(conds$others, subconds$others)
    return(conds)
  }
  conds <- reduce(subconds, append_error_elements)

  saveRDS(conds, file.path(output_folder, "caught_conditions.rds"))

}

message(Sys.time(), " >> finished script.")
