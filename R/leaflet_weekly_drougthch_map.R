#' leaflet_weekly_droughtch_map
#' 
#' A leaflet map with polygons, popup plots and weekly aggregated colours.
#' 
#' Creates a leaflet map. Plots a polygon with colours representing the most frequent class in the 
#'  aggregated weeks. The tint of the colour represents the frequency of this class.
#'  Plot as many weeks as there are layers in the data. 
#'  The levels in the classified data determines the classes plotted.
#'  On click of polygons there are a popup plots, which 
#'  need to be provided as input. 
#'  Requires the packages leaflet, sf, purrr, dplyr, tidyr, ggplot2 and lubridate to be installed.
#'  
#' @param polygon_shape an object representing a polygon shape file, e.g. a Simple feature or a Spatial object. 
#'  Will be projected to WGS84.
#' @param forecast_df a data frame which contains a column "class" with forecast categories as factors, 
#'  and a "date" column of class date. 
#' @param ID_col <data-masking> one column which contains the catchment IDs, 
#'  needs to be present in both forecast_df and polygon_shape.
#' @param popup_plot_list a list of ggplots used as popup plot. Needs to be of the same length (i.e. number of catchment) as polygon_shape, 
#'  or a multiple of it to be used as different plot for each weekly layer.
#'  An example plot can be created with class_dist_plot.R 
#' @param class_colours a vector with class colours, of the same length as levels in forecast_df$class.
#' @param ncolours_tint an integer which is used to generate class_colours with tint to reflect the class_distribution.
#' @param fillOpacity an integer from 0 to 1 (default) indicating the opacity of the colours to fill the polygons
#' @param pop_width an integer of the width of the pop up plots in pixel (default is 300)
#' @param pop_height an integer of the height of the pop up plots in pixel (default is 300)
#' 
#' @return a leaflet map
#' 
#' @author Michael Schirmer, michael.schirmer[at]wsl.ch, Konrad Bogner, konrad.bogner[at]wsl.ch 
#' @note todo:
#'  1) make the maximum lead time week_lead_max_ind as input
#'
#' @export
#' @examples
#' # Qsim_long is a df with the forecast data for all catchments, mach_ID is the catchment ID column in this df
#' 
#' # create colours
#' batlow <- khroma::colour("batlow")
#' class_colours <- batlow(Qsim_long$class %>% levels %>% length, range = c(0, .83)) %>% rev()
#' 
#' popup_plot_list <- Qsim_long  %>% group_by(mach_ID) %>%  group_map(~ class_dist_plot(Qsim_long))
#' 
#' # function call
#' leaflet_map <- leaflet_weekly_droughtch_map(shape_lv03, Qsim_long, mach_ID, popup_plot_list, class_colours)
#' 
#' # save html                                             
#' htmlwidgets::saveWidget(widgetframe::frameableWidget(leaflet_map), leaflet.html", selfcontained = FALSE)                                            
leaflet_weekly_droughtch_map <- function(polygon_shape, forecast_df, ID_col, popup_plot_list, 
                                         class_colours, ncolours_tint = 10, fillOpacity = 1,
                                         pop_width = 300, pop_height = 300) {
  
  library(leaflet)
  library(sf)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  
  # internal functions -------------------
  
  # estimate the probability of the most frequent class for each catchment and for defined week with lead index week_lead_ind, 
  # returns a df with columns ID, class, n and prob and as many rows as catchments
  calc_max_class_prob_per_week <- function(forecast_df, week_lead_ind) {
    # determine week with lead index
    forecast_date1 <- forecast_df$date %>% min
    week_lead <- seq(forecast_date1 + weeks(week_lead_ind - 1), forecast_date1 + weeks(week_lead_ind) - days(1), by = 1)
    class_prob_df <- forecast_df %>% 
      filter(date %in% week_lead) %>% 
      group_by({{ ID_col }}) %>% count(class) %>% 
      mutate(prob = n/sum(n)) %>% slice_max(n, with_ties = FALSE) %>% ungroup()
    
    return(class_prob_df)
    
  }
  
  # create background colour shading to white, i.e. tint
  find_tint <- function(colour, prob, ncolours_tint) {
    tint_pal_fn <- colorRampPalette(c("white", colour))
    tint_pal <- tint_pal_fn(ncolours_tint)
    tint_colour <- tint_pal[findInterval(prob, seq(1, ncolours_tint, by = 1)/ncolours_tint)]
    return(tint_colour)
  }
  
  # this function can be looped over all classes, in order to avoid hardcoding the class levels
  find_class_colour_with_tint_class <- function(class_label, class_prob_df, class_colour_value, ncolours_tint) {
    # if class_colour not present, add this column with NAs
    if (!"class_colour" %in% colnames(class_prob_df)) {
      class_prob_df$class_colour <- NA
    }
    # override NAs with class colour for this specific class
    colour_vec <- class_prob_df %>% 
      mutate(class_colour = if_else(
        class == class_label, find_tint(class_colour_value, prob, ncolours_tint), class_colour))
    
    return(colour_vec)
  }
  
  # loop above function over all levels
  find_class_colour_with_tint <- function(levels, class_prob_df, class_colours, ncolours_tint) {
    colour_df <- class_prob_df
    for (i in seq_along(levels)) {
      # subsequently override test
      colour_df <- find_class_colour_with_tint_class(levels[i], colour_df, class_colours[i], ncolours_tint)
    }
    
    return(colour_df)
  }
  
  # add a polygon to a leaflet map
  leaflet_add_polygon <- function(leaflet_map, polygon_sf, popup_plot_list, group_name, 
                                  pop_width, pop_height, fillOpacity) {
    leaflet_map <- leaflet_map %>% 
      addPolygons(
        data = polygon_sf,
        opacity = .2,
        color = "black",
        weight = 2,
        # todo: make width and height an input 
        popup = leafpop::popupGraph(popup_plot_list, width = pop_width, height = pop_height),
        options = list(clickable = TRUE),
        fill = TRUE,
        fillColor = polygon_sf$class_colour,
        fillOpacity = fillOpacity,
        # fillColor = "white",
        # fillOpacity = 0.5,
        layerId = seq_along(polygon_sf),
        label = polygon_sf$owner_ID,
        group = group_name
      )
    
    return(leaflet_map)
  }
  
  
  # main code -----------------
  # get levels in forecast_df
  levels <- forecast_df$class %>% levels
  
  # transform polygon_shape to sf
  polygon_sf <- polygon_shape %>% sf::st_as_sf() #%>% 
  # project to WGS84 lat lon, does not work currently on vilan myR
  #  st_transform("+proj=longlat +datum=WGS84") 
  
  # arrange by area to plot smaller on top of large ones
  if ("Shape_Area" %in% polygon_sf) {
    polygon_sf <- arrange(polygon_sf, desc(Shape_Area))
  }
  
  # determine forecast weeks for weekly plots
  forecast_start_date <- forecast_df$date %>% min
  forecast_end_date <- forecast_df$date %>% max
  week_lead_max_ind <- difftime(forecast_end_date + 1, forecast_start_date, unit = "weeks") %>% as.double() %>% ceiling()
  week_lead_ind <- 1:week_lead_max_ind
  names(week_lead_ind) <- paste("week", week_lead_ind)
  
  # check the length of popup_plot_list
  if (length(popup_plot_list) %% nrow(polygon_shape) != 0) stop("popup_plot_list is not a multiple of the nrow(polygon_shape).")
  
  if (length(popup_plot_list) %/% nrow(polygon_shape) == week_lead_max_ind) {
    
    # reorder long list as many as there are weeks in the forecast
    popup_plot_list_weekly <- list()
    for (i in seq_along(week_lead_ind)) {
      j <-  i-1
      popup_plot_list_weekly[[i]] <- popup_plot_list[(j*length(catchment_IDs)+1):((j+1)*length(catchment_IDs))]
    }
    rm(popup_plot_list)
    popup_plot_list <- popup_plot_list_weekly
    
    identical_plots <- FALSE
    
  } 
  else if (length(popup_plot_list) == nrow(polygon_shape)) identical_plots <- TRUE
  # here one assumes in second hierarchy a length of nrow(polygon_shape), which is not tested
  else if (length(popup_plot_list) == week_lead_max_ind) identical_plots <- FALSE
  else stop("popup_plot_list is not of length of 1 (or the anticipated weely mean time) times the nrow(polygon_shape).")

  # todo: think about doing the next two commands with group_map
  # estimate the probability of the most frequent class per catchment and store the resultant df per week in a list
  class_prob_list <- map(week_lead_ind, \(x) calc_max_class_prob_per_week(forecast_df, x))
  
  # calculate the colour with tint for each week and catchment
  # this returns a df of dfs with names(class_prob_list) in first hierarchy
  colour_df <- map(class_prob_list, \(x) find_class_colour_with_tint(levels, x, class_colours, ncolours_tint = 10)) %>% 
    list_cbind()
  
  # join colour with polygons via catchment ID to ensure that each colour is at the right place
  polygon_sf_weekly_list <- imap(week_lead_ind, 
                                 ~ full_join(polygon_sf, colour_df %>% pull(.y), by = join_by({{ ID_col }})))
  
  
  # build leaflet map
  leaf_map <- leaflet() %>% 
    addProviderTiles(providers$CartoDB.Positron) 
  
  # determine if a plot list for each week separately is provided
  
  
  # loop over all forecast weeks and add polygons
  for (i in week_lead_ind) {
    if (identical_plots) {
      leaf_map <- leaf_map %>% leaflet_add_polygon(polygon_sf_weekly_list[[i]], popup_plot_list, names(week_lead_ind[i]), 
                                                   pop_width = pop_width, pop_height = pop_height, fillOpacity = fillOpacity)
    } else {
      leaf_map <- leaf_map %>% leaflet_add_polygon(polygon_sf_weekly_list[[i]], popup_plot_list[[i]], names(week_lead_ind[i]), 
                                                   pop_width = pop_width, pop_height = pop_height, fillOpacity = fillOpacity)
    }
  }
  
  # add buttons for each week
  leaf_map <- leaf_map %>% addLayersControl(
    baseGroups = polygon_sf_weekly_list %>% names(),
    options = layersControlOptions(collapsed = FALSE))
  
  # Add common legend
  leaf_map <- leaf_map %>%
    addLegend(colors = class_colours,
              opacity = 1,
              labels = forecast_df$class %>% levels)
  
  return(leaf_map)
}



