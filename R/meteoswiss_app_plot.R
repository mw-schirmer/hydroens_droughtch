#' A reprogramming exercise of the meteoswiss app plot on daily values
#'
#' @param mach_ID catchment ID shown in the title
#' @param owner_ID secondary catchment ID shown in the title
#' @param forecast_df a data frame which contains a column "class" with forecast categories as factors,
#'  and a "date" column of class date.
#' @param precip_breaks breaks for the precipitation legend
#' @param precip_labels labels for the precipitation legend
#' @param precip_colours colours for the precipitation legend
#' @param transfer_add additive transfer function between precipitation and temperature axis
#' @param transfer_mult multiplicative transfer function between precipitation and temperature axis
#'
#' @return a ggplot2 object
#'
#' @author Michael Schirmer, michael.schirmer[at]wsl.ch
#'
#' @examples meteoswiss_app_plot("CHFO-1234", "2589", forecast_extract)
meteoswiss_app_plot <- function(mach_ID, owner_ID,
                                forecast_df,
                                precip_breaks = c(-Inf, 2, 4, 10, 20, 50, 100, Inf),
                                precip_labels = c("<2", "<4", "<10", "<20", "<50", "<100", ">100"),
                                precip_colours = c("lightpink4", "blue", "darkgreen", "green", "yellow", "orange", "red"),
                                transfer_add = 10, transfer_mult = 1
) {

  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(colorspace)
  library(patchwork)

  # title
  init_date <- min(forecast_df$date)
  title_str <- paste(owner_ID, init_date, sep = " - ")

  # data for one catchment
  is_relevant <- forecast_df$mach_ID == mach_ID
  forecast_catchment_df <- forecast_df[is_relevant,] %>% ungroup()

  # calculate quantiles
  probs <- c(quantile10 = .1, quantile50 = .5, quantile90 = .9)
  plot_data_perc <- forecast_catchment_df %>% group_by(date) %>%
    group_modify(~tibble(values = quantile(.x$P, probs = probs), probs = names(probs))) %>%
    pivot_wider(names_from = probs, values_from = values)

  # add to forecast_catchment_df
  forecast_catchment_df <- left_join(forecast_catchment_df, plot_data_perc, by = join_by(date)) %>%
    rename(P_ucl = quantile10, P_median = quantile50, P_ucu = quantile90)

  # same for T
  plot_data_perc <- forecast_catchment_df %>% group_by(date) %>%
    group_modify(~tibble(values = quantile(.x$T, probs = probs), probs = names(probs))) %>%
    pivot_wider(names_from = probs, values_from = values)

  # add to forecast_catchment_df
  forecast_catchment_df <- left_join(forecast_catchment_df, plot_data_perc, by = join_by(date)) %>%
    rename(T_ucl = quantile10, T_median = quantile50, T_ucu = quantile90) %>%
    filter(member == "ctr")

  # Identify days where the value exceeds the maximum for each category
  exceeds_max <- map(seq(2, length(precip_breaks) - 1), \(cat) forecast_catchment_df$P_median > precip_breaks[cat]) %>%
    reduce(cbind)

  # add days to data frame with exceeding breaks, to achieve multiple coloured bars
  df_breaks_added <- forecast_catchment_df
  for (i in seq(ncol(exceeds_max), 1)) {
    df_breaks_added <- df_breaks_added %>% bind_rows(filter(forecast_catchment_df, exceeds_max[, i]) %>% mutate(P_median = precip_breaks[i + 1]))
  }

  # categorize all real data and added data
  df_breaks_added <- df_breaks_added %>%
    mutate(P_cat = cut(P_median, breaks = precip_breaks, labels = precip_labels))

  # Create the plot
  gg <- ggplot(df_breaks_added) +

    # mark 0 lines of both y axes
    geom_hline(yintercept = 0, linewidth = 1, colour = "grey") +
    geom_hline(yintercept = (0 + transfer_add) * transfer_mult, linewidth = 1, colour = "grey") +

    # bars with height value and colours of categories
    geom_col(aes(x = date, y = P_median, fill = P_cat), position = "identity", width = 1) +
    scale_fill_manual(name = "P categories", values = precip_colours) +

    # precip lower uc
    geom_col(data = forecast_catchment_df, aes(x = date, y = P_ucl),
             position = "identity", fill = "grey90", alpha = 0.7, width = 1) +

    # precip upper uc
    # if I plot this last, all bars are nicely separated with a grey line
    geom_col(data = forecast_catchment_df, aes(x = date, y = P_ucu),
             position = "identity", fill = NA, colour = "grey", width = 1) +

    # temperature
    geom_line(aes(x = date, y = (T_median + transfer_add) * transfer_mult), color = "red", linewidth = 1) +
    geom_ribbon(aes(x = date, ymin = (T_ucl + transfer_add) * transfer_mult, ymax = (T_ucu + transfer_add) * transfer_mult),
                fill = "red", alpha = 0.2) +

    # as I need to plot precip in first axis because of 0 start values for bars, I put first axis right, second axis left
    scale_y_continuous(expand = expansion(),
                       name = "mm", position = "right",
                       sec.axis = sec_axis(name = "°C", ~ . / transfer_mult - transfer_add,
                                           guide = guide_axis(position = "left"))) +
    # not expand beyond the data
    scale_x_date(
      date_labels = "%d.%m",
      breaks = scales::breaks_width(width = "week"),
      expand = expansion(mult = 0, add = 0)) +

    theme_minimal() +
    theme(axis.title.y.left = element_text(
      angle = 0,
      face = "bold",
      vjust = 1,
      colour = "red",
      # to center it above labels
      margin = margin(r = -15, l = 15)
    ),
    # larger text size
    # text = element_text(size = 30),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(colour = "grey60"),
    # longer ticks, needs to be tested
    # axis.ticks.length.x.bottom = unit(.2, "cm"),
    axis.text = element_text(colour = "grey60"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    )

  # data for the coloured precip y-axis
  df_ybar <- tibble(
    ymin = c(0, precip_breaks[seq(2, length(precip_breaks) - 1)]),
    ymax = c(precip_breaks[seq(2, length(precip_breaks))]),
    b = factor(seq_along(ymin)))

  # coloured precip y-axis
  precip_ybar <- ggplot(df_ybar, aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = b)) +
    geom_rect() +
    scale_fill_manual(values = precip_colours) +
    scale_y_continuous(expand = expansion(),
                       name = "mm", position = "right") +
    coord_cartesian(ylim = layer_scales(gg)$y$get_limits()) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.y.right = element_text(
            angle = 0,
            face = "bold",
            vjust = 1,
            colour = "#5b9bd5",
            margin = margin(l = -15, r = 15)
          ),
          # text = element_text(size = 30),
          axis.text = element_text(colour = "grey60"),
          axis.text.x = element_blank(),
          panel.grid = element_blank())


  # patch the plots
  gg_out <- gg + ggtitle(title_str) + precip_ybar + plot_layout(widths = c(1, .03))

  return(gg_out)

}


