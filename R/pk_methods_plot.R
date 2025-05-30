#' Plot a [pk()] object.
#'
#' Plot data and model fits from a [pk()] object.
#'
#' If the [pk()] object has not been fitted, then only the data will be plotted
#' (because no curve fits exist).
#'
#' @param x A [pk()] object. In this case `x` is used to align with generic method.
#' @param newdata Optional: A `data.frame` containing new data to plot. Must
#'  contain at least variables `Chemical`, `Species`, `Route`, `Media`, `Dose`,
#'  `Time`, `Time.Units`, `Conc`, `Detect`, `Conc_SD`.  Default `NULL`, to use
#'  the data in `obj$data`.
#' @param model Character: One or more of the models fitted. Curve fits will be
#'  plotted for these models. Default `NULL` to plot fits for all models in
#'  `x$stat_model`.
#' @param method Character: One or more of the [optimx::optimx()] methods used.
#'  Default `NULL` to plot fits for all methods in `x$settings_optimx$method`.
#' @param use_scale_conc Possible values: `TRUE`, `FALSE`, or a named list with
#'  elements `dose_norm` and `log10_trans` which themselves should be either
#'  `TRUE` or `FALSE`.  If `use_scale_conc = FALSE` (the default for this
#'  function), then the data and fits will be plotted without any
#'  dose-normalization or log-transformation. If `use_scale_conc = TRUE` , then
#'  the concentration scaling/transformations in `x` will be applied to the
#'  y-axis (concentration axis). If `use_scale_conc = list(dose_norm = ...,
#'  log10_trans = ...)`, then the specified dose normalization and/or
#'  log10-transformation will be applied to the y-axis (concentration axis) of
#'  the plots.
#' @param time_trans Default `FALSE`. Determines whether time values will be transformed.
#' @param log10_C Default `NULL`. Determines whether y-axis (concentration) should
#'  be log10 transformed. Takes `TRUE` or `FALSE` values. Otherwise it defaults
#'  to the value determined from `use_scale_conc`.
#' @param plot_data_aes Optional: Aesthetic mapping for the plot layer that
#'  visualizes the data. Default `NULL`, in which case a default mapping will be
#'  used based on the value of `use_scale_conc`.
#' @param plot_point_aes Optional: Aesthetic mappings for geom_point layer
#'  that determines the fill of the points. Defaults to `NULL`.
#' @param plot_fit_aes Optional: Aesthetic mapping for the plot layer that
#'  visualizes the fitted curves. Default `NULL`, in which case a default mapping will be
#'  used based on the value of `use_scale_conc`.
#' @param facet_fun  Default `"facet_grid"`. Optional: The name of the `ggplot2` faceting function to use:
#'  [ggplot2::facet_grid()], [ggplot2::facet_wrap()], or `'none'` to do no
#'  faceting. Default `NULL`, in which case a default faceting will be
#'  applied based on the value of `use_scale_conc`.
#' @param facet_fun_args A named list of arguments to the faceting function in
#'  `facet_fun` (if any). Default:
#' ```
#' list(rows = ggplot2::vars(Route),
#' cols= ggplot2::vars(Media),
#' scales = "free_y",
#' labeller = "label_both")
#' ```
#' @param drop_nonDetect Default `FALSE`. Whether to eliminate observations below
#'  the level of quantification (LOQ).
#' @param n_interp For plotting: the number of time points to interpolate between
#'  each observed time point. Default 10.
#' @param fit_limits Default `NULL`. c(Upper Bound, Lower Bound).
#'  Supply a numeric vector. These values filter the predicted
#'  values for fits to not exceed 2.25x of the maximum observed concentration values
#'  for each `data_group` in the `pk` object. When there is a log10 transformation
#'  of concentration values, it limits predicted values to 1/20th of the minimum
#'  observed concentration values and 5 times the maximum value.
#' @param print_out For plotting: whether the output of the function should be
#'  the list of plots. Default `FALSE`.
#' @param best_fit Default FALSE. Determines whether fit plot outputs only the
#'  best fit from `get_winning_model()`
#' @param ... Additional arguments not in use.
#' @return A [ggplot2::ggplot()]-class plot object.
#' @import ggplot2
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
plot.pk <- function(x,
                    newdata = NULL,
                    model = NULL,
                    method = NULL,
                    use_scale_conc = FALSE,
                    # Plotting arguments
                    time_trans = FALSE,
                    log10_C = NULL,
                    plot_data_aes = NULL,
                    plot_point_aes = NULL,
                    facet_fun = NULL,
                    facet_fun_args = NULL,
                    drop_nonDetect = FALSE,
                    # Predict/interpolation arguments
                    plot_fit_aes = NULL,
                    n_interp = 10,
                    fit_limits = NULL,
                    print_out = FALSE,
                    best_fit = FALSE,
                    ...) {
  # ensure that the model has at least been preprocessed
  check <- check_required_status(obj = x, required_status = status_preprocess)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(model))
    model <- names(x$stat_model)
  if (is.null(method))
    method <- x$settings_optimx$method
  if (is.null(newdata))
    newdata <- x$data

  newdata_ok <- check_newdata(
    newdata = newdata,
    olddata = x$data,
    req_vars = c("Time", "Time.Units", "Conc.Units",
                 "Dose", "Route", "Media"),
    exclude = FALSE
  )

  # check method and model
  method_ok <- check_method(obj = x, method = method)
  model_ok <- check_model(obj = x, model = model)

  # apply transformations if so specified
  conc_scale <- conc_scale_use(obj = x, use_scale_conc = use_scale_conc)

  if (drop_nonDetect %in% TRUE) {
    newdata <- subset(newdata, Detect %in% TRUE)
  }

  newdata <- subset(newdata, exclude %in% FALSE)

  if (is.null(log10_C))
    log10_C <- conc_scale$log10_trans


  if (is.null(fit_limits)) {
    fit_limits <- c(2.25, 0.05)
    if (log10_C) {
      fit_limits[1] <- 5
    }
  }

  if (is.numeric(fit_limits) & length(fit_limits) <= 2) {
    limit_predicted <- TRUE
    if (log10_C) {
      if (length(fit_limits) == 1) {
        fit_limits[2] <- 0.05
        warning("Lower Bound not defined, setting it to the default.")
      } else {
        message("Lower and Upper Bounds supplied will be used.")
      }
    } else {
      fit_limits <- fit_limits[1]
      warning("No lower bound needed, Upper Bounds set to first value supplied.")
    }
  } else {
    if (is.character(fit_limits) & (fit_limits == "none")) {
      limit_predicted <- FALSE
    } else {
      stop(
        "fit_limits must be numeric vector of length 1 or 2. Set fit_limits = 'none' to prevent possible predicted data filtering."
      )
    }
  }


  time_var <- ifelse(time_trans, as.name("Time_trans"), as.name("Time"))

  time_units_var <- ifelse(time_trans,
                           as.name("Time_trans.Units"),
                           as.name("Time.Units"))

  common_vars <- ggplot2::vars(Time,
                               Time.Units,
                               Dose,
                               Route,
                               Media)

  if (time_trans) {
    common_vars <- union(common_vars,
                         ggplot2::vars(Time_trans, Time_trans.Units))
  }

  obs_vars <- setdiff(
    ggplot2::vars(
    Conc,
    Conc_SD,
    Conc.Units,
    Value,
    Value.Units,
    Detect,
    exclude,
    Conc_trans,
    Conc_trans.Units,
    Reference
  ),
  x$data_group)


  # Default arguments and functions
  plot_data_aes_default <- ggplot2::aes(
    x = !!time_var,
    y = Conc_set,
    ymin = Conc_set - Conc_set_SD,
    ymax = Conc_set + Conc_set_SD,
    color = factor(Dose)
  )

  plot_point_aes_default <- ggplot2::aes(fill = factor(Dose),
                                         shape = Reference,
                                         alpha = Detect)

  plot_fit_aes_default <- ggplot2::aes(
    x = !!time_var,
    y = Conc_est,
    linetype = interaction(model, method)
  )

  if (conc_scale$dose_norm) {
    facet_fun_default <- "facet_grid"
    facet_fun_args_default <- list(
      rows = ggplot2::vars(Route),
      cols = ggplot2::vars(Media),
      scales = "free_y",
      labeller = "label_both"
    )
  } else {
    facet_fun_default <- "facet_wrap"
    facet_fun_args_default <- list(ggplot2::vars(Route, Media, Dose),
                                   scales = "free_y",
                                   labeller = "label_both")
  }

  if (is.null(plot_data_aes))
    plot_data_aes <- plot_data_aes_default
  if (is.null(plot_point_aes))
    plot_point_aes <- plot_point_aes_default
  if (is.null(plot_fit_aes))
    plot_fit_aes <- plot_fit_aes_default
  if (is.null(facet_fun))
    facet_fun <- facet_fun_default
  if (is.null(facet_fun_args))
    facet_fun_args <- facet_fun_args_default

  if (facet_fun %in% "none")
    facet_fun <- NULL

  # I think ideally there should be required variables and
  # a way to ensure all the aes() variables get added
  newdata <- newdata %>%
    dplyr::select(!!!union(union(x$data_group, common_vars), obs_vars)) %>%
    dplyr::rowwise() %>%
    # apply dose-norm if specified, but do not apply log10-trans even if specified
    # (if log10-trans specified, the y axis will be log10-scaled)
    dplyr::mutate(
      Conc_set = ifelse(conc_scale$dose_norm, Conc / Dose, Conc),
      Conc_set_SD = ifelse(conc_scale$dose_norm, Conc_SD / Dose, Conc_SD),
      Detect = ifelse(Detect, "Detect", "Non-Detect"),
      Dose = Dose
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!x$data_group) %>%
    tidyr::nest(.key = "observations")


  newdata <- newdata %>%
    dplyr::mutate(observation_plot =
                    purrr::map(observations, \(x) {
                      # Need to write this into a mutate + map pattern function
                      # initialize plot
                      p <- ggplot(data = x, mapping = plot_data_aes) +
                        geom_point(mapping = plot_point_aes, stroke = 1) + # plots points
                        geom_errorbar(mapping = plot_data_aes)

                      t_units <- x %>%
                        dplyr::pull(!!time_units_var) %>%
                        unique()

                      # if alpha = Detect then we have to do some trickery to implement that
                      # that is the objective of plot_point_aes
                      if (rlang::as_label(plot_point_aes$alpha) %in% "Detect") {
                        # plot the data
                        p <- p +
                          geom_point()
                        # plots here are un-filled points (only inherits color variable by default)

                        if ("shape" %in% c(names(plot_data_aes), names(plot_point_aes))) {
                          # Ensure shapes can take both fill and color if alpha is set
                          p <- p +
                            scale_shape_manual(values = 21:25)
                        }

                        # set an alpha scale
                        p <- p + scale_alpha_manual(
                          values = c("Detect" = 1, "Non-Detect" = 0),
                          breaks = c("Detect", "Non-Detect"),
                          drop = FALSE,
                          name = NULL
                        )

                        # Need to check whether there are two values in Detect
                        # some might be all Non-Detect or all Detect
                        if (length(unique(x$Detect)) > 1) {
                          alpha_fill <- c("black", NA)
                        } else {
                          alpha_fill <- ifelse(unique(x$Detect) == "Detect", "black", NA)
                        }
                        p <- p +
                          guides(alpha = guide_legend(
                            override.aes = list(
                              shape = 21,
                              color = "black",
                              stroke = 1,
                              fill = alpha_fill,
                              alpha = 1
                            ),
                            order = 2
                          ))

                      }

                      p <- p +
                        do.call(facet_fun, args = facet_fun_args)

                      if (log10_C) {
                        p <- p + scale_y_continuous(trans = "log10", labels = scales::label_log())
                      }

                      p +
                        labs(
                          #title = paste(Chemical, Species),
                          title = paste(!!!x$data_group),
                          x = paste0("Time (", t_units, ")"),
                          y = ifelse(
                            conc_scale$dose_norm,
                            "Concentration/Dose",
                            "Concentration"
                          )
                        ) +
                        theme_bw() +
                        theme(
                          panel.border = element_rect(
                            color = "black",
                            fill = NA,
                            linewidth = 1
                          ),
                          plot.title = element_text(hjust = 0.5, face = "bold"),
                          strip.background = element_rect(
                            fill = "white",
                            color = "black",
                            size = 1
                          )
                        )



                    }))

  # For predictions, interpolate time
  if (get_status(obj = x) == 5) {
    interp_data <- newdata

    interp_data <- interp_data %>%
      dplyr::mutate(interpolated = purrr::map(observations, \(x) {
        # if needed and not present, add Time_trans.Units column
        if (time_trans %in% TRUE) {
          if (!("Time_trans.Units" %in% names(x))) {
            x <- x %>%
              dplyr::mutate(Time_trans.Units = Time.Units)
          }

          x %>%
            dplyr::select(!!!union(common_vars, obs_vars)) %>%
            dplyr::group_by(Dose, Route, Media) %>%
            dplyr::reframe(Time = max(Time), # Change to Time
                           Time.Units, Time_trans.Units,
                           Conc.Units) %>%
            dplyr::mutate(
              maxTime = max(Time),
              Time.Units = unique(Time.Units),
              Time_trans.Units = unique(Time_trans.Units)
            ) %>%
            tidyr::uncount(n_interp) %>%
            dplyr::group_by(Dose, Route, Media) %>%
            dplyr::mutate(Time = (maxTime / (dplyr::n() - 1)) *
                            (dplyr::row_number() - 1))
        } else if (time_trans %in% FALSE) {
          x %>%
            dplyr::select(!!!union(common_vars, obs_vars)) %>%
            dplyr::group_by(Dose, Route, Media) %>%
            dplyr::reframe(Time = max(Time), # Change to Time
                           Time.Units, Conc.Units) %>%
            dplyr::mutate(maxTime = max(Time),
                          Time.Units = unique(Time.Units)) %>%
            tidyr::uncount(n_interp) %>%
            dplyr::group_by(Dose, Route, Media) %>%
            dplyr::mutate(Time = (maxTime / (dplyr::n() - 1)) *
                            (dplyr::row_number() - 1))
        }
      }))

    interp_data <- interp_data %>%
      dplyr::select(!!!x$data_group, interpolated) %>%
      tidyr::unnest(cols = c(interpolated)) %>%
      dplyr::ungroup()


    # Make predictions --
    # apply dose-normalization if it has been specified --
    # but do not use log10-transformation for now, even if it has been specified
    conc_scale_tmp <- conc_scale
    conc_scale_tmp$log10_trans <- FALSE

    interp_data <- predict.pk(
      object = x,
      newdata = interp_data,
      use_scale_conc = conc_scale_tmp,
      model = model,
      method = method,
      exclude = FALSE,
      include_NAs = TRUE
    )

    if (best_fit %in% TRUE) {
      interp_data <- dplyr::left_join(get_winning_model(obj = x), interp_data)
    }

    # if plotting transformed time, then transform interpolated time points
    if (time_trans %in% TRUE) {
      conversion_table <- time_conversions %>%
        dplyr::filter(
          TimeFrom %in% interp_data$Time.Units,
          TimeTo %in% interp_data$Time_trans.Units
        ) %>%
        dplyr::rename(Time.Units = "TimeFrom", Time_trans.Units = "TimeTo")

      interp_data <- interp_data %>%
        dplyr::left_join(conversion_table,
                         by = dplyr::join_by(Time.Units, Time_trans.Units)) %>%
        dplyr::mutate(Time_trans = Time * conversion) %>%
        dplyr::select(!conversion)
    }




    interp_data <- interp_data %>%
      dplyr::group_by(!!!x$data_group) %>%
      tidyr::nest(.key = "predicted")

    newdata <- dplyr::left_join(newdata, interp_data)

    # Need to check if predicted is NULL before filtering below
    # NULL predicted values will occur when there was no adequate fit for the model

    if (limit_predicted %in% TRUE) {
      if (log10_C %in% TRUE) {
        newdata <- newdata %>%
          dplyr::mutate(predicted = purrr::map2(observations, predicted, \(x, y) {
            if (!is.null(y)) {
              dplyr::filter(y,
                            Conc_est <= (fit_limits[1] * max(x$Conc_set)),
                            Conc_est >= (fit_limits[2] * min(x$Conc_set)))
            } else {
              message(paste(
                Chemical, Species, "did not have adequate fit"
              ))
            }

          }))
      } else {
        newdata <- newdata %>%
          dplyr::mutate(predicted = purrr::map2(observations, predicted, \(x, y) {
            if (!is.null(y)) {
              dplyr::filter(y, Conc_est <= (fit_limits[1] * max(x$Conc_set)))
            } else {
              message(paste(
                Chemical, Species, "did not have adequate fit"
              ))
            }
          }))
      }
    }

    newdata <- newdata %>%
      dplyr::mutate(predicted_plot = purrr::map(predicted, \(x) {
        if (!is.null(x)) {
        ggplot2::geom_line(data = x,
                           mapping = plot_fit_aes,
                           inherit.aes = FALSE)
        } else {
          NULL
        }
      })) %>%
      dplyr::mutate(final_plot = purrr::map2(
        observation_plot,
        predicted_plot,
        \(x, y) x + y +
          guides(
            color = guide_legend(title = "Dose", order = 1),
            fill = "none",
            linetype = guide_legend(title = "Model & Method", order = 3),
            shape = guide_legend(order = 4)
          )
      ))

  } else {
    newdata <- newdata %>%
      dplyr::rename(final_plot = "observation_plot")

    message("Note that the final plots do not contain any fits")
  }

  if (print_out)
    return(newdata$final_plot)

  return(newdata)
}
