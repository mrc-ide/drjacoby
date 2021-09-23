
#------------------------------------------------
#' @title Plot loglikelihood 95\% credible intervals
#'   
#' @description Plot loglikelihood 95\% credible intervals.
#'   
#' @param x an object of class \code{drjacoby_output}
#' @param chain which chain to plot.
#' @param phase which phase to plot. Must be either "burnin" or "sampling".
#' @param x_axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
#'   the thermodynamic power.
#' @param y_axis_type how to format the y-axis. 1 = raw values, 2 = truncated at
#'   auto-chosen lower limit. 3 = double-log scale.
#'
#' @export

plot_rung_loglike <- function(x, chain = 1, phase = "sampling", x_axis_type = 1, y_axis_type = 1) {

  # check inputs
  assert_class(x, "drjacoby_output")
  assert_single_pos_int(chain)
  assert_leq(chain, length(x))
  assert_in(phase, c("burnin", "sampling"))
  assert_single_pos_int(x_axis_type)
  assert_in(x_axis_type, 1:2)
  assert_single_pos_int(y_axis_type)
  assert_in(y_axis_type, 1:3)
  
  # get useful quantities
  thermo_power <- x$diagnostics$rung_details$thermodynamic_power
  rungs <- length(thermo_power)
  
  # define x-axis type
  if (x_axis_type == 1) {
    x_vec <- 1:rungs
    x_lab <- "rung"
  } else {
    x_vec <- thermo_power
    x_lab <- "thermodynamic power"
  }
  
  # get plotting values (loglikelihoods)
  phase_get <- phase
  chain_get <- chain
  data <- dplyr::filter(x$rung, .data$chain == chain_get, .data$phase == phase_get)
  y_lab <- "log-likelihood"
  
  # move to plotting deviance if specified
  if (y_axis_type == 3) {
    data$loglikelihood <- -2 * data$loglikelihood
    y_lab <- "deviance"
    
    # if needed, scale by adding/subtracting a power of ten until all values are
    # positive
    if (min(data$loglikelihood) < 0) {
      dev_scale_power <- ceiling( log(abs(min(data$loglikelihood))) / log(10) )
      dev_scale_sign <- -sign(min(data$loglikelihood))
      data$loglikelihood <- data$loglikelihood + dev_scale_sign*10^dev_scale_power
      
      dev_scale_base <- ifelse(dev_scale_power == 0, 1, 10)
      dev_scale_power_char <- ifelse(dev_scale_power <= 1, "", paste("^", dev_scale_power))
      dev_scale_sign_char <- ifelse(dev_scale_sign < 0, "-", "+")
      y_lab <- parse(text = paste("deviance", dev_scale_sign_char, dev_scale_base, dev_scale_power_char))
    }
  }
  
  # get 95% credible intervals over plotting values
  y_intervals <- data %>%
    dplyr::group_by(.data$rung) %>%
    dplyr::summarise(Q2.5 = quantile(.data$loglikelihood, 0.025),
                     Q50 =  quantile(.data$loglikelihood, 0.5),
                     Q97.5 = quantile(.data$loglikelihood, 0.975))
  
  # get data into ggplot format and define temperature colours
  df <- y_intervals
  df$col <- thermo_power
  df$x <- x_vec
  
  # produce plot
  plot1 <- df %>% 
    ggplot() + 
    geom_vline(aes(xintercept = x_vec), col = grey(0.9)) +
    geom_pointrange(aes_(x = ~x, ymin = ~Q2.5, y = ~Q50, ymax = ~Q97.5, col = ~col)) + 
    xlab(x_lab) + 
    ylab(y_lab) + 
    scale_colour_gradientn(colours = c("red", "blue"), name = "thermodynamic\npower", limits = c(0,1)) +
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  # define y-axis
  if (y_axis_type == 2) {
    y_min <- quantile(df$Q2.5, probs = 0.5)
    y_max <- max(df$Q97.5)
    plot1 <- plot1 + coord_cartesian(ylim = c(y_min, y_max))
  } else if (y_axis_type == 3) {
    plot1 <- plot1 + scale_y_continuous(trans = "log10")
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot Metropolis coupling acceptance rates
#'
#' @description Plot Metropolis coupling acceptance rates between all rungs.
#'
#' @inheritParams plot_rung_loglike
#' @param chain which chain to plot. If \code{NULL} then plot all chains.
#'
#' @import ggplot2 dplyr
#' @importFrom grDevices grey
#' @export

plot_mc_acceptance <- function(x, chain = NULL, phase = "sampling", x_axis_type = 1) {
  
  # check inputs
  assert_class(x, "drjacoby_output")
  if (!is.null(chain)) {
    assert_single_pos_int(chain, zero_allowed = FALSE)
    assert_in(chain, unique(x$output$chain))
  }
  assert_in(phase, c("burnin", "sampling"))
  assert_single_pos_int(x_axis_type)
  assert_in(x_axis_type, 1:2)
  
  # get useful quantities
  thermo_power <- x$diagnostics$rung_details$thermodynamic_power
  thermo_power_mid <- thermo_power[-1] - diff(thermo_power)/2
  rungs <- length(thermo_power)
  
  # exit if rungs = 1
  if (rungs == 1) {
    stop("no metropolis coupling when rungs = 1")
  }
  
  # define x-axis type
  if (x_axis_type == 1) {
    breaks_vec <- 1:rungs
    x_vec <- (2:rungs) - 0.5
    x_lab <- "rung"
  } else {
    breaks_vec <- thermo_power
    x_vec <- thermo_power_mid
    x_lab <- "thermodynamic power"
  }
  
  # get chain properties
  phase_get <- phase
  if (is.null(chain)) {
    mc_accept <- dplyr::filter(x$diagnostics$mc_accept, .data$phase == phase_get) %>%
      dplyr::pull(.data$value)
  } else {
    chain_get <- chain
    mc_accept <- dplyr::filter(x$diagnostics$mc_accept, .data$phase == phase_get, .data$chain == chain_get) %>%
      dplyr::pull(.data$value)
  }
  
  # get data into ggplot format and define temperature colours
  df <- data.frame(x_vec = x_vec, mc_accept = mc_accept, col = thermo_power_mid)
  
  # produce plot
  plot1 <- ggplot(df) + 
    geom_vline(aes(xintercept = x), col = grey(0.9), data = data.frame(x = breaks_vec)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
    geom_point(aes(x = x_vec, y = mc_accept, color = col)) + 
    xlab(x_lab) + ylab("coupling acceptance rate") + 
    scale_colour_gradientn(colours = c("red", "blue"), name = "thermodynamic\npower", limits = c(0,1)) +
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot autocorrelation
#'
#' @description Plot autocorrelation for specified parameters
#'
#' @inheritParams plot_rung_loglike
#' @param lag maximum lag. Must be an integer between 1 and 500.
#' @param par vector of parameter names. If \code{NULL} all parameters are
#'   plotted.
#' 
#' @export

plot_autocorrelation <- function(x, lag = 20, par = NULL, chain = 1, phase = "sampling") {
  
  # check inputs
  assert_class(x, "drjacoby_output")
  assert_single_bounded(lag, 1, 500)
  if (is.null(par)) {
    par <- setdiff(names(x$output), c("chain", "iteration", "phase",
                                      "logprior", "loglikelihood"))
  }
  assert_vector_string(par)
  assert_in(par, names(x$output))
  assert_single_pos_int(chain)
  assert_leq(chain, length(x))
  assert_in(phase, c("burnin", "sampling"))

  # get output for the chosen chain, phase
  chain_get <- chain
  phase_get <- phase
  data <- dplyr::filter(x$output, .data$chain == chain_get, .data$phase == phase_get) %>%
    dplyr::select(-.data$chain, -.data$iteration, -.data$phase, -.data$logprior, -.data$loglikelihood) %>%
    as.data.frame()
  
  # select parameters
  data <- data[, par, drop = FALSE]
  
  # estimate autocorrelation
  out <- as.data.frame(apply(data, 2, acf_data, lag = lag))
  
  # format data for plotting
  out$lag <- 0:lag
  out <- do.call(rbind, mapply(function(i) {
    data.frame(lag = out$lag,
               parameter = names(out)[i],
               autocorrelation = out[,i])
  }, seq_len(ncol(data)), SIMPLIFY = FALSE))
  
  # produce plot
  ggplot2::ggplot(data = out,
                  ggplot2::aes(x = .data$lag, y = 0, xend = .data$lag, yend =.data$autocorrelation)) + 
    ggplot2::geom_hline(yintercept = 0, lty = 2, col = "red") + 
    ggplot2::geom_segment(size = 1.5) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Autocorrelation") +
    ggplot2::xlab("Lag") +
    ggplot2::ylim(min(0, min(out$autocorrelation)), 1) +
    ggplot2::facet_wrap(~ .data$parameter)
}

#------------------------------------------------
#' @title Plot parameter estimates
#'
#' @description Produce a series of plots corresponding to each parameter,
#'   including the raw trace, the posterior histogram and an autocorrelation
#'   plot. Plotting objects can be cycled through interactively, or can be
#'   returned as an object allowing them to be viewed/edited by the user.
#'
#' @inheritParams plot_rung_loglike
#' @param show optional vector of parameter names to plot. Parameters matching
#'   show will be included.
#' @param hide optional vector of parameter names to filter out. Parameters
#'   matching hide will be hidden.
#' @param lag maximum lag. Must be an integer between 1 and 500.
#' @param downsample boolean. Whether to downsample chain to make plotting more
#'   efficient.
#' @param rung which temperature rung to plot. If \code{NULL} then defaults to
#'   the cold rung.
#' @param display boolean. Whether to show plots, if \code{FALSE} then plotting
#'   objects are returned without displaying.
#'
#' @export

plot_par <- function(x, show = NULL, hide = NULL, lag = 20,
                     downsample = TRUE, phase = "sampling", rung = NULL,
                     chain = NULL, display = TRUE) {
  
  # check inputs and define defaults
  assert_class(x, "drjacoby_output")
  assert_non_null(x$output)
  param_names <- setdiff(names(x$output), c("chain", "rung", "iteration", "phase", "logprior", "loglikelihood"))
  if (!is.null(show)) {
    assert_vector_string(show)
    assert_in(show, param_names)
    param_names <- show
  }
  if (!is.null(hide)) {
    assert_vector_string(hide)
    assert_in(hide, param_names)
    param_names <- setdiff(param_names, hide)
  }
  assert_gr(length(param_names), 0, message = "at least one parameter must be specified")
  if (length(param_names) > 10){
    message("More than 10 parameters to summarise, consider using the show or hide arguments 
            to select parameters and reduce computation time.")
  }
  assert_single_bounded(lag, 1, 500)
  assert_single_logical(downsample)
  assert_in(phase, c("burnin", "sampling", "both"))
  max_rungs <- max(x$diagnostics$rung_details$rung)
  if (is.null(rung)) {
    rung <- max_rungs
  }
  assert_single_pos_int(rung)
  assert_leq(rung, max_rungs)
  if (is.null(chain)) {
    chain <- unique(x$output$chain)
  }
  assert_pos_int(chain, zero_allowed = FALSE)
  assert_single_logical(display)
  
  # deal with phase = "both" situation
  if (phase == "both") {
    phase <- c("burnin", "sampling")
  }
  
  # subset based on rung, chain and phase
  rung_get <- rung
  chain_get <- chain
  phase_get <- phase
  data <- dplyr::filter(x$output, rung == rung_get, phase %in% phase_get, chain %in% chain_get) 
  
  # Downsample
  if (downsample & nrow(data) > 2000) {
    data <- data[round(seq(1, nrow(data), length.out = 2000)),]
  }
  
  # Autocorrealtion (on downsample)
  ac_data <- as.data.frame(apply(data[, param_names, drop = FALSE], 2, acf_data, lag = lag))
  ac_data$lag <- 0:lag
  
  # Set minimum bin number
  b <- min(nrow(data) / 4, 40)
  
  # produce plots over all parameters
  plot_list <- c()
  for(j in 1:length(param_names)){
    
    # create plotting data
    pd <- data[, c("chain", "iteration", param_names[j])]
    names(pd) <- c("chain", "iteration", "y")
    
    pd2 <- ac_data[, c("lag", param_names[j])]
    names(pd2) <- c("lag", "Autocorrelation")
    
    # Histogram
    p1 <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$y)) + 
      ggplot2::geom_histogram(bins = b, fill = "deepskyblue3", col = "darkblue") + 
      ggplot2::ylab("Count") + 
      ggplot2::xlab(param_names[j]) +
      ggplot2::theme_bw()
    
    # Trace plots
    p2 <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$iteration, y = .data$y, col = as.factor(.data$chain))) + 
      ggplot2::geom_line() +
      scale_color_discrete(name = "Chain") +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(param_names[j]) +
      ggplot2::theme_bw()
    
    # Autocorrealtion
    p3 <- ggplot2::ggplot(data = pd2,
                          ggplot2::aes(x = .data$lag, y = 0, xend = .data$lag, yend =.data$Autocorrelation)) + 
      ggplot2::geom_hline(yintercept = 0, lty = 2, col = "red") + 
      ggplot2::geom_segment(size = 1.5) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Autocorrelation") +
      ggplot2::xlab("Lag") +
      ggplot2::ylim(min(0, min(pd2$Autocorrelation)), 1)
    
    # Arrange
    pc1 <- cowplot::plot_grid(p1, p3, ncol = 2)
    pc2 <- cowplot::plot_grid(p2, pc1, nrow = 2)
    plot_list[[j]] <- list(trace = p2,
                           hist = p1,
                           acf = p3,
                           combined = pc2)
    
    # Display plots, asking user for next page if multiple parameters
    if(display){
      graphics::plot(plot_list[[j]]$combined)
      if (j < length(param_names)) {
        z <- readline("Press n for next plot, f to return the list of all plots or any other key to exit ")
        if(z == "f"){
          display <- FALSE
        } 
        if(!z %in% c("n", "f")){
          return(invisible())            
        }
      }
    }
  }
  names(plot_list) <- paste0("Plot_", param_names)
  
  return(invisible(plot_list))
}

#------------------------------------------------
#' @title Plot parameter correlation
#'
#' @description Plots correlation between two parameters
#'
#' @inheritParams plot_rung_loglike
#' @param parameter1 Name of parameter first parameter.
#' @param parameter2 Name of parameter second parameter.
#' @param downsample Downsample chain for efficiency
#' @param rung Rung
#'
#' @export

plot_cor <- function(x, parameter1, parameter2,
                     downsample = TRUE, phase = "sampling",
                     chain = NULL, rung = NULL) {

  # check inputs
  assert_class(x, "drjacoby_output")
  assert_single_string(parameter1)
  assert_single_string(parameter2)
  assert_in(parameter1, names(x$output))
  assert_in(parameter2, names(x$output))
  assert_single_logical(downsample)
  assert_in(phase, c("burnin", "sampling", "both"))
  if (is.null(chain)) {
    chain <- unique(x$output$chain)
  }
  assert_pos_int(chain)
  max_rungs <- max(x$diagnostics$rung_details$rung)
  if (is.null(rung)) {
    rung <- max_rungs
  }
  assert_single_pos_int(rung)
  assert_leq(rung, max_rungs)
  
  # deal with phase = "both" situation
  if (phase == "both") {
    phase <- c("burnin", "sampling")
  }
  
  # get basic quantities
  rung_get <- rung
  chain_get <- chain
  phase_get <- phase
  data <- dplyr::filter(x$output, rung == rung_get, phase %in% phase_get, chain %in% chain_get)
  
  # subset to corr params
  data <- data[,c("chain", parameter1, parameter2)]  
  colnames(data) <- c("chain", "x", "y")
  
  # downsample
  if (downsample & nrow(data) > 2000) {
    data <- data[round(seq(1, nrow(data), length.out = 2000)),]
  }
  
  # produce plot
  ggplot2::ggplot(data = data,
                  ggplot2::aes(x = .data$x, y = .data$y, col = as.factor(.data$chain))) + 
    ggplot2::geom_point(alpha = 0.5) + 
    ggplot2::xlab(parameter1) +
    ggplot2::ylab(parameter2) +
    scale_color_discrete(name = "Chain") +
    ggplot2::theme_bw()
  
}

#------------------------------------------------
#' @title Plot 95\% credible intervals
#'
#' @description Plots posterior 95\% credible intervals over specified set of
#'   parameters (defauls to all parameters).
#' 
#' @inheritParams plot_rung_loglike
#' @param show Vector of parameter names to plot.
#' @param rung Which rung to plot.
#' @param param_names Optional vector of names to replace the default parameter names.
#'
#' @export

plot_credible <- function(x, show = NULL, phase = "sampling", rung = NULL, param_names = NULL) {
  
  # check inputs
  assert_class(x, "drjacoby_output")
  if (!is.null(show)) {
    assert_string(show)
    assert_in(show, names(x$output))
  }
  assert_in(phase, c("burnin", "sampling", "both"))
  max_rungs <- max(x$diagnostics$rung_details$rung)
  if (is.null(rung)) {
    rung <- max_rungs
  }
  assert_single_pos_int(rung)
  assert_leq(rung, max_rungs)
  
  # define defaults
  if (is.null(show)) {
    show <- setdiff(names(x$output), c("chain", "rung", "iteration", "phase", "logprior", "loglikelihood"))
  }
  if (is.null(param_names)) {
    param_names <- show
  }
  
  # deal with phase = "both" situation
  if (phase == "both") {
    phase <- c("burnin", "sampling")
  }
  
  # subset based on phase and rung
  rung_get <- rung
  phase_get <- phase
  data <- dplyr::filter(x$output, rung == rung_get, phase %in% phase_get)
  data <- data[, show, drop = FALSE]
  
  # get quantiles
  df_plot <- as.data.frame(t(apply(data, 2, quantile_95)))
  df_plot$param <- factor(param_names, levels = param_names)
  
  # produce plot
  ggplot2::ggplot(data = df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(x = .data$param, y = .data$Q50)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$param, y = .data$Q2.5, xend = .data$param, yend = .data$Q97.5)) +
    ggplot2::xlab("") +
    ggplot2::ylab("95% CrI")
  
}

#------------------------------------------------
#' @title Plot posterior correlation matrix
#'
#' @description Produces a matrix showing the correlation between all parameters
#'   from posterior draws.
#' 
#' @inheritParams plot_rung_loglike
#' @param show Vector of parameter names to plot.
#' @param param_names Optional vector of names to replace the default parameter names.
#'
#' @importFrom stats cor
#' @export

plot_cor_mat <- function(x, show = NULL, phase = "sampling", param_names = NULL) {

  # check inputs
  assert_class(x, "drjacoby_output")
  if (!is.null(show)) {
    assert_string(show)
    assert_in(show, names(x$output))
    assert_gr(length(show), 1, message = "must show at least two parameters")
  }
  assert_in(phase, c("burnin", "sampling", "both"))

  
  # define defaults
  if (is.null(show)) {
    show <- setdiff(names(x$output), c("chain", "iteration", "phase", "logprior", "loglikelihood"))
  }
  if (is.null(param_names)) {
    param_names <- show
  }
  
  # deal with phase = "both" situation
  if (phase == "both") {
    phase <- c("burnin", "sampling")
  }
  
  # subset based on phase
  phase_get <- phase
  data <- dplyr::filter(x$output, phase %in% phase_get)
  data <- data[, show, drop = FALSE]
  n <- ncol(data)
  
  # get correlation matrix into dataframe for ggplot
  m <- cor(data)
  df_plot <- data.frame(x = rep(names(data), each = n),
                        xi = rep(1:n, each = n),
                        y = names(data),
                        yi = 1:n,
                        z = as.vector(m))
  df_plot <- subset(df_plot, df_plot$xi < df_plot$yi)
  df_plot$x <- factor(df_plot$x, levels = names(data))
  df_plot$y <- factor(df_plot$y, levels = names(data))
  
  # get colour range
  max_range <- max(abs(range(df_plot$z)))
  max_plot <- ceiling(max_range * 10) / 10
  
  # produce plot
  ggplot2::ggplot(df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_raster(ggplot2::aes_(x = ~x, y = ~y, fill = ~z)) +
    ggplot2::scale_fill_gradientn(colours = c("red", "white", "blue"),
                                  values = c(0, 0.5, 1),
                                  limits = c(-max_plot, max_plot),
                                  name = "correlation") +
    ggplot2::xlab("") + ggplot2::ylab("")
  
}
