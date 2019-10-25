
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
  assert_custom_class(x, "drjacoby_output")
  assert_single_pos_int(chain)
  assert_leq(chain, length(x))
  assert_in(phase, c("burnin", "sampling"))
  assert_single_pos_int(x_axis_type)
  assert_in(x_axis_type, 1:2)
  assert_single_pos_int(y_axis_type)
  assert_in(y_axis_type, 1:3)
  
  # get useful quantities
  beta_raised <- x[[chain]]$diagnostics$beta_raised
  rungs <- length(beta_raised)
  
  # define x-axis type
  if (x_axis_type == 1) {
    x_vec <- rungs:1
    x_lab <- "rung"
  } else {
    x_vec <- beta_raised
    x_lab <- "thermodynamic power"
  }
  
  # get loglikelihoods
  if (phase == "burnin") {
    loglike <- x[[chain]]$loglike_burnin
  } else {
    loglike <- x[[chain]]$loglike_sampling
  }
  
  # take double-logs if needed
  y_lab <- "log-likelihood"
  if (y_axis_type == 3) {
    loglike <- -2*loglike
    y_lab <- "deviance"
  }
  
  # get 95% credible intervals over sampling loglikelihoods
  loglike_intervals <- t(apply(loglike, 2, quantile_95))
  
  # get data into ggplot format and define temperature colours
  df <- as.data.frame(loglike_intervals)
  df$col <- beta_raised
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw() + theme(panel.grid.minor.x = element_blank(),
                                         panel.grid.major.x = element_blank())
  plot1 <- plot1 + geom_vline(aes(xintercept = x_vec), col = grey(0.9))
  plot1 <- plot1 + geom_segment(aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = ~x_vec, y = ~Q50, color = ~col))
  plot1 <- plot1 + xlab(x_lab) + ylab(y_lab)
  plot1 <- plot1 + scale_colour_gradientn(colours = c("red", "blue"), name = "thermodynamic\npower", limits = c(0,1))
  
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
#'
#' @import ggplot2
#' @importFrom grDevices grey
#' @export

plot_mc_acceptance <- function(x, chain = 1, phase = "sampling", x_axis_type = 1) {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_single_pos_int(chain)
  assert_leq(chain, length(x))
  assert_in(phase, c("burnin", "sampling"))
  assert_single_pos_int(x_axis_type)
  assert_in(x_axis_type, 1:2)
  
  # get useful quantities
  beta_raised <- x[[chain]]$diagnostics$beta_raised
  beta_raised_mid <- beta_raised[-1] - diff(beta_raised)/2
  rungs <- length(beta_raised)
  
  # define x-axis type
  if (x_axis_type == 1) {
    breaks_vec <- rungs:2
    x_vec <- (rungs:2) - 0.5
    x_lab <- "rung"
  } else {
    breaks_vec <- beta_raised
    x_vec <- beta_raised_mid
    x_lab <- "thermodynamic power"
  }
  
  # get acceptance rates
  if (phase == "burnin") {
    mc_accept <- x[[chain]]$diagnostics$mc_accept_burnin
  } else {
    mc_accept <- x[[chain]]$diagnostics$mc_accept_sampling
  }
  
  # get data into ggplot format and define temperature colours
  df <- as.data.frame(mc_accept)
  df$col <- beta_raised_mid
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw() + theme(panel.grid.minor.x = element_blank(),
                                         panel.grid.major.x = element_blank())
  plot1 <- plot1 + geom_vline(aes(xintercept = breaks_vec), col = grey(0.9))
  plot1 <- plot1 + scale_y_continuous(limits = c(0,1), expand = c(0,0))
  plot1 <- plot1 + geom_point(aes(x = x_vec, y = mc_accept, color = col))
  plot1 <- plot1 + xlab(x_lab) + ylab("coupling acceptance rate")
  plot1 <- plot1 + scale_colour_gradientn(colours = c("red", "blue"), name = "thermodynamic\npower", limits = c(0,1))
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot autocorrelation
#'
#' @description Plot autocorrelation for specified parameters
#'
#' @inheritParams plot_rung_loglike
#' @param lag Maximum lag. Must be an integer between 20 and 500.
#' @param par Vector of parameter names. If NULL all parameters are plotted
#'
#' @export
 
plot_autocorrelation <- function(x, lag = 20, par = NULL, chain = 1, phase = "sampling") {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_single_pos_int(chain)
  assert_leq(chain, length(x))
  assert_in(phase, c("burnin", "sampling"))
  assert_single_bounded(lag, 1, 500)
  
  # get values
  if (phase == "burnin") {
    data <- x[[chain]]$theta_burnin$rung1
  } else {
    data <- x[[chain]]$theta_sampling$rung1
  }
  
  # Select parameters
  if(!is.null(par)){
    data <- data[, colnames(data) %in% par, drop = FALSE]
  }
  
  # Estimate autocorrelation
  out <- as.data.frame(apply(data, 2, acf_data, lag = lag))
  
  # Format data for plotting
  out$lag <- 0:lag
  out <- tidyr::gather(out, "parameter", "Autocorrelation", -lag)
  
  ggplot2::ggplot(data = out,
                  ggplot2::aes(x = .data$lag, y = 0, xend = .data$lag, yend =.data$Autocorrelation)) + 
    ggplot2::geom_hline(yintercept = 0, lty = 2, col = "red") + 
    ggplot2::geom_segment(size = 1.5) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Autocorrelation") +
    ggplot2::xlab("Lag") +
    ggplot2::ylim(min(0, min(out$Autocorrelation)), 1) +
    ggplot2::facet_wrap(~ parameter)
}

#------------------------------------------------
#' @title Plot parameter estimates
#'
#' @description Plot parameter estimates
#'
#' @inheritParams plot_rung_loglike
#' @param show Optional character (or vector of characters) to filter parameters by.
#'  Parameters matching show will be included.
#' @param hide Optional character (or vector of characters) to filter parameters by.
#'  Parameters matching show will be hidden.
#' @param lag Maximum lag. Must be an integer between 20 and 500
#' @param downsample Downsample chain for efficiency
#'
#' @export

plot_par <- function(x, show = NULL, hide = NULL, lag = 20,
                     downsample = TRUE){
  
  assert_custom_class(x, "drjacoby_output")
  assert_single_bounded(lag, 1, 500)
  
  parameter <- names(x$chain1$theta_sampling$rung1)
  if(!is.null(show)){
    stopifnot(is.character(show))
    parameter <- parameter[grepl(paste(show, collapse = "|"), parameter)]
  }
  if(!is.null(hide)){
    stopifnot(is.character(hide))
    parameter <- parameter[!grepl(paste(hide, collapse = "|"), parameter)]
  }
  if(length(parameter) > 10){
    message("More than 10 parameters to summarise, consider using the show or hide arguments 
            to select parameters and reduce computation time.")
  }
  
  all_chains <- dplyr::bind_rows(lapply(x, function(y){
    y$theta_sampling$rung1
  }))
  # Add chain
  all_chains$chain <- factor(rep(1:length(x), each = nrow(x[[1]]$theta_sampling$rung1)))
  # Downsample
  if(nrow(all_chains) > 2000){
    all_chains <- all_chains[seq.int(1, nrow(all_chains), length.out = 2000),]
  }
  all_chains <- dplyr::group_by(all_chains, chain)
  all_chains <- dplyr::mutate(all_chains, x = 1:dplyr::n())
  all_chains <- dplyr::ungroup(all_chains)
  
  # Autocorrealtion (on downsample)
  ac_data <- as.data.frame(apply(dplyr::select(all_chains, -chain), 2, acf_data, lag = lag))
  ac_data$lag <- 0:lag

  # Set minimum bin number
  b <- min(nrow(all_chains) / 4, 40)
  
  plot_list <- c()
  for(j in 1:length(parameter)){
    pd <- all_chains[, c("chain", "x", parameter[j])]
    names(pd) <- c("chain", "x", "y")
    pd2 <- ac_data[, c("lag", parameter[j])]
    names(pd2) <- c("lag", "Autocorrelation")
    # Histogram
    p1 <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$y)) + 
      ggplot2::geom_histogram(bins = b, fill = "deepskyblue3", col = "darkblue") + 
      ggplot2::ylab("Count") + 
      ggplot2::xlab(parameter[j]) +
      ggplot2::theme_bw()
    # Chains
    p2 <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$x, y = .data$y, col = .data$chain)) + 
      ggplot2::geom_line() +
      scale_color_discrete(name = "Chain") +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(parameter[j]) +
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
    
    pc1 <- cowplot::plot_grid(p1, p3, ncol = 2)
    # Side by side
    plot_list[[j]] <- cowplot::plot_grid(p2, pc1, nrow = 2)
  }
  names(plot_list) <- paste0("Plot_", parameter)
  
  # Display plots, asking user for next page if multiple parameters
  for(j in 1:length(parameter)){
    graphics::plot(plot_list[[j]])
    if(j == 1){
      default_ask <- grDevices::devAskNewPage()
      on.exit(grDevices::devAskNewPage(default_ask))
      grDevices::devAskNewPage(ask = TRUE)
    }
  }
  
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
#'
#' @export

plot_cor <- function(x, parameter1, parameter2,
                     downsample = TRUE){
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_string(parameter1)
  assert_string(parameter2)
  if(!parameter1 %in% names(x$chain1$theta_sampling$rung1)){
    stop("Parameter1 name not recognised")
  }
  if(!parameter2 %in% names(x$chain1$theta_sampling$rung1)){
    stop("Parameter2 name not recognised")
  }
  
  data <- x$chain1$theta_sampling$rung1
  data <- data[ ,c(parameter1, parameter2)]
  colnames(data) <- c("x", "y")
  
  if(downsample & nrow(data) > 5000){
    set.seed(1)
    data <- data[sample(nrow(data),5000),]
  }
  
  ggplot2::ggplot(data = data,
                  ggplot2::aes(x = .data$x, y = .data$y)) + 
    ggplot2::geom_point(alpha = 0.5, col = "darkblue") + 
    ggplot2::xlab(parameter1) +
    ggplot2::ylab(parameter2) +
    ggplot2::theme_bw()
  
}

