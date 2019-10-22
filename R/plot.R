
#------------------------------------------------
#' @title Plot Metropolis Coupling Acceptance
#'
#' @description Plot Metropolis Coupling Acceptance.
#'
#' @param x an object of class \code{drjacoby_output}
#' @param chain which chain to plot.
#' @param phase which phase to plot. Must be either "burnin" or "sampling".
#'
#' @import ggplot2
#' @importFrom grDevices grey
#' @export

plot_mc_acceptance <- function(x, chain = 1, phase = "sampling") {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_single_pos_int(chain)
  assert_leq(chain, length(x))
  assert_in(phase, c("burnin", "sampling"))
  
  # get values
  beta_raised <- x[[chain]]$diagnostics$beta_raised
  beta_raised_mid <- beta_raised[-1] - diff(beta_raised)/2
  if (phase == "burnin") {
    mc_accept <- x[[chain]]$diagnostics$mc_accept_burnin
  } else {
    mc_accept <- x[[chain]]$diagnostics$mc_accept_sampling
  }
  
  # produce plot
  plot1 <- ggplot() + theme_bw() + theme(panel.grid.minor.x = element_blank(),
                                         panel.grid.major.x = element_blank())
  plot1 <- plot1 + geom_vline(aes(xintercept = beta_raised), col = grey(0.9))
  plot1 <- plot1 + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1), expand = c(0,0))
  plot1 <- plot1 + geom_point(aes(x = beta_raised_mid, y = mc_accept), col = "red")
  plot1 <- plot1 + xlab("thermodynamic power") + ylab("coupling acceptance rate")
  
  return(plot1)
}

#' @title Plot autocorrelation
#'
#' @description Plot autocorrelation for specified parameters
#'
#' @inheritParams plot_mc_acceptance
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

#' @title Plot parameter estimates
#'
#' @description Plot parameter estimates
#'
#' @inheritParams plot_mc_acceptance
#' @param show Optional character (or vector of characters) to filter parameters by.
#'  Parameters matching show will be included.
#' @param hide Optional character (or vector of characters) to filter parameters by.
#'  Parameters matching show will be hidden.
#'  @param lag Maximum lag. Must be an integer between 20 and 500
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
  all_chains <- all_chains %>%
    dplyr::group_by(chain) %>%
    dplyr::mutate(x = 1:dplyr::n()) %>%
    dplyr::ungroup()
  
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

#' @title Plot parameter correlation
#'
#' @description Plots correlation between two parameters
#'
#' @inheritParams plot_mc_acceptance
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

