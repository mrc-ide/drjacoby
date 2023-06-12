# return 95% quantile
#' @importFrom stats quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs = c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#' Subset data for plotting
#' @param output_df MCMC output data.frame
#' @param chain Select chain(s)
#' @param phase Select phases can be from: all, tune, burn, sample
#' @noRd
plot_data_subset <- function(output_df, chain, phase){
  data <- output_df
  if(!is.data.frame(data)){
    data <- do.call("rbind", data)
  }
  
  # Phase subset
  if(!is.null(phase)){
    data <- data[data$phase %in% phase, ]
  }
  
  # Chain subset
  if(!is.null(chain)){
    data <- data[data$chain %in% chain, ]
  }
  return(data)
}


create_par_plot <- function(
    par,
    output_df,
    lag,
    downsample,
    phase,
    chain,
    return_elements) {
  
  data <- plot_data_subset(
    output_df = output_df,
    chain = chain,
    phase = phase
  )
  data <- data[,c("chain", "iteration", par)]
  
  # get autocorrelation (on full data, before downsampling)
  ac_data <- as.data.frame(apply(data[, par, drop = FALSE], 2, acf_data, lag = lag))
  ac_data$lag <- 0:lag
  ac_data[,par] <- pmax(ac_data[,par], 0)
  names(ac_data) <- c("lag", "Autocorrelation")
  
  # downsample
  if (downsample & nrow(data) > 2000) {
    data <- data[round(seq(1, nrow(data), length.out = 2000)),]
  }
  
  # set minimum bin number
  b <- min(nrow(data) / 4, 40)
  
  colnames(data) <- c("chain", "iteration", "y")
  
  # Histogram
  p1 <- ggplot2::ggplot(data, ggplot2::aes(x = .data$y)) + 
    ggplot2::geom_histogram(bins = b, fill = "deepskyblue3", col = "darkblue") + 
    ggplot2::ylab("Count") + 
    ggplot2::xlab(par) +
    ggplot2::theme_bw()
  
  # Trace plots
  p2 <- ggplot2::ggplot(data, ggplot2::aes(x = .data$iteration, y = .data$y, col = as.factor(.data$chain))) + 
    ggplot2::geom_line() +
    ggplot2::scale_color_discrete(name = "Chain") +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(par) +
    ggplot2::theme_bw()
  
  # Autocorrealtion
  p3 <- ggplot2::ggplot(data = ac_data,
                        ggplot2::aes(x = .data$Autocorrelation, y = .data$lag)) +
    ggplot2::geom_bar(stat = "identity") + 
    ggplot2::geom_hline(yintercept = 0, lty = 2, col = "red") + 
    ggplot2::theme_bw() +
    ggplot2::ylab("Autocorrelation") +
    ggplot2::xlab("Lag") +
    ggplot2::ylim(0, 1)
  
  if(return_elements){
    out <- list(
      histogram = p1,
      trace = p2,
      autocorrelation = p3
    )
    return(out)
  }
  
  # Arrange
  pc1 <- cowplot::plot_grid(p1, p3, ncol = 2)
  pc2 <- cowplot::plot_grid(p2, pc1, nrow = 2)
  
  return(pc2)
}

create_cor_plot <- function(
    parx,
    pary,
    output_df,
    downsample,
    phase,
    chain)  
{
  
  data <- plot_data_subset(
    output_df = output_df,
    chain = chain,
    phase = phase
  )
  data <- data[,c("chain", "iteration", parx, pary)]
  colnames(data) <- c("chain", "iteration", "x", "y")
  
  # downsample
  if (downsample & nrow(data) > 2000) {
    data <- data[round(seq(1, nrow(data), length.out = 2000)),]
  }
  
  # produce plot
  ggplot2::ggplot(data = data,
                  ggplot2::aes(x = .data$x, y = .data$y, col = as.factor(.data$chain))) + 
    ggplot2::geom_point(alpha = 0.5) + 
    ggplot2::xlab(parx) +
    ggplot2::ylab(pary) +
    ggplot2::scale_color_discrete(name = "Chain") +
    ggplot2::theme_bw()
  
}

create_credible_plot <- function(output_df, pars, chain, phase, param_names) {
  
  data <- plot_data_subset(
    output_df = output_df,
    chain = chain,
    phase = phase
  )
  data <- data[, pars, drop = FALSE]
  
  # get quantiles
  df_plot <- as.data.frame(t(apply(data, 2, quantile_95)))
  df_plot$param <- factor(param_names, levels = param_names)
  
  # produce plot
  ggplot2::ggplot(data = df_plot) +
    ggplot2::geom_point(ggplot2::aes(x = .data$param, y = .data$Q50)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$param, y = .data$Q2.5,
                                       xend = .data$param, yend = .data$Q97.5)) +
    ggplot2::xlab("") +
    ggplot2::ylab("95% CrI") +
    ggplot2::theme_bw() +
    ggplot2::coord_flip()
  
}

create_cor_mat_plot <- function(output_df, pars, chain, phase, param_names) {
  
  data <- plot_data_subset(
    output_df = output_df,
    chain = chain,
    phase = phase
  )
  data <- data[, pars, drop = FALSE]
  n <- ncol(data)
  
  # get correlation matrix into dataframe for ggplot
  m <- stats::cor(data)
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
  ggplot2::ggplot(data = df_plot) + 
    ggplot2::geom_raster(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$z)) +
    ggplot2::scale_fill_gradientn(colours = c("red", "white", "blue"),
                                  values = c(0, 0.5, 1),
                                  limits = c(-max_plot, max_plot),
                                  name = "correlation") +
    ggplot2::xlab("") + 
    ggplot2::ylab("") +
    ggplot2::theme_bw()
  
}

create_mc_acceptance_plot <- function(rungs, beta, ar, x_axis_type) {
  beta_mid <- beta_mid(beta)
  
  # define x-axis type
  if (x_axis_type == 1) {
    breaks_vec <- 1:rungs
    x_vec <- (2:rungs) - 0.5
    x_lab <- "rung"
  } else {
    breaks_vec <- beta
    x_vec <- beta_mid
    x_lab <- "beta"
  }
  
  # get data into ggplot format and define temperature colours
  df <- data.frame(x_vec = x_vec, mc_accept = ar, col = beta_mid)
  
  # produce plot
  mc_ar_plot <- ggplot2::ggplot(data = df) + 
    ggplot2::geom_vline(data = data.frame(x = breaks_vec), ggplot2::aes(xintercept = .data$x), col = "grey90") +
    ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
    ggplot2::geom_point(ggplot2::aes(x = .data$x_vec, y = .data$mc_accept, color = .data$col)) + 
    ggplot2::xlab(x_lab) + 
    ggplot2::ylab("coupling acceptance rate") + 
    ggplot2::scale_colour_gradientn(colours = c("red", "blue"), name = "beta", limits = c(0,1)) +
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())
  
  return(mc_ar_plot)
}

create_rejection_rate_plot <- function(beta, rejection_rate){
  beta_mid <- beta_mid(beta)
  pd <- data.frame(
    beta = rev(beta_mid),
    r = cumsum(rev(rejection_rate))
  )
  ggplot2::ggplot(data = pd, ggplot2::aes(x = .data$beta, y = .data$r, col = .data$beta)) +
    ggplot2::geom_vline(xintercept = beta, col = "grey90") +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::xlab(expression("Schedule" ~ (beta))) +
    ggplot2::ylab(expression("Cumulative rejection rate" ~ (hat("r")))) + 
    ggplot2::scale_colour_gradientn(colours = c("red", "blue"), name = "beta", limits = c(0,1)) +
    ggplot2::ylim(0, NA) +
    ggplot2::theme_bw()
}

create_local_communication_barrier_plot <- function(beta, rejection_rate){
  beta_mid <- beta_mid(beta)
  pd <- data.frame(
    beta = rev(beta_mid),
    r = rev(rejection_rate) / sum(rejection_rate)
  )
  ggplot2::ggplot(data = pd, ggplot2::aes(x = .data$beta, y = .data$r, col = .data$beta)) +
    ggplot2::geom_vline(xintercept = beta, col = "grey90") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::xlab(expression("Schedule" ~ (beta))) +
    ggplot2::ylab(expression("Local communication barrier " ~ over(lambda(beta), Lambda))) +
    ggplot2::scale_colour_gradientn(colours = c("red", "blue"), name = "beta", limits = c(0,1)) +
    ggplot2::ylim(0, NA) +
    ggplot2::theme_bw()
}
