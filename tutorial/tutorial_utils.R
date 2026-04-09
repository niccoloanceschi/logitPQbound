
accuracy <- function(sim, mu, sigma2) {
  n <- length(sim)
  K <- 201
  idx <- seq(from=floor(n/2), to=n, by=5)
  kde <- density(sim[idx])
  lower <- min(kde$x)
  upper <- max(kde$x)
  x <- seq(from=lower, to=upper, length=K)
  f1 <- approxfun(kde$x, kde$y, rule=2)(x)
  f2 <- dnorm(x, mean=mu, sd=sqrt(sigma2))
  dx <- diff(x)
  df <- abs(f1 - f2)
  err <- 0.5*sum((df[-1]+df[-K])*dx)
  acc <- 1 - 0.5*err
  return(acc)
}

panel_cor <- function(x, y, ...){
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=2)[1]
  cex <- 0.3/strwidth(txt)
  text(0.5, 0.5, txt, cex=cex*r, col=2)
  text(0.5, 0.65, "Corr", cex=cex*r/2, col=1)
}

panel_linear <- function(x, y, ...) {
  # usr <- par("usr")
  # on.exit(par(usr))
  points(x, y, pch=20, col=1)
  abline(a=0, b=1, col=2)
}

plot_1run_loglik <- function(fit_list, fit_summary, COLORS, MARKERS){
  
  layout(matrix(1:3, nrow = 1), widths = c(2, 1, 1))
  # Log-likelihood
  niter <- max(sapply(fit_list, \(.) .$niter))
  objval <- matrix(NA, nrow=niter, ncol=length(fit_list))
  colnames(objval) <- names(fit_list)
  for (k in 1:length(fit_list)) {
    objval[1:fit_list[[k]]$niter, k] <- fit_list[[k]]$trace$objval
  }
  matplot(objval[2:min(50, niter),], type="b", lty=1, col=COLORS, pch=MARKERS, xlab="", ylab="")
  title(xlab="Iteration", ylab="Log-Likelihood", main="Penalized Log-Likelihood")
  legend("bottomright", col=COLORS, pch=MARKERS, legend=colnames(objval))
  # Iterations
  with(df_1run_summary, {
    barplt <- barplot(niter, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), niter)
    title(ylab="Iterations", main="Number of Iterations")
  })
  # Exetime
  with(df_1run_summary, {
    barplt <- barplot(exetime, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
  })
}

plot_1run_pairs <- function(fit_list){
  betas <- sapply(fit_list, \(.) .$beta)
  pairs(betas[-1,], lower.panel=panel_cor, upper.panel=panel_linear)
}

plot_path_timegain <- function(path_summary, COLORS, MARKERS){
  with(path_summary, {
    par(mfrow=c(1,3))
    # Total number of iterations
    barplt <- barplot(niter, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), round(niter, 2))
    title(ylab="Iterations", main="Total number of iterations")
    # Total execution time
    barplt <- barplot(exetime, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Total execution Time")
    # Relative time gain
    reltime <- 100*(1-exetime[3]/exetime)
    barplt <- barplot(reltime, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, reltime-0.1*sign(reltime)*max(reltime), paste(floor(reltime), "%"))
    title(ylab="Time Gain", main="Time Gain")
    par(mfrow=c(1,1))
  })
}


plot_fit_path <- function (fit_list, field, D=NULL, log=FALSE, main="", position="bottomright", ...) {
  
  # Check the input field
  fields <- c("dev", "loglik", "logdet", "reml", "norm", "pnorm", "niter", "exetime")
  field  <- match.arg(field, fields)
  
  # Get the lambda grid
  loglambdas <- log10(fit_list[[1]]$lambdas)
  
  # Extract the specified field from each list 
  mat <- fit_list %>% lapply(\(.) .[[field]]) %>% do.call("cbind", .)
  
  # Set the penalty matrix
  if (is.null(D)) D <- diag(p)
  
  # If a "norm" field is required, compute it using the estimated coefficients
  if (field== "norm") {
    mat <- fit_list %>% lapply(\(.) {
      colMeans(.[["beta"]]^2)
    }) %>% do.call("cbind", .)
  }
  
  if (field=="pnorm") {
    mat <- fit_list %>% lapply(\(.) {
      alpha <- .[["alpha"]]
      r <- abs(D %*% .[["beta"]])
      pL1 <- colMeans(r)
      pL2 <- colMeans(r*r)/2
      return(alpha*pL1+(1-alpha)*pL2)
    }) %>% do.call("cbind", .)
  }
  
  # If specified, transform the results in the logarithmic scale
  if (log) mat <- log10(mat)
  
  # Plot the results
  matplot(loglambdas, mat, type="o", xlab="", ylab="", main="", ...)
  title(xlab=expression(log[10](lambda)), ylab=field, main=main)
  legend(position, legend=names(fit_list), ...)
  
  # If required, plot vertical lines identifying the REML/GCV/AIC/BIC solution
  if (field=="reml") {
    abline(v=loglambdas[apply(mat, 2, which.max)], lty=2, col=8)
  }
  if (field %in% c("gcv", "aic", "bic")) {
    abline(v=loglambdas[apply(mat, 2, which.min)], lty=2, col=8)
  }
}


plot_mesh <- function (mesh, incol="grey70", bndcol="grey30", ...) {
  # Set the x-y limits
  xrng <- range(mesh$nodes[, 1])
  yrng <- range(mesh$nodes[, 2])
  
  # Plot the mesh nodes
  plot(mesh$nodes, cex = 0.001, col=col, xlim=xrng, ylim=yrng, 
       xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
  
  # Plot the mesh internal segments
  x_in_start <- mesh$nodes[mesh$edges[, 1], 1]
  y_in_start <- mesh$nodes[mesh$edges[, 2], 1]
  x_in_end <- mesh$nodes[mesh$edges[, 1], 2]
  y_in_end <- mesh$nodes[mesh$edges[, 2], 2]
  segments(x_in_start, x_in_end, y_in_start, y_in_end, col=incol)
  
  # Plot the mesh boundary segments
  x_bnd_start <- mesh$nodes[mesh$segments[, 1], 1]
  y_bnd_start <- mesh$nodes[mesh$segments[, 2], 1]
  x_bnd_end <- mesh$nodes[mesh$segments[, 1], 2]
  y_bnd_end <- mesh$nodes[mesh$segments[, 2], 2]
  segments(x_bnd_start, x_bnd_end, y_bnd_start, y_bnd_end, col=bndcol, lwd=1.5)
}

plot_field <- function(lon, lat, field, ...) {
  plot3D::image2D(x=lon, y=lat, z=field, colvar=field, 
                  col=viridis::inferno(100), xlab="", ylab="", ...)
}


ggplot_field <- function(lon, lat, field, ngrid, locs=NULL, 
                              fun=NULL, palette="viridis", ...) {
  H <- ncol(field)
  
  # Set the x-y limits
  xrng <- range(lon)
  yrng <- range(lat)
  
  # Set the x-y grid
  xs <- lon
  ys <- lat
  
  # Set the x-y lattice
  xx <- rep(xs, ngrid)
  yy <- rep(ys, rep(ngrid, ngrid))
  
  # Set the data-frames
  df_fems <- data.frame(
    g = rep(colnames(field), each=ngrid^2),
    x = rep(xx, times=H),
    y = rep(yy, times=H),
    z = as.vector(field))
  
  df_txt <- data.frame(
    g = colnames(field),
    x = rep(xrng[2]-.225*diff(xrng), H),
    y = rep(yrng[2]-.025*diff(yrng), H),
    t = paste0("TV = ", round(colMeans(field, na.rm=TRUE), 4)))
  
  df_locs <- as.data.frame(locs)
  df_locs$z <- factor(df_locs$z, labels=c("safe","risky"))
  
  # Set the scaling transofrmation
  if (is.null(fun)) fun <- function(t) t
  df_fems$z <- with(df_fems, fun(z)*max(z, na.rm=TRUE)/max(fun(z), na.rm=TRUE))
  
  # Create the data-frame
  plt <- ggplot() + theme_bw() +
    geom_raster(data=df_fems, map=aes(x=x, y=y, fill=z)) + 
    geom_label(data=df_txt, map=aes(x=x, y=y, label=t), size=6) + 
    facet_grid(cols=vars(g)) +
    scale_shape_manual(values=c(1,19)) +
    scale_fill_viridis_c(option=palette, limits=c(0,max(field)), na.value="transparent") +
    labs(x="Longitude", y="Latitude", shape="Response", fill=expression(TV(q[VB],p[MC])))
  
  # Plot the field
  return(plt)
}





