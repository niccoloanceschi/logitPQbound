
get_1run_summary <- function(fitlist, alpha){
  
  exetime_PQ <- fitlist$PQ$exetime
  
  fitsum <- data.frame(
    alpha = rep(alpha, length(fitlist)),
    method = names(fitlist),
    niter = sapply(fitlist, \(.) sum(.$niter)),
    exetime = sapply(fitlist, \(.) round(.$exetime, 4)),
    timeratio = sapply(fitlist, \(.) round(.$exetime/exetime_PQ, 4)),
    loglik = sapply(fitlist, \(.) round(.$loglik, 4)),
    row.names = seq(length(fitlist)))
  
  return(fitsum)
}

get_1run_coeff <- function(fitlist, lambdas){
  K <- length(fitlist)
  fitcoeff <- matrix(NA, nrow=p, ncol=K)
  for(k in 1:K){
    fitcoeff[,k] <- as.vector(fitlist[[k]]$beta)
  }
  dimnames(fitcoeff) <- list(beta = 1:p, method = names(fitlist))
  return(fitcoeff)
}

get_path_summary <- function(fitlist, alpha){
  
  exetime_PQ <- fitlist$PQ$tottime
  
  fitsum <- data.frame(
    alpha = rep(alpha, length(fitlist)),
    method = names(fitlist),
    niter = sapply(fitlist, \(.) sum(.$niter)),
    exetime = sapply(fitlist, \(.) round(.$tottime, 4)),
    timeratio = sapply(fitlist, \(.) round(.$tottime/exetime_PQ, 4)),
    loglik = sapply(fitlist, \(.) round(mean(.$loglik), 4)),
    row.names = seq(length(fitlist)))
  
  return(fitsum)
}

get_path_coeff <- function(fitlist, lambdas){
  L <- length(lambdas)
  K <- length(fitlist)
  fitcoeff <- array(NA, dim = c(p, L, K))
  for(k in 1:K){
    fitcoeff[,,k] <- fitlist[[k]]$beta
  }
  dimnames(fitcoeff) <- list(beta = 1:p,
                             lambda = 1:L,
                             method = names(fitlist))
  return(fitcoeff)
}

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
  
  # Set the scaling transformation
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





