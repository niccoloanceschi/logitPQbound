
## PACKAGE IMPORT ----

library(ggplot2)
library(dplyr)
library(sf)
library(osmdata)

## GLOBAL VARIABLES ----

SHOW <- TRUE
SAVE <- FALSE
DATAPATH <- "data/Portland"
IMGPATH <- "img/portland_DATA"

## PORTLAND DATA ----

### Map and districts ----

# Set the latitude-longitude box and the coordinate reference system
crs <- sf::st_crs(4326)
lims <- c(xmin = -122.85, xmax = -122.4, ymin = 45.44, ymax = 45.7)
bbox <- sf::st_bbox(lims, crs = crs)

# Portland police district data path
zippath <- paste(DATAPATH, "portland-police-districts.zip", sep="/")
shapefile <- "Portland_Police_Districts.shp"
dns <- paste("/vsizip", zippath, shapefile, sep="/")
print(unzip(zippath, list = TRUE))

# Portland police district data import
# IMPORTANT: district- and crime-shape data have different coordinate 
#            systems which must be uniformed by setting a common crs!!
portland_districts <- sf::st_read(dsn=dns) %>% sf::st_transform(crs=crs)
str(portland_districts)

if (FALSE) {
  par(mfrow = c(3,3))
  plot(portland_districts["OBJECTID_1"])
  plot(portland_districts["OBJECTID"])
  plot(portland_districts["AGENCY"])
  plot(portland_districts["DISTRICT"])
  plot(portland_districts["Shape_Leng"])
  plot(portland_districts["Shape_STAr"])
  plot(portland_districts["Shape_STLe"])
  plot(portland_districts["Shape_Area"])
  plot(portland_districts["geometry"])
  par(mfrow = c(1,1))
}

### Crime data in 2015 ----

# Portland crime data path
zippath <- "data/Portland/010115_123115_Data.zip"
shapefile <- "NIJ2015_JAN01_DEC31.shp"
dns <- paste("/vsizip", zippath, shapefile, sep="/")
print(unzip(zippath, list = TRUE))

# Portland crime data import
# IMPORTANT: district- and crime-shape data have different coordinate 
#            systems which must be uniformed by setting a common crs!!
portland_crimes <- sf::st_read(dsn=dns) %>% sf::st_transform(crs=crs)
str(portland_crimes)

if (FALSE) {
  par(mfrow = c(3,3))
  plot(portland_crimes["CATEGORY"])
  plot(portland_crimes["CALL_GROUP"])
  plot(portland_crimes["final_case"])
  plot(portland_crimes["CASE_DESC"])
  plot(portland_crimes["occ_date"])
  plot(portland_crimes["x_coordina"])
  plot(portland_crimes["y_coordina"])
  plot(portland_crimes["census_tra"])
  par(mfrow = c(1,1))
}

# Portland motor vehicle theft data
portland_thefts <- portland_crimes %>% dplyr::filter(CATEGORY=="MOTOR VEHICLE THEFT")
str(portland_thefts)

if (SHOW) {
  plot(portland_districts["geometry"])
  plot(portland_thefts["CATEGORY"], pch = 19, cex = 0.5, add=TRUE)
}

### Grid aggregation ----

# Given an interval obtained with the cut() function, 
# center() finds the central point of the interval
center <- function(interval) {
  limits <- as.numeric(unlist(strsplit(gsub("\\[|\\)|\\]", "", interval), ",")))
  return(mean(limits))
}

# Discretize the Portland region into a regular squared grid
nbin <- 50
coord <- sf::st_coordinates(portland_thefts)
xlim <- c(-122.82, -122.45)
ylim <- c(45.44, 45.625)
xbreaks <- seq(xlim[1], xlim[2], length=nbin+1)
ybreaks <- seq(ylim[1], ylim[2], length=nbin+1)
xbins <- cut(coord[,1], breaks = xbreaks, include.lowest=TRUE, right=FALSE)
ybins <- cut(coord[,2], breaks = ybreaks, include.lowest=TRUE, right=FALSE)
xcenter <- sapply(as.character(xbins), center)
ycenter <- sapply(as.character(ybins), center)

if (SHOW) {
  plot(portland_districts["geometry"])
  points(coord, xlim=xlim, ylim=ylim, col=4, pch=19, cex=0.5)
  abline(v=xbreaks, h=ybreaks, col=8)
  axis(1)
  axis(2)
  box()
}

# Count how many thefts fall into each square of the grid
counts <- as.data.frame(table(xbins, ybins))
counts$xcenter <- sapply(as.character(counts$xbins), center)
counts$ycenter <- sapply(as.character(counts$ybins), center)
counts <- counts %>% dplyr::select(xbins, ybins, xcenter, ycenter, Freq)
str(counts)

if (SAVE) {
  filename <- "portland_count_distribution.pdf"
  height <- 6; width <- 10; zoom <- 2
  pdf(file=paste(IMGPATH, filename, sep="/"), height=zoom*height, width=zoom*width)
  par(mfrow=c(2,3))
  # Count maps
  table(xbins, ybins) %>%             plot3D::image2D(main="Original scale")
  table(xbins, ybins) %>%  sqrt() %>% plot3D::image2D(main="Sqrt scale")
  table(xbins, ybins) %>% log1p() %>% plot3D::image2D(main="Log scale")
  # Count histograms
  table(xbins, ybins) %>%             hist(nclass=10, main="Original scale")
  table(xbins, ybins) %>%  sqrt() %>% hist(nclass=10, main="Sqrt scale")
  table(xbins, ybins) %>% log1p() %>% hist(nclass=10, main="Log scale")
  par(mfrow=c(1,1))
  dev.off()
}

if (SAVE) {
  filename <- "portland_data_distribution.pdf"
  height <- 4; width <- 10; zoom <- 2
  pdf(file=paste(IMGPATH, filename, sep="/"), height=zoom*height, width=zoom*width)
  par(mfrow=c(1,2))
  # Spatial distribution
  plot(0, 0, type="n", axes=FALSE, xlim=lims[1:2], ylim=lims[3:4], xlab="", ylab="")
  plot(portland_districts["geometry"], add=TRUE)
  counts %>% dplyr::filter(Freq > 0) %>% dplyr::select(xcenter, ycenter) %>% points(pch=19, col="white")
  counts %>% dplyr::filter(Freq > 1) %>% dplyr::select(xcenter, ycenter) %>% points(pch=19, col="grey80")
  counts %>% dplyr::filter(Freq > 2) %>% dplyr::select(xcenter, ycenter) %>% points(pch=19, col="grey50")
  counts %>% dplyr::filter(Freq > 3) %>% dplyr::select(xcenter, ycenter) %>% points(pch=19, col="grey10")
  counts %>% dplyr::filter(Freq > 0) %>% dplyr::select(xcenter, ycenter) %>% points(pch=1, col="grey1")
  legend("topright", pch=c(1,19,19,19), col=paste0("grey", c(1,80,50,10)), legend=c(">0", ">1", ">2", ">3"))
  title(main="Spatial distribution")
  # Marginal distribution
  counts %>% dplyr::filter(Freq > 0) %>% dplyr::select(Freq) %>% unlist() %>% hist(nclass=20, main="")
  counts %>% dplyr::filter(Freq > 0) %>% dplyr::select(Freq) %>% unlist() %>% quantile(prob=c(.5,.75,.85,.9,.95,.975,.99,1.)) %>% print()
  title(main="Marginal distribution")
  par(mfrow=c(1,1))
  dev.off()
}

# Filter out the zero-count locations and dichotomize 
# the counts greater than the third quartile
counts_non_zero <- counts %>% 
  dplyr::filter(Freq > 0) %>% 
  dplyr::mutate(x = xcenter, y = ycenter) %>% 
  dplyr::mutate(risk = ifelse(Freq>quantile(Freq, .75), 1, 0)) %>% 
  dplyr::select(x, y, risk)

if (SAVE) {
  filename <- "portland_data_map.pdf"
  height <- 3; width <- 6; zoom <- 2.5
  pdf(file=paste(IMGPATH, filename, sep="/"), height=zoom*height, width=zoom*width)
  with(counts_non_zero, {
    par(mfrow=c(1,2))
    # Original locations
    plot(portland_districts["geometry"], reset=FALSE, border="grey60")
    points(coord, pch=18, cex=0.75)
    legend("topright", pch=18, legend=c("thefts"))
    title(main="Portland theft locations")
    # Dicothomized data
    plot(portland_districts["geometry"], reset=FALSE, border="grey60")
    points(x, y, pch=19, col="white")
    points(x[risk==1], y[risk==1], pch=19, cex=0.75)
    points(x[risk==0], y[risk==0], pch= 1, cex=0.75)
    # plot3D::scatter2D(x=x, y=y, colvar=risk, colkey=FALSE, pch=19, cex=0.75, col=c(4,2), add=TRUE)
    legend("topright", pch=c(1,19), col=1, legend=c("safe", "risky"))
    title(main="Dicothomized thefts data")
  })
  dev.off()
}

### Boundary simplification ----

# Join all the districts
postland_simplified <- portland_districts %>% 
  sf::st_make_valid() %>% sf::st_union() %>% sf::st_transform(crs)

# Get the joined Portland boundary 
border <- postland_simplified[[1]][[2]][[1]]
hole <- postland_simplified[[1]][[2]][[8]]

# Set the boundary nodes to keep for boundary simplification
idx <- c(2,3,10,18,35,38,44,48,51,55,74,77,82,97,99,105,107,166,168,172,
         205,209,263,288,290,300,320,325,326,327,344,346,382,459,468,469,
         472,490,500,570,578,579,598,686,698,704,710,714,722,728,736,772,
         778,811,813,815,870,890,893,929,930,932,934,940,946,953,955,956,
         963,973,988,990,991,1004,1009,1011,1012,1023,1026,1036,1041,1048,
         1049,1056,1058,1075,1077,1078,1079,1080,1082,1195,1233, 1243,1248,
         1259,1264, # 1271,1325,1338,1344,1357,1410,1421,
         1425,1430, 1435,1454,1468,1482,1489,1490,1496,1532,1533,1535,1550,
         1560,1566,1570,1573,1660,1661,1672,1701,1739,1775,1779,1785,1786,
         1794,1799,1816,1818,1819,1820,1821,1822,1824,1825,1826,1827,1832,
         1834,1837,1840,1843,1847,1851,1855,1862,1863,1865,1867,1870,1871,
         1872,1874,1876,1877,1938,1945,1950,1952,1955,1956,1958,1959,1964,
         1967,1968,1976,1977,1978,1979,1980,1981,1983,1998,2009,2015,2032,
         2037,2045,2055,2065,2079,2084,2091,2092,2093,2095)

if (SHOW) {
  plot(border, type="l", xlab="", ylab="")
  lines(hole)
  lines(border[idx,], col=2, lwd=2)
  title(main="Simplified Portland boundary")
}

# plot(border, type="l", xlim=c(-122.52,-122.45), ylim=c(45.45,45.56))
# lines(hole)
# lines(border[idx,], col=2)
# points(border[idx,][115+(1:9),], col=2, pch=19)

## FEM MATRICES CONSTRUCTION ----

# Import the fdaPDE package for mesh construction and FEM utilities
if (!require(fdaPDE, quietly=TRUE)) {
  pkgpath = "C:/Users/crist/Downloads/fdaPDE-master/fdaPDE-master"
  install.packages(pkgpath, type="source", repos=NULL)
}

library(fdaPDE)

### Mesh building ----

plot_mesh <- function (mesh, bnd, incol="grey70", bndcol="grey30", add=FALSE) {
  # Set the x-y limits
  xrng <- range(mesh$nodes[, 1])
  yrng <- range(mesh$nodes[, 2])
  
  # Plot the mesh nodes
  if (add) {
    points(mesh$nodes, cex=0.002, col=incol)
  } else {
    plot(mesh$nodes, cex=0.002, col=incol, 
         xlim=xrng, ylim=yrng, xlab="", ylab="", 
         main="", xaxt="n", yaxt="n", bty="n")
  }

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
  segments(x_bnd_start, x_bnd_end, y_bnd_start, y_bnd_end, col=bndcol, lwd=2)
}

plot_field <- function (coefs, basis, ngrid, ...) {
  # Set the x-y limits
  xrng <- range(basis$mesh$nodes[,1])
  yrng <- range(basis$mesh$nodes[,2])
  
  # Set the x-y grid
  xs <- seq(from=xrng[1], to=xrng[2], length=ngrid)
  ys <- seq(from=yrng[1], to=yrng[2], length=ngrid)
  
  # Set the x-y lattice
  xx <- rep(xvec, ngrid)
  yy <- rep(yvec, rep(ngrid, ngrid))
  
  # Set the FEM object
  fem <- fdaPDE::FEM(coefs, basis)
  
  # Compute the FEM surface
  field <- fdaPDE::eval.FEM(fem, cbind(xx, yy))
  field <- matrix(field, nrow=ngrid, ncol=ngrid)
  
  # Plot the field
  plot3D::image2D(x=xs, y=ys, z=field, colvar=field, xlab="", ylab="", ...)
}

# Set the (simplified) boundary nodes and segments
nodes <- border[idx,]
segments <- cbind(1:length(idx), c(2:length(idx),1))

# Correct the nodes inducing some points to lie outside of the boundary
nodes[118,] <- c(-122.48212, 45.52001)
nodes[119,] <- c(-122.49580, 45.52001)

# plot(border, type="l", xlim=c(-122.52,-122.45), ylim=c(45.45,45.56))
# lines(border[idx,], col=2)
# points(border[idx,][116+(1:4),], col=2, pch=19)
# points(counts_non_zero[,1:2], col=4, pch=19)
# lines(nodes, col=3)
# points(nodes[116+(1:4),], col=3, pch=19)

# Set the spatial locations and data
locs <- counts_non_zero %>% dplyr::select(x, y) %>% as.matrix()
obs <- as.vector(counts_non_zero$risk)*2-1

# Create and refine the mesh
mesh <- fdaPDE::create.mesh.2D(nodes=nodes, segments=segments)
mesh <- fdaPDE::refine.mesh.2D(mesh, maximum_area=0.0000125, delaunay=TRUE)

if (SAVE) {
  filename <- "portland_data_mesh.pdf"
  height <- 3; width <- 6; zoom <- 2.5
  pdf(file=paste(IMGPATH, filename, sep="/"), height=zoom*height, width=zoom*width)
  par(mfrow=c(1,2))
  plot(portland_districts["geometry"], border="grey60", reset=FALSE)
  points(locs, col="grey99", pch=19, cex=0.8)
  points(locs, col="grey15", pch=ifelse(obs==1, 19, 1), cex=0.8)
  title("Portland thefts data")
  plot(portland_districts["geometry"], reset=FALSE)
  plot_mesh(mesh, border, incol="grey70", bndcol="grey50", add=TRUE)
  title("Triangular mesh")
  par(mfrow=c(1,1))  
  dev.off()
}


### FEM basis ----

# Create the FEM basis object
basis <- fdaPDE::create.FEM.basis(mesh)

# Compute the mass and stiffness matrices 
mass <- fdaPDE:::CPP_get.FEM.Mass.Matrix(basis)
stiff <- fdaPDE:::CPP_get.FEM.Stiff.Matrix(basis)

# Compute the FEM basis matrix
coeff <- diag(dim(mesh$nodes)[1])
fun <- fdaPDE::FEM(coeff, basis)
psi <- fdaPDE::eval.FEM(fun, locs)
psi <- Matrix::Matrix(psi, sparse=TRUE)

## DATA SAVING ----

if (SAVE) {
  filename <- "PortlandData.RData"
  save(locs, obs, idx, border, hole, nodes, 
       segments, mesh, basis, mass, stiff, psi,
       file=paste(DATAPATH, filename, sep="/"))
}

## END OF FILE ----
