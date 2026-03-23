
## PACKAGE IMPORT ----

library(dplyr)
library(ggplot2)
library(Matrix)
library(sparseinv)

## GLOBAL VARIABLES ----

SHOW <- FALSE
SAVE <- FALSE

# DATAPATH <- "data/Portland"
# SAVEPATH <- "tutorial/results"
# IMGPATH  <- "img/portland_RIDGE"

DATAPATH <- "data/Portland"
SAVEPATH <- "tutorial/results/Portland"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img/portland_RIDGE", sep="/")


pch <- c(15:18)
col <- c(2:4,7)

## ----

list.files(CSVPATH) %>% matrix(ncol=1)

df_path_file <- rbind(
  data.frame(setting="ridge", alpha=0.00, file="portland_ridge_path_summary.csv"),
  data.frame(setting="enet",  alpha=0.05, file="portland_lasso_path_alpha05_summary.csv"),
  data.frame(setting="enet",  alpha=0.20, file="portland_lasso_path_alpha20_summary.csv"),
  data.frame(setting="enet",  alpha=0.35, file="portland_lasso_path_alpha35_summary.csv"),
  data.frame(setting="enet",  alpha=0.50, file="portland_lasso_path_alpha50_summary.csv"),
  data.frame(setting="enet",  alpha=0.65, file="portland_lasso_path_alpha65_summary.csv"),
  data.frame(setting="enet",  alpha=0.80, file="portland_lasso_path_alpha80_summary.csv"),
  data.frame(setting="lasso", alpha=1.00, file="portland_lasso_path_alpha95_summary.csv"))

df_path_list <- lapply(df_path_file$file, function(file) {
  df <- read.csv2(paste(CSVPATH, file, sep="/"))
  df$niter_avg <- df$niter
  df$niter_tot <- df$niter*29
  df$niter     <- NULL
  df$time_gain <- 1 - df$exetime[df$method=="PQ"] / df$exetime
  if (is.null(df$alpha)) df$alpha <- 0.0
  df <- select(df, alpha, method, niter_avg, niter_tot, exetime, time_gain, loglik)
  return(df)
})

df_path_flat <- do.call(rbind, df_path_list) %>% 
  mutate(setting = case_when(
    alpha == 0 ~ "ridge",
    alpha == 1 ~ "lasso",
    TRUE ~ "enet"
  ), .before=alpha)


df_path_flat %>% filter(method != "NR") %>% 
  ggplot(data = ., mapping = aes(x=alpha, y=niter_tot, color=method)) + 
  geom_point() + geom_line() + theme_bw()

df_path_flat %>% filter(method != "NR") %>% 
  ggplot(data = ., mapping = aes(x=alpha, y=exetime, color=method)) + 
  geom_point() + geom_line() + theme_bw()

df_path_flat %>% filter(method != "NR") %>% 
  ggplot(data = ., mapping = aes(x=alpha, y=time_gain, color=method)) + 
  geom_point() + geom_line() + theme_bw()

