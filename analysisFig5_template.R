#' This script was translated from a MATLAB script. CRAN has a translation file
#' that helped: https://cran.r-project.org/doc/contrib/Hiebeler-matlabR.pdf
#'
#' Running the script in a REPL will print the figure to screen, you need to run
#' it in batch mode to get the figure saved.
#'
#' $ Rscript analysisFigNew-aez.R
#'
library(ggplot2)

#' Geometric mean
#'
#' Compute the numeric mean in log-space to avoid potential numerical issues.
#'
#' @param x numeric vector
#'
geometric_mean <- function(x) {
  n <- length(x)
  return(exp(sum(log(x)) / n))
}

output_csv <- "fig-data.csv"
output_png <- "fig.png"


## ============================================================================
##  Compute geometric mean based metrics for COVID-19
## ============================================================================

# Death delay distribution (Irons 2021)
mdel <- 21
r <- 1 / (1 + 1.1)
p <- mdel / (r + mdel)
# Time scale and weekly bins
TT <- 42 * 4
ids <- seq(from = 6, to = TT, by = 7)
# Infection to death delay CDF at bins and mean
H <- pnbinom(q = ids, size = r, prob = 1 - p)
Geo_H <- geometric_mean(H)

# Reporting delay distribution (Huisman 2020)
FF <- pgeom(q = ids, prob = 1 / (1+10.8))
Geo_FF <- FF[1]

# Under-reporting of deaths (CDC 2021)
sigbnds <- c(1 / 1.34, 1 / 1.29)
M <- 10000
# IFR (Meyerowitz-Katz 2020, Sorensen 2022)
ifrbnds <- c(0.53, 0.82) / 100

# Sample from M trajectories of size TT
Geo_sigma <- numeric(M)
Geo_ifr <- numeric(M)
for (i in 1:M) {
    # Samples of reporting death rates
    psigma = sigbnds[1] + diff(sigbnds) * runif(n = TT)
    Geo_sigma[i] <- geometric_mean(psigma)
    # Uncertainty on ifr
    pifr <- ifrbnds[1] + diff(ifrbnds) * runif(n = TT)
    Geo_ifr[i] <- geometric_mean(pifr)
}

# Under-reporting of cases
Geo_rho <- numeric(M)
Geo_rhoErr <- numeric(M)
# Reporting fractions (Pullano 2021)
rhobnds <- c(0.06, 0.08)

# Fractions from CDC and Russell 2020
for (i in 1:M) {
    # Samples of reporting case rates
    prho <- rhobnds[1] + diff(rhobnds) * runif(n = TT)
    Geo_rho[i] <- geometric_mean(prho)
}

# Derive ordering on cases versus deaths
caseInfoCOVID <- Geo_rho * Geo_FF;
deathInfoCOVID <- Geo_sigma * Geo_ifr * Geo_H;


## ============================================================================
##  Compute geometric mean based metrics for Ebola virus
## ============================================================================

# Death delay distribution (11.4 to onset + 10 to death) 
mdel <- 21.4
r <- 1.5
p <- mdel / (r + mdel)
# Time scale and weekly bins
TT <- 42 * 4
ids <- seq(from = 6, to = TT, by = 7)
# Infection to death delay CDF at bins
H <- pnbinom(q = ids, size = r, prob = 1 - p)
# Geometric mean from delay
Geo_H <- geometric_mean(H)

# Reporting delay distribution (Huisman 2020) but maximised
FF <- pgeom(q = ids, prob = 1 / 10.8)
Geo_FF <- 1

# Under-reporting of deaths (maximised)
sigbnds <- c(1, 1)
M <- 10000
# IFR which is also a reporting rate (WHO 2014)
ifrbnds <- c(0.69, 0.73);
# Sample from M trajectories of size TT
Geo_sigma <- numeric(M)
Geo_ifr <- numeric(M)
for (i in 1:M) {
    # Samples of reporting death rates
    psigma <- sigbnds[1] + diff(sigbnds) * runif(n = TT)
    Geo_sigma[i] <- geometric_mean(psigma)
    # Uncertainty on ifr
    pifr <- ifrbnds[1] + diff(ifrbnds) * runif(n = TT)
    Geo_ifr[i] <- geometric_mean(pifr)
}

# Under-reporting of cases (Dalziel 2018) 
rhobnds <- c(0.33, 0.83)
Geo_rho <- numeric(M)
for (i in 1:M) {
    # Samples of reporting case rates
    prho <- rhobnds[1] + diff(rhobnds) * runif(n = TT)
    Geo_rho[i] <- geometric_mean(prho)
}

# Derive ordering on cases versus deaths
caseInfoEVD <- Geo_rho * Geo_FF
deathInfoEVD <- Geo_sigma * Geo_ifr * Geo_H
# Delays vs reporting
Geo_RHS <- Geo_H / Geo_FF
Geo_LHS <- Geo_rho / (Geo_ifr * Geo_sigma)


#' ============================================================================
#' Make a final figure for presentation (unlike Matlab version no logs taken).
#' ============================================================================

plot_df <- data.frame(
  info_ratio = c(caseInfoCOVID / deathInfoCOVID, caseInfoEVD / deathInfoEVD),
  disease = rep(c("COVID-19", "EVD"), each = M)
)

## We should also write the data used for the figure to a file so that it is
## easy to tweak styling in the future if the need arises.
write.table(x = plot_df,
            file = output_csv,
            sep = ",",
            row.names = FALSE)

## There is a complete list of build-in themes available here:
## https://ggplot2.tidyverse.org/reference/ggtheme.html
fig <- ggplot(data = plot_df,
              mapping = aes(x = info_ratio, y = ..density..)
              ) +
  geom_histogram(
    bins = 20
  ) +
  facet_wrap(~disease, scales = "free") +
  labs(
    x = "Information ratio (confirmed cases vs. reported deaths)",
    y = NULL
  ) +
  theme_classic()

if (interactive()) {
  print(fig)
} else {
  ggsave(
    filename = output_png,
    plot = fig,
    height = 10.5, width = 14.8,
    ## A5 height = 14.8, width = 21.0,
    ## A6 height = 10.5, width = 14.8,
    ## A7 height = 7.4, width = 10.5,
    units = "cm"
  )
}
