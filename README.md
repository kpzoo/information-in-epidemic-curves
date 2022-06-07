# information-in-epidemic-curves
Code to reproduce the simulated analyses and figures from the paper "Quantifying the information in noisy epidemic curves" at https://www.medrxiv.org/content/10.1101/2022.05.16.22275147v1.

- Defines metrics for comparing the influences of under-reporting and reporting delays.
- Provides insight into the relative reliability (theoretically) of case reports vs death counts.
- Based on Fisher information theory and experimental design.
- Template for performing generic analyses also provided in R.
- Slides: https://www.researchgate.net/publication/360845304_Quantifying_the_relative_information_in_noisy_epidemic_time_series


System Requirements
- Should work with Matlab v2020a and above. Tested on macOS v11.6.4.
- Slight dependence on the linespecer package (license included in main folder).
- R file tested on R version 3.6.2 using RStudio version 1.2.1335.

Instructions and installation
- Should work with any standard Matlab installation and generate figures/results.
- Run analyisFigX.m where X is the figure in the manuscript to be reproduced.
- File analysisFig5_template.m can be easily modified for generic analyses and also given in R.
- Code is self contained so no external installations required.
- No empirical data used, only self contained simulations and computations performed.
- Run times of all scripts are of the order of minutes or faster.
