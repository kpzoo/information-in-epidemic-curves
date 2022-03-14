# information-in-epidemic-curves
Code to reproduce analyses and figures from the paper "Quantifying the information in noisy epidemic curves" (will be in preprint form soon).
- Defines metrics for comparing the influences of under-reporting and reporting delays.
- Provides insight into the relative reliability (theoretically) of case reports vs death counts.
- Based on Fisher information theory and experimental design.


System Requirements
- Should work with Matlab v2020a and above. Tested on macOS v11.6.4.
- Slight dependence on the linespecer package (license included in main folder).

Instructions and installation
- Should work with any standard Matlab installation and generate figures/results.
- Run analyisFigX.m where X is the figure in the manuscript to be reproduced.
- Code is self contained so no external installations required.
- Run times of all scripts are of the order of minutes or faster.
