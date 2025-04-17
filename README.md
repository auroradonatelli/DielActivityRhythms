# DielActivityRhythms
Modeling diel activity rhythms, based on hourly movement rates, of brown bear (*Ursus arctos*) populations in relation to anthropogenic, environmental, and climatic variables.

To account for the cyclical and autocorrelated nature of diel activity, the correlation between activity values at different hours of the day was considered in the covariance matrix. Specifically, we assumed that longer time periods between two observations would imply lower correlation between the corresponding activity values. This structure resolves a key challenge in estimation of activity patterns, which rarely accounts for autocorrelation of activity metrics 

To run the model you can use data available at https://doi.org/10.5061/dryad.d51c5b0fd, associated with "Donatelli, A., Ćirović, D., Haroldson, M. A., Huber, D., Kindberg, J., Kojola, I., Kusak, J., Mastrantonio, G., Ordiz, A., Reljić, S., Santini, L., van Manen, F. T. & Ciucci, P. The diel niche of brown bears: Constraints on adaptive capacity in human-modified landscapes. Ecography". Please refer to the aforementioned publication for details and context.

**Software and Package Requirements:**

The versions listed here are the versions used at publication, but future or past versions may be usable.

- R version 4.2.3

- Stan version 2.32.2

- rstan version 2.32.5

- tidyverse version 2.0.0
