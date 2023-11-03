# missing-data-approaches
## Investigation of the Performance of Missing Data Approaches Across Propensity Score Matching Methods

### Background
There have been extensive studies on how to handle missing data in propensity scores. Many discussed the use of inverse probability-weighted estimators while the propensity score matching methods were left unexplored.
### Objective
The goal of the study is to explore how different missing data imputation approaches perform with different propensity score matching methods.
### Methods
A 2x2x2 factorial study design was considered for data simulation scenarios (small vs. large sample size, correlation, and missingness). For each simulation, four missing data imputation approaches (complete case analysis, missing pattern, Across approach, and Within Approach) and five variations of propensity score matching methods (optimal matching, nearest neighbor matching with or without replacement and with or without caliper) were included. We used a series of 50 Monte Carlo Simulations to estimate the average treatment effect among the treated (ATT). The estimand was summarized using bias, standard error (SE), 95% confidence interval, Monte Carlo standard deviation (MCSD), and standardized mean differences (SMD).
