# Design, Identification, and Sensitivity Analysis for Patient Preference Trials

## Author Contributions Checklist Form

All data, analyses, and results for the partisan political news experiment
(described in Sections 2 and 7) are in the "analysis" folder.

The "data" subdirectory contains survey results ("Spring2014omnibus_raw.csv"),
along with a codebook ("Spring2014omnibus_codebook.txt").

The "functions" subdirectory contains two scripts that provide functions to
implement the proposed bounds and sensitivity analyses for arbitrary outcomes
("bounds_frechet_function.R") and binary outcomes ("bounds_binary_function.R").
The provided functions take a data frame with the following columns, defined as
in the paper: S, C, A, D, and Y. Users also specify the desired effects to
calculate, specifics of the posterior simulation procedure (e.g., number of
draws), and desired rho values for use in the sensitivity analysis. Further
details of the implementation are given in comments.

A single script ("replicate.R") contains all preprocessing, analysis, and output
of tables and figures. This script does the following:
- reads in the survey data and recodes media outlets from their raw values (Fox,
  MSNBC, Food Network, Discovery Channel) to "pro-attitudinal",
  "counter-attitudinal", and "entertainment" using self-reported partisanship as
  described in the paper
- constructs the two outcome variables, "media sentiment index" and "likely to
  discuss clip with friend"
- reproduces all summary statistics and tables that are reported in paper
- loads the bounds/sensitivity functions, conducts the analysis, and stores the
  resulting objects
- produces all figures in the paper
- reproduces all interpretation/discussion of the results (e.g., the values of
  rho at which the estimated sensitivity interval no longer includes zero).

Bounds/sensitivity results and figures are output to the "results" subdirectory.

Separately, the "simulations" top-level folder contains all data and code for
the results described in Section 8.

The "run_simulations.R" script loads and preprocesses the raw survey data as in
the main analysis. A simulation population is created based on the "choice
divergence" and "outcome divergence" parameters, according to the procedure
described in the paper. For each simulation population of interest, this script
draws 500 samples and conducts analyses based on three separate methods: (1) the
proposed Bayesian method; (2) a variant of the proposed method in which the
estimated bounds are identical, but a nonparametric bootstrap is used for
uncertainty; and (3) an extension of the parametric model of Long, Little, and
Lin (2008), defined in "/functions/lll.R". For each sample and method, estimates
and confidence intervals are stored in the "results" subdirectory. A second
script, "analyze_simulations.R", loads these results and then computes bias and
coverage rates.
