# Copula-Marginal Variational Inference for ODE-Based Survival Models

### Project Overview
This repository contains the code and analysis for a summer research project exploring the application of Variational Inference (VI) to ODE-based survival models. Our work provides a more computationally efficient alternative to MCMC.

### Methodology
We propose a Copula-Marginal Variational Inference (CMVI) framework. This method has a very flexible compositional structure, which allows for the combination of any valid copula with any marginal family.

* **Model Components**: We use a Gaussian Copula to capture the dependence structure between parameters and employ Two-piece Normal marginals to model the individual parameter distributions. This choice allows for asymmetry and heavy tails, which are not captured by simpler approximations.
* **Optimisation**: We optimise an importance-weighted Evidence Lower Bound (ELBO) to obtain a tighter bound on the marginal likelihood.
* **Post-optimisation**: After optimisation, we apply a self-normalized importance sampling (SNIS) step to improve the quality of the posterior approximation, particularly in the tails. This step provides more accurate estimates for posterior expectations.

### Repository Structure
- `Copula_Marginal_VI.R`: The core R script for implementing the CMVI method.
- `MCMC_vs_NormalApp.R`: Script for comparing MCMC results with a normal approximation.
- `Least_Square_HR.R`: Contains code for least square density fitting
- `VI.indep.R`: Contains code for a variational inference approach with independent two-piece normals.
- `VI.Rmd`: R Markdown file that generates the report (report.html).
- `routines.R`: ODE-based survival model components and functions.
- `tpapp_MLE.R`: Maximum likelihood estimation (MLE) for the two-piece parametric model using MCMC subsample

### Acknowledgements
This project was funded by the Summer Research Studentship from UCL Department of Mathematics. I am deeply thankful to my supervisor, Dr. F. Javier Rubio, for his generous support, constructive feedback, and the amount of time he dedicated throughout the project. Working on this project has been a valuable opportunity for intellectual and personal growth.
