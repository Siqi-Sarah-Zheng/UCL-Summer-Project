# Copula-Marginal Variational Inference for ODE-Based Survival Models

### Project Overview
This repository contains the code and analysis for a summer research project exploring the application of **Variational Inference (VI)** to **ODE-based survival models**. Our work provides a more computationally efficient alternative to MCMC，which is a highly scalable and accurate.

### Background
Bayesian inference for complex models, such as ODE-based survival models, often results in an intractable posterior distribution. While **Markov Chain Monte Carlo (MCMC)** provides a robust solution by drawing samples, it suffers from significant computational costs. This project explores **Variational Inference (VI)** as a computationally efficient alternative that reframes the problem as an optimisation task.

### Methodology
We propose a novel **Copula-Marginal VI (CMVI)** framework. The key innovation of this method is its flexible compositional structure, which allows for the combination of any valid copula with any marginal family.

* **Model Components**: We use a **Gaussian Copula** to capture the dependence structure between parameters and employ **Two-piece Normal marginals** to model the individual parameter distributions. This choice allows for asymmetry and heavy tails, which are not captured by simpler approximations.
* **Optimization**: We optimize an **importance-weighted Evidence Lower Bound (ELBO)** to obtain a tighter bound on the marginal likelihood.
* **Post-processing**: After optimization, we apply a **self-normalized importance sampling (SNIS)** step to improve the quality of the posterior approximation, particularly in the tails. This step helps to correct for residual bias and provides more accurate estimates for posterior expectations.

### Repository Structure
- `Copula_Marginal_VI.R`: The core R script for implementing the CMVI method.
- `MCMC_vs_NormalApp.R`: Script for comparing MCMC results with a simple normal approximation.
- `VI.indep.R`: Contains code for a variational inference approach with independent parameters.
- `VI.Rmd`: R Markdown file that generates the report (report.html).
- `Rotterdam.R`: Data processing or specific analysis script related to the Rotterdam dataset.
- Other `.R` files: Contain auxiliary functions and model definitions.
