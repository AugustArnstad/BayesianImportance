# BayesianImportance
**Bayesian Variable Importance for GLMMs using INLA**
This ia a package developed for my project thesis at NTNU. The project is a 15ECTS project that is followed up by a master thesis the next semester, and therefore lays out the foundation. Tentatively my thesis will be on analyzing heritability traits from sparrows on Helgelandskysten using bayesian statistics.

`BayesianImportance` is an R package designed to compute Bayesian variable importance metrics for Generalized Linear Mixed Models (GLMMs) utilizing the Integrated Nested Laplace Approximation (INLA) methodology.

## Features
- **Bayesian Variable Importance Computation**: Allows for the quantification of the importance of predictors in GLMMs in a Bayesian framework. Currently it only works with LMM's but this is thought to be extended in the near future.
- **INLA Integration**: Leverages the computational advantages of INLA, a popular method for Bayesian inference for latent Gaussian models.
- **Support for Various GLMMs**: Compatible with a wide range of generalized linear mixed models.
- **Extensible**: Designed with the modern R user in mind, offering a range of utilities to further expand upon the base functionality.

## Installation
To install the latest version of `BayesianImportance` from GitHub, use the following command:
```R
# If not already installed, install the 'devtools' package
if(!require(devtools)) install.packages("devtools")

# Install BayesianImportance
devtools::install_github("AugustArnstad/BayesianImportance")
``` 

## Usage
To compute the Bayesian variable importance for your GLMMs, follow the basic structure:

```{r}
set.seed(1234)
model <- run_bayesian_imp(data_bayes, Y ~ V2 + V3 + (1 | gamma) + (1 | eta))

plot_model = plot_posteriors(model, importance=FALSE)
plot_model$posterior_plot

plot_model = plot_posteriors(model, importance=TRUE)
plot_model$posterior_plot

gelman_r2 = gelman_r2_metrics(model, s=1000, plot=TRUE)
gelman_r2$plot
summary(gelman_r2$conditional_gelman_r2)
```
where all functions are found in the R folder.
Detailed examples and tutorials will be made available in future releases.

## Simulation study
In the folder simulation study, we have four files, that contribute to a simulation study where the Bayesian Importance method is compared to other, frequentist and more established, methods in the field of mathematics. Simulation study preparation.Rmd and Simulation study.Rmd are drafts mostly made for how one should do the simulation study and can be viewed as redundant. Simulation run.Rmd contains the code that runs the simulation study and writes the results to the attached csv files. Simulation study analysis.Rmd contains analysis of the resulting files, which is done by violin plots and tables to compare the Bayesian Importance package with other methods that are established in the mathematical field.

## Documentation
Further documentation and function references can be found within the package. Use the standard R help and documentation commands to access detailed information about each function.

## Contributing
Contributions to `BayesianImportance` are welcome. Please ensure that you adhere to standard coding practices and include tests for any new features. Open a pull request with details about your changes.

## License
This project is licensed under the MIT License - see the LICENSE.txt file for details

## Acknowledgements
INLA team for their outstanding work on the INLA methodology and R package.
My counsellor Stefanie Muff, Associate Professor, Department of Mathematical Sciences, NTNU Trondheim

